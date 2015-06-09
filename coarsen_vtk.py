#!/usr/bin/python

from os import path, listdir
from dolfin import Mesh, FunctionSpace, CellSize, project, plot
from argparse import ArgumentParser
from subprocess import check_output, STDOUT
import sys
import re
import numpy as np


def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()
   
    parser.add_argument('--case_path', type=str, default=".", 
                        help="Path to where the  cases are")
    parser.add_argument('--length',  type=float, default=0.001,
                        help="The edge length of the new mesh, this only does" + \
                             "one mesh")
    parser.add_argument('--multiple', default=[1, 1.25, 1.5, 2, 2.5, 3, 4], nargs='*',
                        help="Creates and coursen the files for the original" + \
                        " edge length and edge/l for l in the arguments")
    parser.add_argument('--boundary_layer', default=False)

    args = parser.parse_args()
    return args.case_path, args.length, args.multiple, args.boundary_layer


def success(text):
    """Make sure that there was no errors using vmtk"""
    if not "error: " in text.lower():
        return True, ""
    else:
        error_message = re.search(r'error: (.*)',
        text.lower()).groups()[0]
        return False, error_message


def get_mesh(vtu_path, mesh_path, recompute=False):
    if not path.exists(mesh_path+".gz") or recompute:
        print "Compute dolfin mesh"
        a = check_output(("vmtkmeshreader -ifile %s -ofile %s") % (vtu_path, mesh_path),
                    stderr=STDOUT, shell=True)
        
        status, msg = success(a)
        if not status:
            print "Something went wrong when converting vtu to dolfin xml:\n%s" % msg
            sys.exit(0)

    return Mesh(mesh_path+".gz")


def create_surface(vtu_path, surface_path, recompute=False):
    if not path.exists(surface_path) or recompute:
        a = check_output(("vmtkmeshtosurface -ifile %s -ofile %s") %
                        (vtu_path, surface_path), stderr=STDOUT, shell=True)
        status, msg = success(a)

        if not status:
            print "Something went wrong when converting vtu to vtp:\n%s" % msg
            sys.exit(0)


def create_mesh(ifile, ofile, length, boundary_layer, recompute=False):
    """Runs vmtk commands to create a new mesh"""
    t = ""
    if boundary_layer:
        t = " -boundarylayer 1 -sublayers 2 -boundarylayeroncaps 0" + \
            " -thicknessfactor .95 -sublayerratio .75 -tetrahedralize 1"

    # Create a new uniform courser mesh with/without boundary layers
    if not path.exists(ofile) or recompute:
        a = check_output(("vmtkmeshgenerator -ifile %s -edgelength %s %s" + \
                " -ofile %s") % (ifile, length, t, ofile), stderr=STDOUT, shell=True)

        status, msg = success(a)

        if not status:
            print "Something went wrong making the mesh:\n%s" % msg
            sys.exit(0)

    #return ReadPolyData(ofile)


def get_length(vtu_path, mesh):
    """Get the edgelength from the original mesh"""
    DG = FunctionSpace(mesh, "DG", 0)

    # Compute local edgelength
    h = project(CellSize(mesh), DG)
    
    # Compute median, more robust than mean
    mean = np.median(h.vector().array())

    return mean


def project_mesh(ofile, ifile, mesh_new_path, recompute=False):
    if not path.exists(ofile) or recompute:
        a = check_output("vmtkmeshprojection -ifile %s -rfile %s -ofile %s" % \
                        (mesh_new_path, ifile, ofile), shell=True, stderr=STDOUT)

        status, msg = success(a)
        if not status:
            print "Something went wrong projecting the mesh:"
            print a
            sys.exit(0)


def main(dirpath, factor, bl):
    # Output names
    mesh_old_path = path.join(dirpath, "mesh_old.xml")
    surface_path = path.join(dirpath, "surface.vtp")
    vtu_path = [i for i in listdir(dirpath) if "001.vtu" in i][0]
    vtu_path = path.join(dirpath, vtu_path)
    
    # Get length from old mesh
    if isinstance(factor, type([])):
        mesh_old = get_mesh(vtu_path, mesh_old_path)
        base_length = get_length(vtu_path, mesh_old)
        length = [base_length * i for i in factor]
    else:
        length = [factor]

    # Get the surface
    create_surface(vtu_path, surface_path)
    for l in length:
        # Compute the new mesh
        str_length = ("%.04f" % l).replace(".", "_")
        files_path = path.join(dirpath, "length_"+str_length+"_%s.vtu")
        mesh_new_path = path.join(dirpath, "mesh_coarsen_l_" + str_length +".vtu")
        create_mesh(surface_path, mesh_new_path, l, bl)

        # Coursen and rename
        for file in listdir(dirpath):
            if ".vtu" in file and file.startswith("u_0"):
                print file
                ofile = files_path % file[2:-4]
                ifile = path.join(dirpath, file)
                project_mesh(ofile, ifile, mesh_new_path)


if  __name__ == "__main__":
    case_path, length, multiple, bl = read_command_line()
    
    for folder in listdir(case_path):
        casedir = path.join(case_path, folder)
        if path.isdir(path.join(case_path, folder)):
            print "Working on case", folder, "..."
            if multiple != []:
                for factor in multiple:
                    main(casedir, multiple, bl)
            else:
                main(casedir, length, bl)
