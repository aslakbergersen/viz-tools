from os import path, listdir
from dolfin import *
from argparse import ArgumentParser
from fenicstools import *
from subprocess import check_output, STDOUT
from vtkread import VTKToDOLFIN
import sys

def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()
   
    parser.add_argument('--case_path', type=str, default=".", 
                        help="Path to where the  cases are")
    parser.add_argument('--length',  type=float, default=0.001,
                        help="The edge length of the new mesh, this only does" + \
                             "one mesh")
    parser.add_argument('--multiple', default=[], nargs='*',
                        help="Creates and coursen the files for the original" + \
                        " edge length and edge/l for l in the arguments")

    args = parser.parse_args()
    return args.case_path, args.length, args.multiple


def success(text):
    """Make sure that there was no errors using vmtk"""
    if not "error: " in text.lower():
        return True, ""
    else:
        error_message = re.search(r'error: (.*)',
        text.lower()).groups()[0]
        return False, error_message


def make_mesh(ifile, ofile, length, boundary_layer):
    """Runs vmtk commands to create a new mesh"""
    t = ""
    if boundary_layer:
        t = " -boundarylayer 1 -sublayers 2 -boundarylayeroncaps 0" + \
            " -thicknessfactor .95 -sublayerratio .75 -tetrahedralize 1"

    # Create a new uniform courser mesh with/without boundary layers
    a = check_output(("vmtkmeshgenerator -ifile %s -edgelength %s %s" + \
            " -ofile %s") % (ifile, length, t, ofile), stderr=STDOUT, shell=True)

    status, msg = success(a)

    if not status:
        print "Something went wrong making the mesh:\n%s" % msg
        sys.exit(0)


def get_length(casedir):
    """Get the edgelength from the original mesh"""
    mesh = Mesh(path.join(casedir, "mesh.vtp")
    DG = FunctionSpace(mesh, "DG", 0)

    # Compute local edgelength
    h = CellSize(mesh)
    dl = project(h, DG) 
    
    #TODO: Check if this corresonds to vmtk's edgelength

    # Compute median, more robust than mean
    mean = np.median(dl.vector().array())

    return mean


def main(dirpath, length, bl):
    # Input names
    geometry_path = path.join(dirpath, "geo_name.vtp")
    
    # Output names
    mesh_path = path.join(dirpath, "mesh_coursen_%s.vtp" % length)
    files = path.join(dirpath, "length_"+length+"_%s")

    # Create new mesh is necessary
    if not path.exsists(mesh_path):
        make_mesh(geometry_path, mesh_path, bl)

    mesh = Mesh(mesh_path)

    for file in listdir(dirpath):
        newfile = HDF5File(mpi_comm_world(), path.join(dirpath, file), 'r')

        

if  __name__ == "__main__":
    case_path, length, multiple, bl = read_command_line()

    for folder in listdir(case_path):
        if path.isdir(path.join(case_path, folder)):
            print "Working on case", folder, "..."
            if multiple != []:
                edge_length = get_length(casedir)
                for factor in multiple:
                    main(case_path, edge_length/factor, bl)
            else:
                main(case_path, length)
