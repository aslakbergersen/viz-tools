from dolfin import *
from os import listdir, path
import re
from argparse import ArgumentParser

# General improvements suggestions:
# TODO: Add support for multiple cases
# TODO: Add possibilty to project Q or u down in to a courser mesh

def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()
    parser.add_argument('--d', '--dir_path', type=str, default=".",
                        help="Relative path to the case you want to convert")
    args = parser.parse_args()

    return args.d


def main(dirpath):
    # Get files to convert
    files = listdir(dirpath)
    files = [f for f in files if f.startswith("pipe_ipcs")]
    
    # Read mesh
    case = re.findall("_case(.*?)_", files[0])[0]

    # Read mesh from xdmf if exists, else convert the XML file
    if not path.exists(path.join(dirpath, "case%s_mesh.xdmf") % case):
        mesh = Mesh(path.join(dirpath, "challenge_case%s_zero_pressure.xml.gz" % case))
        file_mesh = XDMFFile(mpi_comm_world(), path.join(dirpath, "case%s_mesh.xdmf") % case)
        file_mesh << mesh
    else:
        mesh = Mesh(path.join(dirpath, "case%s_mesh.xdmf") % case)
    
    # Find order and create mesh
    order = int(re.findall(r"uOrder(\d)", files[0])[0])
    V_seg = FunctionSpace(mesh, "CG", order)
    V_vector = VectorFunctionSpace(mesh, "CG", order)
    u_vector = Function(V_vector)
    Q = Function(V_seg)
    u_vector.rename("u", "velocity")
    Q.rename("Q", "Q-criterion")

    # Get time
    timesteps = []
    times = []
    for text in files:
        times.append(float(re.findall(r"_t_(.*?)_u", text)[0]))
        timesteps.append(int(re.findall(r"ts=(\d*)", text)[0]))
    timesteps.sort()
    times.sort()
    timesteps = timesteps[::3]
    times = times[::3]

    # Get dimention
    dim = mesh.geometry().dim()

    # XDMFFile for storage
    file_u = XDMFFile(mpi_comm_world(), path.join(dirpath, "case%s_u.xdmf") % case)
    file_Q = XDMFFile(mpi_comm_world(), path.join(dirpath, "case%s_Q.xdmf") % case)
    file_u.parameters["rewrite_function_mesh"] = False
    file_u.parameters["flush_output"] = True
    file_Q.parameters["rewrite_function_mesh"] = False
    file_Q.parameters["flush_output"] = True

    # Naming convention
    filename = "pipe_ipcs_ab_cn_challenge_case1_zero_pressure_constant" + \
               "_ts10000_cycles4_uOrder1_curcyc_3_t_%s_u%s_ts=%s.xml.gz"
    filename = path.join(dirpath, filename)

    for j in range(len(times)):
        # Current time
        t = times[j]
        step = timesteps[j]
        
        # Insert segregated velocity fields into the vector space
        u_seg = [Function(V_seg, filename % (t, i, step)) for i in range(dim)]
        [assign(u_vector.sub(i), u_seg[i]) for i in range(dim)]

        # Compute Q criteria
        S = (grad(u_vector) + grad(u_vector).T)/2
        Omega = (grad(u_vector) - grad(u_vector).T)/2
        expr = 0.5*Omega**2 - S**2
        Q.assign(project(expr, V_seg))
        t1 = time()

        # Store in xdmf format
        file_u << (u_vector, float(t))
        file_Q << (Q, float(t))


if __name__ == "__main__":
    dirpath = read_command_line()
    main(dirpath)
    t1 = time()
