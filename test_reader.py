from dolfin import Mesh

mesh = Mesh("../Dropbox/course_test/Pair8a_cl_ex_cl_dist.xml.gz")
reader = VTKToDOLFIN("../Dropbox/course_test/u_000000.vtu", mesh)
