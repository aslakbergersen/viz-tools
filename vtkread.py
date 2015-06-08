from xml.etree.ElementTree import ElementTree
import numpy as np
import os
import vtk
import array_handler

from dolfin import Mesh, MeshEditor, VectorFunctionSpace, FunctionSpace, Function

def vtk_ug_to_dolfin_mesh(ug):
    """
    Create a DOLFIN Mesh from a vtkUnstructuredGrid object
    """
    if not isinstance(ug, vtk.vtkUnstructuredGrid):
        raise TypeError("Expected a 'vtkUnstructuredGrid'")
    
    # Get mesh data
    num_cells = int(ug.GetNumberOfCells())
    num_vertices = int(ug.GetNumberOfPoints())
    
    # Get topological and geometrical dimensions
    cell = ug.GetCell(0)
    gdim = int(cell.GetCellDimension())
    cell_type = cell.GetCellType()                                                                                                                                          
    if cell_type not in [vtk.VTK_TETRA, vtk.VTK_TRIANGLE]:                                                                                                                  
        raise TypeError("DOLFIN only support meshes of triangles " + \
                        "and tetrahedrons.")
    
    tdim = 3 if cell_type == vtk.VTK_TETRA else 2
    
    # Create empty DOLFIN Mesh
    mesh = Mesh()
    editor = MeshEditor()
    editor.open(mesh, tdim, gdim)
    editor.init_cells(num_cells)
    editor.init_vertices(num_vertices)
    editor.close()
    
    # Assign the cell and vertex informations directly from the vtk data
    cells_array = array_handler.vtk2array(ug.GetCells().GetData())
    
    # Get the assumed fixed size of indices and create an index array
    cell_size = cell.GetPointIds().GetNumberOfIds()
    cellinds = np.arange(len(cells_array))
    
    # Each cell_ids_size:th data point need to be deleted from the
    # index array
    ind_delete = slice(0, len(cells_array), cell_size+1)
    
    # Check that the removed value all have the same value which should
    # be the size of the data
    if not np.all(cells_array[ind_delete]==cell_size):
        raise ValueError("Expected all cells to be of the same size")
    
    cellinds = np.delete(cellinds, ind_delete)
    
    # Get cell data from mesh and make it writeable (cell data is non
    # writeable by default) and update the values
    mesh_cells = mesh.cells()
    mesh_cells.flags.writeable = True
    mesh_cells[:] = np.reshape(cells_array[cellinds], \
                              (num_cells , cell_size))
    
    # Set coordinates from vtk data
    vertex_array = array_handler.vtk2array(ug.GetPoints().GetData())
    if vertex_array.shape[1] != gdim:
        vertex_array = vertex_array[:,:gdim]
    mesh.coordinates()[:] = vertex_array
    return mesh
    


class VTKToDOLFIN(object):
    """
    A wrapper around vtk to simplify handling of VTK files
    generated from DOLFIN.

    The class handles reading of data into DOLFIN objects for further processing
    
    """
    def __init__(self, filename, mesh=None, deepcopy=False):
        """
        Initialize a the reader with a pvd or a vtu filename
        """
        if not os.path.isfile(filename):
            raise IOError("File '%s' does not excist"%filename)
        filetype = filename.split(".")[-1]
        self._name = ".".join(filename.split(".")[0:-1])
        if filetype not in ["pvd", "vtu"]:
            raise TypeError("Expected a 'pvd' or a 'vtu' file")

        # Get dirname
        dirname = os.path.dirname(filename)
        
        # Check mesh argument
        if mesh is not None and not isinstance(mesh, Mesh):
            raise TypeError, "Expected a 'Mesh' for the mesh arguments"

        # Store deepcopy argument
        self._deepcopy = deepcopy
        
        # Store mesh
        self._mesh = mesh
        
        # Initialize the filename cache
        self._filenames = []
        if filetype == "vtu":
            self._filenames.append(filename)
            self._times = np.array([])
        else:
            # Parse pvd file
            tree = ElementTree(file=filename)
            times = []
            for item in tree.iter():
                if item.tag == "DataSet":
                    self._filenames.append(os.path.join(\
                        dirname,item.attrib["file"]))
                    times.append(float(item.attrib["timestep"]))
            
            times = np.array(times, dtype='d')

            # If there are no time data stored in the file use an empty array
            if np.all(np.diff(times)==1):
                times = np.array([], "")

            # Store time data
            self._times = times

        # Construct file reader
        self.reader = vtk.vtkXMLUnstructuredGridReader()
        
        # Read in data from file
        self._update_vtk_data()

        # Init dolfin structures (Function, FunctionSpace)
        self._init_dolfin_data()

    def _update_vtk_data(self, index=0):
        "Set a new data file"

        # Update file name
        print "Reading '%s'"%self._filenames[index]
        self.reader.SetFileName(self._filenames[index])
        
        # Read data
        self.reader.Update()
        
        # Set data type (scalar or vector)
        # FIXME: Include Tensors when that is supported by DOLFIN
        self.scalar = self.reader.GetOutput().GetPointData().GetScalars() is not None

        print "Scalar data set" if self.scalar else "Vector data set"
        
    def _init_dolfin_data(self):
        "Update DOLFIN function from vtk data"
        
        if self.reader.GetNumberOfPointArrays() != 1:
            raise ValueError("Expected the vtk file to include one data "\
                             "set per vertex.")

        # Initilize FunctionSpace and Function if not initialized
        if self.scalar:
            self._V = FunctionSpace(self.mesh(), "CG", 1)
        else:
            self._V = VectorFunctionSpace(self.mesh(), "CG", 1)
        
        self._u = Function(self._V)
        
    def _update_dolfin_data(self):
        "Update dolfin data from present VTK file"
        
        # Get VTK point data
        point_data = self.reader.GetOutput().GetPointData()

        # Get data and update Function
        if self.scalar:
            self._u.vector()[:] = array_handler.vtk2array(point_data.GetScalars())
        else:
            values = array_handler.vtk2array(point_data.GetVectors()).transpose()
            self._u.vector()[:] = np.reshape(values, (np.prod(values.shape),))
    
    def functions_space(self):
        "Return the FunctionSpace"
        return self._V
    
    def mesh(self):
        "Return the dolfin mesh"

        # If no mesh is stored read in from UnstructuredGridData
        if self._mesh is None:
            self._mesh = vtk_ug_to_dolfin_mesh(self.reader.GetOutput())

        # Small sanity check
        assert(self._mesh.num_vertices() == \
               self.reader.GetOutput().GetNumberOfPoints() and \
               self._mesh.num_cells() == \
               self.reader.GetOutput().GetNumberOfCells())
        
        return self._mesh
    
    def name(self):
        "Return the name"
        return self._name

    def __getitem__(self, index):
        "x.__getitem__(y) <==> x[y]"
        # Update data structures to next index if not out of files
        if not isinstance(index, int):
            raise TypeError("Expected an int for the index argument")
        if index < 0 or index >= len(self):
            raise IndexError("index need to be smaller than ")

        # Update the stored data
        self._update_vtk_data(index)
        self._update_dolfin_data()

        # Should we return a copy of the stored data?
        u = self._u.copy() if self._deepcopy else self._u

        # If time is registered return with this information
        if len(self._times):
            return self._times[index], u
        
        return u

    def __len__(self):
        "x.__len__() <==> len(x)"
        return len(self._filenames)

    def __iter__(self):
        "x.__iter__() <==> iter(x)"
        for i in xrange(len(self)):
            yield self[i]
