import gmsh
from mpi4py import MPI                                  
from dolfinx.io import gmsh as gmshio 

import sys
sys.path.append('/workspaces/ns-pipe-expansion/src')
from params import Parameters


def create_geometry(d1, d2, pipe_length):
    yinlet_bot = -(d1 / 2)
    youtlet_bot = -(d2 / 2) 
    inlet_rectangle = gmsh.model.occ.addRectangle(0, yinlet_bot , 0, (pipe_length / 2), d1)
    outlet_rectangle = gmsh.model.occ.addRectangle((pipe_length / 2), youtlet_bot, 0, (pipe_length / 2), d2)
    gmsh.model.occ.fragment([(2, inlet_rectangle)], [(2, outlet_rectangle)])
    gmsh.model.occ.synchronize()
    return 

def setting_groups(spacing):
    gmsh.model.mesh.setSize(gmsh.model.getEntities(dim=0), spacing)
    gmsh.model.addPhysicalGroup(1, [4], tag=1)                  # inlet (x=0)
    gmsh.model.addPhysicalGroup(1, [6], tag=2)                  # outlet (x=100)
    gmsh.model.addPhysicalGroup(1, [1, 3, 5, 7, 8, 9], tag=3)   # walls
    gmsh.model.addPhysicalGroup(2, [1, 2], tag=10)              # full domain
    return

def create_mesh(params: Parameters):

    gmsh.initialize()
    gmsh.model.add('Pipe Mesh')

    # using rectangles and constructing joint shape 
    create_geometry(
        d1=params.d1,
        d2=params.d2,
        pipe_length=params.pipe_length)
    setting_groups(spacing=params.spacing)
    gmsh.model.mesh.generate(2)
    mesh_data = gmshio.model_to_mesh(gmsh.model, MPI.COMM_WORLD, rank=0, gdim=2)
    gmsh.finalize()
    
    mesh = mesh_data.mesh
    cell_tags = mesh_data.cell_tags
    facet_tags = mesh_data.facet_tags

    mesh = mesh_data.mesh
    assert mesh_data.facet_tags is not None
    
    return  mesh, facet_tags, cell_tags

def main():
    params = Parameters()
    mesh, facet_tags, cell_tags = create_mesh(params)

if __name__ == "__main__":
    main()
