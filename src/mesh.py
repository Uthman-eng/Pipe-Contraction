import gmsh
from mpi4py import MPI
import gmsh                                   
from dolfinx.io import VTXWriter, gmsh as gmshio   


def create_geometry(d1, d2, yinlet_bot, youtlet_bot, pipe_length):
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


def create_mesh(d1, d2, pipe_length, spacing):

    gmsh.initialize()
    gmsh.model.add('Pipe Mesh')

    yinlet_bot = -(d1 / 2)
    youtlet_bot = -(d2 / 2)

    # using rectangles and constructing joint shape 
    create_geometry(
        d1=d1,
        d2=d2,
        yinlet_bot=yinlet_bot,
        youtlet_bot=youtlet_bot,
        pipe_length=pipe_length)
    setting_groups(spacing=spacing)
    gmsh.model.mesh.generate(2)
    mesh_data = gmshio.model_to_mesh(gmsh.model, MPI.COMM_WORLD, rank=0, gdim=2)
    gmsh.finalize()
    
    mesh = mesh_data.mesh
    cell_tags = mesh_data.cell_tags
    facet_tags = mesh_data.facet_tags

    mesh = mesh_data.mesh
    assert mesh_data.facet_tags is not None
    facet_tags.name = "Facet markers"
    
    return  mesh, facet_tags, cell_tags,  facet_tags.name


if __name__ == "__main__":
    d1 = (7.2 / 1000)
    d2 = (17.2 / 1000)
    pipe_length = (100 / 1000)
    spacing = (1 / 1000)

    mesh, ft, ct , _ = create_mesh(
        d1=d1,
        d2=d2, 
        pipe_length=pipe_length, 
        spacing=spacing)
