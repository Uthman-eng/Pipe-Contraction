from dolfinx.io import XDMFFile
from dolfinx.mesh import meshtags_from_entities
from dolfinx.cpp.mesh import cell_entity_type
from dolfinx.io import distribute_entity_data
from dolfinx.graph import adjacencylist
from dolfinx.mesh import create_mesh
from dolfinx.cpp.mesh import to_type
from dolfinx.cpp.io import perm_gmsh
from mpi4py import MPI
from dolfinx.io.gmsh import (
    model_to_mesh,
    extract_geometry,
    extract_topology_and_markers,
    ufl_mesh,
    read_from_msh,
)

import numpy as np
import gmsh
import warnings

warnings.filterwarnings("ignore")

gmsh.initialize()
gmsh.model.add('3d_model')
L, B, H, r = 2.5, 0.41, 0.41, 0.5
channel = gmsh.model.occ.addBox(0,0,0, L, B, H)
cylinder = gmsh.model.occ.addCylinder(0.5, 0, 0.2, 0, B, 0, r)
fluid = gmsh.model.occ.cut([(3, channel)], [(3, cylinder)])

gmsh.model.occ.synchronize()
volumes = gmsh.model.getEntities(dim=3)

assert volumes == fluid[0]
fluid_marker = 11
gmsh.model.addPhysicalGroup(volumes[0][0], [volumes[0][1]], fluid_marker)
gmsh.model.setPhysicalName(volumes[0][0], fluid_marker, 'Fluid volume')

surfaces = gmsh.model.occ.getEntities(dim=2)
inlet_marker, outlet_marker, wall_marker, obstacle_marker = 1, 3, 5, 7
walls = []
obstacle = []

for surface in surfaces:
    com = gmsh.model.occ.getCenterOfMass(surface[0], surface[1])
    if np.allclose(com, [0, B / 2, H / 2]):
        gmsh.model.addPhysicalGroup(surface[0], [surface[1]], inlet_marker)
        inlet = surface[1]
        gmsh.model.setPhysicalName(surface[0], inlet_marker, "Fluid inlet")
    elif np.allclose(com, [L, B / 2, H / 2]):
        gmsh.model.addPhysicalGroup(surface[0], [surface[1]], outlet_marker)
        gmsh.model.setPhysicalName(surface[0], outlet_marker, "Fluid outlet")
    elif (
        np.isclose(com[2], 0)
        or np.isclose(com[1], B)
        or np.isclose(com[2], H)
        or np.isclose(com[1], 0)
    ):
        walls.append(surface[1])
    else:
        obstacle.append(surface[1])
gmsh.model.addPhysicalGroup(2, walls, wall_marker)
gmsh.model.setPhysicalName(2, wall_marker, 'Walls')
gmsh.model.addPhysicalGroup(2, obstacle, obstacle_marker)
gmsh.model.setPhysicalName(2, obstacle_marker, 'Obstacle')

distance = gmsh.model.mesh.field.add('Distance')
gmsh.model.mesh.field.setNumber(distance, "FacesList", obstacle)

resolution = r / 10 
threshold = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(threshold, "IField", distance)
gmsh.model.mesh.field.setNumber(threshold, "IField", distance)
gmsh.model.mesh.field.setNumber(threshold, "LcMin", resolution)
gmsh.model.mesh.field.setNumber(threshold, "LcMax", 20 * resolution)
gmsh.model.mesh.field.setNumber(threshold, "DistMin", 0.5 * r)
gmsh.model.mesh.field.setNumber(threshold, "DistMax", r)

