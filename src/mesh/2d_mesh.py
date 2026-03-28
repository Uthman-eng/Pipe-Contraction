import gmsh
import matplotlib.pyplot as plt

gmsh.initialize()
gmsh.model.add('2d_Pipe_Contraction')

d1 = 7.2 # mm 
d2 = 17.2 # mm
pipe_length = 100

spacing = 5

# inlet points 7.2mm 
yinlet_bot = -(d1 / 2)
# outlet points 17.2mm
youtlet_bot = -(d2 / 2)
# using rectangles and constructing joint shape 
inlet_rectangle = gmsh.model.occ.addRectangle(0, yinlet_bot , 0, (pipe_length / 2), d1)
outlet_rectangle = gmsh.model.occ.addRectangle((pipe_length / 2), youtlet_bot, 0, (pipe_length / 2), d2)
gmsh.model.occ.fragment([(2, inlet_rectangle)], [(2, outlet_rectangle)])    # correct

gmsh.model.occ.synchronize()

"""# to see the values of the dim 1 and dim 2 objects 
for dim, tag in gmsh.model.getEntities(dim=1):
    bb = gmsh.model.getBoundingBox(dim, tag)
    print(f"tag: {tag},  x={bb[0]:.2f} to {bb[3]:.2f}, y={bb[1]:.2f} to {bb[4]:.2f}")

for dim, tag in gmsh.model.getEntities(dim=2):
    print(f"surface {tag}")"""

gmsh.model.mesh.setSize(gmsh.model.getEntities(dim=0), spacing)


gmsh.model.addPhysicalGroup(1, [4], tag=1)              # inlet (x=0)
gmsh.model.addPhysicalGroup(1, [6], tag=2)              # outlet (x=100)
gmsh.model.addPhysicalGroup(1, [1, 3, 5, 7, 8, 9], tag=3)  # walls
gmsh.model.addPhysicalGroup(2, [1, 2], tag=10)          # full domain


gmsh.model.mesh.generate(2)


# extracting points / nodes
nodes_tags, coords, _ = gmsh.model.mesh.getNodes()
coords = coords.reshape(-1, 3)
x , y = coords[:,0], coords[:, 1]


# extracting lines 
_, _, elem_node_tags = gmsh.model.mesh.getElements(dim=1)
lines = elem_node_tags[0].reshape(-1, 2) - 1

gmsh.finalize()


# visualisation 
fig, ax = plt.subplots(figsize=(10, 8))
for line in lines:
    ax.plot(x[line], y[line], color="steelblue", linewidth=0.5)
ax.scatter(x, y, s=10, color="tomato", zorder=3)
ax.set_aspect("equal")
ax.set_xlabel('pipe length (mm)')
ax.set_ylabel('diameter')
plt.savefig("plotting/2d_mesh")