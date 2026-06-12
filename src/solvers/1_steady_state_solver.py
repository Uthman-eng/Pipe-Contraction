import numpy as np
import tqdm.autonotebook
import pyvista
from petsc4py import PETSc
import gmsh

from basix.ufl import element

from dolfinx.fem import (
    Constant,
    Function,
    functionspace,
    dirichletbc,
    extract_function_spaces,
    form,
    locate_dofs_topological,
    set_bc,
)
from dolfinx.plot import vtk_mesh
from dolfinx.fem.petsc import (
    apply_lifting,
    assemble_matrix,
    assemble_vector,
    create_vector,
    create_matrix,
)

from ufl import (
    TestFunction,
    TrialFunction,
    div,
    dot,
    dx,
    inner,
    lhs,
    grad,
    nabla_grad,
    rhs,
)

from functools import partial
from src.mesh.mesh import create_mesh

import sys
sys.path.append('/workspaces/ns-pipe-expansion/src')
from params import Parameters


def inlet_velocity(x, U_max, d1):
    """
    Parabolic (Poiseuille) inlet velocity profile, centred on y=0.

    Parameters:
    x : np.ndarray, shape (3, N)
        Coordinates of evaluation points, as supplied by `interpolate`.
    U_max : float
        Maximum (centreline) velocity of the parabolic profile.
    d1 : float
        Diameter/height of the inlet channel.

    Returns:
    np.ndarray, shape (2, N)
        Velocity vector (vx, vy) at each point in x.
    """
    yinlet_bot = -(d1 / 2)
    values = np.zeros((2, x.shape[1]))
    values[0] = (4 * U_max * (x[1] - yinlet_bot) * ((yinlet_bot + d1) - x[1])) / d1**2
    return values

def setting_functionspace(mesh):
    """
    Create Taylor-Hood-style function spaces for IPCS (P2 velocity, P2 pressure).

    Parameters:
    mesh : dolfinx.mesh.Mesh

    Returns:
    V : dolfinx.fem.FunctionSpace
        Vector P2 space for velocity.
    Q : dolfinx.fem.FunctionSpace
        Scalar P2 space for pressure.
    fdim : int
        Facet dimension (mesh.topology.dim - 1), for use in locate_dofs_topological.
    """
    v_cg2 = element("Lagrange", mesh.basix_cell(), 2, shape=(mesh.geometry.dim,))
    s_cg2 = element("Lagrange", mesh.basix_cell(), 2)
    V = functionspace(mesh, v_cg2)
    Q = functionspace(mesh, s_cg2)

    fdim = mesh.topology.dim - 1
    return V, Q, fdim

def set_dirichletbc(mesh, V, Q, fdim, facet_tags, d1, U_max):
    """
    Build Dirichlet boundary conditions for pipe flow: no-slip walls,
    parabolic inlet velocity, and zero outlet pressure.

    Parameters
    ----------
    mesh : dolfinx.mesh.Mesh
    V, Q : dolfinx.fem.FunctionSpace
        Velocity and pressure function spaces.
    fdim : int
        Facet dimension.
    facet_tags : dolfinx.mesh.MeshTags
        Tags identifying boundary facets (1=inlet, 2=outlet, 3=walls).
    d1 : float
        Inlet diameter, used in the parabolic velocity profile.
    U_max : float
        Centreline velocity for the inlet profile.

    Returns
    -------
    bcu : list of dolfinx.fem.DirichletBC
        Velocity BCs (inlet + walls).
    bcp : list of dolfinx.fem.DirichletBC
        Pressure BCs (outlet).
    """
    # walls
    u_zero = np.zeros(mesh.geometry.dim, dtype=PETSc.ScalarType)
    bcu_walls = dirichletbc(u_zero, locate_dofs_topological(V, fdim, facet_tags.find(3)), V)

    # inlet
    u_inlet = Function(V)
    u_inlet.interpolate(partial(inlet_velocity, U_max=U_max, d1=d1))
    bcu_inlet = dirichletbc(u_inlet, locate_dofs_topological(V, fdim, facet_tags.find(1)))

    bcp_outlet = dirichletbc(PETSc.ScalarType(0), locate_dofs_topological(Q, fdim, facet_tags.find(2)), Q)

    bcu = [bcu_inlet, bcu_walls]
    bcp = [bcp_outlet]
    return bcu, bcp

def newton_solver(V,Q, mesh, rho, k, mu, T, dt, bcu, bcp):
    """
    NOT YET DONE BTW.

    """

def plotting_figure(mesh, u_n, V, v_avg, Re):
    topology, cell_types, geometry = vtk_mesh(V)
    values = np.zeros((geometry.shape[0], 3), dtype=np.float64)
    values[:, : len(u_n)] = u_n.x.array.real.reshape((geometry.shape[0], len(u_n)))

    function_grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
    function_grid["u"] = values
    glyphs = function_grid.glyph(orient="u", factor=0.001)

    tdim = mesh.topology.dim
    mesh.topology.create_connectivity(tdim, tdim)
    grid = pyvista.UnstructuredGrid(*vtk_mesh(mesh, tdim))

    pyvista.start_xvfb()
    plotter = pyvista.Plotter(off_screen=True)
    plotter.add_text(f"Fluid flow through a Pipe\n\nInlet Average velocity: {v_avg:.4f} (m/s) | Re: {Re:.4f}",
                    position="upper_edge", font_size=12)
    plotter.add_mesh(grid, style="wireframe", color="k", line_width=0.5)
    plotter.add_mesh(glyphs, scalars="GlyphScale", cmap="viridis")
    plotter.view_xy()
    plotter.screenshot("output/time_dependant_solver.png")
    plotter.close()

def main():
    params = Parameters()
    d1 = params.d1
    d2 = params.d2
    pipe_length = params.pipe_length
    spacing = params.pipe_length
    T = params.T
    dt = params.dt

    v_avg = params.v_avg
    U_max = 3 / (2 * v_avg)

    mesh, facet_tags, cell_tags = create_mesh(params)

    k = Constant(mesh, PETSc.ScalarType(dt))
    mu = Constant(mesh, PETSc.ScalarType(0.001))  # Dynamic viscosity
    rho = Constant(mesh, PETSc.ScalarType(1000))
    Re = (rho.value * v_avg * d1) / mu.value  # use d1 as characteristic length

    V, Q, fdim = setting_functionspace(mesh)
    bcu, bcp = set_dirichletbc(mesh, V, Q, fdim, facet_tags, d1, U_max)
    u_n, p_ = newton_solver(V,Q, mesh, rho, k, mu, T, dt, bcu, bcp)

    Re = (rho.value * v_avg * d1) / mu.value  

    plotting_figure(mesh, u_n, V, v_avg, Re)


if __name__ == "__main__":
    main()
