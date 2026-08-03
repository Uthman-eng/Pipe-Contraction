import numpy as np
import tqdm.autonotebook
import pyvista
from petsc4py import PETSc
import gmsh

from basix.ufl import element

from dolfinx.fem import (Constant, Function, functionspace, dirichletbc, 
    extract_function_spaces, form, locate_dofs_topological, set_bc,
)

from dolfinx.plot import vtk_mesh
from dolfinx.fem.petsc import (apply_lifting, assemble_matrix, assemble_vector,
                             create_vector, create_matrix,
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

    Parameters
    ----------
    x : np.ndarray, shape (3, N)
        Coordinates of evaluation points, as supplied by `interpolate`.
    U_max : float
        Maximum (centreline) velocity of the parabolic profile.
    d1 : float
        Diameter/height of the inlet channel.

    Returns
    -------
    np.ndarray, shape (2, N)
        Velocity vector (vx, vy) at each point in x.
    """
    yinlet_bot = -(d1 / 2)
    values = np.zeros((2, x.shape[1]))
    values[0] = (4 * U_max * (x[1] - yinlet_bot) * ((yinlet_bot + d1) - x[1])) / d1**2
    return values

def setting_functionspace(mesh):
    """
    Create Taylor-Hood-style function spaces for IPCS (P2 velocity, P1 pressure).

    Parameters
    ----------
    mesh : dolfinx.mesh.Mesh

    Returns
    -------
    V : dolfinx.fem.FunctionSpace
        Vector P2 space for velocity.
    Q : dolfinx.fem.FunctionSpace
        Scalar P1 space for pressure.
    fdim : int
        Facet dimension (mesh.topology.dim - 1), for use in locate_dofs_topological.
    """
    v_cg2 = element("Lagrange", mesh.basix_cell(), 2, shape=(mesh.geometry.dim,))
    s_cg1 = element("Lagrange", mesh.basix_cell(), 1)
    V = functionspace(mesh, v_cg2)
    Q = functionspace(mesh, s_cg1)

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

def ipcs_solver(V,Q, mesh, rho, k, mu, T, dt, bcu, bcp):
    """
    Solve transient incompressible Navier-Stokes via the IPCS
    (Incremental Pressure Correction Scheme) splitting method.

    Parameters
    ----------
    V, Q : dolfinx.fem.FunctionSpace
        Velocity and pressure function spaces.
    mesh : dolfinx.mesh.Mesh
    rho, mu : dolfinx.fem.Constant
        Density and dynamic viscosity.
    k : dolfinx.fem.Constant
        Timestep size (dt), wrapped as a Constant for use in UFL forms.
    num_steps : int
        Number of timesteps to run.
    bcu : list of dolfinx.fem.DirichletBC
        Velocity boundary conditions.
    bcp : list of dolfinx.fem.DirichletBC
        Pressure boundary conditions.

    Returns
    -------
    u_n : dolfinx.fem.Function
        Velocity field at the final timestep.
    p_ : dolfinx.fem.Function
        Pressure field at the final timestep.
    """
    num_steps = int(T / dt)
    
     
    u = TrialFunction(V)
    v = TestFunction(V)
    u_ = Function(V, name="u")
    u_s = Function(V, name="u_tentative")
    u_n = Function(V)
    u_n1 = Function(V)
    p = TrialFunction(Q)
    q = TestFunction(Q)
    p_ = Function(Q, name="p")
    phi = Function(Q, name="phi")

    f = Constant(mesh, PETSc.ScalarType((0, 0)))
    F1 = rho / k * dot(u - u_n, v) * dx
    F1 += inner(dot(1.5 * u_n - 0.5 * u_n1, 0.5 * nabla_grad(u + u_n)), v) * dx
    F1 += 0.5 * mu * inner(grad(u + u_n), grad(v)) * dx - dot(p_, div(v)) * dx
    F1 += dot(f, v) * dx
    a1 = form(lhs(F1))
    L1 = form(rhs(F1))
    A1 = create_matrix(a1)
    b1 = create_vector(extract_function_spaces(L1))

    # Step 2
    a2 = form(dot(grad(p), grad(q)) * dx)
    L2 = form(-rho / k * dot(div(u_s), q) * dx)
    A2 = assemble_matrix(a2, bcs=bcp)
    A2.assemble()
    b2 = create_vector(extract_function_spaces(L2))

    a3 = form(rho * dot(u, v) * dx)
    L3 = form(rho * dot(u_s, v) * dx - k * dot(nabla_grad(phi), v) * dx)
    A3 = assemble_matrix(a3)
    A3.assemble()
    b3 = create_vector(extract_function_spaces(L3))

    # backend PETSC
    # Solver for step 1
    solver1 = PETSc.KSP().create(mesh.comm)
    solver1.setOperators(A1)
    solver1.setType(PETSc.KSP.Type.BCGS)
    pc1 = solver1.getPC()
    pc1.setType(PETSc.PC.Type.JACOBI)

    # Solver for step 2
    solver2 = PETSc.KSP().create(mesh.comm)
    solver2.setOperators(A2)
    solver2.setType(PETSc.KSP.Type.MINRES)
    pc2 = solver2.getPC()
    pc2.setType(PETSc.PC.Type.HYPRE)
    pc2.setHYPREType("boomeramg")

    # Solver for step 3
    solver3 = PETSc.KSP().create(mesh.comm)
    solver3.setOperators(A3)
    solver3.setType(PETSc.KSP.Type.CG)
    pc3 = solver3.getPC()
    pc3.setType(PETSc.PC.Type.SOR)

    progress = tqdm.autonotebook.tqdm(desc="Solving PDE", total=num_steps)

    for i in range(num_steps):
        progress.update(1)
        # Update current time step


        # Step 1: Tentative velocity step
        A1.zeroEntries()
        assemble_matrix(A1, a1, bcs=bcu)
        A1.assemble()
        with b1.localForm() as loc:
            loc.set(0)
        assemble_vector(b1, L1)
        apply_lifting(b1, [a1], [bcu])
        b1.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        set_bc(b1, bcu)
        solver1.solve(b1, u_s.x.petsc_vec)
        u_s.x.scatter_forward()

        # Step 2: Pressure corrrection step
        with b2.localForm() as loc:
            loc.set(0)
        assemble_vector(b2, L2)
        apply_lifting(b2, [a2], [bcp])
        b2.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        set_bc(b2, bcp)
        solver2.solve(b2, phi.x.petsc_vec)
        phi.x.scatter_forward()

        p_.x.petsc_vec.axpy(1, phi.x.petsc_vec)
        p_.x.scatter_forward()

        # Step 3: Velocity correction step
        with b3.localForm() as loc:
            loc.set(0)
        assemble_vector(b3, L3)
        b3.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        solver3.solve(b3, u_.x.petsc_vec)
        u_.x.scatter_forward()

        # Update variable with solution form this time step
        with (
            u_.x.petsc_vec.localForm() as loc_,
            u_n.x.petsc_vec.localForm() as loc_n,
            u_n1.x.petsc_vec.localForm() as loc_n1,
        ):
            loc_n.copy(loc_n1)
            loc_.copy(loc_n)
    progress.close()


    A1.destroy()
    A2.destroy()
    A3.destroy()
    b1.destroy()
    b2.destroy()
    b3.destroy()
    solver1.destroy()
    solver2.destroy()
    solver3.destroy()
    return u_n, p_  # last timestep velocity (u_n) and pressure (p_)

def plotting_figure(mesh, u_n, V, v_avg, Re):
    topology, cell_types, geometry = vtk_mesh(V)
    values = np.zeros((geometry.shape[0], 3), dtype=np.float64)
    values[:, : len(u_n)] = u_n.x.array.real.reshape((geometry.shape[0], len(u_n)))

    function_grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
    function_grid["u"] = values
    # tolerance thins out glyphs so individual arrows stay legible on the refined mesh
    glyphs = function_grid.glyph(orient="u", factor=0.02, tolerance=0.015)

    tdim = mesh.topology.dim
    mesh.topology.create_connectivity(tdim, tdim)
    grid = pyvista.UnstructuredGrid(*vtk_mesh(mesh, tdim))

    pyvista.start_xvfb()
    plotter = pyvista.Plotter(off_screen=True)
    plotter.add_text(f"Fluid flow through a Pipe\n\nInlet Average velocity: {v_avg:.4f} (m/s) | Re: {Re:.4f}",
                    position="upper_edge", font_size=12)
    plotter.add_mesh(grid, style="wireframe", color="lightgray", line_width=0.4, opacity=0.6)
    plotter.add_mesh(glyphs, scalars="GlyphScale", cmap="viridis")
    plotter.view_xy()
    plotter.enable_anti_aliasing("ssaa")
    plotter.screenshot("output/time_dependent_solver.png", scale=2)
    plotter.close()

def main():
    params = Parameters()
    d1 = params.d1
    T = params.T
    dt = params.dt

    v_avg = params.v_avg
    U_max = 1.5 * v_avg  # parabolic profile: v_avg = (2/3) * U_max

    mesh, facet_tags, cell_tags = create_mesh(params)

    k = Constant(mesh, PETSc.ScalarType(dt))
    mu = Constant(mesh, PETSc.ScalarType(0.001))  # Dynamic viscosity
    rho = Constant(mesh, PETSc.ScalarType(1000))
    Re = (rho.value * v_avg * d1) / mu.value  # use d1 as characteristic length

    V, Q, fdim = setting_functionspace(mesh)
    bcu, bcp = set_dirichletbc(mesh, V, Q, fdim, facet_tags, d1, U_max)
    u_n, p_ = ipcs_solver(V,Q, mesh, rho, k, mu, T, dt, bcu, bcp)

    plotting_figure(mesh, u_n, V, v_avg, Re)


if __name__ == "__main__":
    main()
