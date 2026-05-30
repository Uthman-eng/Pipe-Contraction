# Navier-Stokes Modelling - Technical Log

## 1. Problem 

Solving incompressible Navier-Stokes through an expansion pipe from $d_1 = 7.2\,\text{mm}$ to $d_2 = 17.2\,\text{mm}$. Goal: back-calculate the loss coefficient $K$ and compare to experiment data.


![Fluid through pipe](/assets/pipe_flow.png)


## 2. Governing Equations

**1. Navier-Stokes (Momentum Equation):**

$$\rho\left(\frac{\partial u}{\partial t} + u \cdot \nabla u\right) = \nabla \cdot \sigma(u,\, p) + f \tag{1}$$

where $u$ = velocity, $p$ = pressure, $f$ = force per unit volume.

**1.1 Stress Tensor**

$$\sigma(u,\, p) = 2\mu\,\varepsilon(u) - pI$$

where $\mu$ is the dynamic viscosity and $\varepsilon(u) = \tfrac{1}{2}\!\left(\nabla u + (\nabla u)^T\right)$ is the strain rate tensor.


**2. Incompressibility Condition and Continuity Equations:**


$$\frac{\partial p}{\partial t} + \nabla \cdot (ρu) = 0  $$

Since $ρ$ is kept constant the terms vanishes and becomes:

$$\nabla \cdot u = 0\tag2$$



## 3. Weak Formulation

Inner product notation $\text{(Lebesgue square-integrable)}$ over the domain $\Omega$ and its boundary $\partial\Omega$:

$$\langle v,\, w \rangle = \int_\Omega v \cdot w\ \mathrm{d}x,
\qquad 
\langle v,\, w \rangle_{\partial\Omega} = \int_{\partial\Omega} v \cdot w\, \mathrm{d}s$$



Integrating the stress divergence term by parts gives the boundary term:

$$\langle -\nabla \cdot \sigma,\, v \rangle = \langle \sigma,\, \varepsilon(v) \rangle - \langle T,\, v \rangle_{\partial\Omega}$$

The boundary term $\langle T,\, v \rangle_{\partial\Omega}$ vanishes at inlet and walls; at the outlet it gives the do-nothing (natural outflow) condition.

---

## 4. IPCS Scheme

Monolithic is expensive — splitting reduces to three cheaper sub-problems. Deeper understanding of IPCS schemes vs. monolithic isn't too well-known right now.

**Step 1 — Tentative velocity** (weak momentum, stress term integrated by parts):

$$\langle -\nabla \cdot \sigma,\, v \rangle = \langle \sigma,\, \varepsilon(v) \rangle - \langle T,\, v \rangle_{\partial\Omega}$$

**Step 2 — Pressure correction** (Poisson equation):

$$-\frac{\rho\,\nabla \cdot u^*}{\Delta t} + \nabla^2 p^{n+1} - \nabla^2 p^{n} = 0$$

**Step 3 — Velocity correction**:

$$\rho\,\langle u^{n+1} - u^*,\, v \rangle = -\Delta t\,\langle \nabla(p^{n+1} - p^{n}),\, v \rangle$$

**Crank-Nicolson / Adams-Bashforth time discretisation**

$$u^{n+\frac{1}{2}} \approx \frac{u^{n+1} + u^n}{2}$$

Crank-Nicolson for the time-dependent term — this converges better than implicit Euler for the given mesh and is less stable, so we know where an error may actually come from. Adams-Bashforth approximation for the nonlinear term.

**Taylor-Hood elements P2/P1**

Zero idea about Taylor-Hood elements right now, reading through this.

---

## 5. Implementation

Defined the mesh (2D) in Gmsh following through examples and tagged wall, inlet, outlet physical groups. Kept everything else the same and followed on exactly from FEniCSx tutorials [1, 2].

1. Set up Gmsh and correctly set boundary tags — walls, inlet, outlet. The mesh was made to take input and output parameters so I can easily change the diameters, spacing, and length of the pipe.
2. Set up parameters and function spaces and input variational form.
3. Set up variational UFL forms and PETSc solvers for each of the three IPCS steps — mostly copied from examples.
4. Time loop: stepped through time to get pressure and velocity fields.
5. PyVista plotting of mesh, velocity arrows, and magnitude.

Inlet uses a parabolic (Poiseuille) profile matching zero at the walls and maximum at the centre:

$$u(y) = \frac{4\,U_{\max}(y - y_{\mathrm{bot}})(y_{\mathrm{top}} - y)}{d_1^2}$$

```python
def inlet_velocity(x):
    values = np.zeros((2, x.shape[1]))
    values[0] = (4 * U_max * (x[1] - yinlet_bot) * ((yinlet_bot + d2) - x[1])) / d1**2
    return values
```

Physical parameters set to water; Reynolds number kept low to prevent turbulence — otherwise the model just becomes wrong. I haven't done this yet but I need to — maybe add onto Streamlit, or at the very least make a dashboard that shows me max $u$ and average inlet $\bar{U}$, etc. My experiment data uses average $u$, so I can then modify my $u$ to get the average $u$ matching that of the experiment, and then the back-calculation will be easier. At minimum the dashboard should give me average velocity at inlet and outlet, and the Reynolds number.

---

## 6. Results

*[Placeholder: results once runs are complete]*

---

## 7. What I Still Need to Understand

- Inf-sup condition properly
- Monolithic approach vs. splitting
- Long Chen notes (see references [5, 6])

Still need to read up on inf-sup conditions and play around with methods, and take time to understand.

---

## 8. References

1. Dokken, J. S. — *FEniCSx Tutorial: Navier-Stokes equations.* <https://jsdokken.com/dolfinx-tutorial/chapter2/navierstokes.html>
2. Dokken, J. S. — *Test problem 1: Channel flow.* <https://jsdokken.com/dolfinx-tutorial/chapter2/ns_code1.html>
3. Dokken, J. S. — *Test problem 2: Flow past a cylinder.* <https://jsdokken.com/dolfinx-tutorial/chapter2/ns_code2.html>
4. Dokken, J. S. — *GMSH tutorial.* <https://jsdokken.com/src/tutorial_gmsh.html>
5. Chen, L. — *Finite Element Methods for Stokes Equations.* (to read)
6. Chen, L. — *Inf-Sup Conditions for Operator Equations.* (to read)
