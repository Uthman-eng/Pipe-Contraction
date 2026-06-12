# ns-pipe-expansion

2D incompressible Navier–Stokes simulation of laminar flow through a sudden
pipe expansion (d₁ = 7.2 mm → d₂ = 17.2 mm), implemented in FEniCSx using
the IPCS splitting scheme with Taylor–Hood-style (P2/P2) elements.

The loss coefficient K is back-calculated from the simulated pressure drop
and compared against the analytical Borda–Carnot prediction for a sudden expansion:

$$K = \left(1 - \frac{A_1}{A_2}\right)^2$$

## Project Structure

```
ns-pipe-expansion/
├── src/
│   ├── params.py                      # geometry and fluid parameters
│   ├── mesh/
│   │   └── mesh.py                    # Gmsh mesh generation and boundary tagging
│   ├── solvers/
│   │   ├── 2_time_dependant_sovler.py # IPCS time-dependent Navier-Stokes solver
│   │   └── 1_steady_state_solver.py   # Newton steady-state solver (not yet implemented)
│   └── postprocessing/
│       └── postprocess.py             # pressure extraction and K back-calculation
├── output/                            # solver outputs (figures, data)
├── technical_log/
│   ├── README_time_dependent.md       # IPCS solver: derivations and implementation notes
│   └── README_steady_state.md         # Newton solver: placeholder
├── TODO.md                            # personal task tracking
├── .devcontainer/
│   └── devcontainer.json
└── requirements.txt
```

## Documentation

Technical derivations, weak formulations, and implementation notes:

- [Time-dependent Navier-Stokes (IPCS)](technical_log/README_time_dependent.md)
- [Steady-state Navier-Stokes (Newton)](technical_log/README_steady_state.md) — placeholder, not yet implemented

Personal task tracking: [TODO.md](TODO.md)
