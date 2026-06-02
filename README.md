# ns-pipe-expansion

2D incompressible Navier–Stokes simulation of laminar flow through a sudden 
pipe expansion (d₁ = 7.2 mm → d₂ = 17.2 mm), implemented in FEniCSx using 
the IPCS splitting scheme with Taylor–Hood (P2/P1) elements.

The loss coefficient K is back-calculated from the simulated pressure drop 
and compared against the analytical Borda–Carnot prediction for a sudden expansion:

$$K = \left(1 - \frac{A_1}{A_2}\right)^2$$

## Project Structure

```
ns-pipe-expansion/
├── src/
│   ├── main.py          # IPCS solver and Gmsh mesh generation
│   ├── postprocess.py   # pressure extraction and K back-calculation
│   └── dashboard.py     # Streamlit visualisation
├── assets/
│   └── flow.png         # sample output
├── output/              # solver outputs
├── .devcontainer/
│   └── devcontainer.json
├── requirements.txt
└── Technical_Log.md     # derivations, implementation notes, references
```

## Planned

- Grid convergence study across three mesh refinements
- Reynolds number sweep and K vs Re curve
- Quantified % deviation against Borda-Carnot analytical value

## Documentation

Technical derivations, weak formulation, IPCS scheme, and references: 
[Technical Log](Technical_Log.md)