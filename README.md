# Pipe Contraction Flow — FEniCSx

2D incompressible Navier Stokes simulation of flow through a sudden pipe contraction (17.2mm to 7.2mm), implemented in FEniCSx.

Numerically back calculate the loss coefficient $K$ from the simulated pressure field and compare against experimental measurements from a hydraulic bench lab experiment.

$$K = \frac{\Delta p / \rho g}{v^2 / 2g}$$

## References

Dokken, J.S. (2023). FEniCSx Tutorial. https://jsdokken.com/dolfinx-tutorial