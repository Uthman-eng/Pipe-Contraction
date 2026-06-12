from dataclasses import dataclass


@dataclass(frozen=True)
class Parameters:
    """
    Geometry and fluid properties for the pipe flow problem.
    All lengths in metres (m), density in kg/m^3, viscosity in Pa·s, velocity in m/s.
    """
    d1: float = 7.2 / 1000              # inlet diameter [m]
    d2: float = 17.2 / 1000             # outlet diameter [m]
    pipe_length: float = 100 / 1000     # pipe length [m]
    spacing: float = 1 / 1000           # mesh spacing (circular) [m]
    rho: float = 1000.0                 # water density [kg/m^3]
    mu: float = 0.001                   # dynamic viscosity [Pa·s]
    v_avg: float = 0.2                  # average inlet velocity [m/s]
    T = 60*5                            # Final time
    dt = 1 / 20                       # Time step size