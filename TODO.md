# TODO

Personal working notes — not user-facing documentation.

## Code cleanup / refactoring

- `src/postprocessing/postprocess.py` is currently just a header comment — implement
  pressure extraction and the $K$ back-calculation it describes.
- `src/solvers/1_steady_state_solver.py` (`newton_solver`) is a stub — see
  [Steady-state solver](technical_log/README_steady_state.md). Filenames/numbering
  (`1_steady_state_solver.py`, `2_time_dependant_sovler.py`) and the inconsistencies
  inside that file are known and fine for now.
- `src/solvers/2_time_dependant_sovler.py`: pressure space is degree-2 Lagrange
  (`s_cg1` variable name suggests P1 was intended) — decide whether to switch to a
  true Taylor-Hood P2/P1 pairing.
- `src/solvers/2_time_dependant_sovler.py::main`: `U_max = 3 / (2 * v_avg)` looks
  inverted — `v_avg = (2/3) * U_max` implies `U_max = (3/2) * v_avg`, not
  `3 / (2 * v_avg)`. Double-check.
- Re-add a dashboard (Streamlit/PyVista) reporting average inlet/outlet velocity and
  Reynolds number (`src/dashboard.py` was removed in cleanup).

## Steady-state solver

- Implement the Newton-based steady-state Navier-Stokes solver
  (`src/solvers/1_steady_state_solver.py`).
- Write up the corresponding technical log
  ([technical_log/README_steady_state.md](technical_log/README_steady_state.md)).

## Major future work

- Grid convergence study across three mesh refinements.
- Reynolds number sweep and $K$ vs. $Re$ curve.
- Quantify % deviation of simulated $K$ against the Borda-Carnot analytical value.
