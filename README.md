# Fortran-sampler
This repository holds some basic CFD / PDE solving scripts encountered in a MSc CFD course.\
All scripts are written in Fortran90 and result in .dat format intended for Tecplot.\

## Indvidual Script Description
1. 1D-compressible-euler.f90
   - Base program from Dr, Titarev & Dr. P. Tsoutsanis.
   - Contribution in implementation of following fluxes:
     - HLL & HLLC

2. Parallel_solve_PDE.f90
   - Sole contributor
   - Parallel solving of 1D inviscid Burger's equation with MPI.
   - Measures run time of each processor.
   - Following numerical schemes implemented:
     - Forward Time Centered Space (FTCS)
     - Upwind
     - Lax-Friedrichs
     - Lax-Wendroff
     - MacCormack

3. SIMPLE-TDMA-LidDriven.f90
   - Base program from Dr. Konozsy.
   - SIMPLE solver with Tri-Diagonal Matrix Algo. solver on structured grid lid driven flow.
   - Contribution in modifying program from channel flow case to lid driven flow case.
