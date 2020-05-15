# CLAWPACK-solver-with-CO2-table
Couple CLAWPACK solver with CO2 look-up table

[200~This is the CLAWPACK solver. Here we did not use the latest version of CLAWPACK (written by fortran 97), but it is compatible with the CO2 table written in Fortran90.

Basically, the main codes are in claw/clawpack/1d_CO2 and claw/clawpack/2d_CO2. The CO2 table is in claw/clawpack/CO2. In this version of CO2 table, we added some variables which are used for debugging and boundary condition.

1. set path in your computer for CLAWPACK solver.
2. compile CO2 table
3. compile solver

There are some test cases in 1d and 2d that you can find in the solver fichers. 
