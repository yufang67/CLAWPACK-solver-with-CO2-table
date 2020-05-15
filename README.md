# CLAWPACK solver coupled with CO2 look-up table
This is the CLAWPACK solver. Here we did not use the latest version of CLAWPACK (written by fortran 97), but it is compatible with the CO2 table written in Fortran90.

The original solver is in claw/clawpack. The CO2 table is coupled with HLLC solver in each cases (in rp_HLLC). In this version of CO2 table, we added some variables which are used for debugging and boundary condition.

---
1. set path in your computer for CLAWPACK solver.
2. compile CO2 table
3. compile original solver
4. compile cases

Two test cases in 1d and 2d are presented which you can play with. 
