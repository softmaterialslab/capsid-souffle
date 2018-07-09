# virus souffle

## Testing instructions

* load the necessary modules; module load gsl

* make the code; make clean

* if testing on a separate folder, copy 41part and/or 41part_c

* run the code for the following set of parameters for nose-hoover controlled md (engine selection, capsomere conc (microM), salt conc (mM), stretching constant (kBT), bending constant (kBT))
```echo m 200 0 100 20 | ./simulate_spinach_souffle```

* test the results by comparing energies in gary.traj.out (column1 is kinetic, col7 is total, col8 is potential)
