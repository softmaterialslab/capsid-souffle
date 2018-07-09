# virus souffle

## Testing instructions

* load the necessary modules; module load gsl

* make the code; make clean

* if testing on a separate folder, copy 41part and/or 41part_c

* Enter the following set of parameters for nose-hoover controlled md (engine selection, capsomere conc (microM), salt conc (mM), stretching constant (kBT), bending constant (kBT))
```echo m 200 0 100 20 | ./simulate_spinach_souffle```
