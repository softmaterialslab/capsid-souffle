# virus souffle

## Testing instructions

* load the necessary modules; module load gsl

* make dataclean; make install; make clean

* if testing on a separate folder, copy 41part and/or 41part_c

* run the code for the following set of parameters for nose-hoover controlled md (engine selection, filename, capsomere conc (microM), salt conc (mM), stretching constant (kBT), bending constant (kBT), total time (MD units), timestep (MD units)
'echo m 41part 931 0 100 20 100 0.002 | ./simulate_capsid_souffle'
## echo m 41part 200 0 100 20 100 0.002 | ./simulate_spinach_souffle

* test the results by comparing energies in outfiles/energy.out (column1 is kinetic, col7 is total, col8 is potential)

* To run with electrostatics w/o salt screening, use:
echo m 41part_c 931 9 100 20 100 0.002 | ./simulate_capsid_souffle

* To run with electrostatics w/ moderate salt screening, use:
echo m 41part_c 931 230 100 20 100 0.002 | ./simulate_capsid_souffle

* To run with brownian dynamics w/ moderate salt screening, use:
echo b 41part_c 931 230 100 20 100 0.002 1 | ./simulate_capsid_souffle
