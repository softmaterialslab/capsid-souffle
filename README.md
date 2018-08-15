# virus souffle

## Testing instructions

* load the necessary modules; module load gsl

* make dataclean; make install; make clean

* if testing on a separate folder, copy 41part and/or 41part_c

* run the code for the following set of parameters for nose-hoover controlled md (engine selection, filename, capsomere conc (microM), salt conc (mM), stretching constant (kBT), bending constant (kBT), total time (MD units), timestep (MD units)
```./simulate_capsid_souffle -D m -f 41part -C 931 -c 0 -s 100 -b 20 -T 100 -t 0.002```

* test the results by comparing energies in outfiles/energy.out (column1 is kinetic, col7 is total, col8 is potential)

* To run with electrostatics w/o salt screening, use:
```./simulate_capsid_souffle -D m -f 41part_c -C 931 -c 9 -s 100 -b 20 -T 100 -t 0.002```

* To run with electrostatics w/ moderate salt screening, use:
```./simulate_capsid_souffle -D m -f 41part_c -C 931 -c 230 -s 100 -b 20 -T 100 -t 0.002```

* To run with brownian dynamics w/ moderate salt screening, use:
```./simulate_capsid_souffle -D m -f 41part_c -C 931 -c 230 -s 100 -b 20 -T 100 -t 0.002 -r 1```

* check the outfiles/ folder for information; model.parameters.out contains info on the system
