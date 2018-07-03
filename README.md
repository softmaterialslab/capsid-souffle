# virus souffle

## Testing instructions

* Enter the following set of parameters for nose-hoover controlled md (engine selection, total time (units), LJ energy cutoff (units), salt valency, capsomere conc (microM), salt conc, number of capsomeres, stretching constant, bending constant, input filename, timestep, nose-hoover chain length, nose-hoover mass, temperature)
```echo m 100 2.5 1 200 0 8 100 10 41part 0.002 5 1 1 | ./simulate_spinach_souffle```