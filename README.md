# virus souffle

## Install and run instructions on BigRed2
* First, git clone the project:
```git clone https://github.com/softmaterialslab/capsid-souffle.git```
* Then, load the required modules using following command:
```module swap PrgEnv-cray PrgEnv-gnu && module load boost/1.65.0 && module load gsl```
* Next, go to the root directory:
 ```cd capsid-souffle```
* Then, install the project:
```make cluster-install```
* Next, submit a test job:
```make cluster-test-submit```
* Then, clean the datafiles from the test job:
```make dataclean```
* Fianlly, submit the job:
```make cluster-submit```
* All outputs from the simulation will be stored in the bin folder when the simulation is completed.
* Check and compare files (ex: energy.out) inside the ```bin/outfiles``` directory; model.parameters.out contains info on the system.
* If you want to clean everything and create a new build, use:
```make clean```

## Install and run instructions on Local computer
* Load the necessary modules:
```module load gsl && module load openmpi/3.0.1 && module load boost/1_67_0```
* Also make sure to export OMP_NUM_THREADS environment variable with maximum threads available in your CPU:
```export OMP_NUM_THREADS=16```
* Next, go to the root directory:
 ```cd capsid-souffle```
* Then, install the project:
```make install```
* Next, go to the bin directory:
 ```cd bin```
* Next, run a test job:
``` time mpirun -np 2 -N 16 ./capsid-souffle -f 41part_c -T 10 -t 0.004 -e 12 -C 200 -c 500 -S 64 ```
* All outputs from the simulation will be stored in the bin folder when the simulation is completed.
* Check and compare files (ex: energy.out) inside the ```bin/outfiles``` directory; model.parameters.out contains info on the system.
* If you want to clean everything and create a new build, use:
```make clean```

## Aditional information about different input parameter settings

* if testing on a separate folder, copy 41part and/or 41part_c and/or 41part_cu

##### OPTIONS,                           **FLAG**  DEFAULT   
* engine selection,                   **-D**  m         
    * set to 'b' for brownian, 'm' for molecular dynamics
* filename,                           **-f**  41part
* capsomere conc (microM),            **-C**  75.0
* salt conc (mM),                     **-c**  200.0
* stretching constant (kBT),          **-s**  50.0
* bending constant (kBT),             **-S**  20.0
* total time (MD units),              **-T**  100       
    * # of computational steps = T/t
* timestep (MD units),                **-t**  0.001
* number subunits,                    **-S**  8         
* temperature (K),                    **-K**  298.0
* nose-hoover chain length,           **-q**  5         
    * to turn off thermostat set to 1
* electrostatics cut-off coefficient, **-e**  20.0      
    * to turn off electrostatic cut-off set to 0
* friction coefficient,               **-r**  1.0
* lennard jones attractive E_lj,	      **-E**  2.0
* Restart file bool,                   **-R**  false     
    * set to 'true' if restarting
* Neighbor list build frequency,       **-B**  20
* Neighbor list cutoff,                **-L**  4.0

* verbose                             **-v**
* help                                **-h**

