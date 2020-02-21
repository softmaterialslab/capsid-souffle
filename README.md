# About
SOUFFLE is a molecular dynamics software for use in testing virus self-assemly with elastic capsomeres.


## Install and run instructions on Cluster
* First, git clone the project:
```git clone https://github.com/softmaterialslab/capsid-souffle.git```
* Then, load the required modules: boost, gsl and gnu computing environment.
* * These modules are installed for you in src/Makefile, currently configured for BigRed3 supercomputer
```module swap PrgEnv-intel/6.0.5 PrgEnv-gnu && module load boost/gnu && module load gsl```
* Next, go to the root directory:
 ```cd capsid-souffle```
* Then, install the project:
```make cluster-install```
* Submit an individual job:
```make cluster-submit```
* Use parameter_sweep.sh to submit multiple jobs sweeping over a parameter
```./parameter_sweep```
* All outputs from the simulation will be stored in the bin folder when the simulation is completed.
* Check and compare files (ex: energy.out) inside the ```bin/outfiles``` directory; model.parameters.out contains info on the system.
* If you want to clean everything and create a new build, use:
```make clean```
* If you want to clean the output data for a new run, use:
```make dataclean```

## Install and run instructions on Local computer (Linux)
* Load the necessary modules: boost, gsl and openmpi
```module load gsl && module load openmpi/3.0.1 && module load boost/1_67_0```
* Also make sure to export OMP_NUM_THREADS environment variable with maximum threads available in your CPU:
```export OMP_NUM_THREADS=16```
* Next, go to the root directory:
 ```cd capsid-souffle```
* Then, install the project:
```make install```
* Next, go to the bin directory:
 ```cd bin```
* Next, run the job:
``` time mpirun -np 1 -N 16 ./capsid-souffle -f trimer_MVM -S 20 -T 1000000 -s 500 -b 500 -C 520 -c 1000 -E 1.9 ```
* All outputs from the simulation will be stored in the bin folder when the simulation is completed.
* Check and compare files (ex: energy.out) inside the ```bin/outfiles``` directory; model.parameters.out contains info on the system.
* If you want to clean everything and create a new build, use:
```make clean```
* If you want to clean the output data for a new run, use:
```make dataclean```

## Aditional information about different input parameter settings

* if testing on a separate folder, copy bin/infiles in addition to executable and job script

#### OPTIONS                           FLAG  DEFAULT   
* engine selection,                   **-D**  m         
    * set to 'b' for brownian, 'm' for molecular dynamics
* filename,                           **-f**  41part_c
* capsomere conc (microM),            **-C**  200.0
* salt conc (mM),                     **-c**  500.0
* stretching constant (kBT),          **-s**  50.0
* bending constant (kBT),             **-b**  20.0
* total time (Computational Steps),   **-T**  2500       
    * number of computational steps = T/t
* timestep (MD units),                **-t**  0.004
* number subunits,                    **-S**  64         
* temperature (K),                    **-K**  298.0
* nose-hoover chain length,           **-q**  5         
    * to turn off thermostat set to 1
* electrostatics cut-off coefficient, **-e**  12.0      
    * to turn off electrostatic cut-off set to 0
* friction coefficient,               **-r**  1.0
* lennard jones attractive E_lj,	  **-E**  2
* Restart file bool,                  **-R**  true     
    * set to 'false' to bypass all restart files
* Neighbor list build frequency,      **-B**  20
* Neighbor list cutoff,               **-L**  10.0

* verbose                             **-v**
* help                                **-h**


#### Premade input files
* 41part_c -- this is a rhombic dodecahedron subunit (41part has no charge,  41part_cu has a uniform charge pattern)
* 43part_* -- these files are test files for dimeric T=1 and T3/T4 HBV systems
* 51part_c -- this is an HBV-like subunit that assembles hexamer sheets
* trimer_MVM -- this is an MVM subunit
* xyz_T4 -- this is a cluster file for HBV T3 assembly