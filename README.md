# About
SOUFFLE is a molecular dynamics software for use in testing virus self-assemly with elastic capsomeres. Check out the [Wiki](https://github.com/softmaterialslab/capsid-souffle/wiki) for more detailed information!

Below are brief installation and run instructions. 


## NanoHUB app page:
* https://nanohub.org/tools/capsidsouffle


## Install and run instructions on Cluster(updated to bigred200)(2023/10/24)
* First, git clone the project:
```git clone https://github.com/softmaterialslab/capsid-souffle.git```
* Then, load the required modules: boost, gsl and gnu computing environment.
* * These modules are installed for you in src/Makefile, currently configured for BigRed200 supercomputer
```module load gsl/2.7 && module load boost/1.78.0```
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

## Install and run instructions on Local computer (Linux)(Not updated)
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
