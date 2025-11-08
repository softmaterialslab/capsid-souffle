# About
SOUFFLE is a molecular dynamics software for use in testing virus self-assemly with elastic capsomeres. Check out the [Wiki](https://github.com/softmaterialslab/capsid-souffle/wiki) for more detailed information!

Below are brief installation and run instructions. 


## NanoHUB app page:
* https://nanohub.org/tools/capsidsouffle


## Install and run instructions on Cluster
* First, git clone the project:
```git clone https://github.com/softmaterialslab/capsid-souffle.git```
* Then, load the required modules: boost, gsl and gnu computing environment.
* * These modules are installed for you in src/Makefile, currently configured for BigRed200 supercomputer
```module load gsl/2.8 && module load boost/1.86.0```
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
* Submit a batch script or run with these boost option parameters: ```./capsid-souffle -D b -f trimer_MVM -S 20 -d 100 -T 10000 -s 5000 -b 5000 -C 800 -c 100 -E 2.1 -t 0.005 -W 1000 -M 1000 -K 298```

## Install and run instructions on Local computer (Linux)
* Load the necessary modules: boost, gsl and openmpi
```module load gsl/2.8 && module load openmpi/3.0.1 && module load boost/1.86.0```
* Also make sure to export OMP_NUM_THREADS environment variable with maximum threads available in your CPU:
```export OMP_NUM_THREADS=16```
* Next, go to the root directory:
 ```cd capsid-souffle```
* Then, install the project:
```make install```
* Next, go to the bin directory:
 ```cd bin```
* Next, run the job:
``` time mpirun -np 1 -N 16 ./capsid-souffle -f trimer_MVM -S 20 -d 100 -T 10000 -s 5000 -b 5000 -C 800 -c 100 -E 2.1 -t 0.005 -W 1000 -M 1000 -K 298 ```
* All outputs from the simulation will be stored in the bin folder when the simulation is completed.
* Check and compare files (ex: energy.out) inside the ```bin/outfiles``` directory; model.parameters.out contains info on the system.
* If you want to clean everything and create a new build, use:
```make clean```
* If you want to clean the output data for a new run, use:
```make dataclean```

## Aditional information about different input parameter settings

* if testing on a separate folder, copy bin/infiles in addition to executable and job script
