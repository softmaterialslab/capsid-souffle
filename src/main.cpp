#include <iostream>
#include "md.h"


using namespace std;


int main(int argc, char *argv[]) {

	run_simulation(argc,argv);
}













/*  QUANTITY          SI UNITS           OTHER UNITS      CONVERSION               CODE/REDUCED UNIT
 *
 * mass             6.18e-24 kg           3722 amu          m* => m                      1 m*
 *
 * diameter           1.67e-9 m            16.7 A        L* = L(m)/d(m)                  1 L*
 *
 * energy           3.45E-21 J (Kb*~T)                    E* = E/KbT                     1 KbT
 *
 * ks                                                     k = E/L^2                      1 KbT
 *
 * time             7.06e-11 s            74.7 ps      t* = sqrt(m*d^2/E)                1 steps
 *
 * temperature      ~250K                                T* = (T*Kb)/E                   1 T*
 *
 * charge            1 e-                                  q* => q                       1 q*
 *
 * vacuum        8.854E-12 C^2s^2/(Nm^3)    C^2/J     E0*=E0*m*L^3/q^2*t^2            1.99*10^-3 E0*
 * permittivity
 *
 * Boltzmann: 1.3806E-23 J/K
 * e:         1.6022E-19 C
 * Er(water): 80.1
 * 1M = 0.60022 nm^-3
 *
 * 4.69E-5 m
 */
