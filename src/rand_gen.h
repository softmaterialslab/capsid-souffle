// This is header file for the MY_GENERATE_RAND class.
// This class sets up the needed structure to generate random numbers

#ifndef LEMONSOUFFLE_RAND_GEN_H
#define LEMONSOUFFLE_RAND_GEN_H


#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

class RAND_GEN
{
public:
    unsigned long int seed;
    const gsl_rng_type * T;
    gsl_rng * r;

    RAND_GEN()
    {
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);

        gsl_rng_env_setup();

        // NOTE: if you want to start with the same initial configuration comment three lines below
        // else uncomment them. uncommenting will generate a new randomized distribution of particles every run

//        srand((time_t(0)));                	// srand & time are built-in
//        unsigned long int s = random();  	// gsl_rng_uniform will eventually
//        gsl_rng_set(r,s); 		            // seed the random number generator;
    }

    ~RAND_GEN()
    {
        gsl_rng_free(r);
    }
};

#endif //LEMONSOUFFLE_RAND_GEN_H
