// rand.h
#ifndef RAND_H
#define RAND_H
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define URNG generator


class generator{
    const gsl_rng_type * T;
    public:
    gsl_rng * r;
    generator();
    ~generator();
    void seed(unsigned long int s);
};

#include "types.h"


class r_lognorm{
    double zeta;
    double sigma;
    public:
    r_lognorm(REAL mean,REAL sigma);
    ~r_lognorm();
    REAL operator()(URNG& g);
};

REAL lognorm_func(URNG& g,REAL mean,REAL sig);

class r_hypgeo{
    unsigned int n1;
    unsigned int n2;
    unsigned int t;
    public:
    r_hypgeo(unsigned int population1,unsigned int population2,unsigned int samples);
    ~r_hypgeo();
    unsigned int  operator()(URNG& g);
};

unsigned int hypgeo_func(URNG& g,unsigned int population1,unsigned int population2,unsigned int samples);

class r_uni{
    double a;
    double b;
    public:
    r_uni(REAL min,REAL max);
    ~r_uni();
    REAL operator()(URNG& g);
};

REAL unif_func(URNG& g,REAL min,REAL max);

REAL exp_func(URNG& g,REAL mu);

int ber_func(URNG& g ,REAL p);


#endif


