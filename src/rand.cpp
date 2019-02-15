//rand.c
#include "rand.h"
#include <cmath>

generator::generator(){
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T); 
}
generator::~generator(){
    gsl_rng_free (r);
}
void generator::seed(unsigned long int s){
    gsl_rng_set (r, s);
}

r_lognorm::r_lognorm(REAL mean,REAL sig){
    sigma =(double) sqrt(log((sig*sig)/(mean*mean)+1));
    zeta = (double) log(mean)-sigma*sigma/2;
}
r_lognorm::~r_lognorm(){
}

REAL r_lognorm::operator()(URNG& g){
    return (REAL)gsl_ran_lognormal(g.r,zeta,sigma);
}

REAL lognorm_func(URNG& g,REAL mean,REAL sig){
    double sigma =(double) sqrt(log((sig*sig)/(mean*mean)+1));
    double zeta = (double) log(mean)-sigma*sigma/2;
    return (REAL)gsl_ran_lognormal(g.r,zeta,sigma);
}

r_hypgeo::r_hypgeo(unsigned int population1,unsigned int population2,unsigned int samples):
    n1(population1),
    n2(population2),
    t(samples){}

r_hypgeo::~r_hypgeo(){}

unsigned int r_hypgeo::operator()(URNG& g){
    return gsl_ran_hypergeometric(g.r,n1,n2,t);
}

unsigned int hypgeo_func(URNG& g,unsigned int population1,unsigned int population2,unsigned int samples){
    return gsl_ran_hypergeometric(g.r,population1,population2,samples);
}
r_uni::r_uni(REAL min,REAL max):
    a(min),
    b(max)
{
}
r_uni::~r_uni(){};

REAL r_uni::operator()(URNG& g){
    return (REAL)gsl_ran_flat (g.r,a,b);
}
REAL unif_func(URNG& g,REAL min,REAL max){
    return (REAL)gsl_ran_flat (g.r,min,max);
}
REAL exp_func(URNG& g ,REAL mu){
return (REAL) gsl_ran_exponential (g.r,mu);
}

int ber_func(URNG& g ,REAL p){
    return gsl_ran_bernoulli(g.r,p);
}