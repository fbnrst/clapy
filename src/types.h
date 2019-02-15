// types.h
#ifndef TYPES_H
#define TYPES_H
#include <vector>
#include <iostream>

typedef float REAL;

#include "rand.h"

typedef struct t_InitParameter{
    int nAnimals;
    int nCells;
    int mCells;
    REAL GF;
    REAL sPopulation;
    REAL sAnimal;
    REAL Tc;
    REAL fG1;
    REAL fS;
    REAL fG2M;
    URNG *urng;
    std::vector<REAL> KTimes;
}t_InitParameter;


std::ostream& operator<<(std::ostream& os, const t_InitParameter& obj);

#endif
