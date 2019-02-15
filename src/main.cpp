#include <stdio.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <exception>
#include <random>

#include "types.h"
#include "animal.h"
#include "cell.h"
#include "rand.h"


int main(){
 t_InitParameter InitParameter;
  std::vector<REAL> Result;
  int seed=10; 
  // parse arguments
                                    InitParameter.nAnimals  = 1000;
                                    InitParameter.nCells = 100;
                                    InitParameter.GF = 1;
                                    InitParameter.Tc = 1;
                                    InitParameter.fG1 =(REAL) 0.2;
                                    InitParameter.fS =(REAL) 0.3;
                                    InitParameter.fG2M =(REAL) 0.5;
                                    InitParameter.sAnimal = (REAL)0.3;
                                    InitParameter.sPopulation = (REAL)0.3;

 
  InitParameter.KTimes.push_back((REAL)0.0);
  InitParameter.KTimes.push_back((REAL)0.2);
  InitParameter.KTimes.push_back((REAL)0.4);
  InitParameter.KTimes.push_back((REAL)0.6);
  InitParameter.KTimes.push_back((REAL)0.8);
  InitParameter.KTimes.push_back((REAL)1.0);
  InitParameter.KTimes.push_back((REAL)1.2);
  static URNG generator ;
  generator.seed(seed);
  InitParameter.urng = &generator;
  r_lognorm rr(1,InitParameter.sPopulation);
  for (int aa=0;aa<2;aa++){
      for (int pp=0;pp<2;pp++){
          InitParameter.sAnimal = (REAL)(aa+0.0001)/100;
          InitParameter.sPopulation = (REAL)(pp+0.0001)/100;
  for(unsigned i=0;i<InitParameter.KTimes.size();i++){
      for (int j=0;j<InitParameter.nAnimals;j++){
        animal<cell_sym> a(InitParameter,i, rr(generator));
        a.run();
        Result.push_back(a.get_result());
      }
  }
  }
  }
  return 0;


}
