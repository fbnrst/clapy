//animal.cpp
#include "animal.h"
#include <string>
#include <cmath>
#include "rand.h"

template <typename cell_type>
animal<cell_type>::animal(t_InitParameter &rInit, unsigned nKTime, REAL ATc)
    : RandGenerator(*rInit.urng), rInitParameter(rInit),
      KillTime(rInitParameter.KTimes[nKTime]), AnimalTc(ATc),to_many_cells(0) {
  create_cells(1);
}


template <typename cell_type>
animal<cell_type>::animal(t_InitParameter &rInit, unsigned nKTime, REAL ATc,int mode)
    : RandGenerator(*rInit.urng), rInitParameter(rInit),
      KillTime(rInitParameter.KTimes[nKTime]), AnimalTc(ATc),to_many_cells(0) {
  create_cells(mode);
}

template <typename cell_type> animal<cell_type>::~animal(void) {
  cells.clear();
}

template <typename cell_type> REAL animal<cell_type>::get_killtime() {
  return KillTime;
}

template <typename cell_type> void animal<cell_type>::create_cells(int mode) {
  long int ExpGenerataions = (long int)std::ceil(KillTime / AnimalTc) + 1;
  if (ExpGenerataions > 10){
    //std::cout << "to many cells\n" << "cellls:   " << rInitParameter.nCells * (1 << ExpGenerataions) <<"\n";  
    to_many_cells = 1;
    return ;
    }
  
    cells.reserve(rInitParameter.nCells * (1 << ExpGenerataions)); // ad some for grows

   //std::cout<< "Ktime: " << KillTime << " TC" << AnimalTc << "  ExpCells: " <<  rInitParameter.nCells*(1<<ExpGenerataions) << "\n";
  
  r_lognorm rr(AnimalTc, rInitParameter.sAnimal);
  for (int i = 0; i < rInitParameter.nCells; i++) {
    REAL Tc = rr(RandGenerator);
    REAL inittime = 0;
if (mode == 1){
    inittime = create_inittime(Tc);
}else{
    inittime = Tc*(-(REAL)i/(REAL)(rInitParameter.nCells-1));
}
    int GF_b = ber_func(RandGenerator, rInitParameter.GF);
    cell_type cell(i, -1, rInitParameter, AnimalTc, inittime, -1, KillTime, Tc,
                   GF_b, rInitParameter.fG1, rInitParameter.fS,
                   rInitParameter.fG2M);
    cells.push_back(cell);
  }
}

template <typename cell_type> REAL animal<cell_type>::create_inittime(REAL Tc) {
  return -unif_func(RandGenerator, 0, Tc);
}


template <> REAL animal<cell_sym>::create_inittime(REAL Tc) {
  return -Tc + fmod(exp_func(RandGenerator, Tc * 1.4426950408889634), Tc);
}

template <typename cell_type> void animal<cell_type>::run() {
  if (to_many_cells) return;
  size_t i = 0;
  do {
    cells[i].check_label();                  // label
    cells[i].divideordie(cells, LabelIndex); // divideorharvest
    i++;
  } while (i < cells.size());
  //std::cout << i <<"\n";
}

template <typename cell_type> std::string animal<cell_type>::get_cell(void) {
  std::string s;
  s.append("index,M_Index,Tc ,i_Age,LI,GF \n");
  for (unsigned i = 0; i < cells.size(); i++) {
    s.append(cells[i].get_info());
  }
  s.append("\n");
  return s;
}

template <typename cell_type>
std::string animal<cell_type>::get_result_str(void) {
  std::string s;
  s.append("Unlabeld cells: ");
  s.append(std::to_string(LabelIndex[0]));
  s.append("Labeld cells: ");
  s.append(std::to_string(LabelIndex[1]));
  s.append("LI: ");
  s.append(std::to_string((REAL)LabelIndex[1] /
                          (REAL)(LabelIndex[1] + LabelIndex[0])));
  s.append("\n");
  return s;
}

template <typename cell_type> REAL animal<cell_type>::get_result(void) {
  if (to_many_cells) return 1.0;
  return (REAL)hypgeo_func(RandGenerator, (unsigned int)LabelIndex[1],
                           (unsigned int)LabelIndex[0],
                           (unsigned int)rInitParameter.mCells) /
         (REAL)rInitParameter.mCells;
}

template <> REAL animal<cell_sym>::get_result(void) {
  return (REAL)hypgeo_func(RandGenerator, (unsigned int)LabelIndex[1],
                           (unsigned int)LabelIndex[0],
                           (unsigned int)rInitParameter.mCells) /
         (REAL)rInitParameter.mCells;
}


