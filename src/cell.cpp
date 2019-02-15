//cell.c
#include "cell.h"
#include <string>
#include "rand.h"

cell_base::cell_base(int Index, int Mother, t_InitParameter &rInit,
                     REAL AnimalTc, REAL IAge, int Label, REAL KTime, REAL Tc,
                     int GF, REAL fG1, REAL fS, REAL fG2M)
    : Index(Index), Mother(Mother), rInitParameter(rInit), InitAge(IAge),
      Label(Label), KillTime(KTime), Tc(Tc), GF(GF), fG1(fG1), fS(fS),
      fG2M(fG2M), AnimalTc(AnimalTc) {
  // std::cout << "create" <<  Index  << " "<<InitAge<<"\n";
}

cell_base::cell_base(const cell_base &obj)
    : Index(obj.Index), Mother(obj.Mother), rInitParameter(obj.rInitParameter),
      InitAge(obj.InitAge), Label(obj.Label), KillTime(obj.KillTime),
      Tc(obj.Tc), GF(obj.GF), fG1(obj.fG1), fS(obj.fS), fG2M(obj.fG2M),
      AnimalTc(obj.AnimalTc) {
  // std::cout << "copy" << Index  << "\n";

}
cell_base::~cell_base() {}

std::string cell_base::get_info() {
  std::string s;
  s.append(std::to_string(Index));
  s.append(";");
  s.append(std::to_string(Mother));
  s.append(";");
  s.append(std::to_string(Tc));
  s.append(";");
  s.append(std::to_string(InitAge));
  s.append(";");
  s.append(std::to_string(GF));
  s.append(";");
  s.append(std::to_string(Label));
  s.append("\n");
  return s;
}

cell_as::cell_as(int Index, int Mother, t_InitParameter &rInit, REAL AnimalTc,
                 REAL IAge, int Label, REAL KTime, REAL Tc, int GF, REAL fG1,
                 REAL fS, REAL fG2M)
    : cell_base(Index, Mother, rInit, AnimalTc, IAge, Label, KTime, Tc, GF, fG1,
                fS, fG2M) {}

cell_as::cell_as(const cell_as &obj) : cell_base(obj) {}

cell_as::~cell_as() {}

void cell_as::check_label(void) {
  if (Label == -1) {
    Label = 0;
    if (InitAge > Tc * (-1 + fG1 + fS))
      return;
  }
  if (Label == 0 && GF) {
    if (InitAge + KillTime > Tc * (-1 + fG1))
      Label = 1;
    return;
  }
  return;
}

void cell_as::divideordie(std::vector<cell_as> &cells,
                          std::array<int, 2> &LabelIndex) {
  if (InitAge + KillTime > 0 && GF) {
    int LastCell = cells.size();
    cell_as a(LastCell, Index, rInitParameter, AnimalTc, InitAge - Tc, Label,
              KillTime, Tc, GF, fG1, fS, fG2M);
    if (cells.size()> 10000 ){
         LabelIndex[Label]++;
    }else{
    cells.push_back(a);
     }
  } else {
    LabelIndex[Label]++;
  }
}

cell_as_newd::cell_as_newd(int Index, int Mother, t_InitParameter &rInit, REAL AnimalTc,
                 REAL IAge, int Label, REAL KTime, REAL Tc, int GF, REAL fG1,
                 REAL fS, REAL fG2M)
    : cell_as(Index, Mother, rInit, AnimalTc, IAge, Label, KTime, Tc, GF, fG1,
                fS, fG2M) {}

cell_as_newd::cell_as_newd(const cell_as_newd &obj) : cell_as(obj) {}

cell_as_newd::~cell_as_newd() {}

void cell_as_newd::divideordie(std::vector<cell_as_newd> &cells,
                          std::array<int, 2> &LabelIndex) {
  if (InitAge + KillTime > 0 && GF) {
    int LastCell = cells.size();
    double nTc =
        lognorm_func(*rInitParameter.urng, AnimalTc, rInitParameter.sAnimal);
    //std::cout << "AnimalTC: " << AnimalTc << "   std: "<<  rInitParameter.sAnimal << "\n";
    cell_as_newd a(LastCell, Index, rInitParameter, AnimalTc, InitAge - nTc, Label,
              KillTime, nTc, GF, fG1, fS, fG2M);
    if (cells.size()> 10000 ){
         LabelIndex[Label]++;
    }else{
    cells.push_back(a);
    }
  } else {
    LabelIndex[Label]++;
  }
}


cell_sym::cell_sym(int Index, int Mother, t_InitParameter &rInit, REAL AnimalTc,
                   REAL IAge, int Label, REAL KTime, REAL Tc, int GF, REAL fG1,
                   REAL fS, REAL fG2M)
    : cell_base(Index, Mother, rInit, AnimalTc, IAge, Label, KTime, Tc, GF, fG1,
                fS, fG2M) {}

cell_sym::cell_sym(const cell_as &obj) : cell_base(obj) {}

cell_sym::~cell_sym() {}

void cell_sym::check_label(void) {
  if (Label == -1) {
    Label = 0;
    if (InitAge > Tc * (-1 + fG1 + fS))
      return;
  }
  if (Label == 0 && GF) {
    if (InitAge + KillTime > Tc * (-1 + fG1))
      Label = 1;
    return;
  }
  return;
}

void cell_sym::divideordie(std::vector<cell_sym> &cells,
                           std::array<int, 2> &LabelIndex) {
  if (InitAge + KillTime > 0 && GF) {
    int LastCell = cells.size();
    // r_lognorm rr(AnimalTc,rInitParameter.sAnimal);
    // REAL nTc = rr(*rInitParameter.urng);
    double nTc =
        lognorm_func(*rInitParameter.urng, AnimalTc, rInitParameter.sAnimal);
    cell_sym a(LastCell, Index, rInitParameter, AnimalTc, InitAge - nTc, Label,
               KillTime, nTc, GF, fG1, fS, fG2M);
    cell_sym b(LastCell + 1, Index, rInitParameter, AnimalTc, InitAge - nTc,
               Label, KillTime, nTc, GF, fG1, fS, fG2M);
    if (cells.size()> 10000 ){
         LabelIndex[Label]++;
    }else{
    cells.push_back(a);
    cells.push_back(b);
    }
  } else {
    LabelIndex[Label]++;
  }
}
