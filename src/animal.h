// animal.h
#ifndef ANIMAL_H
#define ANIMAL_H
#include "types.h"
#include <vector>
#include <array>
#include <random>
#include "cell.h"


template < typename cell_type>
class animal
{
    URNG& RandGenerator;
    t_InitParameter& rInitParameter;
    REAL KillTime;
    REAL AnimalTc;
    std::vector<cell_type> cells;
    std::array<int,2> LabelIndex{{0,0}};
    void create_cells(int mode);
    REAL create_inittime(REAL Tc);
    int to_many_cells;
public:
    animal(t_InitParameter& rInit,unsigned nKTime,REAL ATc);
    animal(t_InitParameter& rInit,unsigned nKTime,REAL ATc,int mode);
    ~animal(void);
    REAL get_killtime();
    std::string get_cell(void);
    std::string get_result_str(void);
    REAL get_result(void);
    void run();
};


template <>
REAL animal<cell_sym>::create_inittime(REAL Tc);
template <>
REAL animal<cell_sym>::get_result(void);


#include "animal.cpp"
#endif
