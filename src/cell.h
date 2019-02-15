// cell.h
#ifndef CELL_H
#define CELL_H
#include "types.h"
#include <vector>
#include <array>


class cell_base
{
protected:
    int Index;
    int Mother;
    t_InitParameter& rInitParameter;
    REAL InitAge;
    int  Label;
    REAL KillTime;
    REAL Tc;
    int GF;
    REAL fG1;
    REAL fS;
    REAL fG2M;
    REAL AnimalTc;
    cell_base(int Index,int Mother, t_InitParameter& rInit, REAL AnimalTc,REAL IAge , int Label, REAL KTime ,REAL Tc, int GF,REAL fG1,REAL fS,REAL fG2M);
    cell_base( const cell_base &obj );
    ~cell_base();
public:
    std::string get_info(void);
};

class cell_as : public  cell_base
{
public:
    cell_as(int Index,int Mother, t_InitParameter& rInit, REAL AnimalTc, REAL IAge , int Label, REAL KTime ,REAL Tc, int GF,REAL fG1,REAL fS,REAL fG2M);
    cell_as( const cell_as &obj );
    ~cell_as();
    void check_label(void);
    void divideordie(std::vector<cell_as>& cells, std::array<int,2>& LabelIndex);
};

class cell_as_newd : public  cell_as
{
public:
    cell_as_newd(int Index,int Mother, t_InitParameter& rInit, REAL AnimalTc, REAL IAge , int Label, REAL KTime ,REAL Tc, int GF,REAL fG1,REAL fS,REAL fG2M);
    cell_as_newd( const cell_as_newd &obj );
    ~cell_as_newd();
    void divideordie(std::vector<cell_as_newd>& cells, std::array<int,2>& LabelIndex);
};

class cell_sym : public  cell_base
{
public:
    cell_sym(int Index,int Mother, t_InitParameter& rInit, REAL AnimalTc, REAL IAge , int Label, REAL KTime ,REAL Tc, int GF,REAL fG1,REAL fS,REAL fG2M);
    cell_sym( const cell_as &obj );
    ~cell_sym();
    void check_label(void);
    void divideordie(std::vector<cell_sym>& cells, std::array<int,2>& LabelIndex);
};


#endif
