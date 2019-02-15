#include "types.h"

std::ostream& operator<<(std::ostream& os, const t_InitParameter& obj)
{
os << "nAnimals:\t"  << obj.nAnimals << "\n";
os << "nCells:\t"  << obj.nCells << "\n";
os << "mCells:\t" << obj.mCells << "\n";
os << "MTimes:\t";
for (auto el : obj.KTimes){os << el <<", "  ;}
os << "\n";
os << "GF:\t"  << obj.GF << "\n";
os << "sPopulation:\t"  << obj.sPopulation << "\n";
os << "sAnimal:\t"  << obj.sAnimal << "\n";
os << "Tc:\t"  << obj.Tc << "\n";
os << "fG1:\t"  << obj.fG1 << "\n";
os << "fS:\t"  << obj.fS << "\n";
os << "fG2M:\t"  << obj.fG2M << "\n";
    return os;
}


