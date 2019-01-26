//distroIO1_1.h
//by Mike Williams 04-15-05
//load and save star distrobutions
//simplere than 1_0 which tried to save units and such

#ifndef distroIO_h
#define distroIO_h
 
#include"star6_0.h"
#include<fstream>

void saveDistro(star *stars[],std::ofstream &disOut, double a_dTime, double a_dTimeStep);
void loadDistro(star *stars[],std::ifstream &disIn);

void saveDistroRealUnits(star *stars[],std::ofstream &disOut);
void loadDistroRealUnits(star *stars[],std::ifstream &disIn);

void saveConst(std::ofstream &constOut,char sAdapt[]);

#endif
