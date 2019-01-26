//energy2_0.h Mike Williams 2/05
#ifndef _energy_h
#define _energy_h

#include"star6_0.h"
#include"constants.h"

double energyNotVerlet(star *stars[],double Rs[g_c_iStars][g_c_iStars],const int &n);
double energyNotVerlet(star *stars[],double Rs[g_c_iStars][g_c_iStars],const int &n,double &PE, double &KE);
double energyNotVerlet(star *stars[],const int &n,double &PE, double &KE);
double energyP(star *stars[]);
double energy(star *stars[]);
double energy(star *stars[],double &dPE,double &dKE);
//#include"energy1_0.cpp"

#endif
