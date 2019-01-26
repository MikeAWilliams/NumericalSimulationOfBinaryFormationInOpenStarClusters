//distrobution.h
//by Mike Williams 2-2-05

#ifndef distrobution_h
#define distrobution_h
 
#include<math.h>
#include"constants.h"
#include"star6_0.h"

//spatial distrobution constants
const int g_c_iNBins=15;


//spatial functions
double ro(double a_dR,double a_dR0);
double rMid(int a_iI,double a_dR0);
int nStarsInBinR(double a_dR,double a_dR0);
double rOne(double a_dR0);

//more appropriate boltzmann shaped spatial functions
double roR2(double a_dR,double a_dR0);
double rMidR2(int a_iI,double a_dR0,double a_dROne);
int nStarsInBinRR2(double a_dR,double a_dR0,double a_dROne);
double rOneR2(double a_dR0,double a_dRInit);
double tailR2p(double a_dR,double a_dR0);
double tailR2(double a_dR,double a_dR0);

//energy functions
double boltzmann(double a_dE,double a_dT);
double eMid(int a_iI,double a_dT);
int nStarsInBinE(double a_dE,double a_dT);
double eOne(double a_dT);

//energy functions for v^2Exp[-E/kt]
double boltzmannV2(double a_dE,double a_dT);
double eMid(int a_iI,double a_dT,double a_dE1);
int nStarsInBinEV2(double a_dE,double ad_T,double a_dE1);
double eOneV2(double a_dT,double a_dE0);
double tailV2p(double a_dE, double a_dT);
double tailV2(double a_dE,double a_dT);

//mass functions
double massPower(double a_dM,double a_dMMin, double a_dMMax);
double mMid(int a_iI,double a_dMMin, double a_dMMax);
int nStarsInBinM(double a_dM,double a_dMMin, double a_dMMax);

void uniformMassDistribute(star *stars[]);
void powerMassDistribute(star *stars[],double a_dMMin, double a_dMMax);

void spatialyDistributeMid(star *stars[],const double a_dR0);
void spatialyDistributeMidR2(star *stars[],const double a_dR0, double &a_dROne);
void energyDistribute(star *stars[],double &a_dT);
void energyDistributeV2(star *stars[],double &a_dT,double &a_dE1);
void energyDistribute(star *stars[]);
void energyDistributeV2(star *stars[]);

#endif
