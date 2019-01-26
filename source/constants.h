/*constants.h
programed by Mike Williams 10-02-02
Contains constants that are needed in multiple files of the simulation.
*/
#ifndef constants_h
#define constants_h
#include<math.h>


//physicly meaningfull constants
const float g_c_fG = (float)6.67e-11;//Gravitational constant
const double g_c_dG = 6.67e-11;//Nm^2/kg^2 Gravitational consant
const double g_c_dMSun=1.99E30;//kg
const double g_c_dRSun=6.955E8;//m

//simulation scale units, in standard units.
const double g_c_dR0=(9.46E15)*1.5;//m is my default measure of lenght actualy it is about a light year


const double g_c_dL=g_c_dR0;
const double g_c_dM=g_c_dMSun;
const double g_c_dT=sqrt(g_c_dL*g_c_dL*g_c_dL/(g_c_dG*g_c_dM));

//const double g_c_dT=31536000;//s = 1 year

//Units to help me think
//devide a number by one of these to get that number of (name) in scaled units
//or multiply a scaled quantity by these numbers to get them in days or years or whatever
const double g_c_dDayConvert=g_c_dT/(3600*24);
const double g_c_dYearConvert=g_c_dT/(3600*24*365.25);


//scaled gravitational constant
const double g_c_dSG=g_c_dG*g_c_dM*g_c_dT*g_c_dT/(g_c_dL*g_c_dL*g_c_dL);//G*M*T^2/L^3 should be Force*lenght^2/Mas^2
const double g_c_dEnergyConverter=g_c_dT*g_c_dT/(g_c_dL*g_c_dL*g_c_dM);//T^2/(L^2*M) Jules = Energy/EnergyConverter

const int g_c_iStars = 100;
const double g_c_dPairs = (g_c_iStars*g_c_iStars-g_c_iStars)/2.0;
const double g_c_dDtMin=1E-7;
const double g_c_dDtMax=.001;
const double g_c_dClusterTimes=3;//1.0;//0.1;
const double g_c_dRSystem=1000;
//const double g_c_dRMinM=9.46E10;//m
const double g_c_dRMinM=g_c_dRSun*10;
const double g_c_dRMin=g_c_dRMinM/g_c_dL;
const double g_c_dRMin2=g_c_dRMin*g_c_dRMin;
const int g_c_iRMultipy=2;
const int g_c_iSmallTimeSteps=10;
const int g_c_iDtDevider=1000;
const double g_c_dDtSmall=.1;//how much smaller must he the local dt before I care
const int g_c_iNTimeStepsToBeBinary=0;
const double g_c_dTSystem =0;//= g_c_iNTimeStepsToBeBinary*g_c_dDt;//minimum time to be considered a binary in scaled units
const double g_c_dAdaptPower=0.2;
const double g_c_dError0m=9.46E2;//m 
const double g_c_dError0=g_c_dError0m/g_c_dL;//ensure that my error is unitless

//other constants
const double g_c_dPi=3.14159265;
const double g_c_dSK=1;//Scailed bolzman's constant. I don't think it matters so I am setting it to one

#endif
