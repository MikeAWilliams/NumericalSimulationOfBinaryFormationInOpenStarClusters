//distroStats.h
//by Mike Williams 04-25-05
//Calculate thins ling Rrms Vrms DelRrms

#ifndef distroStats_h
#define distroStats_h

#include"star6_0.h"
#include"Vector3d.h"

void RVDRrms(star *stars[],double &dRrms,double &dVrms,double &dDRrms);
void RVDRrmsL(star *stars[],double &dRrms,double &dVrms,double &dDRrms,Vector3d &vL);

double AveMass(star *stars[]);

#endif
