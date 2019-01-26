//rk2_1,cpp Mike Williams 2-05
#ifndef _rk_h
#define _rk_h

#include"star6_0.h"
#include"constants.h"
#include"analytical.h"

void rk4SoftNC(star *stars[], double Rs[g_c_iStars][g_c_iStars], const int &n,  const double &dt, double &Effort, const double &c_dRMin);
void rk42StepSoft(star *stars[],star *starsp[],star *stars1[], double dRs[g_c_iStars][g_c_iStars],  const double &dDt,double &dDtNew, double &dEffort,const double &c_dError0, const double &c_dRMin);//just calls rk4SoftNC in a 2 step way

#endif
