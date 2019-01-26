//compare2_1.h Mike Williams 3/05
#ifndef _compare_h
#define _compare_h

#include"star6_0.h"
#include"constants.h"


bool compare3Simple(star *star1, star *star2, Vector3d &vF,const double &rMin);
bool compare3Simple2(star *star1, star *star2, Vector3d &vF,double &r, const double &rMin);
bool compareSoftSimple(star *star1, star *star2, Vector3d &vF,double &r,const double &c_dRMin);
#endif
