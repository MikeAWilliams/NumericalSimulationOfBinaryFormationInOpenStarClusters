//binaryCount.h
//written by Mike Williams 02-25-05

#ifndef binaryCount_h
#define binaryCount_h

#include"constants.h"
#include"star6_0.h"

int binaryCount(star *stars[], double Rs[g_c_iStars][g_c_iStars]);
int binaryCount(star *stars[]);
int binaryCountIndex(star *stars[],int iFirst[],int iSecond[], double dRBinarys[]);
#endif
