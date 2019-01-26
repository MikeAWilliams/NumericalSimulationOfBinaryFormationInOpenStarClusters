// rand.cpp
// System Includes
#include <stdlib.h>
#include "rand.h"

// Constants
static float g_fRandMaxReciprocal=1.0f/(float) RAND_MAX;
static float g_dRandMaxReciprocal=1.0/(double) RAND_MAX;

// Random (0 <= i <= RAND_MAX)
int Random()
{
	return rand();
}

// Random (0 <= i <= nRange)
int Random(int nRange)
{
	return rand()%(nRange+1);
}

// Random Range (nStart <= i <= nEnd)
int RandomRange(int nStart,int nEnd)
{
	return ((nStart<nEnd)?(nStart+rand()%(nEnd-nStart+1)):(nEnd+rand()%(nStart-nEnd+1)));
}

// fRandom (0.0f <= f <= 1.0f)
float fRandom()
{
	return g_fRandMaxReciprocal*rand();
}

// fRandom (0.0f <= f <= fRange)
float fRandom(float fRange)
{
	return fRange*g_fRandMaxReciprocal*rand();
}

// fRandom (fStart <= f <= fEnd)
float fRandom(float fStart,float fEnd)
{
	return fStart+(fEnd-fStart)*g_fRandMaxReciprocal*rand();
}

double dRandom(double dRange)
{
	return dRange*g_dRandMaxReciprocal*rand();
}
