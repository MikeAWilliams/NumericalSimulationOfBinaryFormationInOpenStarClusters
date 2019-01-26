//adaptFunctions.cpp Mike Williams 4/05
#include"adaptFunctions.h"
#include"constants.h"
#include"star6_0.h"

void copyStars(star *starData[],star *starTarget[])
{
	for(int i=0;i<g_c_iStars;i++)
	{
		starTarget[i]->m_vPosition=starData[i]->m_vPosition;
		starTarget[i]->m_vVelocity=starData[i]->m_vVelocity;
		starTarget[i]->m_dMass=starData[i]->m_dMass;
	}
}
double errorPosition(star *stars1[],star *starsh[])
{
	double dAnswer;
	
	dAnswer=0;
	for(int i=0;i<g_c_iStars;i++)
	{
		dAnswer+=(starsh[i]->m_vPosition - stars1[i]->m_vPosition).Size();
	}
	return dAnswer;
}
