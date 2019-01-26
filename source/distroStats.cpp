//distroStats.cpp
//by Mike Williams 04-25-05

#include"distroStats.h"
#include"star6_0.h"
#include"constants.h"

void RVDRrms(star *stars[],double &dRrms,double &dVrms,double &dDRrms)
{
	dRrms=0;
	dVrms=0;
	dDRrms=0;
	
	for(int i=0;i<g_c_iStars;i++)
	{
		for(int j=i+1;j<g_c_iStars;j++)
		{
			dDRrms+=(stars[i]->m_vPosition-stars[j]->m_vPosition).SizeSquared();
		}
		dRrms+=stars[i]->m_vPosition.SizeSquared();
		dVrms+=stars[i]->m_vVelocity.SizeSquared();
	}
	
	dRrms=sqrt(dRrms/(double)g_c_iStars);
	dVrms=sqrt(dVrms/(double)g_c_iStars);
	dDRrms=sqrt(dDRrms/g_c_dPairs);
}
void RVDRrmsL(star *stars[],double &dRrms,double &dVrms,double &dDRrms,Vector3d &vL)
{
	dRrms=0;
	dVrms=0;
	dDRrms=0;
	vL.SetZero();
	
	for(int i=0;i<g_c_iStars;i++)
	{
		for(int j=i+1;j<g_c_iStars;j++)
		{
			dDRrms+=(stars[i]->m_vPosition-stars[j]->m_vPosition).SizeSquared();
		}
		dRrms+=stars[i]->m_vPosition.SizeSquared();
		dVrms+=stars[i]->m_vVelocity.SizeSquared();
		vL+=stars[i]->m_vPosition.Cross(stars[i]->m_vVelocity)*stars[i]->m_dMass;
	}
	
	dRrms=sqrt(dRrms/(double)g_c_iStars);
	dVrms=sqrt(dVrms/(double)g_c_iStars);
	dDRrms=sqrt(dDRrms/g_c_dPairs);
}
double AveMass(star *stars[])
{
	double dM;
	dM=0;
	for(int i=0;i<g_c_iStars;i++)
	{
		dM+=stars[i]->m_dMass;
	}
	return dM/g_c_iStars;
}
