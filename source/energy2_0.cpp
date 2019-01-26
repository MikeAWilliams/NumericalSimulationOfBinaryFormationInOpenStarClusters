//energy2_0.cpp Mike Williams 2/05
#include"star6_0.h"
#include"constants.h"
#include"Vector3d.h"
#include<math.h>
#include"energy1_0.h"

double energyNotVerlet(star *stars[],const int &n,double &PE, double &KE)
{
	PE=0;
	KE=0;
	for(int i = 0; i<n;i++)
	{
		for(int j =i+1; j<n; j++)
		{
			PE+=-stars[i]->m_dMass*stars[j]->m_dMass/(stars[i]->m_vPosition-stars[j]->m_vPosition).Size();
		}//end iner for
		KE+=.5*stars[i]->m_dMass*stars[i]->m_vVelocity.SizeSquared();
	}//end outer for
	return KE+PE;
}

double energyNotVerlet(star *stars[],double Rs[g_c_iStars][g_c_iStars],const int &n)
//computes the energy for a non verlet integrator (assumes the velocity is current)
{
	double dPE=0;
	double dKE=0;
	for(int i = 0; i<n;i++)
	{
		for(int j =i+1; j<n; j++)
		{
			dPE+=-stars[i]->m_dMass*stars[j]->m_dMass/sqrt(Rs[i][j]);
		}//end iner for
		dKE+=.5*stars[i]->m_dMass*stars[i]->m_vVelocity.SizeSquared();
	}//end outer for
	return dKE+dPE;
}

double energyNotVerlet(star *stars[],double Rs[g_c_iStars][g_c_iStars],const int &n,double &PE, double &KE)
//computes the energy for a non verlet integrator (assumes the velocity is current)
{
	PE=0;
	KE=0;
	for(int i = 0; i<n;i++)
	{
		for(int j =i+1; j<n; j++)
		{
			PE+=-stars[i]->m_dMass*stars[j]->m_dMass/sqrt(Rs[i][j]);
		}//end iner for
		KE+=.5*stars[i]->m_dMass*stars[i]->m_vVelocity.SizeSquared();
	}//end outer for
	return KE+PE;
}

double energyP(star *stars[])
{
	double dPE;
	Vector3d vR;
	dPE=0;
	
	for(int i = 0; i<g_c_iStars;i++)
	{
		for(int j =i+1; j<g_c_iStars; j++)
		{
			vR=stars[i]->m_vPosition - stars[j]->m_vPosition;
			dPE+=-g_c_dSG*stars[i]->m_dMass*stars[j]->m_dMass/vR.Size();
		}//end iner for
	}
	return dPE;
}
double energy(star *stars[])
{
	double dPE,dKE;
	Vector3d vR;
	dPE=0;
	dKE=0;
	
	for(int i = 0; i<g_c_iStars;i++)
	{
		for(int j =i+1; j<g_c_iStars; j++)
		{
			vR=stars[i]->m_vPosition - stars[j]->m_vPosition;
			dPE+=-g_c_dSG*stars[i]->m_dMass*stars[j]->m_dMass/vR.Size();
		}//end iner for
		dKE+=.5*stars[i]->m_dMass*stars[i]->m_vVelocity.SizeSquared();
	}
	return dPE+dKE;
}
double energy(star *stars[],double &dPE,double &dKE)
{
	Vector3d vR;
	dPE=0;
	dKE=0;
	
	for(int i = 0; i<g_c_iStars;i++)
	{
		for(int j =i+1; j<g_c_iStars; j++)
		{
			vR=stars[i]->m_vPosition - stars[j]->m_vPosition;
			dPE+=-g_c_dSG*stars[i]->m_dMass*stars[j]->m_dMass/vR.Size();
		}//end iner for
		dKE+=.5*stars[i]->m_dMass*stars[i]->m_vVelocity.SizeSquared();
	}
	return dPE+dKE;
}
