//rk2_1,cpp Mike Williams 2-05
#include"star6_0.h"
#include"constants.h"
#include"compare2_1.h"
#include"rk2_1.h"
#include<iostream>
#include<fstream>
#include"util.h"
#include"adaptFunctions.h"



void rk42StepSoft(star *stars[],star *starsp[],star *stars1[], double dRs[g_c_iStars][g_c_iStars],  const double &dDt,double &dDtNew,double &dEffort,const double &c_dError0, const double &c_dRMin)//just calls rk4SoftNC in a 2 step way
{
	double dErrorP;
	copyStars(stars,starsp);
	
	rk4SoftNC(stars,dRs,g_c_iStars,dDt,dEffort,c_dRMin);
	copyStars(stars,stars1);//stor where we are after 1 whole step
	copyStars(starsp,stars);//reset to initial conditions
		
	for(int c=0; c<g_c_iStars;c++)
	{
		//reset the stars for the next step
		stars[c]->m_vAcceleration.SetZero();
		stars[c]->m_vTmpAcceleration.SetZero();
		stars[c]->m_vTmpAcceleration2.SetZero();
		stars[c]->m_vTmpAcceleration3.SetZero();
		stars[c]->m_bClose=false;
	}
		
	rk4SoftNC(stars,dRs,g_c_iStars,dDt/2.0,dEffort,c_dRMin);
	for(int c=0; c<g_c_iStars;c++)
	{
		//reset the stars for the next step
		stars[c]->m_vAcceleration.SetZero();
		stars[c]->m_vTmpAcceleration.SetZero();
		stars[c]->m_vTmpAcceleration2.SetZero();
		stars[c]->m_vTmpAcceleration3.SetZero();
		stars[c]->m_bClose=false;
	}
		
	rk4SoftNC(stars,dRs,g_c_iStars,dDt/2.0,dEffort,c_dRMin);
		
	dErrorP=errorPosition(stars1,stars);
	//copyStars(stars1,stars);//really we want to use the time step suggested last time
	dDtNew=dDt/pow(dErrorP/c_dError0,g_c_dAdaptPower);//set the dt for next time step		
	if(dDtNew<g_c_dDtMin)
	{
		dDtNew=g_c_dDtMin;
	}
	else if(dDtNew>g_c_dDtMax)
	{
		dDtNew=g_c_dDtMax;
	}
}


void rk4SoftNC(star *stars[], double Rs[g_c_iStars][g_c_iStars], const int &n,  const double &dt, double &Effort,const double &c_dRMin)
{	
	double dR;
	Vector3d vF;

	//compute first r
	for(int i = 0; i<n;i++)
	{
		for(int j =i+1; j<n; j++)
		{
			compareSoftSimple(stars[j],stars[i], vF,dR,c_dRMin);
			stars[j]->m_vAcceleration+=vF/stars[j]->m_dMass;
			stars[i]->m_vAcceleration-=vF/stars[i]->m_dMass;
				

			Rs[i][j]=dR;
			Effort++;

			//system counting Calculations
			//trackSystems.setTime(i,j,dR,dt);

		}//end inner for
/*********************************************
			Perform actual integration
********************************************/
		stars[i]->m_vPositionT=stars[i]->m_vPosition;
		
		stars[i]->m_vPosition1=stars[i]->m_vVelocity*dt;
		stars[i]->m_vVelocity1=stars[i]->m_vAcceleration*dt;
		
		stars[i]->m_vPosition2 = (stars[i]->m_vVelocity + stars[i]->m_vVelocity1*0.5)*dt;
		stars[i]->m_vPosition=stars[i]->m_vPositionT + stars[i]->m_vPosition1*0.5;
	}//end outer for
	
	for(int i = 0; i<n;i++)
	{
		for(int j =i+1; j<n; j++)
		{
			compareSoftSimple(stars[j],stars[i], vF,dR,c_dRMin);
			//PE+=dPEtemp;
			stars[j]->m_vTmpAcceleration+=vF/stars[j]->m_dMass;//use TmpAcceleration to save processor time of clearing acceleration
			stars[i]->m_vTmpAcceleration-=vF/stars[i]->m_dMass;
			
			Effort++;
		}//end inner for
/*********************************************
		Perform actual integration
********************************************/
		stars[i]->m_vVelocity2 = stars[i]->m_vTmpAcceleration*dt;
		stars[i]->m_vPosition3 = (stars[i]->m_vVelocity + stars[i]->m_vVelocity2*0.5)*dt;
		
		stars[i]->m_vPosition = stars[i]->m_vPositionT + stars[i]->m_vPosition2*0.5;
	}//end outer for
	
	for(int i = 0; i<n;i++)
	{
		for(int j =i+1; j<n; j++)
		{
			compareSoftSimple(stars[j],stars[i], vF,dR,c_dRMin);
			//PE+=dPEtemp;
			stars[j]->m_vTmpAcceleration2+=vF/stars[j]->m_dMass;
			stars[i]->m_vTmpAcceleration2-=vF/stars[i]->m_dMass;
				
			Effort++;

		}//end inner for
/*********************************************
		Perform actual integration
********************************************/
		stars[i]->m_vVelocity3=stars[i]->m_vTmpAcceleration2*dt;
		stars[i]->m_vPosition4=(stars[i]->m_vVelocity + stars[i]->m_vVelocity3)*dt;
		
		stars[i]->m_vPosition = stars[i]->m_vPositionT + stars[i]->m_vPosition3;
		
	}//end outer for
	
		for(int i = 0; i<n;i++)
		{
			for(int j =i+1; j<n; j++)
			{
				compareSoftSimple(stars[j],stars[i], vF,dR,c_dRMin);
				//PE+=dPEtemp;
				stars[j]->m_vTmpAcceleration3+=vF/stars[j]->m_dMass;
				stars[i]->m_vTmpAcceleration3-=vF/stars[i]->m_dMass;
				
				Effort++;
			}//end inner for
/*********************************************
		Perform actual integration
********************************************/
		stars[i]->m_vVelocity4=stars[i]->m_vTmpAcceleration3*dt;
		
		stars[i]->m_vPosition=stars[i]->m_vPositionT + (stars[i]->m_vPosition1 + (stars[i]->m_vPosition2 + stars[i]->m_vPosition3)*2.0 + stars[i]->m_vPosition4)/6.0;
		stars[i]->m_vVelocity=stars[i]->m_vVelocity + (stars[i]->m_vVelocity1 + (stars[i]->m_vVelocity2 + stars[i]->m_vVelocity3)*2.0 + stars[i]->m_vVelocity4)/6.0;
		
	}//end outer for
	
}
