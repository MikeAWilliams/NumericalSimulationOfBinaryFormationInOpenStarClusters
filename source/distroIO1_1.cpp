//distroIO.cpp
//Mike Williams 02-05-05

#include"distroIO1_1.h"
#include"star6_0.h"
#include"constants.h"
#include<fstream>
#include<iostream>

void saveDistro(star *stars[],std::ofstream &disOut,double a_dTime,double a_dTimeStep)
{
	double dMt,dXt,dYt,dZt,dVxt,dVyt,dVzt;
	disOut<<g_c_iStars<<std::endl;
	//save star data
	for(int i = 0;i<g_c_iStars;i++)
	{
		dMt=stars[i]->m_dMass;
		dXt=stars[i]->m_vPosition.X;
		dYt=stars[i]->m_vPosition.Y;
		dZt=stars[i]->m_vPosition.Z;

		dVxt=stars[i]->m_vVelocity.X;
		dVyt=stars[i]->m_vVelocity.Y;
		dVzt=stars[i]->m_vVelocity.Z;
		disOut<<dMt<<" "<<dXt<<" "<<dYt<<" "<<dZt<<" "<<dVxt<<" "<<dVyt<<" "<<dVzt<<std::endl;
	}	
}

void loadDistro(star *stars[],std::ifstream &disIn)
{
	double dMt,dXt,dYt,dZt,dVxt,dVyt,dVzt;
	int iNStars;
	
	disIn>>iNStars;
	if(iNStars!=g_c_iStars)
	{
		std::cout<<"Error stars in file "<<iNStars<<" g_c_iStars "<<g_c_iStars<<std::endl;
		exit(1);
	}
	//load star data
	for(int i = 0;i<g_c_iStars;i++)
	{
		disIn>>dMt>>dXt>>dYt>>dZt>>dVxt>>dVyt>>dVzt;
		stars[i]->initialize(Vector3d(dXt,dYt,dZt),Vector3d(dVxt,dVyt,dVzt),dMt);
	}
}

void saveDistroRealUnits(star *stars[],std::ofstream &disOut)
{
	double dMt,dXt,dYt,dZt,dVxt,dVyt,dVzt;
	
	disOut<<g_c_iStars<<std::endl;
	//save star data
	for(int i = 0;i<g_c_iStars;i++)
	{
		dMt=stars[i]->m_dMass*g_c_dM;
		dXt=stars[i]->m_vPosition.X*g_c_dL;
		dYt=stars[i]->m_vPosition.Y*g_c_dL;
		dZt=stars[i]->m_vPosition.Z*g_c_dL;

		dVxt=stars[i]->m_vVelocity.X*g_c_dL/g_c_dT;
		dVyt=stars[i]->m_vVelocity.Y*g_c_dL/g_c_dT;
		dVzt=stars[i]->m_vVelocity.Z*g_c_dL/g_c_dT;
		disOut<<dMt<<" "<<dXt<<" "<<dYt<<" "<<dZt<<" "<<dVxt<<" "<<dVyt<<" "<<dVzt<<std::endl;
	}
}
void loadDistroRealUnits(star *stars[],std::ifstream &disIn)
{
	double dMt,dXt,dYt,dZt,dVxt,dVyt,dVzt;
	int iNStars;
	
	disIn>>iNStars;
	if(iNStars!=g_c_iStars)
	{
		std::cout<<"Error stars in file "<<iNStars<<" g_c_iStars "<<g_c_iStars<<std::endl;
		exit(1);
	}
	//load star data
	for(int i = 0;i<g_c_iStars;i++)
	{
		disIn>>dMt>>dXt>>dYt>>dZt>>dVxt>>dVyt>>dVzt;
		dMt/=g_c_dM;
		dXt/=g_c_dL;
		dYt/=g_c_dL;
		dZt/=g_c_dL;
		dVxt*=g_c_dT/g_c_dL;
		dVyt*=g_c_dT/g_c_dL;
		dVzt*=g_c_dT/g_c_dL;
		stars[i]->initialize(Vector3d(dXt,dYt,dZt),Vector3d(dVxt,dVyt,dVzt),dMt);
	}
}

void saveConst(std::ofstream &constOut,char sAdapt[])
{
	//save every possible parameter
	
	//simulation scale units, in standard units.
	constOut<<"Sim Constants"<<std::endl;
	constOut<<"g_c_dR0 =\t\t\t"<<g_c_dR0<<std::endl;

	constOut<<"g_c_dL =\t\t\t"<<g_c_dL<<std::endl;
	constOut<<"g_c_dM =\t\t\t"<<g_c_dM<<std::endl;
	constOut<<"g_c_iStars =\t\t\t"<<g_c_iStars<<std::endl;
	constOut<<"g_c_dDtMin =\t\t\t"<<g_c_dDtMin<<std::endl;
	constOut<<"g_c_dDtMax=\t\t\t"<<g_c_dDtMax<<std::endl;
	constOut<<"g_c_dClusterTimes =\t\t"<<g_c_dClusterTimes<<std::endl;
	constOut<<"g_c_dRSystem =\t\t\t"<<g_c_dRSystem<<std::endl;
	constOut<<"g_c_dRMinM =\t\t\t"<<g_c_dRMinM<<std::endl;
	constOut<<"g_c_iRMultipy =\t\t\t"<<g_c_iRMultipy<<std::endl;
	constOut<<"g_c_iSmallTimeSteps =\t\t"<<g_c_iSmallTimeSteps<<std::endl;
	constOut<<"g_c_iDtDevider =\t\t"<<g_c_iDtDevider<<std::endl;
	constOut<<"g_c_dDtSmall =\t\t\t"<<g_c_dDtSmall<<std::endl;
	constOut<<"g_c_iNTimeStepsToBeBinary =\t"<<g_c_iNTimeStepsToBeBinary<<std::endl;
	constOut<<"g_c_dTSystem =\t\t\t"<<g_c_dTSystem<<std::endl;
	constOut<<"g_c_dAdaptPower =\t\t"<<g_c_dAdaptPower<<std::endl;
	constOut<<"g_c_dError0m =\t\t\t"<<g_c_dError0m<<std::endl;
	
	constOut<<"Adapt description\t\t"<<sAdapt<<std::endl;	
}
