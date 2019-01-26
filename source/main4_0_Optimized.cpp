/*main4_0_RealDistroSwaper.cpp
Written by Mike Williams 05-16-05
Built to exchange real distrobutions with others
Step by steps
	Generate a distro and ouptut it in real units.
	Send it to your friends
	Load the same distro
	convert it to code units
	run simulation
	compare to friends.
Implements classic RK4 adaption
Takes arguments to determin distro size
*/
#include<iostream>
#include<fstream>
#include<math.h>
#include"star6_0.h"
#include"group1_0.h"
#include"constants.h"
#include"Vector3d.h"
#include"rand.h"
#include<sys/time.h>
#include"rk2_1.h"
#include"energy1_1.h"
#include"distrobution.h"
#include"distroIO1_2.h"
#include"util.h"
#include"binaryCount.h"
#include"outputManager1_1.h"
#include"adaptFunctions.h"
#include"distroStats.h"



//so stars is an aray of pointers to stars to be updated, n is the number of stars in the aray, t is the total time to go through, and dt is the time step


int main(int argc, char *argv[])
{
/*
*********************************************************************************************************
								Declare variables
*********************************************************************************************************
*/	
	
	char *cEnd;
	//units
	double dR0,dL,dT,dRMin,dError0;
	
	//pars input	
	dR0=strtod(argv[1],&cEnd);
	std::cout<<dR0<<std::endl;
	
	dL=dR0;
	dT=sqrt(dL*dL*dL/(g_c_dG*g_c_dM));
	dRMin=g_c_dRMinM/dL;
	dError0=g_c_dError0m/dL;
	

	//Sim constants
	//using namespace std;
	const int c_iNDataPtsS=1000;
	const int c_iNDataPtsF=1000;
	const double c_dChangePerOutputS=1/(double)c_iNDataPtsS;
	const double c_dChangePerOutputF=1/(double)c_iNDataPtsF;
	const int c_iNSteps=100000;
	const double c_dSimTime =g_c_dClusterTimes*sqrt(dR0*dR0*dR0/(g_c_dG*g_c_dM))/dT;
	double dLastTPS,dLastTPF;//last time percent screen, file
	double dDt = 1E-5;
	double dDtNew;
	
	char sAdaptDesc[13]="2StepRk4Soft";
	
	double dErrorP;

	
	outputManager allOutput;
	//const double c_dSimTime =2.418;//Total time for cecil example
	//loading related variables
	const bool c_bLoad = false;//to load the 
	const bool c_bSave = true;//to save or not to save

	//Sim variables
	double dTime = 0;
	double dDtMin=dDt;
	double dStopTime = c_dSimTime;
	star *stars[g_c_iStars],*starsp[g_c_iStars],*stars1[g_c_iStars];
	double dEffort=0;//a number representing the effort the code has done so far. So if you are in the main loop it is time steps*Total stars sqyared. however if you are
						//in the small group processing it is time steps*the number of stars in that small group squared

	//Output variables and temps

	//Energy related
	double dPE=0;//total potential engergy
	double dKE=0;//total kenetic energy
	double dEnergy;
	double dMinR;
	double dROne;

	//Noodle related

	//Binary related
	double dRs[g_c_iStars][g_c_iStars];
	int iNGroups,iNMembers;
	int iMembers[g_c_iMaxSystem];
//	double dR;// 0 to c_dRdist during init, in sim radious of an interaction
	int iNaries[g_c_iMaxSystem];
	int iNOrbits=0;
	int iBinary;
	bool bSplody=false;
	
	Vector3d vVCm,vCm;
	double dTotalMass;
	
	
	std::ofstream ruOut,constOut,massOut;
	std::ifstream saveIn;
	
	//Output controles
	int iCounter = 0;//used to controle number of outputs
	int c_iOutDevide = c_dSimTime/dDt/c_iNDataPtsS;//output will only happen every c_iOutDevide time

	c_iOutDevide = c_dSimTime/dDt/c_iNDataPtsS;//output will only happen every c_iOutDevide time
	dTime=0;
	iCounter=0;
	dEffort=0;

/*
*********************************************************************************************************
									Initialize variables
*********************************************************************************************************
*/

	
	
	dTime=0;
	allOutput.openOpt();

	for(int i = 0; i<g_c_iStars; i++)
	{
		stars[i]=new star;
		stars1[i]=new star;
		starsp[i]=new star;
	}

	if(c_bLoad)//get the initial conditions from file
	{
		//saveIn.open("./codeIo/input/clusterDistro.txt");
		saveIn.open("./codeIo/input/50VMass4-22.txt");
		//saveIn.open("./codeIo/input/50VariableMass.txt");
		
		//saveIn.open("./codeIo/input/4BODY.TXT");
		loadDistroRealUnits(stars,saveIn,dL,dT);
		saveIn.close();
	}
	else//generate the initial conditions yourself
	{
   		srand( (unsigned)time( NULL ) );
		//set up stars
		//uniformMassDistribute(stars);
		powerMassDistribute(stars,.5,10);
		spatialyDistributeMidR2(stars,dR0/dL,dROne);
		//energyDistribute(stars);
		double dJunk;
		//energyDistributeV2(stars,dJunk,dJunk);
		energyDistributeV2(stars);
		//store the distrobution for repetition
	}
	//Initialize output files
	ruOut.open("./codeIo/clusterDistro.txt");
	saveDistroRealUnits(stars,ruOut,dL,dT);
	ruOut.close();
	
	constOut.open("./codeIo/const.txt");
	saveConst(constOut,sAdaptDesc,dR0,dL);
	constOut.close();
	
	//allOutput.doCmOut(stars,dTime);
	massOut.open("./codeIo/averageMass.txt");
	massOut<<"Average Mass "<<AveMass(stars)<<std::endl;
	
	
	vVCm.SetZero();
	vCm.SetZero();
	dTotalMass=0;
	for(int i=0;i<g_c_iStars;i++)
	{
		dTotalMass+=stars[i]->m_dMass;
		vVCm+=stars[i]->m_vVelocity*stars[i]->m_dMass;
		vCm+=stars[i]->m_vPosition*stars[i]->m_dMass;
	}
	vVCm/=dTotalMass;
	vCm/=dTotalMass;
	
	for(int i=0;i<g_c_iStars;i++)
	{
		stars[i]->m_vPosition-=vCm;
		stars[i]->m_vVelocity-=vVCm;
	}
	
	allOutput.doDistroOut(stars,dTime,dDt);
	allOutput.doInitialOutOpt(stars,dTime);
	//allOutput.doCmOut(stars,dTime);
	//do 5 integration steps
	for(int k=0;k<5;k++)
	{
		for(int c=0; c<g_c_iStars;c++)
		{
			//reset the stars for the next step
			stars[c]->m_vAcceleration.SetZero();
			stars[c]->m_vTmpAcceleration.SetZero();
			stars[c]->m_vTmpAcceleration2.SetZero();
			stars[c]->m_vTmpAcceleration3.SetZero();
			stars[c]->m_bClose=false;
		}
		
		rk42StepSoft(stars,starsp,stars1,dRs,dDt,dDtNew,dEffort,dError0,dRMin);
		
		dTime+=dDt;
		dDt=dDtNew;
		dLastTPS=dTime/dStopTime;
		dLastTPF=dTime/dStopTime;
		
//		allOutput.doMcNeilOut(stars,dTime,7,47);
		//allOutput.doMcNeilOut(stars,dTime,0,1);
		allOutput.doFileOutOpt(stars,dTime,dEnergy);
		allOutput.doScreenOutOpt(stars,dEnergy,dTime,c_dSimTime,dDt);
		
		
	}
	for(int c=0; c<g_c_iStars;c++)
	{
		//reset the stars for the next step
		stars[c]->m_vAcceleration.SetZero();
		stars[c]->m_vTmpAcceleration.SetZero();
		stars[c]->m_vTmpAcceleration2.SetZero();
		stars[c]->m_vTmpAcceleration3.SetZero();
		stars[c]->m_bClose=false;
	}
		
/*
*********************************************************************************************************
									Main Loop
*********************************************************************************************************
*/
	dKE=dPE=0;
	for(;dTime<dStopTime;dKE=0,dPE=0)//main loop, controles when the sim stops
	//for(int iSteps=0;iSteps<c_iNSteps;iSteps++,dTime+=dDt,dKE=0,dPE=0)
	{
		//eliminate old forces
		for(int c=0; c<g_c_iStars;c++)
		{
			//reset the stars for the next step
			stars[c]->m_vAcceleration.SetZero();
			stars[c]->m_vTmpAcceleration.SetZero();
			stars[c]->m_vTmpAcceleration2.SetZero();
			stars[c]->m_vTmpAcceleration3.SetZero();
			stars[c]->m_bClose=false;
			if(stars[c]->m_vPosition.X !=stars[c]->m_vPosition.X)
			{
				std::cout<<"Warning star "<<c<<" is broken at time "<<dTime<<std::endl;
				for(int i=c+1;i<g_c_iStars;i++)
				{
					if(stars[i]->m_vPosition.X !=stars[c]->m_vPosition.X)
					{	
						std::cout<<"Warning star "<<i<<" is broken"<<std::endl;
					}
				}
				exit(1);
			}
		}
		
		rk42StepSoft(stars,starsp,stars1,dRs,dDt,dDtNew,dEffort,dError0,dRMin);
		dTime+=dDt;
		dDt=dDtNew;
		
		//controle output
		if((dTime/dStopTime-dLastTPF)>=c_dChangePerOutputF)
		{
			dLastTPF=dTime/dStopTime;
			//allOutput.doMcNeilOut(stars,dTime,22,40);
			//allOutput.doMcNeilOut(stars,dTime,7,47);
			//allOutput.doMcNeilOut(stars,dTime,0,1);
			allOutput.doFileOutOpt(stars,dTime,dEnergy);
		}//endif
		if((dTime/dStopTime-dLastTPS)>=c_dChangePerOutputS)
		{
			dLastTPS=dTime/dStopTime;
			//allOutput.doMcNeilOut(stars,dTime,0,1);
			//allOutput.doMainScrenOut1(stars,dRs,dTime,c_dSimTime,dDt);
			allOutput.doScreenOutOpt(stars,dEnergy,dTime,c_dSimTime,dDt);
		}//endif
		iCounter++;
	}//outermost while

	//do one last outputs
	//allOutput.doMainOut(stars,dTime,dRs);
	allOutput.doFileOutOpt(stars,dTime,dEnergy);
	
	
	//allOutput.doCmOut(stars,dTime);
	//allOutput.close();
	allOutput.closeOpt();
	
	allOutput.doFinalScreenOut();	
	for(int i = 0; i<g_c_iStars; i++)
	{
		delete stars[i];
		delete stars1[i];
		delete starsp[i];
	}
	
	system("xmms -p /mnt/ROUTER/public/reaserch/sounds/done.wav&");
	
	return 1;
}
