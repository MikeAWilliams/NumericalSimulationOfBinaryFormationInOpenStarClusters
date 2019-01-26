//distrobution.h Mike Williams 2/05
#include"distrobution.h"
#include"constants.h"
#include"star6_0.h"
#include"energy1_0.h"
#include"rand.h"
#include<math.h>
#include<iostream>
#include<fstream>
#include"util.h"

//spatial functions
double ro(double a_dR,double a_dR0)
{
	return g_c_iStars/a_dR0*exp(-a_dR/a_dR0);
} 

double rMid(int a_iI,double a_dR0)
{
	return (a_iI*rOne(a_dR0)/g_c_iNBins - (a_iI-1)*rOne(a_dR0)/g_c_iNBins)/2 + (a_iI-1)*rOne(a_dR0)/g_c_iNBins;
}

int nStarsInBinR(double a_dR, double a_dR0)
{
	return round(rOne(a_dR0)/g_c_iNBins*ro(a_dR,a_dR0));
}

double rOne(double a_dR0)
{
	return a_dR0*-log(1.0/g_c_iStars);
}

void uniformMassDistribute(star *stars[])
{
	for(int i=0;i<g_c_iStars;i++)
	{
		stars[i]->m_dMass=g_c_dMSun/g_c_dM;
	}
}

double boltzmann(double a_dE,double a_dT)
{
	return g_c_iStars/g_c_dSK/a_dT*exp(-a_dE/g_c_dSK/a_dT);
}

double boltzmannV2(double a_dE,double a_dT)
{
	return g_c_iStars/(2*g_c_dSK*g_c_dSK*g_c_dSK*a_dT*a_dT*a_dT)*a_dE*a_dE*exp(-a_dE/g_c_dSK/a_dT);
}

double eMid(int a_iI,double a_dT)
{
	return (a_iI*eOne(a_dT)/g_c_iNBins - (a_iI-1)*eOne(a_dT)/g_c_iNBins)/2 + (a_iI-1)*eOne(a_dT)/g_c_iNBins;
}

int nStarsInBinE(double a_dE, double a_dT)
{
	return round(eOne(a_dT)/g_c_iNBins*boltzmann(a_dE,a_dT));
}
int nStarsInBinEV2(double a_dE,double a_dT,double a_dE1)
{
	double dDE=a_dE1/g_c_iNBins;
	
	return round(dDE*boltzmannV2(a_dE,a_dT));
}

double eOne(double a_dT)
{
	return -g_c_dSK*a_dT*log(1.0/g_c_iStars);
}

double eMid(int a_iI,double a_dT,double a_dE1)
{
	double dDE=a_dE1/g_c_iNBins;
	return dDE/2 + (a_iI-1)*dDE;
}

double tailV2(double a_dE,double a_dT)
{
	double dK2=g_c_dSK*g_c_dSK;
	double dT2=a_dT*a_dT;
	return exp(-1*a_dE/(g_c_dSK*a_dT))*g_c_iStars*(a_dE*a_dE+2*a_dE*a_dT*g_c_dSK+2*dK2*dT2)/(2*dK2*dT2);
}

double tailV2p(double a_dE, double a_dT)
{
	double dK2=g_c_dSK*g_c_dSK;
	double dT2=a_dT*a_dT;
	
	double dTerm1,dTerm2;
	
	dTerm1=exp(-1*a_dE/(g_c_dSK*a_dT))*g_c_iStars*(2*a_dE+2*g_c_dSK*a_dT)/(2*dK2*dT2);
	
	dTerm2=exp(-1*a_dE/(g_c_dSK*a_dT))*g_c_iStars*(a_dE*a_dE+2*a_dE*a_dT*g_c_dSK+2*dK2*dT2)/(2*dK2*dT2*g_c_dSK*a_dT);
	
	return dTerm1-dTerm2;
	
}

double eOneV2(double a_dT,double a_dE0)
{
	const double c_dErrorTollerence=1E-10;//error convergence tollerance
	const int c_iNIterations = 100;//maximum iterations to try befroe giving up
	double dError=1;
	double dE,dNE;
	
	dE=a_dE0;
	
	for(int i=0;(i<c_iNIterations&&dError>c_dErrorTollerence);i++)
	{
			
		dNE=tailV2(dE,a_dT);
		dError = abs(dNE-1);
		dE=dE - (dNE-1)/tailV2p(dE,a_dT);
	}
	return dE;
}
void spatialyDistributeMid(star *stars[],const double a_dR0)
{
	int iNInit=0;
	int iStarsInBin,iIndex;
	double dR;
	bool bInit[g_c_iStars];
	//std::ofstream spOut;
	
//	spOut.open("spaceOut.csv");
	
	//spOut<<"i,x,y,z"<<std::endl;
	
	for(int i=0;i<g_c_iStars;i++)
	{
		bInit[i]=false;
	}
	
	for(int i=1;i<=g_c_iNBins;i++)
	{
		dR=rMid(i,a_dR0);
		iStarsInBin=nStarsInBinR(dR,a_dR0);
		for(int j=0;j<iStarsInBin;j++)
		{
			iIndex=Random(g_c_iStars-1);
			while(bInit[iIndex])
			{
				iIndex=Random(g_c_iStars-1);
			}
			bInit[iIndex]=true;
			stars[iIndex]->m_vPosition=Vector3d(fRandom(-1,1),fRandom(-1,1),fRandom(-1,1)).Normal()*dR;
			iNInit++;
			//spOut<<iIndex<<","<<stars[iIndex]->m_vPosition.X<<","<<stars[iIndex]->m_vPosition.Y<<","<<stars[iIndex]->m_vPosition.Z<<std::endl;
		}
	}
	if(iNInit<g_c_iStars)
	{
		//put one at rOne
		//hopefully there arn't to many so rather than randomize i will be systematic
		//so the first uninitiated star I find I initalize to be at rOne
		//spOut<<",,,,ROne"<<std::endl;
		for(int i=0;i<g_c_iStars;i++)
		{
			if(!bInit[i])
			{
				stars[i]->m_vPosition=Vector3d(fRandom(-1,1),fRandom(-1,1),fRandom(-1,1)).Normal()*rOne(a_dR0);
				bInit[i]=true;
				iNInit++;
				//spOut<<i<<","<<stars[i]->m_vPosition.X<<","<<stars[i]->m_vPosition.Y<<","<<stars[i]->m_vPosition.Z<<std::endl;
				//kill the for
				break;
			}
		}
		//it could happen that there are sill a couple who didn't get placed...
		//so place any unpaced stars completely randomly
		//spOut<<",,,,Extras"<<std::endl;
		for(int i=0;i<g_c_iStars;i++)
		{
			if(!bInit[i])
			{
				stars[i]->m_vPosition=Vector3d(fRandom(-1,1),fRandom(-1,1),fRandom(-1,1)).Normal()*fRandom(rOne(a_dR0));
				bInit[i]=true;
				iNInit++;
				//spOut<<i<<","<<stars[i]->m_vPosition.X<<","<<stars[i]->m_vPosition.Y<<","<<stars[i]->m_vPosition.Z<<std::endl;
			}
		}
	}
//	spOut.close();
}
void energyDistribute(star *stars[],double &a_dT)
{
	int iNInit=0;
	int iStarsInBin,iIndex;
	double dE,dPE,dKE,dV;
	bool bInit[g_c_iStars];
	
	dPE=energyP(stars);
	dKE=-dPE/2.0;
	a_dT=dKE/g_c_iStars/g_c_dSK;
	
	for(int i=0;i<g_c_iStars;i++)
	{
		bInit[i]=false;
	}
	
	for(int i=1;i<=g_c_iNBins;i++)
	{
		dE=eMid(i,a_dT);
		iStarsInBin=nStarsInBinE(dE,a_dT);
		for(int j=0;j<iStarsInBin;j++)
		{
			iIndex=Random(g_c_iStars-1);
			while(bInit[iIndex])
			{
				iIndex=Random(g_c_iStars-1);
			}
			bInit[iIndex]=true;
			dV=sqrt(2*dE/stars[iIndex]->m_dMass);
			stars[iIndex]->m_vVelocity=Vector3d(fRandom(-1,1),fRandom(-1,1),fRandom(-1,1)).Normal()*dV;
			iNInit++;
		}
	}
	if(iNInit<g_c_iStars)
	{
		//put one at rOne
		//hopefully there arn't to many so rather than randomize i will be systematic
		//so the first uninitiated star I find I initalize to be at rOne
		for(int i=0;i<g_c_iStars;i++)
		{
			if(!bInit[i])
			{
				dE=eOne(a_dT);
				dV=sqrt(2*dE/stars[i]->m_dMass);
				stars[i]->m_vVelocity=Vector3d(fRandom(-1,1),fRandom(-1,1),fRandom(-1,1)).Normal()*dV;
				bInit[i]=true;
				iNInit++;
				//kill the for
				break;
			}
		}
		//it could happen that there are sill a couple who didn't get placed...
		//so place any unpaced stars completely randomly
		for(int i=0;i<g_c_iStars;i++)
		{
			if(!bInit[i])
			{
				dE=fRandom(eOne(a_dT));
				dV=sqrt(2*dE/stars[i]->m_dMass);
				stars[i]->m_vVelocity=Vector3d(fRandom(-1,1),fRandom(-1,1),fRandom(-1,1)).Normal()*dV;
				bInit[i]=true;
				iNInit++;
			}
		}
	}
}

void energyDistributeV2(star *stars[],double &a_dT,double &a_dE1)
{
	int iNInit=0;
	int iStarsInBin,iIndex;
	double dE,dPE,dKE,dV,dE1;
	bool bInit[g_c_iStars];
	
	dPE=energyP(stars);
	dKE=-dPE/2.0;
	a_dT=dKE/(g_c_iStars*g_c_dSK*3);
	
	dE1=eOneV2(a_dT,dKE/g_c_iStars);
	a_dE1=dE1;
	
	for(int i=0;i<g_c_iStars;i++)
	{
		bInit[i]=false;
	}
	
	for(int i=1;i<=g_c_iNBins;i++)
	{
		dE=eMid(i,a_dT,dE1);
		iStarsInBin=nStarsInBinEV2(dE,a_dT,dE1);
		for(int j=0;j<iStarsInBin;j++)
		{
			iIndex=Random(g_c_iStars-1);
			while(bInit[iIndex])
			{
				iIndex=Random(g_c_iStars-1);
			}
			bInit[iIndex]=true;
			dV=sqrt(2*dE/stars[iIndex]->m_dMass);
			stars[iIndex]->m_vVelocity=Vector3d(fRandom(-1,1),fRandom(-1,1),fRandom(-1,1)).Normal()*dV;
			iNInit++;
		}
	}
	if(iNInit<g_c_iStars)
	{
		//put one at eOne
		//hopefully there arn't to many so rather than randomize i will be systematic
		//so the first uninitiated star I find I initalize to be at eOne
		for(int i=0;i<g_c_iStars;i++)
		{
			if(!bInit[i])
			{
				dE=dE1;
				dV=sqrt(2*dE/stars[i]->m_dMass);
				stars[i]->m_vVelocity=Vector3d(fRandom(-1,1),fRandom(-1,1),fRandom(-1,1)).Normal()*dV;
				bInit[i]=true;
				iNInit++;
				//kill the for
				break;
			}
		}
		//it could happen that there are sill a couple who didn't get placed...
		//so place any unpaced stars completely randomly
		for(int i=0;i<g_c_iStars;i++)
		{
			if(!bInit[i])
			{
				dE=fRandom(dE1);
				dV=sqrt(2*dE/stars[i]->m_dMass);
				stars[i]->m_vVelocity=Vector3d(fRandom(-1,1),fRandom(-1,1),fRandom(-1,1)).Normal()*dV;
				bInit[i]=true;
				iNInit++;
			}
		}
	}
	
}

void energyDistributeV2(star *stars[])
{
	int iNInit=0;
	int iStarsInBin,iIndex;
	double dE,dPE,dKE,dV,dE1,dT;
	bool bInit[g_c_iStars];
	
	dPE=energyP(stars);
	dKE=-dPE/2.0;
	dT=dKE/(g_c_iStars*g_c_dSK*3);
	
	dE1=eOneV2(dT,dKE/g_c_iStars);
	
	for(int i=0;i<g_c_iStars;i++)
	{
		bInit[i]=false;
	}
	
	for(int i=1;i<=g_c_iNBins;i++)
	{
		dE=eMid(i,dT,dE1);
		iStarsInBin=nStarsInBinEV2(dE,dT,dE1);
		for(int j=0;j<iStarsInBin;j++)
		{
			iIndex=Random(g_c_iStars-1);
			while(bInit[iIndex])
			{
				iIndex=Random(g_c_iStars-1);
			}
			bInit[iIndex]=true;
			dV=sqrt(2*dE/stars[iIndex]->m_dMass);
			stars[iIndex]->m_vVelocity=Vector3d(fRandom(-1,1),fRandom(-1,1),fRandom(-1,1)).Normal()*dV;
			iNInit++;
		}
	}
	if(iNInit<g_c_iStars)
	{
		//put one at rOne
		//hopefully there arn't to many so rather than randomize i will be systematic
		//so the first uninitiated star I find I initalize to be at rOne
		for(int i=0;i<g_c_iStars;i++)
		{
			if(!bInit[i])
			{
				dE=dE1;
				dV=sqrt(2*dE/stars[i]->m_dMass);
				stars[i]->m_vVelocity=Vector3d(fRandom(-1,1),fRandom(-1,1),fRandom(-1,1)).Normal()*dV;
				bInit[i]=true;
				iNInit++;
				//kill the for
				break;
			}
		}
		//it could happen that there are sill a couple who didn't get placed...
		//so place any unpaced stars completely randomly
		for(int i=0;i<g_c_iStars;i++)
		{
			if(!bInit[i])
			{
				dE=fRandom(dE1);
				dV=sqrt(2*dE/stars[i]->m_dMass);
				stars[i]->m_vVelocity=Vector3d(fRandom(-1,1),fRandom(-1,1),fRandom(-1,1)).Normal()*dV;
				bInit[i]=true;
				iNInit++;
			}
		}
	}
	
}

void energyDistribute(star *stars[])
{
	int iNInit=0;
	int iStarsInBin,iIndex;
	double dE,dPE,dKE,dV,dT;
	bool bInit[g_c_iStars];
	//std::ofstream spOut;
	
//	spOut.open("velocityOut.csv");
	
	//spOut<<"i,Vx,Vy,Vz"<<std::endl;
	
	dPE=energyP(stars);
	dKE=-dPE/2.0;
	dT=dKE/g_c_iStars/g_c_dSK;
	
	for(int i=0;i<g_c_iStars;i++)
	{
		bInit[i]=false;
	}
	
	for(int i=1;i<=g_c_iNBins;i++)
	{
		dE=eMid(i,dT);
		iStarsInBin=nStarsInBinE(dE,dT);
		for(int j=0;j<iStarsInBin;j++)
		{
			iIndex=Random(g_c_iStars-1);
			while(bInit[iIndex])
			{
				iIndex=Random(g_c_iStars-1);
			}
			bInit[iIndex]=true;
			dV=sqrt(2*dE/stars[iIndex]->m_dMass);
			stars[iIndex]->m_vVelocity=Vector3d(fRandom(-1,1),fRandom(-1,1),fRandom(-1,1)).Normal()*dV;
			//spOut<<iIndex<<","<<stars[iIndex]->m_vVelocity.X<<","<<stars[iIndex]->m_vVelocity.Y<<","<<stars[iIndex]->m_vVelocity.Z<<std::endl;
			iNInit++;
		}
	}
	if(iNInit<g_c_iStars)
	{
		//spOut<<",,,,eOne"<<std::endl;
		//put one at rOne
		//hopefully there arn't to many so rather than randomize i will be systematic
		//so the first uninitiated star I find I initalize to be at eOne
		for(int i=0;i<g_c_iStars;i++)
		{
			if(!bInit[i])
			{
				dE=eOne(dT);
				dV=sqrt(2*dE/stars[i]->m_dMass);
				stars[i]->m_vVelocity=Vector3d(fRandom(-1,1),fRandom(-1,1),fRandom(-1,1)).Normal()*dV;
				bInit[i]=true;
				iNInit++;
				//spOut<<i<<","<<stars[i]->m_vVelocity.X<<","<<stars[i]->m_vVelocity.Y<<","<<stars[i]->m_vVelocity.Z<<std::endl;
				//kill the for
				break;
			}
		}
		//spOut<<",,,,Extras"<<std::endl;
		//it could happen that there are sill a couple who didn't get placed...
		//so place any unpaced stars completely randomly
		for(int i=0;i<g_c_iStars;i++)
		{
			if(!bInit[i])
			{
				dE=dRandom(eOne(dT));
				dV=sqrt(2*dE/stars[i]->m_dMass);
				stars[i]->m_vVelocity=Vector3d(fRandom(-1,1),fRandom(-1,1),fRandom(-1,1)).Normal()*dV;
				bInit[i]=true;
				iNInit++;
				//spOut<<i<<","<<stars[i]->m_vVelocity.X<<","<<stars[i]->m_vVelocity.Y<<","<<stars[i]->m_vVelocity.Z<<std::endl;
			}
		}
	}
	
}

//mass realted functions

double massPower(double a_dM,double a_dMMin, double a_dMMax)
{
	return 0.35*pow(a_dM,-1.35)*g_c_iStars/(1.0/pow(a_dMMin,0.35)-1.0/pow(a_dMMax,0.35));
}
double mMid(int a_iI,double a_dMMin, double a_dMMax)
{
	double dDM=(a_dMMax-a_dMMin)/g_c_iNBins;
	return a_dMMin + dDM/2.0 + (a_iI-1)*dDM;
}
int nStarsInBinM(double a_dM,double a_dMMin, double a_dMMax)
{
	double dDM=(a_dMMax-a_dMMin)/g_c_iNBins;
	//return round(dDM*massPower(a_dM,a_dMMin,a_dMMax));
	double dA=a_dM-dDM/2.0;//integrating between a and b
	double dB=a_dM+dDM/2.0;
	
	double dAnswer;
	
	dAnswer=round((1.0/pow(dA,0.35)-1.0/pow(dB,0.35))*g_c_iStars/(1.0/pow(a_dMMin,0.35)-1.0/pow(a_dMMax,0.35)));
	
	return dAnswer;
	
}

void powerMassDistribute(star *stars[],double a_dMMin, double a_dMMax)
{
	int iNInit=0;
	int iStarsInBin,iIndex;
	double dM;
	bool bInit[g_c_iStars];
	double dDM=(a_dMMax-a_dMMin)/g_c_iNBins;
	int iTotalStars=0;
		
	//test massPower
	for(int i=1;i<=(g_c_iNBins);i++)
	{
		dM=mMid(i,a_dMMin,a_dMMax);
		iTotalStars+=nStarsInBinM(dM,a_dMMin,a_dMMax);
	}
	if(iTotalStars>g_c_iStars)
	{
		std::cout<<"Error to course a mass distrobution"<<std::endl;
		std::cout<<"Total stars initialzed "<<iTotalStars<<" Total stars "<<g_c_iStars<<std::endl;
	}
	
	for(int i=0;i<g_c_iStars;i++)
	{
		bInit[i]=false;
	}
	
	for(int i=1;i<=(g_c_iNBins);i++)
	{
		dM=mMid(i,a_dMMin,a_dMMax);
		iStarsInBin=nStarsInBinM(dM,a_dMMin,a_dMMax);
		for(int j=0;j<iStarsInBin;j++)
		{
			iIndex=Random(g_c_iStars-1);
			while(bInit[iIndex])
			{
				iIndex=Random(g_c_iStars-1);
			}
			bInit[iIndex]=true;
			stars[iIndex]->m_dMass=dM;
			iNInit++;
			//spOut<<iIndex<<","<<stars[iIndex]->m_vPosition.X<<","<<stars[iIndex]->m_vPosition.Y<<","<<stars[iIndex]->m_vPosition.Z<<std::endl;
		}
	}
	if(iNInit<g_c_iStars)
	{
		//it could happen that there are sill a couple who didn't get placed...
		//so place any unpaced stars completely randomly
		//spOut<<",,,,Extras"<<std::endl;
		for(int i=0;i<g_c_iStars;i++)
		{
			if(!bInit[i])
			{
				stars[i]->m_dMass=fRandom(a_dMMin,a_dMMax);
				bInit[i]=true;
				iNInit++;
			}
		}
	}
}

//R2 spatial functions
double roR2(double a_dR,double a_dR0)
{
	return g_c_iStars*pow(a_dR,2)/(2*pow(a_dR0,3))*exp(-1*a_dR/a_dR0);
}

double tailR2p(double a_dR,double a_dR0)
{
	double dTerm1,dTerm2;
	
	dTerm1=exp(-1*a_dR/a_dR0)*g_c_iStars*2*(a_dR0+a_dR)/(2*pow(a_dR0,2));
	dTerm2=exp(-1*a_dR/a_dR0)*g_c_iStars*(2*pow(a_dR0,2)+2*a_dR0*a_dR+pow(a_dR,2))/(2*pow(a_dR0,3));
	
	return dTerm1 - dTerm2;
}

double tailR2(double a_dR,double a_dR0)
{
	return exp(-1*a_dR/a_dR0)*g_c_iStars*(2*pow(a_dR0,2) + 2*a_dR0*a_dR + pow(a_dR,2))/(2*pow(a_dR0,2));
}

double rMidR2(int a_iI,double a_dR0,double a_dROne)
{
	double dDR=a_dROne/g_c_iNBins;
	
	return  dDR/2 + (a_iI-1)*dDR;
}

int nStarsInBinRR2(double a_dR,double a_dR0,double a_dROne)
{
	double dDR=a_dROne/g_c_iNBins;
	return round(roR2(a_dR,a_dR0)*dDR);
}

double rOneR2(double a_dR0,double a_dRInit)
{
	const double c_dErrorTollerence=1E-10;//error convergence tollerance
	const int c_iNIterations = 100;//maximum iterations to try befroe giving up
	double dError=1;
	double dR,dNR;
	
	dR=a_dRInit;
	
	for(int i=0;(i<c_iNIterations&&dError>c_dErrorTollerence);i++)
	{
			
		dNR=tailR2(dR,a_dR0);
		dError = abs(dNR-1);
		dR=dR - (dNR-1)/tailR2p(dR,a_dR0);
	}
	return dR;
}

void spatialyDistributeMidR2(star *stars[],const double a_dR0, double &a_dROne)
{
	int iNInit=0;
	int iStarsInBin,iIndex;
	double dR,dR1;
	bool bInit[g_c_iStars];
	
	dR1=rOneR2(a_dR0,a_dR0);
	a_dROne=dR1;
	
	for(int i=0;i<g_c_iStars;i++)
	{
		bInit[i]=false;
	}
	
	for(int i=1;i<=g_c_iNBins;i++)
	{
		dR=rMidR2(i,a_dR0,dR1);
		iStarsInBin=nStarsInBinRR2(dR,a_dR0,dR1);
		for(int j=0;j<iStarsInBin;j++)
		{
			iIndex=Random(g_c_iStars-1);
			while(bInit[iIndex])
			{
				iIndex=Random(g_c_iStars-1);
			}
			bInit[iIndex]=true;
			stars[iIndex]->m_vPosition=Vector3d(fRandom(-1,1),fRandom(-1,1),fRandom(-1,1)).Normal()*dR;
			iNInit++;
		}
	}
	if(iNInit<g_c_iStars)
	{
		//put one at eOne
		//hopefully there arn't to many so rather than randomize i will be systematic
		//so the first uninitiated star I find I initalize to be at eOne
		for(int i=0;i<g_c_iStars;i++)
		{
			if(!bInit[i])
			{
				dR=dR1;
				stars[i]->m_vPosition=Vector3d(fRandom(-1,1),fRandom(-1,1),fRandom(-1,1)).Normal()*dR;
				bInit[i]=true;
				iNInit++;
				//kill the for
				break;
			}
		}
		//it could happen that there are sill a couple who didn't get placed...
		//so place any unpaced stars completely randomly
		for(int i=0;i<g_c_iStars;i++)
		{
			if(!bInit[i])
			{
				dR=fRandom(dR1);
				stars[i]->m_vPosition=Vector3d(fRandom(-1,1),fRandom(-1,1),fRandom(-1,1)).Normal()*dR;
				bInit[i]=true;
				iNInit++;
			}
		}
	}
}
