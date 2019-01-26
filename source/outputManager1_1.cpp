/*outputManager1_1.cpp
Written by Mike Williams 04-13-05
*/

#include"outputManager1_1.h"
#include<fstream>
#include<iostream>
#include"group1_0.h"
#include"Vector3d.h"
#include"distroIO.h"
#include"binaryCount.h"
#include"energy1_1.h"
#include"util.h"
#include"distroStats.h"
#include"binarytracker.h"
#include<sys/time.h>

void outputManager::open()
{
	system("rm ./codeIo/last/*.csv");
	system("rm ./codeIo/last/*.txt");
	system("cp ./codeIo/*.* ./codeIo/last/");
	system("rm ./codeIo/*.csv");
	system("rm ./codeIo/*.txt");
	
	mcOut.open("./codeIo/mcneil.csv");
	
	eOut.open("./codeIo/energy.csv");

	noodleOut.open("./codeIo/noodle.csv");
	
	rOut.open("./codeIo/rs.csv");
	
	binaryOut.open("./codeIo/binary.csv");
	
	cmOut.open("./codeIo/cm.csv");
	
	rvOut.open("./codeIo/rvMike.txt");
	
	abOut.open("./codeIo/averageBinary.csv");
	
	dtOut.open("./codeIo/dt.csv");
	
	bIndexOut.open("./codeIo/binaryIndex.csv");
	
	
	//initialize files
	mcOut<<"T,DE,Nb,xi,yi,zi,xj,yj,zj,rij,cm,vcm,Rrms,Vrms,DRrms,Lx,Ly,Lz"<<std::endl;
	bIndexOut<<"Time,Binaries,i0,j0,r0,i1,j1,r1"<<std::endl;
	dtOut<<"Time,Dt"<<std::endl;
	noodleOut<<"Time,";
	rOut<<"Time,";

	binaryOut<<"Time,Binaries"<<std::endl;
	rvOut<<"T ";
	for(int i = 0; i<g_c_iStars;i++)
	{
		for(int j =i+1; j<g_c_iStars; j++)
		{
			rOut<<"r"<<i<<"-"<<j<<",";
		}
		rvOut<<"Rx"<<i<<" Ry"<<i<<" Rz"<<i<<" Vx"<<i<<" Vy"<<i<<" Vz"<<i<<" ";
		noodleOut<<"X"<<i<<","<<"Y"<<i<<","<<"Z"<<i<<",";
	}
	rvOut<<std::endl;
	rOut<<std::endl;
	noodleOut<<std::endl;
	eOut<<"T,Etot,PE,KE,abs(E-.5PE)/(-Eo),abs(E-Eo)/(-Eo)"<<std::endl;
	
	cmOut<<"t,CMx,CMy,CMz,VCMx,VCMy,VCMz"<<std::endl;
	abOut<<"t,AverageBinary,=Sum(#Binaries*dt)/dTime"<<std::endl;
	m_dAvBinSum=0;
	m_dLastT=0;
	
}
void outputManager::doMainOut(star *stars[],const double &dTime, double dRs[g_c_iStars][g_c_iStars])
{
	int iBinary,iNGroups,iNMembers;
	int iNaries[g_c_iMaxSystem];
	int iMembers[g_c_iMaxSystem];
	double dPE,dKE,dEnergy;
	
	iBinary=binaryCount(stars);
	binaryOut<<dTime<<","<<iBinary<<std::endl;
	
	m_dAvBinSum+=iBinary*(dTime-m_dLastT);
	abOut<<dTime<<","<<m_dAvBinSum/dTime<<std::endl;
	m_dLastT=dTime;
	
	for(int i=0;i<g_c_iMaxSystem;i++)
		iNaries[i]=0;
	
	noodleOut<<dTime<<",";
	rOut<<dTime<<",";
	rvOut<<dTime<<" ";
	for(int i = 0; i<g_c_iStars; i++)
	{
		for(int j =i+1; j<g_c_iStars; j++)
		{
			rOut<<sqrt(dRs[i][j])<<",";
		}
		rvOut<<stars[i]->m_vPosition.X<<" "<<stars[i]->m_vPosition.Y<<" "<<stars[i]->m_vPosition.Z<<" "<<stars[i]->m_vVelocity.X<<" "<<stars[i]->m_vVelocity.Y<<" "<<stars[i]->m_vVelocity.Z<<" ";
		noodleOut<<stars[i]->m_vPosition.X<<","<<stars[i]->m_vPosition.Y<<","<<stars[i]->m_vPosition.Z<<",";
	}
	rvOut<<std::endl;
	rOut<<std::endl;
	noodleOut<<std::endl;
	dEnergy=energyNotVerlet(stars,dRs,g_c_iStars,dPE,dKE);
	eOut<<dTime<<","<<dEnergy<<","<<dPE<<","<<dKE<<","<<abs(dEnergy-.5*dPE)/-m_dE0<<","<<abs(dEnergy-m_dE0)/-m_dE0<<std::endl;
}
void outputManager::doMcNeilOut(star *stars[],const double dTime,int i, int j)
{
	double dEnergy,dCm,dVCm,dTotalMass,dRrms,dVrms,dDRrms;
	Vector3d vVCm,vCm,vL;
	
	dEnergy=energy(stars);
	for(int i=0;i<g_c_iStars;i++)
	{
		dTotalMass+=stars[i]->m_dMass;
		vVCm+=stars[i]->m_vVelocity*stars[i]->m_dMass;
		vCm+=stars[i]->m_vPosition*stars[i]->m_dMass;
	}
	vVCm/=dTotalMass;
	vCm/=dTotalMass;
	dCm=vCm.Size();
	dVCm=vVCm.Size();
	
	RVDRrmsL(stars,dRrms,dVrms,dDRrms,vL);
	
	mcOut<<dTime<<","<<abs(dEnergy-m_dE0)/-m_dE0<<","<<binaryCount(stars)<<","<<stars[i]->m_vPosition.X<<","<<stars[i]->m_vPosition.Y<<","<<stars[i]->m_vPosition.Z<<","<<stars[j]->m_vPosition.X<<","<<stars[j]->m_vPosition.Y<<","<<stars[j]->m_vPosition.Z<<","<<(stars[i]->m_vPosition - stars[j]->m_vPosition).Size()<<","<<dCm<<","<<dVCm<<","<<dRrms<<","<<dVrms<<","<<dDRrms<<","<<vL.X<<","<<vL.Y<<","<<vL.Z<<std::endl;
}
void outputManager::doScreenOutOpt(star *stars[],const double &dE, const double &dTime, const double &dStopTime,const double c_dDt)
{
	double dTimeOut;
	char sTimeOutTag;

	gettimeofday(&m_tEnd, NULL);
	dTimeOut=(m_tEnd.tv_sec+((m_tEnd.tv_usec/1000000.0))-(m_tStart.tv_sec+(m_tStart.tv_usec/1000000.0)))/(dTime/dStopTime)-(m_tEnd.tv_sec+((m_tEnd.tv_usec/1000000.0))-(m_tStart.tv_sec+(m_tStart.tv_usec/1000000.0)));
	//convert to better units
	if(dTimeOut<60)//time best measured in seconds
	{
		sTimeOutTag='s';
	}
	else if(dTimeOut<3600)//time best measured in minutes
	{
		dTimeOut/=60.0;
		sTimeOutTag='m';
	}
	else if(dTimeOut<86400)
	{
		dTimeOut/=3600.0;
		sTimeOutTag='h';
	}
	else
	{
		dTimeOut/=86400.0;
		sTimeOutTag='d';
	}
	std::cout<<dTime/dStopTime<<" "<<dTimeOut<<sTimeOutTag<<" "<<abs(dE-m_dE0)/-m_dE0<<" "<<c_dDt<<std::endl;
}
void outputManager::doFileOutOpt(star *stars[],const double &dTime,double &dEnergy)
{
	double dCm,dVCm,dTotalMass,dRrms,dVrms,dDRrms,dPE,dKE;
	double dMu,dEP,dRMin,dR,dRTemp,dMinR,dSplodyTime;
	Vector3d vVCm,vCm,vL,vR,vV;
	int iBinary,iRight,iLeft;
	char sSplodyTag;
	
	dTotalMass=0;
	dPE=0;
	dKE=0;
	iBinary=0;
	dRrms=0;
	dVrms=0;
	dDRrms=0;
	dMinR=1E100;
	vL.SetZero();
	for(int i=0;i<g_c_iStars;i++)
	{
		//calc cCM and CM sums
		dTotalMass+=stars[i]->m_dMass;
		vVCm+=stars[i]->m_vVelocity*stars[i]->m_dMass;
		vCm+=stars[i]->m_vPosition*stars[i]->m_dMass;
		
		for(int j =i+1; j<g_c_iStars; j++)
		{
			//Binary counting stuff
			dMu=stars[i]->m_dMass*stars[j]->m_dMass/(stars[i]->m_dMass+stars[j]->m_dMass);
			vV=stars[i]->m_vVelocity - stars[j]->m_vVelocity;
			vR=stars[i]->m_vPosition-stars[j]->m_vPosition;
			dR=vR.Size();
			if(dR<dMinR)
			{
				dMinR=dR;
			}
			dEP=0.5*dMu*vV.SizeSquared() - stars[i]->m_dMass*stars[j]->m_dMass/dR;
			if(dEP<0)//could be a binary check for local stars
			{
				dRMin=100*dR;
				for(int k=0;k<g_c_iStars;k++)
				{
					if((k!=i)&&(k!=j))
					{
						dRTemp=(stars[i]->m_vPosition-stars[k]->m_vPosition).Size();
						if(dRTemp<dRMin)
							dRMin=dRTemp;
						
						dRTemp=(stars[j]->m_vPosition-stars[k]->m_vPosition).Size();
						if(dRTemp<dRMin)
							dRMin=dRTemp;
					}
				}
				
				if(dRMin>g_c_iRMultipy*dR)
				{
					iBinary++;
					m_binTracker.addBinary(i,j,dR,stars[i]->m_dMass,stars[j]->m_dMass,dTime);
				}
			}//end energy if
			
			//calc PE
			dPE+=-g_c_dSG*stars[i]->m_dMass*stars[j]->m_dMass/vR.Size();
			
			//Average seperation
			dDRrms+=vR.SizeSquared();
		}//end iner for
		//Calc KE
		dKE+=.5*stars[i]->m_dMass*stars[i]->m_vVelocity.SizeSquared();
		//rms quantities
		dRrms+=stars[i]->m_vPosition.SizeSquared();
		dVrms+=stars[i]->m_vVelocity.SizeSquared();
		vL+=stars[i]->m_vPosition.Cross(stars[i]->m_vVelocity)*stars[i]->m_dMass;
	}
	//calc cCM and CM
	vVCm/=dTotalMass;
	vCm/=dTotalMass;
	dCm=vCm.Size();
	dVCm=vVCm.Size();
	
	//rms finalization
	dRrms=sqrt(dRrms/(double)g_c_iStars);
	dVrms=sqrt(dVrms/(double)g_c_iStars);
	dDRrms=sqrt(dDRrms/g_c_dPairs);
	
	dEnergy=dPE+dKE;
	
	m_dAvBinSum+=iBinary*(dTime-m_dLastT);
	m_dLastT=dTime;
	
	m_binTracker.checkTerm(dTime);
	
	optOut<<dTime<<","<<dEnergy<<","<<abs(dEnergy-m_dE0)/-m_dE0<<","<<abs(dEnergy-.5*dPE)/-m_dE0<<","<<iBinary<<","<<m_dAvBinSum/dTime<<","<<m_binTracker.m_dAveR<<","<<m_binTracker.m_dAveT<<","<<m_binTracker.m_dAvePM<<","<<m_binTracker.m_dAveSM<<","<<dMinR<<","<<dCm<<","<<dVCm<<","<<dRrms<<","<<dVrms<<","<<dDRrms<<","<<vL.Size()<<std::endl;
	if((dEnergy>0&&!m_bSplody)||abs(dEnergy-m_dE0)/-m_dE0>0.01)
	{
		system("xmms -p /mnt/ROUTER/public/reaserch/sounds/splody.wav&");
		m_bSplody=true;
		gettimeofday(&m_tEnd, NULL);
		dSplodyTime=(m_tEnd.tv_sec+((m_tEnd.tv_usec/1000000.0))-(m_tStart.tv_sec+(m_tStart.tv_usec/1000000.0)));
		//convert to better units
		if(dSplodyTime<60)//time best measured in seconds
		{
			sSplodyTag='s';
		}
		else if(dSplodyTime<3600)//time best measured in minutes
		{
			dSplodyTime/=60.0;
			sSplodyTag='m';
		}
		else if(dSplodyTime<86400)
		{
			dSplodyTime/=3600.0;
			sSplodyTag='h';
		}
		else
		{
			dSplodyTime/=86400.0;
			sSplodyTag='d';
		}
		
		splodyOut.open("./codeIo/splody.txt");
		splodyOut<<"Simulation Time of explosion "<<dTime<<std::endl;
		splodyOut<<"Run Time of explosion "<<dSplodyTime<<sSplodyTag<<std::endl;
		splodyOut.close();
		exit(1);
	}
	
}
void outputManager::openOpt()
{
	system("rm ./codeIo/last/*.csv");
	system("rm ./codeIo/last/*.txt");
	system("cp ./codeIo/*.* ./codeIo/last/");
	system("rm ./codeIo/*.csv");
	system("rm ./codeIo/*.txt");
	
	optOut.open("./codeIo/output.csv");
	optOut<<"T,E,DE,DV,Nb,AB,ABR,ABLT,ABPM,ABSM,MinR,cm,vcm,Rrms,Vrms,DRrms,L"<<std::endl;
}
void outputManager::doInitialOutOpt(star *stars[],const double &dTime)
{
	double dEnergy,dCm,dVCm,dTotalMass,dRrms,dVrms,dDRrms,dPE,dKE;
	double dMu,dEP,dRMin,dR,dRTemp,dMinR;
	Vector3d vVCm,vCm,vL,vR,vV;
	int iBinary,iRight,iLeft;
	
	dTotalMass=0;
	dPE=0;
	dKE=0;
	iBinary=0;
	dRrms=0;
	dVrms=0;
	dDRrms=0;
	dMinR=1E100;
	vL.SetZero();
	for(int i=0;i<g_c_iStars;i++)
	{
		//calc cCM and CM sums
		dTotalMass+=stars[i]->m_dMass;
		vVCm+=stars[i]->m_vVelocity*stars[i]->m_dMass;
		vCm+=stars[i]->m_vPosition*stars[i]->m_dMass;
		
		for(int j =i+1; j<g_c_iStars; j++)
		{
			//Binary counting stuff
			dMu=stars[i]->m_dMass*stars[j]->m_dMass/(stars[i]->m_dMass+stars[j]->m_dMass);
			vV=stars[i]->m_vVelocity - stars[j]->m_vVelocity;
			vR=stars[i]->m_vPosition-stars[j]->m_vPosition;
			dR=vR.Size();
			if(dR<dMinR)
			{
				dMinR=dR;
			}
			dEP=0.5*dMu*vV.SizeSquared() - stars[i]->m_dMass*stars[j]->m_dMass/dR;
			if(dEP<0)//could be a binary check for local stars
			{
				dRMin=100*dR;
				for(int k=0;k<g_c_iStars;k++)
				{
					if((k!=i)&&(k!=j))
					{
						dRTemp=(stars[i]->m_vPosition-stars[k]->m_vPosition).Size();
						if(dRTemp<dRMin)
							dRMin=dRTemp;
						
						dRTemp=(stars[j]->m_vPosition-stars[k]->m_vPosition).Size();
						if(dRTemp<dRMin)
							dRMin=dRTemp;
					}
				}
				
				if(dRMin>g_c_iRMultipy*dR)
				{
					iBinary++;
					m_binTracker.addBinary(i,j,dR,stars[i]->m_dMass,stars[j]->m_dMass,dTime);
				}
			}//end energy if
			
			//calc PE
			dPE+=-g_c_dSG*stars[i]->m_dMass*stars[j]->m_dMass/vR.Size();
			
			//Average seperation
			dDRrms+=vR.SizeSquared();
		}//end iner for
		//Calc KE
		dKE+=.5*stars[i]->m_dMass*stars[i]->m_vVelocity.SizeSquared();
		//rms quantities
		dRrms+=stars[i]->m_vPosition.SizeSquared();
		dVrms+=stars[i]->m_vVelocity.SizeSquared();
		vL+=stars[i]->m_vPosition.Cross(stars[i]->m_vVelocity)*stars[i]->m_dMass;
	}
	//calc cCM and CM
	vVCm/=dTotalMass;
	vCm/=dTotalMass;
	dCm=vCm.Size();
	dVCm=vVCm.Size();
	
	//rms finalization
	dRrms=sqrt(dRrms/(double)g_c_iStars);
	dVrms=sqrt(dVrms/(double)g_c_iStars);
	dDRrms=sqrt(dDRrms/g_c_dPairs);
	
	dEnergy=dPE+dKE;
	m_dE0=dEnergy;
	
	gettimeofday(&m_tStart, NULL);
	
	m_dAvBinSum=0;
	m_dLastT=0;
	
	optOut<<dTime<<","<<dEnergy<<","<<abs(dEnergy-m_dE0)/-m_dE0<<","<<abs(dEnergy-.5*dPE)/-m_dE0<<","<<iBinary<<","<<iBinary<<","<<m_binTracker.m_dAveR<<","<<m_binTracker.m_dAveT<<","<<m_binTracker.m_dAvePM<<","<<m_binTracker.m_dAveSM<<","<<dMinR<<","<<dCm<<","<<dVCm<<","<<dRrms<<","<<dVrms<<","<<dDRrms<<","<<vL.Size()<<std::endl;
}
void outputManager::closeOpt()
{
	optOut.close();
}
void outputManager::doMainOut1(star *stars[],const double &dTime, double dRs[g_c_iStars][g_c_iStars],const double c_dDt)
{
	int iBinary,iNGroups,iNMembers;
	int iNaries[g_c_iMaxSystem];
	int iMembers[g_c_iMaxSystem];
	int iFirst[g_c_iStars/2];
	int iSecond[g_c_iStars/2];
	double dRBinary[g_c_iStars/2];
	double dPE,dKE,dEnergy;
	
	dtOut<<dTime<<","<<c_dDt<<std::endl;
	
	iBinary=binaryCountIndex(stars,iFirst,iSecond,dRBinary);
	binaryOut<<dTime<<","<<iBinary<<std::endl;
	if(iBinary>0)
	{
		bIndexOut<<dTime<<","<<iBinary<<",";
	}
	for(int i=0;i<iBinary;i++)
	{
		bIndexOut<<iFirst[i]<<","<<iSecond[i]<<","<<dRBinary[i]<<",";
	}
	if(iBinary>0)
	{
		bIndexOut<<std::endl;
	}

	
	m_dAvBinSum+=iBinary*(dTime-m_dLastT);
	abOut<<dTime<<","<<m_dAvBinSum/dTime<<std::endl;
	m_dLastT=dTime;
	
	noodleOut<<dTime<<",";
	rOut<<dTime<<",";
	rvOut<<dTime<<" ";
	for(int i = 0; i<g_c_iStars; i++)
	{
		if(g_c_iStars<11)
		{
			for(int j =i+1; j<g_c_iStars; j++)
			{
				rOut<<sqrt(dRs[i][j])<<",";
			}
		}
		rvOut<<stars[i]->m_vPosition.X<<" "<<stars[i]->m_vPosition.Y<<" "<<stars[i]->m_vPosition.Z<<" "<<stars[i]->m_vVelocity.X<<" "<<stars[i]->m_vVelocity.Y<<" "<<stars[i]->m_vVelocity.Z<<" ";
		noodleOut<<stars[i]->m_vPosition.X<<","<<stars[i]->m_vPosition.Y<<","<<stars[i]->m_vPosition.Z<<",";
	}
	rvOut<<std::endl;
	rOut<<std::endl;
	noodleOut<<std::endl;
	dEnergy=energy(stars,dPE,dKE);
	eOut<<dTime<<","<<dEnergy<<","<<dPE<<","<<dKE<<","<<abs(dEnergy-.5*dPE)/-m_dE0<<","<<abs(dEnergy-m_dE0)/-m_dE0<<std::endl;
}

void outputManager::doCmOut(star *stars[],const double &dTime)
{
	Vector3d vVCm,vCm;
	double dTotalMass;
	
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
	
	cmOut<<dTime<<","<<vCm.X<<","<<vCm.Y<<","<<vCm.Z<<","<<vVCm.X<<","<<vVCm.Y<<","<<vVCm.Z<<std::endl;
	
}
void outputManager::doDistroOut(star *stars[],const double &dTime,const double &dDt)
{
	saveOut.open("./codeIo/initDistro.txt");
	saveDistro(stars,saveOut,dTime,dDt);
	saveOut.close();
}
void outputManager::doInitialOut(star *stars[],const double &dTime)
{
	int iBinary;
	double dEnergy,dPE,dKE;
	
	iBinary=binaryCount(stars);
	m_iOldBinary=iBinary;
	binaryOut<<dTime<<","<<iBinary<<std::endl;
	noodleOut<<dTime*g_c_dYearConvert<<",";
	rvOut<<dTime<<" ";
	for(int i=0;i<g_c_iStars;i++)
	{
		noodleOut<<stars[i]->m_vPosition.X<<","<<stars[i]->m_vPosition.Y<<","<<stars[i]->m_vPosition.Z<<",";
		rvOut<<stars[i]->m_vPosition.X<<" "<<stars[i]->m_vPosition.Y<<" "<<stars[i]->m_vPosition.Z<<" "<<stars[i]->m_vVelocity.X<<" "<<stars[i]->m_vVelocity.Y<<" "<<stars[i]->m_vVelocity.Z<<" ";
	}
	rvOut<<std::endl;
	noodleOut<<std::endl;
	dEnergy=energyNotVerlet(stars,g_c_iStars,dPE, dKE);
	m_dE0=dEnergy;
	eOut<<dTime<<","<<dEnergy<<","<<dPE<<","<<dKE<<","<<abs(dEnergy-.5*dPE)/-m_dE0<<","<<abs(dEnergy-m_dE0)/-m_dE0<<std::endl;
	gettimeofday(&m_tStart, NULL);
	m_bSplody=false;
}
void outputManager::close()
{
	eOut.close();
	noodleOut.close();
	rOut.close();
	binaryOut.close();
	cmOut.close();
	rvOut.close();
	dtOut.close();
	bIndexOut.close();
	mcOut.close();
}
void outputManager::doMainScrenOut(star *stars[],double dRs[g_c_iStars][g_c_iStars], const double &dTime, const double &dStopTime)
{
	double dTimeOut,dEnergy,dKE,dPE;
	char sTimeOutTag;

	gettimeofday(&m_tEnd, NULL);
	dTimeOut=(m_tEnd.tv_sec+((m_tEnd.tv_usec/1000000.0))-(m_tStart.tv_sec+(m_tStart.tv_usec/1000000.0)))/(dTime/dStopTime)-(m_tEnd.tv_sec+((m_tEnd.tv_usec/1000000.0))-(m_tStart.tv_sec+(m_tStart.tv_usec/1000000.0)));
	//convert to better units
	if(dTimeOut<60)//time best measured in seconds
	{
		sTimeOutTag='s';
	}
	else if(dTimeOut<3600)//time best measured in minutes
	{
		dTimeOut/=60.0;
		sTimeOutTag='m';
	}
	else if(dTimeOut<86400)
	{
		dTimeOut/=3600.0;
		sTimeOutTag='h';
	}
	else
	{
		dTimeOut/=86400.0;
		sTimeOutTag='d';
	}
	dEnergy=energyNotVerlet(stars,dRs,g_c_iStars);
	if(dEnergy>0&&!m_bSplody)
	{
		system("xmms -p /mnt/ROUTER/public/reaserch/sounds/splody.wav&");
		m_bSplody=true;
	}
	dEnergy=energyNotVerlet(stars,dRs,g_c_iStars,dPE,dKE);
	std::cout<<dTime/dStopTime<<" "<<dTimeOut<<sTimeOutTag<<" "<<abs(dEnergy-m_dE0)/-m_dE0<<std::endl;
}
void outputManager::doMainScrenOutClose(star *stars[],double dRs[g_c_iStars][g_c_iStars], const double &dTime, const double &dStopTime)
{
	double dTimeOut,dEnergy,dKE,dPE;
	char sTimeOutTag;

	gettimeofday(&m_tEnd, NULL);
	dTimeOut=(m_tEnd.tv_sec+((m_tEnd.tv_usec/1000000.0))-(m_tStart.tv_sec+(m_tStart.tv_usec/1000000.0)))/(dTime/dStopTime)-(m_tEnd.tv_sec+((m_tEnd.tv_usec/1000000.0))-(m_tStart.tv_sec+(m_tStart.tv_usec/1000000.0)));
	//convert to better units
	if(dTimeOut<60)//time best measured in seconds
	{
		sTimeOutTag='s';
	}
	else if(dTimeOut<3600)//time best measured in minutes
	{
		dTimeOut/=60.0;
		sTimeOutTag='m';
	}
	else if(dTimeOut<86400)
	{
		dTimeOut/=3600.0;
		sTimeOutTag='h';
	}
	else
	{
		dTimeOut/=86400.0;
		sTimeOutTag='d';
	}
	dEnergy=energyNotVerlet(stars,dRs,g_c_iStars);
	if(dEnergy>0&&!m_bSplody)
	{
		system("xmms -p /mnt/ROUTER/public/reaserch/sounds/splody.wav&");
		m_bSplody=true;
	}
	dEnergy=energyNotVerlet(stars,dRs,g_c_iStars,dPE,dKE);
	std::cout<<dTime/dStopTime<<" "<<dTimeOut<<sTimeOutTag<<" "<<abs(dEnergy-m_dE0)/-m_dE0<<" "<<dRs[0][1]<<std::endl;
}
void outputManager::doMainScrenOut1(star *stars[],double dRs[g_c_iStars][g_c_iStars], const double &dTime, const double &dStopTime,const double c_dDt)
{
	double dTimeOut,dEnergy,dKE,dPE;
	char sTimeOutTag;

	gettimeofday(&m_tEnd, NULL);
	dTimeOut=(m_tEnd.tv_sec+((m_tEnd.tv_usec/1000000.0))-(m_tStart.tv_sec+(m_tStart.tv_usec/1000000.0)))/(dTime/dStopTime)-(m_tEnd.tv_sec+((m_tEnd.tv_usec/1000000.0))-(m_tStart.tv_sec+(m_tStart.tv_usec/1000000.0)));
	//convert to better units
	if(dTimeOut<60)//time best measured in seconds
	{
		sTimeOutTag='s';
	}
	else if(dTimeOut<3600)//time best measured in minutes
	{
		dTimeOut/=60.0;
		sTimeOutTag='m';
	}
	else if(dTimeOut<86400)
	{
		dTimeOut/=3600.0;
		sTimeOutTag='h';
	}
	else
	{
		dTimeOut/=86400.0;
		sTimeOutTag='d';
	}
	dEnergy=energyNotVerlet(stars,dRs,g_c_iStars);
	if(dEnergy>0&&!m_bSplody)
	{
		system("xmms -p /mnt/ROUTER/public/reaserch/sounds/splody.wav&");
		m_bSplody=true;
	}
	dEnergy=energy(stars,dPE,dKE);
	std::cout<<dTime/dStopTime<<" "<<dTimeOut<<sTimeOutTag<<" "<<abs(dEnergy-m_dE0)/-m_dE0<<" "<<c_dDt<<std::endl;
}
void outputManager::doFinalScreenOut()
{
	double dTimeOut;
	char sTimeOutTag;
	std::ofstream timeOut;
	
	//errorOut<<g_c_dDt*365<<","<<dMaxErrorE<<","<<dMaxErrorR<<std::endl;
	gettimeofday(&m_tEnd, NULL);
	dTimeOut=(m_tEnd.tv_sec+((m_tEnd.tv_usec/1000000.0))-(m_tStart.tv_sec+(m_tStart.tv_usec/1000000.0)));
	//convert to better units
	if(dTimeOut<60)//time best measured in seconds
	{
		sTimeOutTag='s';
	}
	else if(dTimeOut<3600)//time best measured in minutes
	{
		dTimeOut/=60.0;
		sTimeOutTag='m';
	}
	else if(dTimeOut<86400)
	{
		dTimeOut/=3600.0;
		sTimeOutTag='h';
	}
	else
	{
		dTimeOut/=86400.0;
		sTimeOutTag='d';
	}
	
	std::cout<<"Time to complete = "<<dTimeOut<<sTimeOutTag<<std::endl;
	
	timeOut.open("./codeIo/runtime.txt");
	timeOut<<"Time to complete = "<<dTimeOut<<sTimeOutTag<<std::endl;
	timeOut.close();
}
