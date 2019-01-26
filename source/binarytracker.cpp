//binarytracker.cpp
//written by Mike Williams 05-18-05

#include"binarytracker.h"

binaryTracker::binaryTracker()
{
	m_dAveR=0;
	m_dAveT=0;
	m_iBinaries=0;
	m_dRSum=0;
	m_iRMeasures=0;
	m_dTSum=0;
	m_iTMeasures=0;
	m_dPMSum=0;
	m_dSMSum=0;
	m_iMMeasures=0;
	m_dAvePM=0;
	m_dAveSM=0;
}
binaryTracker::~binaryTracker()
{
	//empty
}

void binaryTracker::addBinary(const int &i,const int &j,const double &dR,const double &dM1,const double &dM2,const double &dT)
{
	bool bFound;
	
	bFound=false;
	
	
	//check if this is a new binary
	for(int iPair=0;iPair<m_iBinaries;iPair++)
	{
		if(((m_iIndicies[iPair][0]==i)&&(m_iIndicies[iPair][1]==j))||((m_iIndicies[iPair][0]==j)&&(m_iIndicies[iPair][1]==i)))
		{
			m_bCurrent[iPair]=bFound=true;
			break;
		}
	}
	if(!bFound)
	{
		m_iIndicies[m_iBinaries][0]=i;
		m_iIndicies[m_iBinaries][1]=j;
		m_dTimes[m_iBinaries]=dT;
		m_bCurrent[m_iBinaries]=true;
		m_iBinaries++;
		//update average r stuff
		m_dRSum+=dR;
		m_iRMeasures++;
		m_dAveR=m_dRSum/m_iRMeasures;
		
		//add masses
		if(dM1>dM2)
		{
			m_dPMSum+=dM1;
			m_dSMSum+=dM2;
		}
		else
		{
			m_dPMSum+=dM2;
			m_dSMSum+=dM1;	
		}
		
		m_iMMeasures++;
		m_dAvePM=m_dPMSum/m_iMMeasures;
		m_dAveSM=m_dSMSum/m_iMMeasures;
	}
}

void binaryTracker::checkTerm(const double &dT)
{
	bool bSwitch;
	bSwitch=true;
	//make sure the last pair is current
	if(m_iBinaries>0)
	{	
		while(bSwitch)
		{
			if(!m_bCurrent[m_iBinaries-1])
			{
				m_dTSum+=(dT-m_dTimes[m_iBinaries-1]);
				m_iTMeasures++;
				m_dAveT=m_dTSum/m_iTMeasures;
				m_iBinaries--;
			}
			else
			{
				bSwitch=false;
			}
		}
	}
	//now that I know the last pair is current i can check the rest
	if(m_iBinaries>1)
	{
		for(int iPair=0;iPair<m_iBinaries;iPair++)
		{
			if(!m_bCurrent[iPair])
			{
				m_dTSum+=(dT-m_dTimes[iPair]);
				m_iTMeasures++;
				m_dAveT=m_dTSum/m_iTMeasures;
				//put the last one here
				m_iIndicies[iPair][0]=m_iIndicies[m_iBinaries-1][0];
				m_iIndicies[iPair][1]=m_iIndicies[m_iBinaries-1][1];
				m_dTimes[iPair]=m_dTimes[m_iBinaries-1];
				
				m_iBinaries--;
			}
			//reset current for next time step
			m_bCurrent[iPair]=false;
		}
	}
	else
	{
		m_bCurrent[0]=false;
	}
}
