//binaryCount.cpp
//written by Mike Williams 02-25-05

#include"binaryCount.h"
#include"util.h"

int binaryCount(star *stars[], double Rs[g_c_iStars][g_c_iStars])
{
	double dMu,dE,dRMin,dR,dRTemp;
	Vector3d vV;
	int iRight,iLeft,iBinary;
	
	
	iBinary=0;
	for(int i=0;i<g_c_iStars;i++)
	{
		for(int j=i+1;j<g_c_iStars;j++)
		{
			dMu=stars[i]->m_dMass*stars[j]->m_dMass/(stars[i]->m_dMass+stars[j]->m_dMass);
			vV=stars[i]->m_vVelocity - stars[j]->m_vVelocity;
			dR=sqrt(Rs[i][j]);
			dE=0.5*dMu*vV.SizeSquared() - stars[i]->m_dMass*stars[j]->m_dMass/dR;
			if(dE<0)//could be a binary check for local stars
			{
				dRMin=100*dR;
				for(int k=0;k<g_c_iStars;k++)
				{
					if((k!=i)&&(k!=j))
					{
						iRight=max(i,k);
						iLeft=min(i,k);
						dRTemp=sqrt(Rs[iLeft][iRight]);
						if(dRTemp<dRMin)
							dRMin=dRTemp;
						
						iRight=max(j,k);
						iLeft=min(j,k);
						dRTemp=sqrt(Rs[iLeft][iRight]);
						if(dRTemp<dRMin)
							dRMin=dRTemp;
					}
				}
				
				if(dRMin>g_c_iRMultipy*dR)
				{
					iBinary++;
				}
			}//end energy if
		}//end for j
	}//end for i
	
	return iBinary;
}//end function

int binaryCount(star *stars[])
{
	double dMu,dE,dRMin,dR,dRTemp;
	Vector3d vV,vR;
	int iRight,iLeft,iBinary;
	
	
	iBinary=0;
	for(int i=0;i<g_c_iStars;i++)
	{
		for(int j=i+1;j<g_c_iStars;j++)
		{
			dMu=stars[i]->m_dMass*stars[j]->m_dMass/(stars[i]->m_dMass+stars[j]->m_dMass);
			vV=stars[i]->m_vVelocity - stars[j]->m_vVelocity;
			vR=stars[i]->m_vPosition-stars[j]->m_vPosition;
			dR=vR.Size();
			dE=0.5*dMu*vV.SizeSquared() - stars[i]->m_dMass*stars[j]->m_dMass/dR;
			if(dE<0)//could be a binary check for local stars
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
				}
			}//end energy if
		}//end for j
	}//end for i
	
	return iBinary;
}//end function

int binaryCountIndex(star *stars[],int iFirst[],int iSecond[], double dRBinarys[])
{
	double dMu,dE,dRMin,dR,dRTemp;
	Vector3d vV,vR;
	int iRight,iLeft,iBinary;
	
	
	iBinary=0;
	for(int i=0;i<g_c_iStars;i++)
	{
		for(int j=i+1;j<g_c_iStars;j++)
		{
			dMu=stars[i]->m_dMass*stars[j]->m_dMass/(stars[i]->m_dMass+stars[j]->m_dMass);
			vV=stars[i]->m_vVelocity - stars[j]->m_vVelocity;
			vR=stars[i]->m_vPosition-stars[j]->m_vPosition;
			dR=vR.Size();
			dE=0.5*dMu*vV.SizeSquared() - stars[i]->m_dMass*stars[j]->m_dMass/dR;
			if(dE<0)//could be a binary check for local stars
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
					iFirst[iBinary]=i;
					iSecond[iBinary]=j;
					dRBinarys[iBinary]=dR;
					iBinary++;
				}
			}//end energy if
		}//end for j
	}//end for i
	
	return iBinary;
}//end function
