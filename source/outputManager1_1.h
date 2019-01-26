/*outputManager1_1.h
Written by Mike Williams 04-13-05
*/

#ifndef outputManager_h
#define outputManager_h

#include<fstream>
#include<sys/time.h>
#include"constants.h"
#include"star6_0.h"
#include"binarytracker.h"

// Vector3d Class
class outputManager
{
public:
	void open();
	void openOpt();
	void doMainOut(star *stars[],const double &dTime, double dRs[g_c_iStars][g_c_iStars]);
	void doMainOut1(star *stars[],const double &dTime, double dRs[g_c_iStars][g_c_iStars],const double c_dDt);
	void doCmOut(star *stars[],const double &dTime);
	void doDistroOut(star *stars[],const double &dTime,const double &dDt);
	void doInitialOut(star *stars[],const double &dTime);
	void doMainScrenOut(star *stars[],double dRs[g_c_iStars][g_c_iStars], const double &dTime, const double &dStopTime);
	void doMainScrenOutClose(star *stars[],double dRs[g_c_iStars][g_c_iStars], const double &dTime, const double &dStopTime);
	void doMainScrenOut1(star *stars[],double dRs[g_c_iStars][g_c_iStars], const double &dTime, const double &dStopTime,const double c_dDt);
	void doFinalScreenOut();
	void doMcNeilOut(star *stars[],const double dTime,int i, int j);
	void doScreenOutOpt(star *stars[], const double &dE, const double &dTime, const double &dStopTime,const double c_dDt);
	void doFileOutOpt(star *stars[],const double &dTime,double &dEnergy);
	void doInitialOutOpt(star *stars[],const double &dTime);
	void close();
	void closeOpt();

private:
	std::ofstream noodleOut,saveOut,eOut,rOut,binaryOut,cmOut,rvOut,abOut,dtOut,bIndexOut,mcOut,optOut,splodyOut;
	double m_dE0,m_dAvBinSum,m_dLastT;
	int m_iOldBinary;
	timeval m_tStart, m_tEnd;
	bool m_bSplody;
	binaryTracker m_binTracker;
};

#endif
