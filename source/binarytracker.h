//binarytracker.h
//written by Mike Williams 05-18-05

#ifndef binarytrack_h
#define binarytrack_h

const int g_c_iMaxBin=10;

class binaryTracker
{
public:
	//constructors, destructor
	binaryTracker();
	~binaryTracker();

	void addBinary(const int &i,const int &j,const double &dR,const double &dM1,const double &dM2,const double &dT);
	void checkTerm(const double &dT);

	//public variables
	double m_dAveR;
	double m_dAveT;
	double m_dAvePM;//primary mass average
	double m_dAveSM;//secondary mass average

private:
	//for binary tracking
	int m_iBinaries;
	int m_iIndicies[g_c_iMaxBin][2];
	double m_dTimes[g_c_iMaxBin];
	bool m_bCurrent[g_c_iMaxBin];
	
	//for average r
	double m_dRSum;
	int m_iRMeasures;

	//for average life time T
	double m_dTSum;
	int m_iTMeasures;

	//mass average
	double m_dPMSum,m_dSMSum;
	int m_iMMeasures;


};

#endif
