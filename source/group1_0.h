/*group1_0.h
 *Written by Mike Williams 07-16-04
 *was system renamed to group since system is a reserved word in this version of c++
 *Added a remove member function
 
 
*/

#ifndef group_h
#define group_h
//#include"star5_1.h"
#include"analytical.h"

const int g_c_iMaxSystem = 10;//maxiumum stars to expect in a system
const int UNUSED_INDEX = -1;//since a negative index is impossible i will use it as a tag to identify that an aray element contains no data

class group
{
public:
	//constructors, destructor
	group();
	~group();

	//functions
	bool isMember(const int &index);//search m_iMembers[] return true if index is present
	void addMember(const int &index);//add index to list
	void removeMember(const int &index);//remove index from the list
	void clear();
	void getMembers(int members[],int &n);//tell me who is in the system
	int numberOfMembers();//return m_iNumber
	void suggestDt(const double &a_dDt);

	//public variables
	double m_dDt;

private:
	//Variables
	int m_iMembers[g_c_iMaxSystem];//aray to store member indicies
	int m_iNumber;//number of stars in system
};


#endif
