/*group1_0.cpp
 *Written by Mike Williams 07-16-03
 has getMembers function
*/

#include"group1_0.h"
#include<iostream>

group::group()
{//Constructor
	clear();
}

group::~group()
{//destructor
	//do nothing
}

bool group::isMember(const int &index)
{//search m_iMembers[] return true if index is present
	for(int i = 0; i<m_iNumber; i++)
	{
		if(m_iMembers[i]==index)
			return true;
	}
	return false;
}

void group::addMember(const int &index)
{//add index to list
	if(!isMember(index))
	{
		m_iMembers[m_iNumber]=index;
		m_iNumber++;
	}
}

void group::removeMember(const int &index)
{
	for(int i = 0; i<m_iNumber-1; i++)
	{
		if(m_iMembers[i]==index)
		{
			m_iMembers[i]=m_iMembers[m_iNumber];
			m_iMembers[m_iNumber]=UNUSED_INDEX;
			m_iNumber--;
		}
	}
	if(m_iMembers[m_iNumber]==index)
	{
		m_iMembers[m_iNumber]=UNUSED_INDEX;
		m_iNumber--;
	}
}

void group::getMembers(int members[],int &n)
{
	for(int i = 0; i<m_iNumber;i++)
	{
		members[i]=m_iMembers[i];
	}
	n=m_iNumber;
}

void group::clear()
{
	m_iNumber = 0;
	for(int i = 0; i<g_c_iMaxSystem; i++)
	{
		m_iMembers[i]=UNUSED_INDEX;
	}
	m_dDt=1E10;//set dt arbitrarily large so that the first suggestion is taken
}

int group::numberOfMembers()
{
	return m_iNumber;
}

void group::suggestDt(const double &a_dDt)
{
	if(a_dDt<m_dDt)//ensure that the smallest dt of the system is used
	{
		m_dDt=a_dDt;
	}
}
