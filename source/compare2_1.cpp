/*	compare2_1,cpp
	written by Mike Williams 03-15-05
*/

#include"star6_0.h"
#include"constants.h"
#include"compare2_1.h"
#include"util.h"

bool compareSoftSimple(star *star1, star *star2, Vector3d &vF,double &r,const double &c_dRMin)
{
	Vector3d vR;//radious between stars
	double dR;

	vR=star2->m_vPosition-star1->m_vPosition;
	r=(vR).SizeSquared();
	dR=sqrt(r);
	vF=vR.Normal()*star1->m_dMass*star2->m_dMass/((dR+c_dRMin)*(dR+c_dRMin));
	return true;
}



bool compare3Simple(star *star1, star *star2, Vector3d &vF,const double &rMin)
{
	Vector3d vR;//radious between stars

	vR=star2->m_vPosition-star1->m_vPosition;
	if(vR.SizeSquared()>rMin)//fix this
	{
		//calc force and acceleration
		//vF=vR.Normal()*g_c_dSG*star1->m_dMass*star2->m_dMass/vR.SizeSquared();
		vF=vR.Normal()*star1->m_dMass*star2->m_dMass/vR.SizeSquared();
		return true;
	}
	return false;
}

bool compare3Simple2(star *star1, star *star2, Vector3d &vF,double &r, const double &rMin)
{
	Vector3d vR;//radious between stars

	vR=star2->m_vPosition-star1->m_vPosition;
	r=vR.SizeSquared();
	if(r>rMin)//fix this
	{
		//calc force and acceleration
		//vF=vR.Normal()*g_c_dSG*star1->m_dMass*star2->m_dMass/vR.SizeSquared();
		vF=vR.Normal()*star1->m_dMass*star2->m_dMass/vR.SizeSquared();
		return true;
	}
	return false;
}
