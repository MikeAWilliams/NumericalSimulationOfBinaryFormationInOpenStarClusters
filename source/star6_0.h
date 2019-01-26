/*star6_0.h 
written by Mike Williams 04-1-05
basic point particle.
Uses Vector3d.h
A diferent aproach then previous versons. uses public variables and limited functions.
Tracks previous position for verlet integration
Tracks jerk for hermete integration
Tracks intermediate force and jerk for hermete integration
Stors 4 position and velocity vectors for rk4
keeps 5 old positions for legrange interpolating polynomials.
*/
#ifndef _star_h
#define _star_h

#include "Vector3d.h"

class star
{
public:
	//Functions
	star();
	~star();
	void initialize(Vector3d p, Vector3d v, double m);
	
	//p=position of the star, v= velocity of the star, m=mass of the star, min=minimum size of radious squared for comparison to be calculated

	//public data
	Vector3d m_vPosition;
	Vector3d m_vPosition1;
	Vector3d m_vPosition2;
	Vector3d m_vPosition3;
	Vector3d m_vPosition4;
	Vector3d m_vPositionT;
	Vector3d m_vVelocity;
	Vector3d m_vVelocity1;
	Vector3d m_vVelocity2;
	Vector3d m_vVelocity3;
	Vector3d m_vVelocity4;
	Vector3d m_vAcceleration;
	Vector3d m_vTmpAcceleration;
	Vector3d m_vTmpAcceleration2;
	Vector3d m_vTmpAcceleration3;
	double m_dMass;
	bool m_bClose;
};

#endif
