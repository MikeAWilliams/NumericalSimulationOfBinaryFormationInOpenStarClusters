/*star6_0.cpp
written by Mike Williams 10-05-04
uses scailed units
Tracks previous position for Verlet type integration
Implements Jerk for Hermete integration
Stors intermediate force and jerk for hermete integration
Stors 4 velocity and position vectors for RK4
*/
#include"star6_0.h"
#include "Vector3d.h"
#include "constants.h"
#include <math.h>

star::star()
{

}

star::~star()
{
	//Null
}

void star::initialize(Vector3d p, Vector3d v,double m)
{
	m_vPosition = p;
	m_vVelocity = v;
	m_dMass = m;
	m_bClose=false;
}
