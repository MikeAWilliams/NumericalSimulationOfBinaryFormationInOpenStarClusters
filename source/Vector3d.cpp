//Vector3d.cpp used throughout the project

#include <math.h>

// Local Includes
#include "Vector3d.h"

// Static Initialization
const Vector3d Vector3d::ms_Null(0.0,0.0,0.0);

//==================================================================================
// Size Functions

// Size
double Vector3d::Size() const 
{
	// Return Euclidean Norm
	return sqrt(SizeSquared());
}

//==================================================================================
// Rotation Functions

// Rotate X
Vector3d Vector3d::RotateX(double dA) const
{
	// Compute Cosine & Sine
	double dC=cos(dA);
	double dS=sin(dA);
	// Compute Rotation
	return Vector3d(X,dC*Y-dS*Z,dS*Y+dC*Z);
}

// Set Rotate X
Vector3d &Vector3d::SetRotateX(double dA)
{
	// Compute Cosine & Sine
	double dC=cos(dA);
	double dS=sin(dA);
	// Compute Rotation
	Set(X,dC*Y-dS*Z,dS*Y+dC*Z);
	return *this;
}

// Rotate Y
Vector3d Vector3d::RotateY(double dA) const
{
	// Compute Cosine & Sine
	double dC=cos(dA);
	double dS=sin(dA);
	// Compute Rotation
	return Vector3d(dC*X+dS*Z,Y,-dS*X+dC*Z);
}

// Set Rotate Y
Vector3d &Vector3d::SetRotateY(double dA)
{
	// Compute Cosine & Sine
	double dC=cos(dA);
	double dS=sin(dA);
	// Compute Rotation
	Set(dC*X+dS*Z,Y,-dS*X+dC*Z);
	return *this;
}

// Rotate Z
Vector3d Vector3d::RotateZ(double dA) const
{
	// Compute Cosine & Sine
	double dC=cos(dA);
	double dS=sin(dA);
	// Compute Rotation
	return Vector3d(dC*X-dS*Y,dS*X+dC*Y,Z);
}

// Set Rotate Z
Vector3d &Vector3d::SetRotateZ(double dA)
{
	// Compute Cosine & Sine
	double dC=cos(dA);
	double dS=sin(dA);
	// Compute Rotation
	Set(dC*X-dS*Y,dS*X+dC*Y,Z);
	return *this;
}

//Rotate R
Vector3d Vector3d::RotateN(const double dAngle,const Vector3d vN) const
{
	double dC=cos(dAngle);
	double dS=sin(dAngle);
	
	return *this*dC+vN*vN.Dot(*this)*(1-dC) + vN.Cross(*this)*dS;
}
