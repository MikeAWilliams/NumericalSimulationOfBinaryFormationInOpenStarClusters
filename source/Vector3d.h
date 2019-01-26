/*Vector3d.h
Written by Mike Williams 10-29-02
re worked to use only double valuses as an attempt to explain wierd behaviour noted in notes today
*/
#ifndef _vector3d_h
#define _vector3d_h

// Vector3d Class
class Vector3d
{
public:
	double X;									// Components
	double Y;
	double Z;
	double W;									// Padding

private:
	const static Vector3d ms_Null;				// Null Instance

public:
	// Constructors
	Vector3d() : X(0.0),Y(0.0),Z(0.0),W(0.0) {}
	Vector3d(double dX,double dY) : X(dX),Y(dY),Z(0.0),W(0.0) {}
	Vector3d(double dX,double dY,double dZ) : X(dX),Y(dY),Z(dZ),W(0.0) {}
	Vector3d(const Vector3d &v) : X(v.X),Y(v.Y),Z(v.Z),W(0.0) {}
	Vector3d(const Vector3d &A,const Vector3d &B) : X(B.X-A.X),Y(B.Y-A.Y),Z(B.Z-A.Z),W(0.0) {}

	// Null Functions
	const static Vector3d &GetNull() {return ms_Null;}

	// Set Functions
	void SetZero() {X=Y=Z=W=0.0;}
	void Set(double dX,double dY,double dZ) {X=dX; Y=dY; Z=dZ;}
	void SetW(double dW) {W=dW;}
	void Negate() {X=-X; Y=-Y; Z=-Z;}

	// Operators
	Vector3d &operator=(const Vector3d &v) {X=v.X; Y=v.Y; Z=v.Z; return *this;}
	Vector3d &operator+=(const Vector3d &v) {X+=v.X; Y+=v.Y; Z+=v.Z; return *this;}
	Vector3d &operator-=(const Vector3d &v) {X-=v.X; Y-=v.Y; Z-=v.Z; return *this;}
	Vector3d &operator*=(const double d) {X*=d; Y*=d; Z*=d; return *this;}
	Vector3d &operator/=(const double d) {return operator*=(1.0/d);}

	Vector3d operator+(const Vector3d &v) const {return Vector3d(X+v.X,Y+v.Y,Z+v.Z);}
	Vector3d operator-(const Vector3d &v) const {return Vector3d(X-v.X,Y-v.Y,Z-v.Z);}
	Vector3d operator*(const double d) const {return Vector3d(X*d,Y*d,Z*d);}
	double operator*(const Vector3d &v) const {return (X*v.X+Y*v.Y+Z*v.Z);}
	Vector3d operator/(const double d) const {return operator*(1.0/d);}

	int operator==(const Vector3d &v) const {return ((X==v.X && Y==v.Y && Z==v.Z)?1:0);}
	int operator!=(const Vector3d &v) const {return ((X!=v.X || Y!=v.Y || Z!=v.Z)?1:0);}

	// Products
	double Dot(const Vector3d &v) const {return (X*v.X+Y*v.Y+Z*v.Z);}
	Vector3d Cross(const Vector3d &v) const {return Vector3d(Y*v.Z-v.Y*Z,-X*v.Z+v.X*Z,X*v.Y-v.X*Y);}

	// Angle Functions
	double Cosine(const Vector3d &v) {return Dot(v)/(Size()*v.Size());}

	// Size Functions
	double Size() const;
	double SizeSquared() const {return (X*X+Y*Y+Z*Z);}
	Vector3d Normal() const {return operator/(Size());}
	Vector3d &Normalize() {return operator/=(Size());}

	// Project Functions
	double ProjectedDistance(const Vector3d &v) const {return Dot(v)/v.Size();}
	Vector3d Project(const Vector3d &v) const {double d=Dot(v)/v.SizeSquared(); return Vector3d(v.X*d,v.Y*d,v.Z*d);}
	void SetProject(const Vector3d &v) {double d=Dot(v)/v.SizeSquared(); X=v.X*d; Y=v.Y*d; Z=v.Z*d;}
	Vector3d ProjectNormal(const Vector3d &v) const {double d=Dot(v)/v.SizeSquared(); return Vector3d(X-v.X*d,Y-v.Y*d,Z-v.Z*d);}
	void SetProjectNormal(const Vector3d &v) {double d=Dot(v)/v.SizeSquared(); X-=v.X*d; Y-=v.Y*d; Z-=v.Z*d;}

	// Rotation Functions
	Vector3d RotateX(const double dAngle) const;
	Vector3d &SetRotateX(const double dAngle);
	Vector3d RotateY(const double dAngle) const;
	Vector3d &SetRotateY(const double dAngle);
	Vector3d RotateZ(const double dAngle) const;
	Vector3d &SetRotateZ(const double dAngle);
	Vector3d RotateN(const double dAngle,const Vector3d vN)const;
};


#endif
