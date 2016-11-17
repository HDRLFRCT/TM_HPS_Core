//-----------------------------------------------------------------------------------
//
//	3D Vector class
//
//	Anton M. Krivtsov
//
//	11.04.2001
//
//-----------------------------------------------------------------------------------

#ifndef ___Vect3D_h___
#define ___Vect3D_h___

#include "def.h"
#include <iostream>
using namespace std;
#define Vect3D_ const Vect3D&

class Box3D;
class Tensor;

//-----------------------------------------------------------------------------------

class Vect3D
{
public:
	struct	 SFLOAT								{ float x;	float y;  float z; };
	struct	 SDOUBLE							{ double x;	double y; double z; };

	Vect3D()									{}
	Vect3D(FLA x, FLA y, FLA z)					{ Set(x, y, z); }
	Vect3D(const SFLOAT&  r)					{ x = r.x; y = r.y; z = r.z; }
	Vect3D(const SDOUBLE& r)					{ x = r.x; y = r.y; z = r.z; }

	void Set(FLA x, FLA y, FLA z)				{ this->x = x; this->y = y; this->z = z; }

	void SetRand(FLA x_min, FLA y_min, FLA z_min, FLA x_max, FLA y_max, FLA z_max);
	void SetRand(const Box3D&);					//	defined in Box3D.cpp
	void SetRand(FLA r_max);

	operator SFLOAT()					const	{ SFLOAT r = {float(x), float(y), float(z)}; return r; }
	operator SDOUBLE()					const	{ SDOUBLE r = {x, y, z}; return r; }

	Vect3D& operator +=(Vect3D_ r)				{ x += r.x;	y += r.y; z += r.z; return *this; }
	Vect3D& operator -=(Vect3D_ r)				{ x -= r.x;	y -= r.y; z -= r.z; return *this; }
	Vect3D& operator *=(FLA k)					{ x *= k;	y *= k;   z *= k;	return *this; }
	Vect3D& operator /=(FLA k)					{ return *this *= 1 / k; }
	Vect3D& operator *=(const Tensor& T);		//	defined in Tensor.cpp

	Vect3D operator -()					const	{ Vect3D r(-x, -y, -z); return r; }

	BOOL operator ==(Vect3D_ r)			const	{ return x==r.x && y==r.y && z==r.z; }
	BOOL operator !=(Vect3D_ r)			const	{ return !operator==(r); }
	BOOL IsNull()						const	{ return !(x || y || z); }

	void FlipX()								{ x = -x; }
	void FlipY()								{ y = -y; }
	void FlipZ()								{ z = -z; }

	void RotateX(FLA angle);
	void RotateY(FLA angle);
	void RotateZ(FLA angle);

	void RotateX(FLA sin_f, FLA cos_f)			{ Rotate(y, z, sin_f, cos_f); }
	void RotateY(FLA sin_f, FLA cos_f)			{ Rotate(z, x, sin_f, cos_f); }
	void RotateZ(FLA sin_f, FLA cos_f)			{ Rotate(x, y, sin_f, cos_f); }

	void RotateZXZ(FLA phi, FLA theta, FLA psi)	{ RotateZ(phi); RotateX(theta); RotateZ(psi); }

	void   Swap		  (BYTE mode = 1);			// mode = 0,1,2
	Vect3D GetSwapped (BYTE mode = 1)	const;	// mode = 0,1,2

	FLT Sqr()							const	{ return x * x + y * y + z * z; }
	FLT Abs()							const;

	#ifdef _DEBUG
	void Dump()							const	{ afxDump << x << "\n"; afxDump << y << "\n"; afxDump << z << "\n"; }
	#endif

	FLT x, y, z;

private:
	static void Rotate(FLT& x, FLT& y, FLA sin_f, FLA cos_f);
};

//-----------------------------------------------------------------------------------

static const Vect3D VECT3D_0(0, 0, 0);		// zero vector

//-----------------------------------------------------------------------------------

static inline Vect3D operator +(Vect3D_ r1, Vect3D_ r2)	{ Vect3D rr(r1);	return rr += r2;	}
static inline Vect3D operator -(Vect3D_ r1, Vect3D_ r2)	{ Vect3D rr(r1);	return rr -= r2;	}
static inline Vect3D operator *(Vect3D_ r,  FLA k)		{ Vect3D rr(r);		return rr *= k;		}
static inline Vect3D operator *(FLA k,		Vect3D_ r)	{					return r * k;		}
static inline Vect3D operator /(Vect3D_ r,	FLA k)		{ Vect3D rr(r);		return rr /= k;		}

static inline Vect3D operator %(Vect3D_ r1, Vect3D_ r2)	// vector product
{ return Vect3D(
	r1.y * r2.z - r1.z * r2.y,
	r1.z * r2.x - r1.x * r2.z,
	r1.x * r2.y - r1.y * r2.x);
}

static inline FLT	 operator *(Vect3D_ r1, Vect3D_ r2)	{ return r1.x * r2.x + r1.y * r2.y + r1.z * r2.z; }

//-----------------------------------------------------------------------------------

static inline ostream& operator <<(ostream& stream, Vect3D_ r)	{ return stream << "(" << r.x << ", " << r.y << ", " << r.z << ")"; }
static inline istream& operator >>(istream& stream, Vect3D& r)	{ return stream >> r.x >> r.y >> r.z; }

//-----------------------------------------------------------------------------------

#endif //___Vect3D_h___

