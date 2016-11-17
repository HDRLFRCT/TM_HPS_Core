//-----------------------------------------------------------------------------------
//
//	2D Vector class
//
//	Anton M. Krivtsov
//
//	27.10.2001
//
//-----------------------------------------------------------------------------------

#ifndef ___Vect2D_h___
#define ___Vect2D_h___

#include "def.h"
#include <iostream>


#define Vect2D_ const Vect2D&

class Box2D;

//-----------------------------------------------------------------------------------

class Vect2D
{
public:
	struct	 SFLOAT								{ float x;	float y; };
	struct	 SDOUBLE							{ double x;	double y; };

	Vect2D()									{}
	Vect2D(FLA x, FLA y)						{ Set(x, y); }
	Vect2D(const SFLOAT&  r)					{ x = r.x; y = r.y; }
	Vect2D(const SDOUBLE& r)					{ x = r.x; y = r.y; }

	void Set(FLA x, FLA y)						{ this->x = x; this->y = y; }

	void SetRand(FLA x_min, FLA y_min, FLA x_max, FLA y_max);
	void SetRand(const Box2D&);					//	defined in Box2D.cpp
	void SetRand(FLA r_max);

	operator SFLOAT()					const	{ SFLOAT r = {float(x), float(y)}; return r; }
	operator SDOUBLE()					const	{ SDOUBLE r = {x, y}; return r; }

	Vect2D& operator +=(Vect2D_ r)				{ x += r.x;	y += r.y; return *this; }
	Vect2D& operator -=(Vect2D_ r)				{ x -= r.x;	y -= r.y; return *this; }
	Vect2D& operator *=(FLA k)					{ x *= k;	y *= k;   return *this; }
	Vect2D& operator /=(FLA k)					{ return *this *= 1 / k; }

	Vect2D operator -()					const	{ Vect2D r(-x, -y); return r; }

	BOOL operator ==(Vect2D_ r)			const	{ return x==r.x && y==r.y; }
	BOOL operator !=(Vect2D_ r)			const	{ return !operator==(r); }
	BOOL IsNull()						const	{ return !(x || y); }

	void FlipX()								{ x = -x; }
	void FlipY()								{ y = -y; }
	void Swap();

	void Rotate(FLA angle);
	void Rotate(FLA sin_f, FLA cos_f);

	FLT Sqr()							const	{ return x * x + y * y; }
	FLT Abs()							const;

	#ifdef _DEBUG
	void Dump()							const; //
//	{ afxDump << x << "\n"; afxDump << y << "\n"; }
	#endif

	FLT x, y;
	~Vect2D(){};
};

//-----------------------------------------------------------------------------------

static const Vect2D VECT2D_0(0, 0);		// zero vector

//-----------------------------------------------------------------------------------

static inline Vect2D operator +(Vect2D_ r1, Vect2D_ r2)	{ Vect2D rr(r1);	return rr += r2;	}
static inline Vect2D operator -(Vect2D_ r1, Vect2D_ r2)	{ Vect2D rr(r1);	return rr -= r2;	}
static inline Vect2D operator *(Vect2D_ r,  FLA k)		{ Vect2D rr(r);		return rr *= k;		}
static inline Vect2D operator *(FLA k,		Vect2D_ r)	{					return r * k;		}
static inline Vect2D operator /(Vect2D_ r,	FLA k)		{ Vect2D rr(r);		return rr /= k;		}

static inline FLT operator %(Vect2D_ r1, Vect2D_ r2)	// vector product
{ return r1.x * r2.y - r1.y * r2.x; }

static inline FLT	 operator *(Vect2D_ r1, Vect2D_ r2)	{ return r1.x * r2.x + r1.y * r2.y; }

//-----------------------------------------------------------------------------------

//static inline ostream& operator <<(ostream& stream, Vect2D_ r)	{ return stream << "(" << r.x << ", " << r.y << ")"; }
//static inline istream& operator >>(istream& stream, Vect2D& r)	{ return stream >> r.x >> r.y; }

//static inline ostream& operator <<(ostream& stream, Vect2D_ r)	{ return stream << "(" << r.x << ", " << r.y << ", " << r.z << ")"; }
//static inline istream& operator >>(istream& stream, Vect2D& r)	{ return stream >> r.x >> r.y >> r.z; }

//-----------------------------------------------------------------------------------

#endif //___Vect2D_h___

