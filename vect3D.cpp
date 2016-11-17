//-----------------------------------------------------------------------------------
//
//	3D Vector class
//
//	Anton M. Krivtsov
//
//	11.04.2001
//
//-----------------------------------------------------------------------------------

#include "vect3D.h"
#include "util.h"
#include <math.h>

//-----------------------------------------------------------------------------------

FLT Vect3D::Abs() const	{ return sqrt(Sqr()); }

//-----------------------------------------------------------------------------------

void Vect3D::Swap(BYTE mode)				//	 mode = 0,1,2;  default mode = 1
{
	if (!mode) return;
	FLT zz = z;
	if (mode == 1)	{ z = x; x = y; y = zz; return; } 
	if (mode == 2)	{ z = y; y = x; x = zz;	return; }
}

Vect3D Vect3D::GetSwapped(BYTE mode) const	//	 mode = 0,1,2;  default mode = 1
{
	Vect3D r = *this; r.Swap(mode); return r;
}

//-----------------------------------------------------------------------------------

void Vect3D::SetRand(FLA x_min, FLA y_min, FLA z_min, FLA x_max, FLA y_max, FLA z_max) 
{
	x = rand(x_min, x_max);
	y = rand(y_min, y_max);
	z = rand(z_min, z_max);
}

//-----------------------------------------------------------------------------------

void Vect3D::SetRand(FLA r_max) 
{
	FLT rr_max = r_max * r_max;
	do SetRand(-r_max, -r_max, -r_max, r_max, r_max, r_max);
	while (Sqr() > rr_max);
}

//-----------------------------------------------------------------------------------

void Vect3D::RotateZ(FLA angle) { Rotate(x, y, sin(angle), cos(angle)); }
void Vect3D::RotateX(FLA angle) { Rotate(y, z, sin(angle), cos(angle)); }
void Vect3D::RotateY(FLA angle) { Rotate(z, x, sin(angle), cos(angle)); }

//-----------------------------------------------------------------------------------
//	private
//-----------------------------------------------------------------------------------

void Vect3D::Rotate(FLT& x, FLT& y, FLA sin_f, FLA cos_f)
{ 
	FLT x1;
	x1 =  x * cos_f - y * sin_f; 
	y  =  x * sin_f + y * cos_f; 
	x  =  x1;
}

//-----------------------------------------------------------------------------------

