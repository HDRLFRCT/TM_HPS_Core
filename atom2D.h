//-----------------------------------------------------------------------------------
//
//	3D Atom Class and utilities
//
//  Anton M. Krivtsov
//  Vadim A. Tsaplin
//
//  28.07.2002
//
//-----------------------------------------------------------------------------------

#include "vect2D.h"		// 3D Vector class
//#include "Tens3D.h"		// 3D tensor class
//#include "box2D.h"		// 3D Box class
#include "colors.h"
#include <vector>
#include "connection.h"

//-----------------------------------------------------------------------------------

#ifndef ___Atom3D_h___
#define ___Atom3D_h___

//-----------------------------------------------------------------------------------





class Atom2D
{
public:
    Atom2D(): f(VECT2D_0) {}
    Atom2D(float x, float y);
    Atom2D(float x, float y,bool _im);
	Atom2D(double x, double y, bool _im);
    Atom2D& operator=(const Atom2D &at2d );
	int numSpringBond;
    int number_x;//номер на мнимую частицу вдоль оси х
	int number_y;//номер на мнимаую частицу вдоль оси у
    int number_xy;//номе на мнимаую частицу вдоль оси xу
    int r_obr_x;//показывает в какую сторону образ вдоль оси х (1 - в положительном направлении образ; -1 в отрицательном)
	int r_obr_y;//показывает в какую сторону образ вдоль оси y (1 - в положительном направлении образ; -1 в отрицательном)
    int r_obr_xy;//показывает в какую сторону по диагонали образ (22 - по хи у положиьельный, 12 по х отрицательный у положительный, 21 по у отрицательный х положительный, 11 оба отрицательных
	Vect2D r, v;	// Position and velocity
	Vect2D f;	    // Force
	std::vector<int> arrNumBond;//массив номеров связей
	LONG c;
	bool pressure;
	bool Spring;//имеет упругих связей
	bool im;//мнимая частица или нет
	bool obr_x;//есть ли мнимая частица вдоль оси х
	bool obr_y;//есть ли мнимая частица вдоль оси у
	bool obr_xy;//есть ли мнимая частица вдоль оси xу
	int numBond;
	//LONG c;			// Reserved to mark atoms. Can be used as a color.
	std::vector <connection> bonds;
	~Atom2D();
};










/*


class Atom3D
{
public:
    Vect3D r, v, v_old;	// Position and velocity
	Vect3D v_flow;	FLT rho; FLT P;
	FLT xplus;
	FLT xminus;
	BOOL fluid;
	FLT R, m;
//	Vect3D nn[4];
	//FLT con;
	//BOOL stop;
	//Vect3D r_stop;
//	Atom3D* nei_T[NT];
//	short int nei_n[4];
 	LONG real;
//	Atom3D* miror;
//	BOOL bound;
//	COLORREF c;
	//LONG edge;
//	BOOL free;
//	BOOL reflect;
  //  Atom3D(): v0(&v)  {}   // конструктор
//	Atom3D(): omega0(&omega)  {}
 //   Atom3D& operator=(Atom3D& at)           { r = at.r; v0 = &at.v; return at; }

	Atom3D& operator=(Atom3D& at)
	{
		r = at.r;  v = at.v;
		v_old = at.v_old;
		R = at.R;
		fluid = at.fluid;
		P = at.P;
		xplus = at.xplus;
		xminus = at.xminus;
		rho = at.rho;
		m = at.m;
		real = at.real;
		return at;
	}
//  Atom3D& operator+=(const Vect3D& vect)  { *v0 += vect;  return *this; }
//  Atom3D& operator-=(const Vect3D& vect)  { *v0 -= vect; return *this; }\
//private:
//	Atom3D *a0;// ссылка на атом
//	Vect3D *v0;     // Used at boundary copied atoms as pointer to original velocities
};

class Atom3D_0
{
public:
    Vect3D r, v, v_old, omega;	// Position and velocity
	FLT m, R;
	BOOL fluid;
	FLT con;
	FLT xplus;
	FLT xminus;
	BOOL stop;
	Vect3D r_stop;
	Vect3D r0; //Vect3D phi0;
	Vect3D nn[4];
    BOOL n[4];
	LONG real;
	LONG edge;
};






class Atom3D_0
{
public:
    Vect3D r, v, phi, omega;	// Position and velocity

	Vect3D r0; Vect3D phi0;

	BOOL real;


    Atom3D(): v0(&v)  {}
    Atom3D& operator=(Atom3D& at)           { r = at.r; v0 = &at.v; return at; }
    Atom3D& operator+=(const Vect3D& vect)  { *v0 += vect; return *this; }
    Atom3D& operator-=(const Vect3D& vect)  { *v0 -= vect; return *this; }

private:
	Vect3D *v0;     // Used at boundary copied atoms as pointer to original velocities
};*/




//-----------------------------------------------------------------------------------
/*
Box3D GetBoundBox		(Atom3D* as, LONG n);
Box3D GetBoundBox_0     (Atom3D_0* as_0, LONG n);

FLT   GetDeviation		(Atom3D* as, LONG n);

void AddRandVelocity    (Atom3D* as, LONG n, FLT v_max, LONG seed = 1);
void CentrateVelocities (Atom3D* as, LONG n);
Vect3D GetCenterVelocity(Atom3D* as, LONG n);
Vect3D GetCenter		(Atom3D* as, LONG n);
*/
//-----------------------------------------------------------------------------------

#endif // ___Atom3D_h___

//-----------------------------------------------------------------------------------

