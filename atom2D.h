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
    int number_x;//����� �� ������ ������� ����� ��� �
	int number_y;//����� �� ������� ������� ����� ��� �
    int number_xy;//���� �� ������� ������� ����� ��� x�
    int r_obr_x;//���������� � ����� ������� ����� ����� ��� � (1 - � ������������� ����������� �����; -1 � �������������)
	int r_obr_y;//���������� � ����� ������� ����� ����� ��� y (1 - � ������������� ����������� �����; -1 � �������������)
    int r_obr_xy;//���������� � ����� ������� �� ��������� ����� (22 - �� �� � �������������, 12 �� � ������������� � �������������, 21 �� � ������������� � �������������, 11 ��� �������������
	Vect2D r, v;	// Position and velocity
	Vect2D f;	    // Force
	std::vector<int> arrNumBond;//������ ������� ������
	LONG c;
	bool pressure;
	bool Spring;//����� ������� ������
	bool im;//������ ������� ��� ���
	bool obr_x;//���� �� ������ ������� ����� ��� �
	bool obr_y;//���� �� ������ ������� ����� ��� �
	bool obr_xy;//���� �� ������ ������� ����� ��� x�
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
  //  Atom3D(): v0(&v)  {}   // �����������
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
//	Atom3D *a0;// ������ �� ����
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

