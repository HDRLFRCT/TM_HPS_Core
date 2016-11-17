#include "atom2D.h"

Atom2D::Atom2D(float x, float y)
{

}

 Atom2D& Atom2D::operator=(const Atom2D &at2d)
 {
     r.x = at2d.r.x;
	 r.y = at2d.r.y;
     v.x = at2d.v.x;
	 v.x = at2d.v.x;
     //f = at2d.f;
	 im = at2d.im;//������ ������� ��� ���
	 obr_x = at2d.obr_x;//���� �� ������ ������� ����� ��� �
	 obr_y = at2d.obr_y;//���� �� ������ ������� ����� ��� �
	 obr_xy = at2d.obr_xy;//���� �� ������ ������� ����� ��� x�
	 number_x = at2d.number_x;//����� �� ������ ������� ����� ��� �
	 number_y = at2d.number_y;//����� �� ������� ������� ����� ��� �
	 number_xy = at2d.number_xy;//���� �� ������� ������� ����� ��� x�
	 r_obr_x = at2d.r_obr_x;//������ � ������ ������� ����� ��� �
	 r_obr_y = at2d.r_obr_y;//������ � ������ ������� ����� ��� y
	 r_obr_xy = at2d.r_obr_xy;
	 arrNumBond.clear();
	 for (int i = 0;i < at2d.arrNumBond.size();i++)
	 {
		 int pr = at2d.arrNumBond.at(i);
		 arrNumBond.push_back(pr);
	 }
	 bonds.clear();
	 for (int i = 0;i < at2d.bonds.size();i++)
	 {
		 connection pp = at2d.bonds.at(i);
		 bonds.push_back(pp);
	 }
	 pressure = at2d.pressure;
	 Spring = at2d.Spring;
     return *this;
 }
Atom2D::Atom2D(float x, float y,bool _im)
{
    r.x = x;
    r.y = y;
    v.x = 0;
    v.y = 0;
    f.x = 0;
    f.y = 0;
    obr_x = false;
    obr_y = false;
    obr_xy = false;
	//r_obr_x = 0;//VECT2D_0;
	//r_obr_xy = 0;//VECT2D_0;
	//r_obr_y = 0;//VECT2D_0;
    im = _im;
    bonds.clear();
	numSpringBond = 0;
	pressure = false;
	Spring = true;
}

Atom2D::Atom2D(double x, double y, bool _im)
{
	r.x = x;
	r.y = y;
	v.x = 0;
	v.y = 0;
	f.x = 0;
	f.y = 0;
	   obr_x = false;
	   obr_y = false;
	   obr_xy = false;
	r_obr_x = 0;//VECT2D_0;
	r_obr_xy = 0;//VECT2D_0;
	r_obr_y = 0;//VECT2D_0;
	im = _im;
	bonds.clear();
	numSpringBond = 0;
	pressure = false;
	Spring = true;
}


Atom2D::~Atom2D()
{

    bonds.clear();
	arrNumBond.clear();

}
