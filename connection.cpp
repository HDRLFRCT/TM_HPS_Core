#include "connection.h"
#include <iostream>

connection::connection(int a, int b, int k, double eq)
{
    i = a;
	j = b;
	//koef = k;
	eq_st_sq = eq;
	//critical_state = cr;
	//ff = 0;
	origin = false;
	pressInBond = 0;
	pulled = false;
	pressure = false;
	cheek = false;
	numLiquidBond.clear();
	Q = 0;
	p_Q = 0;
}

connection& connection::operator=(const connection &at2d)
 {
     i = at2d.i;
	 j = at2d.j;
	 //koef = at2d.koef;
     eq_st_sq = at2d.eq_st_sq;
     //critical_state = at2d.critical_state;
     //ff = at2d.ff;
	 pulled = at2d.pulled;
	 cheek = at2d.cheek;
	 pressure = at2d.pressure;//связь под давлением
	 free = at2d.free;
	 pressInBond = at2d.pressInBond;
	 pulled = at2d.pulled;
	 Q = at2d.Q;
	 p_Q = at2d.p_Q;
	 return *this;
 }

connection::~connection()
{
	numLiquidBond.clear();
	//std::cout<<"dest connect \n";
}

connection::connection(int a, int b, double p)
{
    i = a;
	j = b;
	eq_st_sq = p;
	//critical_state = 1.01*p;//1.01*p;
	pulled = false;
	pressure = false;
	free = false;
	cheek = false;
	origin = false;
	pressInBond = 0;
	Q = 0;
	p_Q = 0;
	numLiquidBond.clear();
}
