#include "cell.h"

Cell::Cell(int N)
{
	maxPart = N;
	numberPart = 0;
	arrNumPart.clear();
	//arrNumPart = new int [maxPart];
}
Cell::Cell()
{

}
void Cell::Add(int numPart)
{
	arrNumPart.push_back(numPart);
	numberPart++;
}
void Cell::Delete()
{
	//delete[]arrNumPart;
	arrNumPart.clear();
	numberPart = 0;
}
Cell& Cell::operator=(const Cell &at2d)
 {
     maxPart = at2d.maxPart;
	 numberPart = at2d.numberPart;
	 //arrNumPart = new int[maxPart];
	 arrNumPart.clear();
	 for (int i = 0;i < at2d.arrNumPart.size();i++)
	 {
		 arrNumPart.push_back(at2d.arrNumPart.at(i));
	 }
     return *this;
 }
Cell::~Cell()
{
	arrNumPart.clear();
}
