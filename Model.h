#ifndef MODEL_H_INCLUDED
#define MODEL_H_INCLUDED
#include <list>
#include "connection.h"
#include "atom2D.h"
#include <vector>
#include "const.h"
#include "a3r.h"

#pragma once

class Model
{
public:
	//int N; //���������� ������ ������
	std::list<connection> Connections;//������ ������ ������
	//std::vector<Atom2D> Particals;//������ ������ � ������
	Atom2D* Particles;//������ ������ � ������

	Model(void);
	~Model(void);
	void SetConfiguration();//������� �������� ������������ ������.
	void Step();
	void Save_configuration(const char* file_name, double r);
	void Run_Modeling(int kol);
};


#endif // MODEL_H_INCLUDED
