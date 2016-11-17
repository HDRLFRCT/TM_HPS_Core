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
	//int N; //количество частиц модели
	std::list<connection> Connections;//список связей модели
	//std::vector<Atom2D> Particals;//массив частиц в моделе
	Atom2D* Particles;//массив частиц в моделе

	Model(void);
	~Model(void);
	void SetConfiguration();//функция создания конфигурации модели.
	void Step();
	void Save_configuration(const char* file_name, double r);
	void Run_Modeling(int kol);
};


#endif // MODEL_H_INCLUDED
