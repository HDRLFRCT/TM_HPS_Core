#include "Model.h"
#include <math.h>
#include <fstream>

Model::Model(void)
{
	//Particles = new Atom2D[N];
}


Model::~Model(void)
{
	delete []Particles;
	Connections.clear();
}

void Model::Save_configuration(const char* file_name, double r)
{
	//Save_A3R(file_name,Particles,N,r);
}

void Model::SetConfiguration()
{
	//Здесь будет код от дениса//
	//----------------ВРЕМЕННАЯ ЗАГЛУШКА ЦЕПОЧКА----------
	//for (int i = 0; i<N;i++)//интегрируем по каждой частице скороксти
    //{
        //Particles[i].obr = false;
        // cout<<" velocity na  "<<i<<"chastice"<<Particles[i].v.x<<"\n";
    //}
	/*double dx = 1;
	std::cout<<"particals\n";
	for (int i = 0; i<N;i++)
	{
		Atom2D a(i*dx,0);
		if (i ==1){a.r.x+=0.1;}
		//if (i ==2){a.r.x+=0.1;}
		Particles[i] = a;
		std::cout<<Particles[i].r.x<<"    "<<Particles[i].r.y<<"\n";
	}
	std::cout<<"=============================================\n";
	for (int j = 0; j<(N-1);j++)
	{

		connection cc(&Particles[j],&Particles[j+1],10.0,0.5,2.0);
		std::cout<<j+1<<"csvayz   first atom"<<cc.atom_A->r.x<<"  csvayz   second atom"<<cc.atom_B->r.x<<"\n";
		Connections.push_back(cc);
		//Connections(0).atom_A = 1
	}
	/*Particles[0].obraz = &Particles[N-1];
	Particles[0].obr = true;
	std::cout<<"obr   first atom"<<Particles[0].obraz->r.x<<"\n";
	Particles[N-1].obraz = &Particles[0];
	Particles[N-1].obr = true;
    std::cout<<"obr   posl atom"<<Particles[N-1].obraz->r.x<<"\n";
    */

	//============================================
}
void Model::Step()
{
    /*std::list<connection>::iterator iter;
    int jjj = 1;
    for (iter = Connections.begin();iter != Connections.end();++iter)//проход по всем связям. Перед циклом силы в частицах 0
    {
        Vect2D dr = (iter->atom_B->r) - (iter->atom_A->r);
        double length = sqrt(dr.Sqr()) - iter->eq_st; //растяжение-сжатие пружины
        if (dr.Sqr() > pow(iter->critical_state,2))//проверка на состоянии пружины
        {
            //случай разрыва связи
        }
        double ff = length*iter->koef;//модуль силы

        iter->atom_B->f -= (ff/sqrt(dr.Sqr()))*dr;
        iter->atom_A->f += (ff/sqrt(dr.Sqr()))*dr;

    }

    for (int i = 0; i<N;i++)//интегрируем по каждой частице скороксти
    {

        Particles[i].v+=Particles[i].f*dt;
/*        if (Particles[i].obr)
        {
            cout<<"chasti "<<i<<"  \n";
            Vect2D fff = Particles[i].obraz->f;
            Particles[i].v+=fff*dt;
        }*/
     /*    cout<<" velocity na  "<<i<<"chastice"<<Particles[i].v.x<<"\n";
    }
    for (int i = 0; i<N;i++)//интегрируем перемещения каждой частице
    {
        Particles[i].r+=Particles[i].v*dt;
        cout<<" posotia na  "<<i<<"chastice"<<Particles[i].r.x<<"\n";
        Particles[i].f=VECT2D_0;
    }
    std::cout<<"  end step\n  ";*/
}

void Model::Run_Modeling(int k)
{
    /*SetConfiguration();
    for (k = 0;k< kol_iter; k++)//цикл интегрирования (может быть изменен на while)
    {
        Step();

        if (k % 10 == 0)
		{
			//cout << k << "\n" << flush;
			char file_name[10];	sprintf(file_name, "a%04u.a3r", k);
			Save_A3R(file_name, Particles, N, a0 / 2);
		}
    }*/
}
