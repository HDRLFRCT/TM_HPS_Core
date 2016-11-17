#include "connection.h"
#include "atom2D.h"
#include <list>
#include "const.h"
#include "a3r.h"
#include "fstream"
#include "math.h"
using namespace std;
/*
Atom2D* particles = new Atom2D [N];//������ ������
std::list<connection> connect_list;

void SetConfiguration(Atom2D* particles);
void Step(Atom2D* particles);
////������� ������
void Run_Modeling()
{
    //������ ������
    SetConfiguration(particles);
    std::list<connection>::iterator it = connect_list.begin();
	for (;it!=connect_list.end();++it)
    {
        std::cout<<it->atom_A;
         std::cout<<it->atom_B;
          std::cout<<it->koef;
    }
    for (int step = 0;step< kol_iter; step++)//���� �������������� (����� ���� ������� �� while)
    {
        Step(particles);
        if (step % 10 == 0)
		{
			cout << step << "\t" << flush;
			char file_name[10];	sprintf(file_name, "a%04u.a3r", step);
			Save_A3R(file_name, particles, N, a0 / 2);
		}
    }
}
//������� �������� ��������� ������������
void SetConfiguration(Atom2D* particles)
{
    double dx = 1;
	std::cout<<"particals\n";
	for (int i = 0; i<N;i++)
	{
		Atom2D a(i*dx,0);
		particles[i] = a;
		std::cout<<particles[i].r.x<<"    "<<particles[i].r.y<<"\n";
	}
	for (int j = 0; j<(N-1);j++)
	{
		connection cc(j,j+1,10.0,0.1,1.5);
		//std::cout<<cc.atom_A<<"   "<<cc.atom_B<<"\n";
		connect_list.push_back(cc);
		//Connections(0).atom_A = 1
	}
	std::list<connection>::iterator it = connect_list.begin();
	for (;it!=connect_list.end();++it)
    {
        std::cout<<it->atom_A;
         std::cout<<it->atom_B;
          std::cout<<it->koef;
    }
}
//������� ������ ���� ��������������
void Step(Atom2D* particles)
{
    std::list<connection>::iterator iter = connect_list.begin();
    //������ ���
    for (iter = connect_list.begin();iter != connect_list.end();++iter)//������ �� ���� ������. ����� ������ ���� � �������� 0
    {
        Vect2D dr = particles[iter->atom_B].r - particles[iter->atom_A].r;
        std::cout<<"po svyazam";
        double length = sqrt(dr.Sqr()) - iter->eq_st; //����������-������ �������
        if (sqrt(dr.Sqr()) > iter->critical_state)//�������� �� ��������� �������
        {
            //������ ������� �����
        }
        iter->ff = length*iter->koef;//������ ����
        particles[iter->atom_B].f -= ((iter->ff)/sqrt(dr.Sqr()))*dr;
        particles[iter->atom_A].f += ((iter->ff)/sqrt(dr.Sqr()))*dr;
    }
    for (int i = 0; i<N;i++)//����������� �� ������ �������
    {
        particles[i].v+=particles[i].f*dt;
        particles[i].r+=particles[i].v*dt;
        particles[i].f=VECT2D_0;
        std::cout<<particles[i].r.x<<"    ";
    }
    std::cout<<"  \n  ";
}
*/
