#include <vector>
#include <stdlib.h>
#include "atom2D.h"
#include "connection.h"
#include "const.h"
#include "cell.h"
#include "a3r.h"
#include <fstream>
#include <math.h>
#include <cmath>
#include <random>
#include <ctime>
#include <algorithm>
#include <iostream>
#include "timer.h"
#include <fstream>
#include <string>
//#include <vld.h>

#define axis_x 1
#define axis_y 2
#define volume_ext 3

using namespace std;
void Modeling();//������� ������� �������������
/****================������� �������� ��������� ������������================***/
void add_new_ball(vector<Atom2D>* P, double x, double y, double w, double h, double bond,
	double r2, double periodic_border);//���������� ������ ����(��������������-����������)

void add_new_ball1(vector<Atom2D>* P, double x, double y, Cell**grid, int numCellsX, int numCellsY,
	double w, double h, double bond, double r2, double periodic_border);//���������� ������ ����(�����������)

void add_new_ball_with_periodic(vector<Atom2D>* P, double x, double y, Cell**grid, int numCellsX, int numCellsY,
	double w, double h, double bond, double r2, double periodic_border);//����������� ���� ��� ������������� ��������� ��������

void start_new_system(vector<Atom2D> *P, vector<connection> *B, double w, double h,
	double bond_len, double r, double periodic_border, int max_bonds, int attempts);//�������� ����� ������������ ������� (���������-��������������)

void start_new_system1(vector<Atom2D> *P, vector<connection> *B, Cell**grid, Cell** gridBonds,
	int numCellsX, int numCellsY, double w, double h, double bond_len, double r, double periodic_border,
	int max_bonds, int attempts);//�������� ����� ������������ ������� (�����������)

void start_new_system_treangle(vector<Atom2D> *P, vector<connection> *B, Cell**grid, Cell**gridBonds,
	int numCellsX, int numCellsY, double w, double h, double bond_len, double r,
	double periodic_border, int max_bonds, int attempts);//�������� ����� ������������ ������� c ���������� ��������

void new_free_fracture(vector<Atom2D>*P, vector<connection>*B, double x, double y, double alf_rad, double len);//�������� ��������� �������
void add_rand_free_fractures(vector<Atom2D>*P, vector<connection>*B, int n, double l, double w, double h);//������� n ������������ ������ ����� l
void CreateLiquidSourse(vector<Atom2D>*as, vector<connection> *l, double x, double y);//�������� ��������� ��������

void CreateFracture(vector<Atom2D>*as, vector<connection> *l, vector<int> *PressBond, Cell **gridBonds,
	double x, double y, double alfa, double lenght, double periodic_border, double bond_len, int numCellsX, int numCellsY);//�������� ������� 

void ConnectLiquidNodes(vector<Atom2D>*P, vector<connection> *B, Cell** gridBonds, int numgridX, 
	int numgridY, int numCellX, int numCellsY, double bond_len);//��������� ���������� �������� � ������ ������



/***====================�����������====================***/
void Step(vector<Atom2D>* as, vector<connection>*l, double w, double h);//����������

void Step1(vector<Atom2D>* as, vector<connection>*l, bool**gridPressBond, int numCellsX, int numCellsY, 
	double w, double h, double bond, double periodic_border, vector<int>*PressBond);//����������

void StepWithViscosity(vector<Atom2D>* as, vector<connection>* l, bool ** gridPressBond, int numCellsX, int numCellsY, 
	double w, double h, double bond, double periodic_border, vector<int>* PressBond, double x_c0, double x_c1, double lenBetweenNodes, 
	int numNodes, double** widthFracture, double** Press, double** Visc, double* v_cur);//���������� c �������� ����������� ������ ��������

void StepWithViscosity(vector<Atom2D>* as, vector<connection>* l, double w, double h, double bond, double periodic_border, 
	vector<int>* PressBond, double x_c0, double x_c1, double lenBetweenNodes, int numNodes, vector<double>* widthFracture, 
	vector<double>* Press, vector<double>* Visc, double* v_cur);//���������� c �������� ����������� ������ ��������

void StepManyFrackPress(vector<Atom2D>* as, vector<connection>*l, bool**gridPressBond, int numCellsX, int numCellsY, double w, double h,
	double bond, double periodic_border, vector<int>*PressBond, vector<int>*FrackPress, double * masPress, int numFrack);//����������

void StepRungeKutta(vector<double>* w, vector<double>* p, vector<double>* v, int numElements, 
	double deltaX, double x_c0, double x_c1, double h, double* v_cur);//��� �������������� ������� ����� �����

void StepWithSimpleModle(vector<Atom2D>*as, vector<connection> *l, vector<int> *PressBond, double w, double h, 
	int numNodes, int numNodesDisplacement, double lenBetweenNodes, vector<double>*Press, vector<double> *widht, 
	vector<double>*q, vector<double>*u, vector<double>*pr_u, vector<int>*num_u);//��� ������� ������

void StepWithLiquidPipe(vector<Atom2D>*as, vector<connection> *l, vector<int> *PressBond, Cell**gridBonds, double w, 
	double h, double periodic_border, double bond_len, int numCellsX, int numCellsY);// int numNodes, int numNodesDisplacement, double lenBetweenNodes, vector<double>*Press, vector<double> *widht, vector<double>*q, vector<double>*u, vector<double>*pr_u, vector<int>*num_u);

void RecalculateVelosityPositionPeriodic(vector<Atom2D>*as, double w, double h);//�������� ��������� � ��������� ������������� ��������� �������
void RecalculateVelosityPositionFixed(vector<Atom2D>*as, double w, double h);//�������� ��������� � ��������� ��������������� ����

/***====================��������������� �������====================***/
void Clear_all(vector<Atom2D>* as);
void SavePressBondDisplacement(vector<Atom2D>*as, vector<connection>*l, vector<int>*PressBond);
bool intersection(double A_x,double A_y,double B_x, double B_y,double C_x, double C_y,double D_x, double D_y);//�������� ������������ �� ������� �� � CD

void SetConf(vector<Atom2D>* as, vector<connection>* l, double w, double h, double bond_len);
void AproximationPress(vector<double>* Press, double deltaX,int numNodes);
double GetDeterminant(double** matrix);

void output_to_file(vector<Atom2D> *P, vector<connection> *B, double r, double w, double h,int num);//������ � ����
void extention_system(vector<Atom2D>*P, double e,int ax,double &w,double &h);//���������� �������, ax - �������� ����������, 
void dublicate_config(vector<Atom2D>*Initial,vector<Atom2D>*Final);//����������� � ������������ vector ������ ��������� ������������
void Save_StressField(vector<Atom2D>*P, vector<connection> *B);

void output_pressured_area(vector<Atom2D>* P, vector<connection> *B, int j, int number_iter, double w, double h, double r);// ������ � ���� ������� �������
double FindFractureSquare(vector<Atom2D>* P, vector<connection> *B, vector<int> *PressBond, double bond_len);
double LenConnection(Atom2D A, Atom2D B);//����� �����
double LenBetweenConnection(Atom2D A1, Atom2D A2, Atom2D B1, Atom2D B2);
void CreateAvtoModelSolution(vector<Atom2D>* P, vector<connection> *B, vector<int> *PressBond, double LenFracture, double*w_avt,  int numNodes, double h, double w);
void CreateLiquidSourse(vector<Atom2D>*as,vector<connection> *l,double x,double y);//�������� ��������� ��������
void CreateFracture(vector<Atom2D>*as,vector<connection> *l, vector<int> *PressBond, Cell **gridBonds,double x, double y,double alfa,double lenght, double periodic_border, double bond_len, int numCellsX, int numCellsY);//�������� ������� 
void ConnectLiquidNodes(vector<Atom2D>*P, vector<connection> *B, Cell** gridBonds, int numgridX, int numgridY, int numCellX, int numCellsY, double bond_len);//��������� ���������� �������� � ������ ������
void SaveConfiguration(vector<Atom2D>*P, vector<connection>*B);//������ � ���� ������������ ������:��������� ������, ����� � ������ ������. ��� ��������, �������� � �������. ������ �����.
void LoadConfiguration(vector<Atom2D>*P, vector<connection>*B,Cell** gridBonds,int numCellsX, int numCellsY, double bond_len,double periodic_border);//�������� ������ �� �����
void SaveLiquidConnection(vector<Atom2D>*P,vector<connection>*B);
void FindPressInFractures(double V_liq, vector<int>*FrackPress, double*masPress, int numFracks,double bond);