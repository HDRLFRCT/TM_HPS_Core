#ifndef CONNECTION_H_INCLUDED
#define CONNECTION_H_INCLUDED
//#include "atom2D.h"
#include <vector>

class connection
{
    public:
        //Atom2D* atom_A;//������ ������� �����, ��� ����� � ������ ������
        //Atom2D* atom_B;//������ ������� �����, ��� ����� � ������ ������
        
        //double koef;//����������
        double eq_st_sq; //���������� ����� � ��������
        //double critical_state; //���������� �������
		int i;//����� ������� �������� � ������� ������
		int j;//����� ������� �������� � ������� ������
        //double ff;//���� ����������� �� ������� �� �����
		double pressInBond;
		bool origin;
        connection(int a, int b, int k, double eq);
        connection(int a, int b, double p);
		connection& operator=(const connection &at2d );
        virtual ~connection();
        bool pulled;
		bool cheek;//���������� ��� ������� �� ������
		bool pressure;//����� ��� ���������
		bool free;//����� �������(��������� �����)
		double Q;//����� �������� � ����
		float p_Q;//����� �� ���������� �������� ��� �������� ������� �������
		std::vector<int>numLiquidBond;
    protected:
    private:
};


#endif // CONNECTION_H_INCLUDED
