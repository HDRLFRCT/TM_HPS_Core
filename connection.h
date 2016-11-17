#ifndef CONNECTION_H_INCLUDED
#define CONNECTION_H_INCLUDED
//#include "atom2D.h"
#include <vector>

class connection
{
    public:
        //Atom2D* atom_A;//первый элемент связи, его номер в списке частиц
        //Atom2D* atom_B;//второй элемент связи, его номер в списке частиц
        
        //double koef;//коэфициент
        double eq_st_sq; //расстояние покоя в квадрате
        //double critical_state; //расстояние разрыва
		int i;//номер первого элемента в векторе частиц
		int j;//номер второго элемента в векторе частиц
        //double ff;//сила действующая на частицы от связи
		double pressInBond;
		bool origin;
        connection(int a, int b, int k, double eq);
        connection(int a, int b, double p);
		connection& operator=(const connection &at2d );
        virtual ~connection();
        bool pulled;
		bool cheek;//переменная для прохода по частям
		bool pressure;//связь под давлением
		bool free;//связь порвана(свободная связь)
		double Q;//объем жидкости в узле
		float p_Q;//Объем на предыдущей итерации для подсчета размера каналов
		std::vector<int>numLiquidBond;
    protected:
    private:
};


#endif // CONNECTION_H_INCLUDED
