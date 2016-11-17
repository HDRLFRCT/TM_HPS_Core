#ifndef CONST_H_INCLUDED
#define CONST_H_INCLUDED

// онстанты модели
//const int N = 10;
const int number_iteraction = 10000000;//10000;//20000000;// оличество итераций в расчете
const int attempts = 40000000; //600000;//20000000; //100000;
const double pi = 3.14159265358979;
const float dt = float(0.1);
const double dt_liq = 0.05;//инкремент времени дл€ жидкости
const double a0 = 0.1;
const double e0 = 1.008*1.008;// оэффициент дл€ критического рассто€ни€ (квадрат)
const double k_bond = 0.5;//0.5
const float visc_k = float(0.3);
const double visc_liq = 0.5*1e1;//в€зкость жидкости
const int kol_iter = 50;//20;
const double E = 1;
const double gamma = 1;
const int numInitialNodes = 8;
const double koef_dq = 0.1*1e3;// 1e4;
const double koef_pp = 1*1e0;//1e-2;//1e-4;//0.5*1e-1;//1e-1;
const double koef_q = 1e7;//1e7;
//const double d_t_q = 1e1;//множитель инкримента времени по жидкости
const double partInitialFluidFracture = 0.8;//часть котора€ наполнена жидкостью в начальный момент времени
const double q_0 = 0.5*1e-2;//поток загон€емый в трещину за 1 итериацию
const double part_fixed = 0.01;//часть от длины которую закрепл€ем
const int save_step =25;//100;  //Ўаг сохранени€ файлоы
const double part_bond_liquid = 0.3;//0.15;//0.1;
const int typeBC = 1;//“ип граничных условий (0 - периодические, 1 - фиксированные)
//const double press = 0.001;
#endif // CONST_H_INCLUDED