#include "modeling.h"

#define USI unsigned short int
float press = 0;//1e-4;
double num_errase = 0;
int numPressBond = 0;
int stepp = 0;
double max_rx = 0;
double min_rx = 0;
double max_ry = 0;
double min_ry = 0;
bool nodes = false;

void Modeling()//Главная функция моделирования
{
	vector<Atom2D> Particles;
    std::vector<connection> connect_list;
    double w = 8.0;                        // ширина поля для генерации частиц
    double h = 8.0;                           // высота поля для генерации частиц
    double r =  0.025;//0.025;                       // минимальное расстояние между частицами
    double bond_len = 2*r;//2*r;               // максимальная длина связи
    double periodic_border = 1.0 * bond_len;// ширина границы, на которой добавляются ghost-частицы
	double part_frack = 0.02;//отношение длины трещины к ширине
	float				 alfa;
	float           extr_x;
	float			 extr_y;
	float			dist_bwn_frcks;
	int				num_fracks;
	int final_count_break_bond;
	float step_press;
	float w_t = 0;
	float h_t = 0;
	float init_frac = 1;
	float q_in_fr = 1e-3;
	FILE *file;
	file = fopen("data.dat", "r");
	string pr;
	char str[80];

	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%g", &(alfa));
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%g", &(extr_x));
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%g", &(extr_y));
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%g", &(dist_bwn_frcks));
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%u", &(num_fracks));
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%u", &(final_count_break_bond));
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%e", &(step_press));
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%g", &(w_t));
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%g", &(h_t));
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%g", &(init_frac));
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%e", &(q_in_fr));
	fclose(file);
	alfa *= pi / 2;
	if ((w_t > 0) && (w_t < 150)) {
		w = w_t;
	}
	if ((h_t > 0) && (h_t < 150)) {
		h = h_t;
	}
	if ((init_frac > 0) && (init_frac < 150)) {
		part_frack = init_frac;
	}
	if ((q_in_fr < 0) || (q_in_fr >100)) {
		q_in_fr = 0;//1e-3;
	}
	/******************************/
	USI n_x = int(bond_len/r)+1;//максимальное количество частиц в ячейке вдоль х при максимаьно плотной упаковке
	USI n_y = int(bond_len/(sqrt(3)*0.5*r))+1;
	USI n_max = n_x*n_y;//max количество возможное частиц в ячейке
	USI numCellsX = USI((w+2*periodic_border)/bond_len)+1;//количество ячеек вдоль х
	USI numCellsY = USI((h+2*periodic_border)/bond_len)+1;//количество ячеек вдоль у
	USI numCellsXLiq = USI(w /bond_len/part_bond_liquid) + 2;
	USI numCellsYLiq = USI(h / bond_len/part_bond_liquid) + 2;
	std::vector<int>FracturesBond;//Каждой связи соответсвует пара чисел-номер связи+номер трещины
	double V_liq = 0;//Объем загнанный в каждую трещину 
	cout<<numCellsX<<" "<<numCellsY;
	Cell **gridBonds = new Cell*[numCellsX];
	for (int g1 = 0;g1 < numCellsX;g1++)
	{
		gridBonds[g1] = new Cell[numCellsY];
	}
	for (USI x = 0;x<numCellsX;x++)
		for (USI y = 0;y<numCellsY;y++)
		{
			gridBonds[x][y] = Cell(n_max*n_max);
		}
	Cell **grid = new Cell*[numCellsX];//сетка области, в ячейки которой храняться номера частиц в соответсвующей ячейке. Обращение [X][Y]
	for (USI i = 0;i<numCellsX;i++)
	{
		grid[i] = new Cell[numCellsY];
	}
	for (USI x = 0;x<numCellsX;x++)
		for(USI y = 0;y<numCellsY;y++)
		{
			grid[x][y] =  Cell(n_max);
		}
	
	/******************СОЗДАЮ ЗАПИСЬ ИЗ ВНЕШНЕГО ФАЙЛА*********************/
	
	//output_to_file(&Particles, &connect_list, r, w, h,1);
	//vector<Atom2D> InitialConfiguration;
	//start_new_system_treangle(&Particles, &connect_list, grid, gridBonds, numCellsX, numCellsY, w, h, bond_len, r, periodic_border, 0, attempts);
	start_new_system1(&Particles, &connect_list, grid, gridBonds, numCellsX, numCellsY, w, h, bond_len, r, periodic_border, 0, attempts);
	//start_new_system(&Particles, &connect_list, w,h, bond_len, r, periodic_border, max_bonds, attempts);
	//start_new_system(&Particles, &connect_list, w,h, bond_len, r, periodic_border, max_bonds, attempts);



	USI k = 0;
	bool ff = false;
	vector<int> PressBond;//номера связей с давлением
	double e = 0.02;//задают отступ от центра где удаляю связи
	double e_tr = 0.000005;//параметр растяжения
	double lenInitialFracture = part_frack * h;//50;//длина начальнйо трещины
	//double x;
	//double y;
	for (USI ii = 0; ii<numCellsX;ii++)
	{
		delete[] grid[ii];
	}
	delete [] grid;
	grid = 0;
	
	int kol_bond = 0;
	//LoadConfiguration(&Particles, &connect_list,gridBonds,numCellsX,numCellsY,bond_len,periodic_border);
	//SaveConfiguration(&Particles, &connect_list);
	
   	bool **gridPressBond = new bool*[numCellsX];//сетка области, в ячейки которой храняться номера частиц в соответсвующей ячейке. Обращение [X][Y]
	for (USI i = 0;i<numCellsX;i++)
	{
		gridPressBond[i] = new bool [numCellsY];
	}
	for (USI x = 0;x<numCellsX;x++)
		for (USI y = 0;y<numCellsY;y++)
		{
			gridPressBond[x][y] = false;
		}
	/**Создаем необходимое количество трещин**/
	int* numArr = new int[num_fracks];
	if (num_fracks % 2 == 1) {
		
		CreateFracture(&Particles, &connect_list, &PressBond, gridBonds, w/2, h / 2, 
			alfa, lenInitialFracture, periodic_border, bond_len, numCellsX, numCellsY);
		numArr[0] = PressBond.size();
		for (int i = 0;i < 0.5*(num_fracks-1);i++) {

			double dx = w*0.5 + (i+1)*lenInitialFracture*dist_bwn_frcks;
			CreateFracture(&Particles, &connect_list, &PressBond, gridBonds, dx, 2.7*(i+1)*lenInitialFracture+(h / 2),
				alfa, lenInitialFracture, periodic_border, bond_len, numCellsX, numCellsY);
			numArr[1 + 2 * i] = PressBond.size();
			dx = w*0.5 - (i + 1)*lenInitialFracture*dist_bwn_frcks;
			CreateFracture(&Particles, &connect_list, &PressBond, gridBonds, dx, -2.7*(i + 1)*lenInitialFracture + h / 2,
				alfa, lenInitialFracture, periodic_border, bond_len, numCellsX, numCellsY);
			numArr[1 + 2 * i+1] = PressBond.size();
		}
	}
	else {
		for (int i = 0;i < 0.5*num_fracks;i++) {
			double dx = w*0.5 - 0.5*lenInitialFracture + (i + 1)*lenInitialFracture*dist_bwn_frcks;
			CreateFracture(&Particles, &connect_list, &PressBond, gridBonds, dx, h / 2,	
				alfa, lenInitialFracture, periodic_border, bond_len, numCellsX, numCellsY);
			if (i == 0) {
				numArr[2 * i] = PressBond.size();
			}
			else {
				numArr[2 * i] = PressBond.size();
			}
			dx = w*0.5 + 0.5*lenInitialFracture - (i + 1)*lenInitialFracture*dist_bwn_frcks;
			CreateFracture(&Particles, &connect_list, &PressBond, gridBonds, dx, h / 2, 
				alfa, lenInitialFracture, periodic_border, bond_len, numCellsX, numCellsY);
			numArr[2 * i + 1] = PressBond.size();
		}
	}
	int cur_int =  0;
	for (int iii = 0;iii < num_fracks;iii++) {
		cout << numArr[iii] << "\n";
	}
	for (int ii = 0;ii < PressBond.size();ii++) {
		if (ii <= numArr[cur_int]) {
			FracturesBond.push_back(PressBond.at(ii));
			FracturesBond.push_back(cur_int);
		}
		else {
			cur_int++;
			FracturesBond.push_back(PressBond.at(ii));
			FracturesBond.push_back(cur_int);
		}
	}
	//cout << cur_int << "\n";
	delete[]numArr;
	
	//CreateFracture(&Particles,&connect_list,&PressBond,gridBonds,w/2,h/2,alfa,lenInitialFracture,periodic_border,bond_len,numCellsX,numCellsY);
	
	/*Задаю источник жидоктси*/
	//CreateLiquidSourse(&Particles,&connect_list,w*0.5,h*0.5);
	//SaveLiquidConnection(&Particles,&connect_list);
	//new_free_fracture(&Particles,&connect_list,w*0.3,h*0.85,pi*0.5,0.1);
	//new_free_fracture(&Particles,&connect_list,w*0.5,h*0.5,-pi*0.25,0.8);
	//add_rand_free_fractures(&Particles,&connect_list,20,0.4,w,h);
	double* pressInFr = new double[num_fracks];
	for (int ii = 0;ii < num_fracks;ii++) {
		pressInFr[ii] = 0;
	}
	string nameLog = "Log.txt";
	ofstream log_txt;
	log_txt.open(nameLog);
	log_txt << "Number particles " << Particles.size() << "\n";
	log_txt << "Number bonds " << connect_list.size() << "\n";
	log_txt.close();
	/*ОСНОВНОЙ ЦИКЛ МОДЕЛИРОВАНИЕ */
	for (int j = 0;j<number_iteraction;j++)//увеличение давления
	{
			for (int number_iter = 0; number_iter<kol_iter+1;number_iter++ )
			{
				if ((j % save_step== 0) && (number_iter == 0))
					{
						cout << "Iteration " << stepp << "\n";
						char file_name[20];	sprintf(file_name, "a%08u%04u.a3r",j+1, number_iter);
						//Save_A3R(file_name, &Particles, Particles.size(), r / 2,w,h);
						output_pressured_area(&Particles,&connect_list,j,number_iter,w,h,r);
						//Save_StressField(&Particles, &connect_list);
						
					}
				
				/*Интегратор выбирается в зависимости от типа задачи (с жидкостью или без)*/

				//StepWithSimpleModle(&Particles, &connect_list, &PressBond, w, h, numElement,numElementsDisplasment, deltaX, &PressDistrib, &widthFracture, &q, &u,&pr_u,&num_u);
 				//StepWithLiquidPipe(&Particles, &connect_list, &PressBond,gridBonds, w, h,periodic_border,bond_len,numCellsX,numCellsY);
				//StepWithViscosity(&Particles, &connect_list, w, h, bond_len, periodic_border, &PressBond, x_c0, x_c1, deltaX, numElement, &widthFracture, &PressDistrib, &velosity,v_cur);
				//Step1(&Particles, &connect_list, gridPressBond, numCellsX, numCellsY, w, h, bond_len, periodic_border, &PressBond);
				StepManyFrackPress(&Particles, &connect_list, gridPressBond, numCellsX, numCellsY, w, h, bond_len, periodic_border, &PressBond, 
					&FracturesBond, pressInFr, num_fracks);
				
			}
			//if (num_errase<final_count_break_bond) { press += step_press; }
			if (j % save_step == 0) {
				SaveConfiguration(&Particles, &connect_list);
				Save_StressField(&Particles, &connect_list);
			}
			if (abs(extr_x)>0) { extention_system(&Particles, extr_x, axis_x, w, h); }
			if (abs(extr_y) > 0) { extention_system(&Particles, extr_y, axis_y, w, h); }
			//V_liq += q_0*dt;
			V_liq += q_in_fr*dt;
			FindPressInFractures(V_liq, &FracturesBond, pressInFr, num_fracks,bond_len);
			//press = V_liq / PressBond.size() / PressBond.size();
			
			//w *= (1 + extr_x);
			//h *= (1 + extr_y);
		
	}
	//delete[] grid;
	delete[]pressInFr;
	for (USI ii = 0; ii<numCellsX;ii++)
	{
	delete[] gridBonds[ii];
	}
	delete[] gridBonds;
	gridBonds = 0;
	for (USI ii = 0; ii<numCellsX;ii++)
	{
		delete[] gridPressBond[ii];
	}
	delete[] gridPressBond;
	gridPressBond = 0;
    Particles.clear();
    connect_list.clear();
}
void add_new_ball(vector<Atom2D>* P, double x, double y, double w, double h, double bond, double r2, double periodic_border) {
    double aCut = periodic_border;
    if (x < 0 || x > w || y < 0 || y > h) return;

    for (int i = 0; i < P->size(); i++) {
        double rx = P->at(i).r.x - x;
        double ry = P->at(i).r.y - y;
        double rLen2 = rx * rx + ry * ry;
        if (rLen2 < r2) return;
    }
    P->push_back(*(new Atom2D(x, y, false)));//истинная частица

    // проверка, стоит ли частица рядом с границей, и если да, то создание ghost частиц и добавление их в массив
    // 4 случая в углах, если нет, то два у стен и два у потолка/пола

    // расстояния до left, rigth, top, bottom
    double ul = abs(0 - x), ur = abs(w - x), ut = abs(0 - y), ub = abs(h - y);
    int nomer;
    if ((ul < aCut) && (ut < aCut)) {
        // в каждом угловом случае нужно сделать 3 ghost частицы
        P->push_back(*(new Atom2D(x + w, y, true)));
        P->push_back(*(new Atom2D(x, y + h, true)));
        P->push_back(*(new Atom2D(x + w, y + h, true)));
        nomer = P->size() - 1;
  //      P->at(nomer-3).obr_x = true;
  //      P->at(nomer-3).number_x = nomer-2;
		//P->at(nomer - 3).r_obr_x = 1;
  //      //P->at(nomer-3).r_obr_x = P->at(nomer-2).r - P->at(nomer-3).r;
  //      P->at(nomer-3).obr_y = true;
  //      P->at(nomer-3).number_y = nomer-1;
		//P->at(nomer - 3).r_obr_y = 1;
  //     // P->at(nomer-3).r_obr_y = P->at(nomer-1).r - P->at(nomer-3).r;
  //      P->at(nomer-3).obr_xy = true;
  //      P->at(nomer-3).number_xy = nomer;
		//P->at(nomer - 3).r_obr_xy = 22;
        //P->at(nomer-3).r_obr_xy = P->at(nomer).r - P->at(nomer-3).r;
    } else if ((ul < aCut) && (ub < aCut)) {
        P->push_back(*(new Atom2D(x + w, y, true)));
        P->push_back(*(new Atom2D(x, y - h, true)));
        P->push_back(*(new Atom2D(x + w, y - h, true)));
    } else if ((ur < aCut) && (ut < aCut)) {
        P->push_back(*(new Atom2D(x - w, y, true)));
        P->push_back(*(new Atom2D(x, y + h, true)));
        P->push_back(*(new Atom2D(x - w, y + h, true)));
        nomer = P->size() - 1;
  
    } else if ((ur < aCut) && (ub < aCut)) {
        P->push_back(*(new Atom2D(x - w, y, true)));
        P->push_back(*(new Atom2D(x, y - h, true)));
        P->push_back(*(new Atom2D(x - w, y - h, true)));
        nomer = P->size() - 1;
  
    } 
}

void add_new_ball1(vector<Atom2D>* P, double x, double y,Cell**grid,int numCellsX, int numCellsY, double w, double h, double bond, double r2, double periodic_border) {
    double aCut = periodic_border;
    if (x < 0 || x > w || y < 0 || y > h) return;
	int numAddPart = 0;
	int numgridX = int((x+periodic_border)/bond);//номер ячейки по х куда попала частица
	int numgridY = int((y+periodic_border)/bond);
	int arrNumCellX[3];//массив номеров ячеек проверяем по х
	int arrNumCellY[3];
	int X_num = 3;//количество ячеек проверяем по х
	int Y_num = 3;//количество ячеек проверяем по н
	if (numgridX == numCellsX-1){
		//arrNumCellX = new int [2];
		arrNumCellX[0] = numgridX-1;
		arrNumCellX[1] = numgridX;
		X_num = 2;
	}
	else if (numgridX == 0){
		//arrNumCellX = new int [2];
		arrNumCellX[0] = 0;
		arrNumCellX[1] = 1;
		X_num = 2;
	}
	else {
		//arrNumCellX = new int [3];
		arrNumCellX[0] = numgridX-1;
		arrNumCellX[1] = numgridX;
		arrNumCellX[2] = numgridX+1;
	}
	if (numgridY == numCellsY-1){
		//arrNumCellY = new int [2];
		arrNumCellY[0] = numgridY-1;
		arrNumCellY[1] = numgridY;
		Y_num = 2;
	}
	else if (numgridY == 0){
		//arrNumCellY = new int [2];
		arrNumCellY[0] = 0;
		arrNumCellY[1] = 1;
		Y_num = 2;
	}
	else {
		//arrNumCellY = new int [3];
		arrNumCellY[0] = numgridY-1;
		arrNumCellY[1] = numgridY;
		arrNumCellY[2] = numgridY+1;
	}
    for (int i = 0; i<X_num;i++)
	{
		for (int j = 0; j<Y_num;j++)
		{
			int ii = arrNumCellX[i];
			int jj = arrNumCellY[j];
			int NN = grid[ii][jj].numberPart;
			for (int k = 0;k<NN;k++)
			{
				int nomer = grid[ii][jj].arrNumPart[k];
				double rx = P->at(nomer).r.x - x;
				double ry = P->at(nomer).r.y - y;
				double rLen2 = rx * rx + ry * ry;
				if (rLen2 < r2) return;
			}
		}
	}
	//delete []arrNumCellX;
	//delete []arrNumCellY;
	Atom2D pr(x, y, false);
	double ul = abs(0 - x), ur = w-x, ut = y, ub = h-y;
	if ((ul<w*part_fixed)||(ut<h*part_fixed)||(ur<w*part_fixed)||(ub<h*part_fixed))
	{
		pr.im = true;
	}
	P->push_back(pr);//истинная частица
	grid[numgridX][numgridY].Add(P->size()-1);
    
}
void add_new_ball_with_periodic(vector<Atom2D>* P, double x, double y, Cell ** grid, int numCellsX, int numCellsY, double w, double h, double bond, double r2, double periodic_border)
{
	double aCut = periodic_border;
	if (x < 0 || x > w || y < 0 || y > h) return;
	int numAddPart = 0;
	int numgridX = int((x + periodic_border) / bond);//номер ячейки по х куда попала частица
	int numgridY = int((y + periodic_border) / bond);
	int arrNumCellX[3];//массив номеров ячеек проверяем по х
	int arrNumCellY[3];
	int X_num = 3;//количество ячеек проверяем по х
	int Y_num = 3;//количество ячеек проверяем по н
	if (numgridX == numCellsX - 1) {
		//arrNumCellX = new int [2];
		arrNumCellX[0] = numgridX - 1;
		arrNumCellX[1] = numgridX;
		X_num = 2;
	}
	else if (numgridX == 0) {
		//arrNumCellX = new int [2];
		arrNumCellX[0] = 0;
		arrNumCellX[1] = 1;
		X_num = 2;
	}
	else {
		//arrNumCellX = new int [3];
		arrNumCellX[0] = numgridX - 1;
		arrNumCellX[1] = numgridX;
		arrNumCellX[2] = numgridX + 1;
	}
	if (numgridY == numCellsY - 1) {
		//arrNumCellY = new int [2];
		arrNumCellY[0] = numgridY - 1;
		arrNumCellY[1] = numgridY;
		Y_num = 2;
	}
	else if (numgridY == 0) {
		//arrNumCellY = new int [2];
		arrNumCellY[0] = 0;
		arrNumCellY[1] = 1;
		Y_num = 2;
	}
	else {
		//arrNumCellY = new int [3];
		arrNumCellY[0] = numgridY - 1;
		arrNumCellY[1] = numgridY;
		arrNumCellY[2] = numgridY + 1;
	}
	for (int i = 0; i<X_num;i++)
	{
		for (int j = 0; j<Y_num;j++)
		{
			int ii = arrNumCellX[i];
			int jj = arrNumCellY[j];
			int NN = grid[ii][jj].numberPart;
			for (int k = 0;k<NN;k++)
			{
				int nomer = grid[ii][jj].arrNumPart[k];
				double rx = P->at(nomer).r.x - x;
				double ry = P->at(nomer).r.y - y;
				double rLen2 = rx * rx + ry * ry;
				if (rLen2 < r2) return;
			}
		}
	}
	//delete []arrNumCellX;
	//delete []arrNumCellY;
	Atom2D pr(x, y, false);
	P->push_back(pr);//истинная частица
	grid[numgridX][numgridY].Add(P->size() - 1);
	// проверка, стоит ли частица рядом с границей, и если да, то создание ghost частиц и добавление их в массив
	// 4 случая в углах, если нет, то два у стен и два у потолка/пола

	// расстояния до left, rigth, top, bottom
	double ul = abs(0 - x), ur = abs(w - x), ut = abs(0 - y), ub = abs(h - y);
	int nomer;
	if ((ul < aCut) && (ut < aCut)) {
		// в каждом угловом случае нужно сделать 3 ghost частицы
		pr = Atom2D(x + w, y, true);
		P->push_back(pr);
		pr = Atom2D(x, y + h, true);
		P->push_back(pr);
		pr = Atom2D(x + w, y + h, true);
		P->push_back(pr);
		nomer = P->size() - 1;
		P->at(nomer - 3).obr_x = true;
		P->at(nomer - 3).number_x = nomer - 2;
		P->at(nomer - 3).r_obr_x = 1;
		//P->at(nomer-3).r_obr_x = P->at(nomer-2).r - P->at(nomer-3).r;
		P->at(nomer - 3).obr_y = true;
		P->at(nomer - 3).number_y = nomer - 1;
		P->at(nomer - 3).r_obr_y = 1;
		// P->at(nomer-3).r_obr_y = P->at(nomer-1).r - P->at(nomer-3).r;
		P->at(nomer - 3).obr_xy = true;
		P->at(nomer - 3).number_xy = nomer;
		P->at(nomer - 3).r_obr_xy = 22;
		//P->at(nomer-3).r_obr_xy = P->at(nomer).r - P->at(nomer-3).r;
	}
	else if ((ul < aCut) && (ub < aCut)) {
		pr = Atom2D(x + w, y, true);
		P->push_back(pr);
		pr = Atom2D(x, y - h, true);
		P->push_back(pr);
		pr = Atom2D(x + w, y - h, true);
		P->push_back(pr);
		nomer = P->size() - 1;
		P->at(nomer - 3).obr_x = true;
		P->at(nomer - 3).number_x = nomer - 2;
		P->at(nomer - 3).r_obr_x = 1;
		// P->at(nomer-3).r_obr_x = P->at(nomer-2).r - P->at(nomer-3).r;
		P->at(nomer - 3).obr_y = true;
		P->at(nomer - 3).number_y = nomer - 1;
		P->at(nomer - 3).r_obr_y = -1;
		//P->at(nomer-3).r_obr_y = P->at(nomer-1).r - P->at(nomer-3).r;
		P->at(nomer - 3).obr_xy = true;
		P->at(nomer - 3).number_xy = nomer;
		P->at(nomer - 3).r_obr_xy = 21;
		// P->at(nomer-3).r_obr_xy = P->at(nomer).r - P->at(nomer-3).r;
	}
	else if ((ur < aCut) && (ut < aCut)) {
		pr = Atom2D(x - w, y, true);
		P->push_back(pr);
		pr = Atom2D(x, y + h, true);
		P->push_back(pr);
		pr = Atom2D(x - w, y + h, true);
		P->push_back(pr);
		nomer = P->size() - 1;
		P->at(nomer - 3).obr_x = true;
		P->at(nomer - 3).number_x = nomer - 2;
		P->at(nomer - 3).r_obr_x = -1;
		// P->at(nomer-3).r_obr_x = P->at(nomer-2).r - P->at(nomer-3).r;
		P->at(nomer - 3).obr_y = true;
		P->at(nomer - 3).number_y = nomer - 1;
		P->at(nomer - 3).r_obr_y = 1;
		//P->at(nomer-3).r_obr_y = P->at(nomer-1).r - P->at(nomer-3).r;
		P->at(nomer - 3).obr_xy = true;
		P->at(nomer - 3).number_xy = nomer;
		P->at(nomer - 3).r_obr_xy = 12;
		//P->at(nomer-3).r_obr_xy = P->at(nomer).r - P->at(nomer-3).r;
	}
	else if ((ur < aCut) && (ub < aCut)) {
		pr = Atom2D(x - w, y, true);
		P->push_back(pr);
		pr = Atom2D(x, y - h, true);
		P->push_back(pr);
		pr = Atom2D(x - w, y - h, true);
		P->push_back(pr);
		nomer = P->size() - 1;
		P->at(nomer - 3).obr_x = true;
		P->at(nomer - 3).number_x = nomer - 2;
		P->at(nomer - 3).r_obr_x = -1;
		// P->at(nomer-3).r_obr_x = P->at(nomer-2).r - P->at(nomer-3).r;
		P->at(nomer - 3).obr_y = true;
		P->at(nomer - 3).number_y = nomer - 1;
		P->at(nomer - 3).r_obr_y = -1;
		//P->at(nomer-3).r_obr_y = P->at(nomer-1).r - P->at(nomer-3).r;
		P->at(nomer - 3).obr_xy = true;
		P->at(nomer - 3).number_xy = nomer;
		P->at(nomer - 3).r_obr_xy = 11;
		//P->at(nomer-3).r_obr_xy = P->at(nomer).r - P->at(nomer-3).r;
	}
	else {
		if (ul < aCut) {
			pr = Atom2D(x + w, y, true);
			P->push_back(pr);
			nomer = P->size() - 1;
			P->at(nomer - 1).obr_x = true;
			P->at(nomer - 1).number_x = nomer;
			P->at(nomer - 1).r_obr_x = 1;
			//P->at(nomer-1).r_obr_x = P->at(nomer).r - P->at(nomer-1).r;
		}
		if (ur < aCut) {
			pr = Atom2D(x - w, y, true);
			P->push_back(pr);
			nomer = P->size() - 1;
			P->at(nomer - 1).obr_x = true;
			P->at(nomer - 1).number_x = nomer;
			P->at(nomer - 1).r_obr_x = -1;
			//P->at(nomer-1).r_obr_x = P->at(nomer).r - P->at(nomer-1).r;
		}
		if (ut < aCut) {
			pr = Atom2D(x, y + h, true);
			P->push_back(pr);
			nomer = P->size() - 1;
			P->at(nomer - 1).obr_y = true;
			P->at(nomer - 1).number_y = nomer;
			P->at(nomer - 1).r_obr_y = 1;
			//P->at(nomer-1).r_obr_y = P->at(nomer).r - P->at(nomer-1).r;
		}
		if (ub < aCut) {
			pr = Atom2D(x, y - h, true);
			P->push_back(pr);
			nomer = P->size() - 1;
			P->at(nomer - 1).obr_y = true;
			P->at(nomer - 1).number_y = nomer;
			P->at(nomer - 1).r_obr_y = -1;
			//P->at(nomer-1).r_obr_y = P->at(nomer).r - P->at(nomer-1).r;
		}
	}

}
bool sortFunc(connection a, connection b) {
    return (a.eq_st_sq < b.eq_st_sq);
}
void output_to_file(vector<Atom2D> *P, vector<connection> *B, double r, double w, double h,int num) {
    ofstream out_P;
	char file_name[20];
	sprintf(file_name, "Data.txt", num);
    out_P.open(file_name);
	out_P << "Coordinates"<<endl;
        for (int i = 0; i < P->size(); i++) {
            out_P << P->at(i).r.x << " " << P->at(i).r.y << " " << P->at(i).im << endl;
        }
        //out_P << "\";" << endl;

        out_P << "Bonds"<<endl;
        for (int i = 0; i < B->size(); i++) {
            out_P << B->at(i).i << " " << B->at(i).j << endl;
        }
        /*out_P << "\";" << endl;
        out_P << "r = " << r << ";" << endl;
        out_P << "w = " << w << ";" << endl;
        out_P << "h = " << h << ";" << endl;*/
    out_P.close();
}
void start_new_system(vector<Atom2D> *P, vector<connection> *B, double w, double h, double bond_len, double r, double periodic_border, int max_bonds, int attempts) {
    uniform_real_distribution<double> unif(0,1);
    // random_device rd;
    // mt19937 mt(rd());

    // time(NULL) - не самый лучший вариант, но для нашей задачи подойдет (пока без распараллеливания)
    // для распараллеливания нужен random_device
    // для варианта с random_device на windows нужен свежий MinGW
    // на Linux работает без проблем
    mt19937 mt(time(NULL));

    for (int i = 0; i < attempts; i++){
        add_new_ball(P, unif(mt) * w, unif(mt) * h, w, h, bond_len, r*r, periodic_border);
    }

    for (int i = 0; i < P->size() - 1; i++) {
        if (1)//(!(P->at(i).im))
            {
                for (int j = i + 1; j < P->size(); j++){
                if (1)//(!(P->at(j).im))
                {
                    double rx = P->at(i).r.x - P->at(j).r.x;
                    double ry = P->at(i).r.y - P->at(j).r.y;
                    double rLen2 = rx * rx + ry * ry;
                    if (rLen2 < bond_len*bond_len) {
                    connection *bond_obj = new connection(i, j, rLen2);
                    P->at(i).bonds.push_back(*bond_obj);

                    //P->at(j).bonds.push_back(*bond_obj);
                    }
                }

                }
            }
    }

    // сортировка по длине связи по возрастанию
    for (int i = 0; i < P->size(); i++) {
        sort (P->at(i).bonds.begin(), P->at(i).bonds.end(), sortFunc);
    }

    /*for (int i = 0; i < P->size(); i++) {
        for (int j = 0; j < P->at(i).bonds.size(); j++) {
            if (0){//(j >= max_bonds) {
                P->at(i).bonds.at(j).pulled = true;
            }
        }
    }*/

    for (int i = 0; i < P->size(); i++) {
        for (int j = 0; j < P->at(i).bonds.size(); j++) {
			P->at(i).arrNumBond.push_back(B->size());
			P->at(i).numSpringBond++;
            if (!P->at(i).bonds.at(j).pulled) B->push_back(P->at(i).bonds.at(j));
        }
		P->at(i).bonds.clear();
    }
    /*for (int i= 0; i<P->size();i++)
    {
        P->at(i).bonds.clear();
    }*/

}
void start_new_system1(vector<Atom2D> *P, vector<connection> *B, Cell**grid, Cell** gridBonds, int numCellsX, int numCellsY, double w, double h, double bond_len, double r, double periodic_border, int max_bonds, int attempts) {
	uniform_real_distribution<double> unif(0, 1);
	// random_device rd;
	// mt19937 mt(rd());

	// time(NULL) - не самый лучший вариант, но для нашей задачи подойдет (пока без распараллеливания)
	// для распараллеливания нужен random_device
	// для варианта с random_device на windows нужен свежий MinGW
	// на Linux работает без проблем
	mt19937 mt(time(NULL));

	for (int i = 0; i < attempts; i++) {
		if (typeBC == 1) {
			add_new_ball1(P, unif(mt) * w, unif(mt) * h, grid, numCellsX, numCellsY, w, h, bond_len, r*r, periodic_border);
		}
		else {
			add_new_ball_with_periodic(P, unif(mt) * w, unif(mt) * h, grid, numCellsX, numCellsY, w, h, bond_len, r*r, periodic_border);
		}
		
	}
	cout << "size  " << P->size() << "\n";


	for (int i = 0; i < P->size() - 1; i++)
	{
		float x = P->at(i).r.x;
		float y = P->at(i).r.y;
		if ((x > w*1.5) || (y>h*1.5))
		{
			system("pause");
		}
		int numgridX = int((x + periodic_border) / bond_len);//номер ячейки по х куда попала частица
		int numgridY = int((y + periodic_border) / bond_len);
		int* arrNumCellX;//массив номеров ячеек проверяем по х
		int* arrNumCellY;
		int X_num = 3;//количество ячеек проверяем по х
		int Y_num = 3;//количество ячеек проверяем по н
		if (numgridX == numCellsX - 1) {
			arrNumCellX = new int[2];
			arrNumCellX[0] = numgridX - 1;
			arrNumCellX[1] = numgridX;
			X_num = 2;
		}
		else if (numgridX == 0) {
			arrNumCellX = new int[2];
			arrNumCellX[0] = 0;
			arrNumCellX[1] = 1;
			X_num = 2;
		}
		else {
			arrNumCellX = new int[3];
			arrNumCellX[0] = numgridX - 1;
			arrNumCellX[1] = numgridX;
			arrNumCellX[2] = numgridX + 1;
		}
		if (numgridY == numCellsY - 1) {
			arrNumCellY = new int[2];
			arrNumCellY[0] = numgridY - 1;
			arrNumCellY[1] = numgridY;
			Y_num = 2;
		}
		else if (numgridY == 0) {
			arrNumCellY = new int[2];
			arrNumCellY[0] = 0;
			arrNumCellY[1] = 1;
			Y_num = 2;
		}
		else {
			arrNumCellY = new int[3];
			arrNumCellY[0] = numgridY - 1;
			arrNumCellY[1] = numgridY;
			arrNumCellY[2] = numgridY + 1;
		}

		if (1)//(!(P->at(i).im))
		{
			for (int i1 = 0; i1 < X_num;i1++)
			{
				for (int j1 = 0; j1 < Y_num;j1++)
				{

					int ii = arrNumCellX[i1];
					int jj = arrNumCellY[j1];
					int NN = grid[ii][jj].numberPart;
					for (int k = 0;k < NN;k++)
					{
						if (1)//(!(P->at(j).im))
						{
							int nomer = grid[ii][jj].arrNumPart[k];
							if (nomer != i) {
								float rx = P->at(i).r.x - P->at(nomer).r.x;
								float ry = P->at(i).r.y - P->at(nomer).r.y;
								float rLen2 = rx * rx + ry * ry;
								if (rLen2 < bond_len*bond_len) {
									//connection *bond_obj = new connection(i, nomer, rLen2);
									connection bond_obj(i, nomer, rLen2);
									P->at(i).bonds.push_back(bond_obj);
									P->at(i).numSpringBond++;
									//P->at(i).bonds.push_back(*bond_obj);
								}
								//P->at(j).bonds.push_back(*bond_obj);
							}
						}
					}
				}
			}


		}
		delete[]arrNumCellX;
		delete[]arrNumCellY;
	}

	// сортировка по длине связи по возрастанию
	for (int i = 0; i < P->size(); i++) {
		sort(P->at(i).bonds.begin(), P->at(i).bonds.end(), sortFunc);
	}

	for (int i = 0; i < P->size(); i++) {
		for (int j = 0; j < P->at(i).bonds.size(); j++) {
			if (0) {//(j >= max_bonds) {
				P->at(i).bonds.at(j).pulled = true;
			}
		}
	}

	/*for (int i = 0; i < P->size(); i++) {
		//P->at(i).CreateArrNumBond(P->at(i).bonds.size());

		for (int j = 0; j < P->at(i).bonds.size(); j++) {
			P->at(i).arrNumBond.push_back(B->size());

			if (!P->at(i).bonds.at(j).pulled) B->push_back(P->at(i).bonds.at(j));
		}
		P->at(i).bonds.clear();
	}
	*/
	for (int i = 0; i < P->size(); i++) {
		//P->at(i).CreateArrNumBond(P->at(i).bonds.size());

		for (int j = 0; j < P->at(i).bonds.size(); j++) {
			P->at(i).arrNumBond.push_back(B->size());

			if (!P->at(i).bonds.at(j).pulled)
			{
				int a2i = P->at(i).bonds.at(j).i;
				int  a2j = P->at(i).bonds.at(j).j;
				double x = (P->at(a2i).r.x + P->at(a2j).r.x) / 2;
				double y = (P->at(a2i).r.y + P->at(a2j).r.y) / 2;
				int numgridX = int((x + periodic_border) / bond_len);//номер ячейки по х куда попала частица
				int numgridY = int((y + periodic_border) / bond_len);
				gridBonds[numgridX][numgridY].Add(B->size());
				B->push_back(P->at(i).bonds.at(j));
			}
		}
		P->at(i).bonds.clear();
	}
}
void start_new_system_treangle(vector<Atom2D>* P, vector<connection>* B, Cell ** grid, Cell**gridBonds, int numCellsX, int numCellsY, double w, double h, double bond_len, double r, double periodic_border, int max_bonds, int attempts)
{
	double step_x = 0;//шаг вдоль х для частиц
	double step_y = 0;//шаг вдоль y для частиц
	step_y = bond_len*0.48;//Чтобы через ряд соединялись частицы
	step_x = step_y * 2 / sqrt(3);
	int numPartX = int(w / step_x);//число частиц вдоль х
	int numPartY = int(h / step_y);
	for (int j = 0;j < numPartY; j++)
	{
		int otkl = j % 2;
		for (int i = 0;i < numPartX-otkl;i++)
		{
			double x_part = i*step_x+0.5*step_x*otkl;
			double y_part = j*step_y;

			add_new_ball1(P,  y_part, x_part, grid, numCellsX, numCellsY, w, h, bond_len, r*r, periodic_border);
		}
	}
	cout << "size  " << P->size() << "\n";


	for (int i = 0; i < P->size() - 1; i++)
	{
		double x = P->at(i).r.x;
		double y = P->at(i).r.y;
		if ((x > w*1.5) || (y>h*1.5))
		{
			system("pause");
		}
		int numgridX = int((x + periodic_border) / bond_len);//номер ячейки по х куда попала частица
		int numgridY = int((y + periodic_border) / bond_len);
		int* arrNumCellX;//массив номеров ячеек проверяем по х
		int* arrNumCellY;
		int X_num = 3;//количество ячеек проверяем по х
		int Y_num = 3;//количество ячеек проверяем по н
		if (numgridX == numCellsX - 1) {
			arrNumCellX = new int[2];
			arrNumCellX[0] = numgridX - 1;
			arrNumCellX[1] = numgridX;
			X_num = 2;
		}
		else if (numgridX == 0) {
			arrNumCellX = new int[2];
			arrNumCellX[0] = 0;
			arrNumCellX[1] = 1;
			X_num = 2;
		}
		else {
			arrNumCellX = new int[3];
			arrNumCellX[0] = numgridX - 1;
			arrNumCellX[1] = numgridX;
			arrNumCellX[2] = numgridX + 1;
		}
		if (numgridY == numCellsY - 1) {
			arrNumCellY = new int[2];
			arrNumCellY[0] = numgridY - 1;
			arrNumCellY[1] = numgridY;
			Y_num = 2;
		}
		else if (numgridY == 0) {
			arrNumCellY = new int[2];
			arrNumCellY[0] = 0;
			arrNumCellY[1] = 1;
			Y_num = 2;
		}
		else {
			arrNumCellY = new int[3];
			arrNumCellY[0] = numgridY - 1;
			arrNumCellY[1] = numgridY;
			arrNumCellY[2] = numgridY + 1;
		}

		if (1)//(!(P->at(i).im))
		{
			for (int i1 = 0; i1<X_num;i1++)
			{
				for (int j1 = 0; j1<Y_num;j1++)
				{

					int ii = arrNumCellX[i1];
					int jj = arrNumCellY[j1];
					int NN = grid[ii][jj].numberPart;
					for (int k = 0;k<NN;k++)
					{
						if (1)//(!(P->at(j).im))
						{
							int nomer = grid[ii][jj].arrNumPart[k];
							if (nomer != i) {
								double rx = P->at(i).r.x - P->at(nomer).r.x;
								double ry = P->at(i).r.y - P->at(nomer).r.y;
								double rLen2 = rx * rx + ry * ry;
								if (rLen2 < pow(0.56*bond_len,2)) {
									//connection *bond_obj = new connection(i, nomer, rLen2);
									connection bond_obj(i, nomer, rLen2);
									P->at(i).bonds.push_back(bond_obj);
									P->at(i).numSpringBond++;
									//P->at(i).bonds.push_back(*bond_obj);
								}
								//P->at(j).bonds.push_back(*bond_obj);
							}
						}
					}
				}
			}


		}
		delete[]arrNumCellX;
		delete[]arrNumCellY;
	}

	// сортировка по длине связи по возрастанию
	for (int i = 0; i < P->size(); i++) {
		sort(P->at(i).bonds.begin(), P->at(i).bonds.end(), sortFunc);
	}

	for (int i = 0; i < P->size(); i++) {
		for (int j = 0; j < P->at(i).bonds.size(); j++) {
			if (0) {//(j >= max_bonds) {
				P->at(i).bonds.at(j).pulled = true;
			}
		}
	}

	for (int i = 0; i < P->size(); i++) {
		//P->at(i).CreateArrNumBond(P->at(i).bonds.size());
		
		for (int j = 0; j < P->at(i).bonds.size(); j++) {
			P->at(i).arrNumBond.push_back(B->size());

			if (!P->at(i).bonds.at(j).pulled)
			{
				int a2i = P->at(i).bonds.at(j).i;
				int  a2j = P->at(i).bonds.at(j).j;
				double x = (P->at(a2i).r.x + P->at(a2j).r.x) / 2;
				double y = (P->at(a2i).r.y + P->at(a2j).r.y) / 2;
				int numgridX = int((x + periodic_border) / bond_len);//номер ячейки по х куда попала частица
				int numgridY = int((y + periodic_border) / bond_len);
				gridBonds[numgridX][numgridY].Add(B->size());
				B->push_back(P->at(i).bonds.at(j));
			}
		}
		P->at(i).bonds.clear();
	}
	/*СОЕДИНЯЕМ СВЯЗИ КАНАЛАМИ С ЖИДКОСТЯМИ*/
	if (nodes)
	{
		bool presence = false;
		//int numggX = int((w + 2 * periodic_border) / bond_len) + 1;//количество ячеек вдоль х
		//int numggY = int((h + 2 * periodic_border) / bond_len) + 1;//количество ячеек вдоль у
		for (int jjCon = 0;jjCon < B->size();jjCon++)//проход по связьям
		{
			int a2i = B->at(jjCon).i;
			int  a2j = B->at(jjCon).j;
			double x = (P->at(a2i).r.x + P->at(a2j).r.x) / 2;
			double y = (P->at(a2i).r.y + P->at(a2j).r.y) / 2;
			int numgridX = int((x + periodic_border) / bond_len);//номер ячейки по х куда попала связь
			int numgridY = int((y + periodic_border) / bond_len);
			int* arrNumCellX;//массив номеров ячеек проверяем по х
			int* arrNumCellY;
			int X_num = 3;//количество ячеек проверяем по х
			int Y_num = 3;//количество ячеек проверяем по н
			if (numgridX == numCellsX - 1) {
				arrNumCellX = new int[2];
				arrNumCellX[0] = numgridX - 1;
				arrNumCellX[1] = numgridX;
				X_num = 2;
			}
			else if (numgridX == 0) {
				arrNumCellX = new int[2];
				arrNumCellX[0] = 0;
				arrNumCellX[1] = 1;
				X_num = 2;
			}
			else {
				arrNumCellX = new int[3];
				arrNumCellX[0] = numgridX - 1;
				arrNumCellX[1] = numgridX;
				arrNumCellX[2] = numgridX + 1;
			}
			if (numgridY == numCellsY - 1) {
				arrNumCellY = new int[2];
				arrNumCellY[0] = numgridY - 1;
				arrNumCellY[1] = numgridY;
				Y_num = 2;
			}
			else if (numgridY == 0) {
				arrNumCellY = new int[2];
				arrNumCellY[0] = 0;
				arrNumCellY[1] = 1;
				Y_num = 2;
			}
			else {
				arrNumCellY = new int[3];
				arrNumCellY[0] = numgridY - 1;
				arrNumCellY[1] = numgridY;
				arrNumCellY[2] = numgridY + 1;
			}
			for (int iiCon1 = 0;iiCon1 < X_num;iiCon1++)//По ячейкам Х
			{
				for (int iiCon2 = 0;iiCon2 < Y_num;iiCon2++)//По ячейкам У
				{
					int ii = arrNumCellX[iiCon1];
					int jj = arrNumCellY[iiCon2];
					int NN = gridBonds[ii][jj].numberPart;
					for (int k = 0;k < NN;k++)
					{
						int nomer = gridBonds[ii][jj].arrNumPart.at(k);
						auto res = std::find(B->at(jjCon).numLiquidBond.begin(),B->at(jjCon).numLiquidBond.end(),nomer);
						if((jjCon==nomer)||((res!=B->at(jjCon).numLiquidBond.end())))
						{
							break;
						}
						int ai1 = B->at(jjCon).i;//итый атом 1 связи
						int ai2 = B->at(nomer).i;//итый атом 2 связи
						int aj1 = B->at(jjCon).j;
						int aj2 = B->at(nomer).j;
						double x1 = 0.5*(P->at(ai1).r.x + P->at(aj1).r.x);
						double x2 = 0.5*(P->at(ai2).r.x + P->at(aj2).r.x);
						double y1 = 0.5*(P->at(ai1).r.y + P->at(aj1).r.y);
						double y2 = 0.5*(P->at(ai2).r.y + P->at(aj2).r.y);
						double dist = sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));
						if (dist < part_bond_liquid* bond_len)
						{
							B->at(nomer).numLiquidBond.push_back(jjCon);
							B->at(jjCon).numLiquidBond.push_back(nomer);
						}
					}
				}


			}
			delete[]arrNumCellX;
			delete[]arrNumCellY;
		}
	
	}
}
void SetConf(vector<Atom2D>* as, vector<connection>* l, double w, double h, double bond_len)
{
	int n = int(w/2*bond_len);
	double x,y;
	for (int j = 0; j<n;j++)
		for (int i = 0;i<n;i++)
		{
			if(j%2==0)
			{
				x = i*2*bond_len;
				y = j*bond_len;
				as->push_back(*(new Atom2D(x,y,true)));
			}
			else
			{
				if (i!=(n-1))
				{
					x = bond_len+i*2*bond_len;
					y = j*bond_len;
					as->push_back(*(new Atom2D(x,y,true)));
				}
			}
		}
		for (int i = 0;i<as->size();i++)
		{


		}

}
void Step(vector <Atom2D>* as, vector<connection>* l,double w,double h)
{
    std::vector<connection>::iterator it;
    double length;
    double ff;
    Vect2D dr;
	Vect2D ww = Vect2D(float(w), 0);
	Vect2D hh = Vect2D(0, float(h));
	Vect2D wh = Vect2D(float(w), float(h));
	for (int ii = (l->size()-1);ii>=0;ii--)//проход по всем связям. Перед циклом силы в частицах 0
    {
		dr = as->at(l->at(ii).j).r - as->at(l->at(ii).i).r;
		if (dr.Abs()==0){
				//system("pause");
		}
		length = sqrt(dr.Sqr()) - sqrt(l->at(ii).eq_st_sq); //растяжение-сжатие пружины
		if (length!=0){
				//system("pause");
		}
		if (l->at(ii).pressure)
		{
			ff = -press;
			if (ff==0){
				system("pause");
			}
			as->at(l->at(ii).j).f -= (float(ff)/sqrt(dr.Sqr()))*dr;
			as->at(l->at(ii).i).f += (float(ff) /sqrt(dr.Sqr()))*dr;
		}
		else
		if (dr.Sqr() > pow(l->at(ii).eq_st_sq*e0,1))//проверка на состоянии пружины
        {//случай разрыва связи
			l->at(ii).pressure = true;
		
			num_errase++;
			//l->erase(l->begin()+ii);
			//system("pause");
		}
		else if (!(l->at(ii).free))
		{
			ff = length*k_bond;//it->koef;//модуль силы
			if (ff!=0){
				//system("pause");
			}
			as->at(l->at(ii).j).f-= (float(ff) /sqrt(dr.Sqr()))*dr;
			as->at(l->at(ii).i).f += (float(ff) /sqrt(dr.Sqr()))*dr;
		}

	}
	for (int i = 0; i<as->size();i++)//интегрируем по каждой частице скороксти
    {
        if (!(as->at(i).im))
        {
			as->at(i).f-=as->at(i).v*visc_k;
        }
    }
    for (int i = 0; i<as->size();i++)//интегрируем по каждой частице скороксти
    {
        if (!(as->at(i).im))
        {
            as->at(i).v+=as->at(i).f*dt;
        }
    }
    for (int i = 0; i<as->size();i++)//интегрируем перемещения каждой частице
    {
        if ((!(as->at(i).im)))//&&(!(as->at(i).pressure)))
        {
            as->at(i).r+=as->at(i).v*dt;
    //        if  (((as->at(i).obr_x)||(as->at(i).obr_y)||(as->at(i).obr_xy))&&(as->at(i).r.x < 0 || as->at(i).r.x > w || as->at(i).r.y < 0 || as->at(i).r.y > h))//проверяем стала ли частица мнимой
    //        {
    //            if (as->at(i).obr_x)//пересчет координат
    //            {
				//	as->at(as->at(i).number_x).r = as->at(i).r + as->at(i).r_obr_x*ww; //+as->at(i).r_obr_x;
    //            }
    //            if (as->at(i).obr_y)
    //            {
				//	as->at(as->at(i).number_y).r = as->at(i).r + as->at(i).r_obr_y*hh;//+as->at(i).r_obr_y;
    //            }
    //            if (as->at(i).obr_xy)
    //            {
				//	int px = 0;
				//	int py = 0;
				//	if (as->at(i).r_obr_xy == 22)
				//	{
				//		px = 1;
				//		py = 1;
				//	}
				//	if (as->at(i).r_obr_xy == 12)
				//	{
				//		px = -1;
				//		py = 1;
				//	}
				//	if (as->at(i).r_obr_xy == 21)
				//	{
				//		px = 1;
				//		py = -1;
				//	}
				//	if (as->at(i).r_obr_xy == 22)
				//	{
				//		px = -1;
				//		py = -1;
				//	}
				//	Vect2D pr = wh;
				//	pr.x *= px;
				//	pr.y *= py;
				//	as->at(as->at(i).number_xy).r = as->at(i).r+pr;
    //            }
    //            as->at(i).im = true;
				//if ((as->at(i).obr_x)&&(!(as->at(as->at(i).number_x).r.x < 0 || as->at(as->at(i).number_x).r.x > w || as->at(as->at(i).number_x).r.y < 0 || as->at(as->at(i).number_x).r.y > h)))//проверяем стал ли образ по оси х ральной частицей
				//{
				//	as->at(as->at(i).number_x).im = false;//меняем образ х который теперь стал реальной с рассматриываемой частицей, которая стала теперь мнимой
				//	as->at(as->at(i).number_x).obr_x = true;
				//	as->at(as->at(i).number_x).number_x = i;
				//	if ((as->at(i).r - as->at(as->at(i).number_x).r).x > 0)
				//	{
				//		as->at(as->at(i).number_x).r_obr_x = 1;
				//	}
				//	else {
				//		as->at(as->at(i).number_x).r_obr_x = -1;
				//	}
				//	//as->at(as->at(i).number_x).r_obr_x = as->at(i).r-as->at(as->at(i).number_x).r;
				//	if (as->at(i).obr_y)//у рассматриваемой частицы 3 образа, оставшиеся меняем местами
				//	{
				//		//if ((abs(as->at(i).number_y)>2000)||(abs(as->at(i).number_xy)>2000)){system("pause");}
				//		int px = 0;
				//		int py = 0;
				//		as->at(as->at(i).number_x).number_xy = as->at(i).number_y;
				//		as->at(as->at(i).number_x).number_y = as->at(i).number_xy;
				//		as->at(as->at(i).number_x).obr_xy = true;
				//		as->at(as->at(i).number_x).obr_y = true;
				//		if ((as->at(as->at(i).number_y).r - as->at(as->at(i).number_x).r).x > 0) {
				//			px = 2;
				//		}
				//		else {
				//			px = 1;
				//		}
				//		if ((as->at(as->at(i).number_y).r - as->at(as->at(i).number_x).r).y > 0) {
				//			py = 2;
				//		}
				//		else {
				//			py = 1;
				//		}
				//		as->at(as->at(i).number_x).r_obr_xy = px * 10 + py;
				//		as->at(as->at(i).number_x).r_obr_y = ((py == 2) ? (1) : (-1));
				//		//as->at(as->at(i).number_x).r_obr_xy = as->at(as->at(i).number_y).r - as->at(as->at(i).number_x).r;
				//		//as->at(as->at(i).number_x).r_obr_y = as->at(as->at(i).number_xy).r - as->at(as->at(i).number_x).r;
				//	}
				//	as->at(i).obr_x = false;
				//	as->at(i).obr_y = false;
				//	as->at(i).obr_xy = false;

				//}
				//else if ((as->at(i).obr_y)&&(!(as->at(as->at(i).number_y).r.x < 0 || as->at(as->at(i).number_y).r.x > w || as->at(as->at(i).number_y).r.y < 0 || as->at(as->at(i).number_y).r.y > h)))//проверяем стал ли образ по оси y ральной частицей
				//{
				//	as->at(as->at(i).number_y).im = false;//меняем образ y который теперь стал реальной с рассматриываемой частицей, которая стала теперь мнимой
				//	as->at(as->at(i).number_y).obr_y = true;
				//	as->at(as->at(i).number_y).number_y = i;
				//	if ((as->at(i).r - as->at(as->at(i).number_y).r).y > 0)
				//	{
				//		as->at(as->at(i).number_y).r_obr_y = 1;
				//	}
				//	else {
				//		as->at(as->at(i).number_y).r_obr_y = -1;
				//	}
				//	//as->at(as->at(i).number_y).r_obr_y = as->at(i).r-as->at(as->at(i).number_y).r;
				//	if (as->at(i).obr_x)//у рассматриваемой частицы 3 образа, оставшиеся меняем местами
				//	{
				//		//if ((abs(as->at(i).number_x)>2000)||(abs(as->at(i).number_xy)>2000)){system("pause");}
				//		as->at(as->at(i).number_y).number_xy = as->at(i).number_x;
				//		as->at(as->at(i).number_y).number_x = as->at(i).number_xy;
				//		as->at(as->at(i).number_y).obr_xy = true;
				//		as->at(as->at(i).number_y).obr_x = true;
				//		int px = 0;
				//		int py = 0;
				//		if ((as->at(as->at(i).number_x).r - as->at(as->at(i).number_y).r).x > 0)
				//		{
				//			px = 2;
				//		}
				//		else {
				//			px = 1;
				//		}
				//		if ((as->at(as->at(i).number_x).r - as->at(as->at(i).number_y).r).y > 0)
				//		{
				//			py = 2;
				//		}
				//		else {
				//			py = 1;
				//		}
				//		as->at(as->at(i).number_y).r_obr_xy = px * 10 + py;//as->at(as->at(i).number_x).r - as->at(as->at(i).number_y).r;
				//		if (px == 2)
				//		{
				//			px = 1;
				//		}
				//		else {
				//			px = -1;
				//		}
				//		as->at(as->at(i).number_y).r_obr_x = px;//as->at(as->at(i).number_xy).r - as->at(as->at(i).number_y).r;
				//	}
				//	as->at(i).obr_x = false;
				//	as->at(i).obr_y = false;
				//	as->at(i).obr_xy = false;
				//}
				//else if (as->at(i).obr_xy)//реальной стал образ по ху
				//{
				//	as->at(as->at(i).number_xy).im = false;//меняем образ х который теперь стал реальной с рассматриываемой частицей, которая стала теперь мнимой
				//	as->at(as->at(i).number_xy).obr_xy = true;
				//	as->at(as->at(i).number_xy).number_xy = i;
				//	if ((as->at(i).r - as->at(as->at(i).number_xy).r).x > 0)
				//	{
				//		as->at(as->at(i).number_xy).r_obr_x = 1;
				//	}
				//	else {
				//		as->at(as->at(i).number_xy).r_obr_x = -1;
				//	}
				//	//as->at(as->at(i).number_xy).r_obr_x = as->at(i).r-as->at(as->at(i).number_xy).r;
				//	if (as->at(i).obr_x)//у рассматриваемой частицы 3 образа, оставшиеся меняем местами
				//	{
				//		//if ((abs(as->at(i).number_xy)>2000)||(abs(as->at(i).number_y)>2000)){system("pause");}
				//		as->at(as->at(i).number_xy).number_y = as->at(i).number_x;
				//		as->at(as->at(i).number_xy).number_x = as->at(i).number_y;
				//		as->at(as->at(i).number_xy).obr_x = true;
				//		as->at(as->at(i).number_xy).obr_y = true;
				//		if ((as->at(as->at(i).number_x).r - as->at(as->at(i).number_xy).r).y > 0)
				//		{
				//			as->at(as->at(i).number_xy).r_obr_y = 1;
				//		}
				//		else {
				//			as->at(as->at(i).number_xy).r_obr_y = -1;
				//		}
				//		if ((as->at(as->at(i).number_y).r - as->at(as->at(i).number_xy).r).x > 0)
				//		{
				//			as->at(as->at(i).number_xy).r_obr_x = 1;
				//		}
				//		else {
				//			as->at(as->at(i).number_xy).r_obr_x = -1;
				//		}
				//		//as->at(as->at(i).number_xy).r_obr_y = as->at(as->at(i).number_x).r - as->at(as->at(i).number_xy).r;
				//		//as->at(as->at(i).number_xy).r_obr_x = as->at(as->at(i).number_y).r - as->at(as->at(i).number_xy).r;
				//	}
				//	as->at(i).obr_x = false;
				//	as->at(i).obr_y = false;
				//	as->at(i).obr_xy = false;
				//}


    //        }

    //        else //частица не пересекала границы
    //        {
    //            if (as->at(i).obr_x)
    //            {
    //               
				//	as->at(as->at(i).number_x).r = as->at(i).r+as->at(i).r_obr_x*ww;
    //            }
    //            if (as->at(i).obr_y)
    //            {
    //                as->at(as->at(i).number_y).r = as->at(i).r+as->at(i).r_obr_y*hh;
    //            }
    //            if (as->at(i).obr_xy)
    //            {
				//	Vect2D pr = wh;
				//	int px = 0; 
				//	int py=0;
				//	if (as->at(i).r_obr_xy == 22)
				//	{
				//		px = 1;
				//		py = 1;
				//	}
				//	if (as->at(i).r_obr_xy == 21)
				//	{
				//		px = 1;
				//		py = -1;
				//	}
				//	if (as->at(i).r_obr_xy == 12)
				//	{
				//		px = -1;
				//		py = 1;
				//	}
				//	if (as->at(i).r_obr_xy == 11)
				//	{
				//		px = -1;
				//		py = -1;
				//	}
				//	pr.x *= px;
				//	pr.y *= py;
				//	as->at(as->at(i).number_xy).r = as->at(i).r+pr;
    //            }
    //        }
        }
        as->at(i).f=VECT2D_0;
    }
//    std::cout<<"step finished \n";
}
void Step1(vector <Atom2D>* as, vector<connection>* l, bool**gridPressBond, int numCellsX, int numCellsY, double w, double h, double bond, double periodic_border, vector<int>*PressBond)
{
	std::vector<connection>::iterator it;
	double length;
	double ff;
	Vect2D dr;
	Vect2D ww = Vect2D(float(w), 0);
	Vect2D hh = Vect2D(0, float(h));
	Vect2D wh = Vect2D(float(w), float(h));
	for (int kk = 0; kk < as->size();kk++)//проходим по частицам
	{
		for (int j1 = 0;j1 < as->at(kk).arrNumBond.size();j1++)//проходим по номерам связей имеющихся у частицы
		{
			int jj = as->at(kk).arrNumBond.at(j1);//номер связи
			if (!(l->at(jj).cheek))
			{
				int numSecPart = l->at(jj).j;
				if (numSecPart == kk)
				{
					numSecPart = l->at(jj).i;
				}
				dr = as->at(kk).r - as->at(numSecPart).r;
				length = sqrt(dr.Sqr()) - sqrt(l->at(jj).eq_st_sq); //растяжение-сжатие пружины
				if (l->at(jj).pressure)
				{
					ff = -press;
					as->at(kk).f -= (float(ff) / sqrt(dr.Sqr()))*dr;
					as->at(numSecPart).f += (float(ff) / sqrt(dr.Sqr()))*dr;
					/*if (length < 0) {
						ff = 1 / dr.Abs();
						as->at(kk).f -= (float(ff) / sqrt(dr.Sqr()))*dr;
						as->at(numSecPart).f += (float(ff) / sqrt(dr.Sqr()))*dr;
					}*/
				}
				else
					if (dr.Sqr() > pow(l->at(jj).eq_st_sq*e0, 1))//проверка на состоянии пружины
					{//случай разрыва связи
						double xCurrBond = (as->at(kk).r.x + as->at(numSecPart).r.x)*0.5;
						double yCurrBond = (as->at(kk).r.y + as->at(numSecPart).r.y)*0.5;
						bool f = false;

						int numgridX = int((xCurrBond + periodic_border) / bond);//номер ячейки по х куда попала частица
						int numgridY = int((yCurrBond + periodic_border) / bond);
						int arrNumCellX[3];//массив номеров ячеек проверяем по х
						int arrNumCellY[3];
						int X_num = 3;//количество ячеек проверяем по х
						int Y_num = 3;//количество ячеек проверяем по н
						if (numgridX == numCellsX - 1) {
							//arrNumCellX = new int [2];
							arrNumCellX[0] = numgridX - 1;
							arrNumCellX[1] = numgridX;
							X_num = 2;
						}
						else if (numgridX == 0) {
							//arrNumCellX = new int [2];
							arrNumCellX[0] = 0;
							arrNumCellX[1] = 1;
							X_num = 2;
						}
						else {
							//arrNumCellX = new int [3];
							arrNumCellX[0] = numgridX - 1;
							arrNumCellX[1] = numgridX;
							arrNumCellX[2] = numgridX + 1;
						}
						if (numgridY == numCellsY - 1) {
							//arrNumCellY = new int [2];
							arrNumCellY[0] = numgridY - 1;
							arrNumCellY[1] = numgridY;
							Y_num = 2;
						}
						else if (numgridY == 0) {
							//arrNumCellY = new int [2];
							arrNumCellY[0] = 0;
							arrNumCellY[1] = 1;
							Y_num = 2;
						}
						else {
							//arrNumCellY = new int [3];
							arrNumCellY[0] = numgridY - 1;
							arrNumCellY[1] = numgridY;
							arrNumCellY[2] = numgridY + 1;
						}
						arrNumCellX[0] = numgridX - 1;
						arrNumCellX[1] = numgridX;
						arrNumCellX[2] = numgridX + 1;
						arrNumCellY[0] = numgridY - 1;
						arrNumCellY[1] = numgridY;
						arrNumCellY[2] = numgridY + 1;
						if (true)//((f))
						{
							l->at(jj).pressure = true;
							int ii1 = l->at(jj).i;
							int jj1 = l->at(jj).j;
							as->at(ii1).numSpringBond--;
							as->at(ii1).pressure = true;
							as->at(jj1).numSpringBond--;
							as->at(jj1).pressure = true;
							PressBond->push_back(jj);
							numPressBond++;
							num_errase++;
							//cout << num_errase << "\n";
						}
						else
						{
							l->at(jj).free = true;
							int ii1 = l->at(jj).i;
							int jj1 = l->at(jj).j;
							as->at(ii1).numSpringBond--;
							as->at(jj1).numSpringBond--;
							
						}
						
					}
					else if (!(l->at(jj).free))
					{
						ff = length*k_bond;//it->koef;//модуль силы
						if (ff != 0) {
							//system("pause");
						}
						as->at(kk).f -= (float(ff) / sqrt(dr.Sqr()))*dr;
						as->at(numSecPart).f += (float(ff) / sqrt(dr.Sqr()))*dr;
					}
				l->at(jj).cheek = true;
			}
			else {
				l->at(jj).cheek = false;
			}
		}
		
		//}
	}
	if (typeBC == 1) {
		RecalculateVelosityPositionFixed(as, w, h);
	}
	else {
		RecalculateVelosityPositionPeriodic(as, w, h);
	}
	//

	stepp++;
//    std::cout<<"step finished \n";
}
//void StepWithViscosity(vector<Atom2D>* as, vector<connection>* l, bool ** gridPressBond, int numCellsX, int numCellsY, double w, double h, double bond, double periodic_border, vector<int>* PressBond,double x_c0,double x_c1, double lenBetweenNodes,int numNodes,double** widthFracture, double** Press,double** Visc,double* v_cur)
void StepWithViscosity(vector<Atom2D>* as, vector<connection>* l, double w, double h, double bond, double periodic_border, vector<int>* PressBond, double x_c0, double x_c1, double lenBetweenNodes, int numNodes, vector<double>* widthFracture, vector<double>* Press, vector<double>* Visc, double* v_cur)
{
	//std::vector<connection>::iterator it;
	double length;
	Vect2D dr;
	Vect2D ww = Vect2D(float(w), 0);
	Vect2D hh = Vect2D(0, float(h));
	Vect2D wh = Vect2D(float(w), float(h));
	int* numBondNodes = new int(numNodes);
	
	int ff = 0;
	double** prW = new double*[2];
	for (int ipr = 0;ipr < 2;ipr++)
	{
		prW[ipr] = new double[numNodes];
	}
	for (int i11 = 0;i11 < 2;i11++)
		for (int ipr = 0;ipr < numNodes;ipr++)
		{
			if (stepp<100)
			{
				prW[i11][ipr] = 0;//widthFracture->at(i11*numNodes + ipr);//[i11][ipr];
			}
			else
			{
				prW[i11][ipr] = widthFracture->at(i11*numNodes + ipr);
			}
		}
	int** arr_number = new int*[2];
	for (int i12 = 0;i12 < 2;i12++)
	{
		arr_number[i12] = new int[numNodes];
	}
	for (int i11 = 0;i11 < 2;i11++)
		for (int ipr = 0;ipr < numNodes;ipr++)
		{
			arr_number[i11][ipr] = 0;//widthFracture->at(i11*numNodes + ipr);//[i11][ipr];
			Press->at(i11*numNodes + ipr) = 0;
		}
	for (int i = 0; i < PressBond->size();i++)//определяю давления в узлах
	{
		int num = PressBond->at(i);
		int cr_i = l->at(num).i;
		int cr_j = l->at(num).j;
		if (as->at(cr_i).r.x < as->at(cr_j).r.x)
		{
			int pr = cr_i;
			cr_i = cr_j;
			cr_j = pr;
		}
		double len = abs((as->at(cr_i).r.y - as->at(cr_j).r.y));
		int curNum = 0;
		double x = 0.5*(as->at(cr_i).r.x + as->at(cr_j).r.x);
		double y = 0.5*(as->at(cr_i).r.y + as->at(cr_j).r.y);
		if (y - h*0.5 > 0)
		{
			curNum = int((y - h*0.5 + lenBetweenNodes*0.5) / lenBetweenNodes);
			ff = 0;
		}
		else
		{
			curNum = int((abs(y - h*0.5) + lenBetweenNodes*0.5) / lenBetweenNodes);
			ff = 1;
		}
		Vect2D f = VECT2D_0;
		double lenSpring = 0;
		for (int jj = 0;jj < as->at(cr_i).arrNumBond.size();jj++)
		{
			int iii = as->at(cr_i).arrNumBond.at(jj);
			if (!(l->at(iii).pressure))
			{
				int at1 = l->at(iii).i;
				int at2 = l->at(iii).j;
				if (cr_i == at2) {
					int pr = at2;
					at2 = at1;
					at1 = pr;
				}
				Vect2D prVect = -(as->at(at2).r - as->at(at1).r);
				prVect /= prVect.Abs();
				lenSpring = LenConnection(as->at(at1), as->at(at2)) - sqrt(l->at(iii).eq_st_sq);
				f += k_bond*lenSpring*prVect;
			}
		}
		if ((curNum > numNodes - 1) || (ff > 1))
		{
			system("pause");
		}
		else
		{
			Press->at(ff*numNodes + curNum) += f.x;
			arr_number[ff][curNum]++;
		}
	}
	ofstream write_press;
	string namePressOut = "PressTest.txt";
	write_press.open(namePressOut, std::ofstream::out | std::ofstream::app);
	ofstream PrWidth;
	string namePrWidth = "PrWithd.txt";
	PrWidth.open(namePrWidth, std::ofstream::out | std::ofstream::app);
	//PrWidth << "Step " << stepp + 1 << "\n";
	if (stepp>100){
		write_press << "111 " << stepp++ << "\n";
		for (int pr_out = 0;pr_out < numNodes;pr_out++)
		{
			write_press << Press->at(pr_out) << "  " << Press->at(numNodes + pr_out) << "\n";
		}
		
	}
	
	
	for (int i11 = 0;i11 < 2;i11++)
	{
		for (int ipr = 0;ipr < numNodes;ipr++)
		{
			//widthFracture->at(i11*numNodes + ipr);//[i11][ipr];
			if (arr_number[i11][ipr] != 0)
			{
				Press->at(i11*numNodes + ipr) /= arr_number[i11][ipr];
			}
			else { Press->at(i11*numNodes + ipr) = Press->at(i11*numNodes + ipr - 1); }
		}
	}
	double pr = 0;
	if (Press->at(0) > Press->at(1))
	{
		pr = Press->at(0);
	}
	else
	{
		pr = Press->at(1);
	}
	if (pr < Press->at(2))
	{
		pr = Press->at(2);
	}
	Press->at(0) = 1.1*pr;
	Press->at(numNodes) = 1.1*pr;
	if (stepp>100) {
		write_press << "2222 " << stepp++ << "\n";
		for (int pr_out = 0;pr_out < numNodes;pr_out++)
		{
			write_press << Press->at(pr_out) << "  " << Press->at(numNodes + pr_out) << "\n";
		}
	}
	for (int i11 = 0;i11 < 2;i11++)
	{
		delete[] arr_number[i11];
	}
	delete[] arr_number;
	int ppp = stepp;
	bool prpp = (stepp > 101);
	
	if (prpp)
	{
		write_press << "3333 " << stepp++ << "\n";
		for (int pr_out = 0;pr_out < numNodes;pr_out++)
		{
			write_press << Press->at(pr_out) << "  " << Press->at(numNodes + pr_out) << "\n";
		}
		AproximationPress(Press, lenBetweenNodes,numNodes);
		write_press << "4444 " << stepp++ << "\n";
		for (int pr_out = 0;pr_out < numNodes;pr_out++)
		{
			write_press << Press->at(pr_out) << "  " << Press->at(numNodes + pr_out) << "\n";
		}
		StepRungeKutta(widthFracture, Press, Visc, numNodes, lenBetweenNodes, x_c0, x_c1, dt, v_cur);
	}
	write_press.close();
	double E_ung = 2 * sqrt(3)*k_bond / 3;
	double time_s = pow(60, gamma);
	double w_n = pow((time_s*E_ung / visc_liq), (1 / (gamma + 2)));
	for (int ipr = 0;ipr < numNodes;ipr++)
	{
		if (stepp <100)
		{
		prW[0][ipr] = widthFracture->at(ipr)/100;//-prW[0][ipr];
		prW[1][ipr] = widthFracture->at(numNodes + ipr)/100;// -prW[1][ipr];
		}
		else if (ipr==0)
		{
			prW[0][ipr] = abs(widthFracture->at(ipr) - prW[0][ipr]);
			prW[1][ipr] = abs(widthFracture->at(numNodes + ipr)-prW[1][ipr]);
			
		}
		else
		{
			prW[0][ipr] = abs(widthFracture->at(ipr) - prW[0][ipr]);
			prW[1][ipr] = abs(widthFracture->at(numNodes + ipr) - prW[1][ipr]);
		}
		if (stepp > 101) { PrWidth << prW[0][ipr] << "  " << prW[0][ipr] << "  "; }
	}
	if (stepp > 101) { PrWidth << "\n"; }
	PrWidth.close();
	for (int i = 0;i < PressBond->size();i++)//задаю перемещения для частиц с давлением
	{
		int num = PressBond->at(i);
		int cr_i = l->at(num).i;
		int cr_j = l->at(num).j;
		double y_i = double(as->at(cr_i).r.y);
		double y_j = double(as->at(cr_j).r.y);
		double x_con = double(as->at(cr_i).r.x + as->at(cr_j).r.x)*0.5;
		int curNum = 0;
		if (y_i-h*0.5 > 0)
		{
			curNum = int((y_i - h*0.5 + lenBetweenNodes*0.5) / lenBetweenNodes);
			ff = 0;
			if (curNum<0)
			{
				system("pause");
			}
			if (curNum>numNodes)
			{
				system("pause");
			}
		}
		else
		{
			curNum = int((abs(y_i - h*0.5) + lenBetweenNodes*0.5) / lenBetweenNodes);
			ff = 1;
			if (curNum<0)
			{
				system("pause");
			}
			if (curNum>numNodes)
			{
				system("pause");
			}
		}
		if (as->at(cr_i).r.x > as->at(cr_j).r.x) {
			as->at(cr_i).r.x += (prW[ff][curNum]*0.5);
		}
		else { as->at(cr_i).r.x -= (prW[ff][curNum]*0.5); }
		if (y_j - h*0.5 > 0)
		{
			curNum = int((y_j - h*0.5 + lenBetweenNodes*0.5) / lenBetweenNodes);
			ff = 0;
			if (curNum<0)
			{
				system("pause");
			}
			if (curNum>numNodes)
			{
				system("pause");
			}
		}
		else
		{
			curNum = abs(int((abs(y_j - h*0.5)+ lenBetweenNodes*0.5) / lenBetweenNodes));
			if (curNum<0)
			{
				system("pause");
			}
			if (curNum>numNodes)
			{
				system("pause");
			}
			ff = 1;
		}
		if (as->at(cr_j).r.x > as->at(cr_i).r.x) {
			as->at(cr_j).r.x += (prW[ff][curNum])*0.5;
		}
		else { as->at(cr_j).r.x -= (prW[ff][curNum])*0.5; }
		
	}
	double absForce = 0;
	for (int kk = 0; kk < as->size();kk++)//проходим по частицам
	{
		//std::cout << as->at(kk).numSpringBond<<"\n";
		if (as->at(kk).numSpringBond>0)
		{

			for (int j1 = 0;j1 < as->at(kk).arrNumBond.size();j1++)//проходим по номерам связей имеющихся у частицы
			{
				int jj = as->at(kk).arrNumBond.at(j1);//номер связи
				if (!(l->at(jj).cheek))
				{
					int numSecPart = l->at(jj).j;
					if (numSecPart == kk)
					{
						numSecPart = l->at(jj).i;
					}
					dr = as->at(kk).r - as->at(numSecPart).r;

					length = sqrt(dr.Sqr()) - sqrt(l->at(jj).eq_st_sq); //растяжение-сжатие пружины
					if ((!l->at(jj).pressure) && (abs(length) > 1e-6)) {
						int pto = 0;//system("pause");
					}
					if (l->at(jj).pressure)
					{
						absForce = -press;
						if (ff != 0) {
							//system("pause");
						}
///////////////////////////as->at(kk).f -= (ff / sqrt(dr.Sqr()))*dr;////////////////////////////////////////////////////////////////////
///////////////////////////as->at(numSecPart).f += (ff / sqrt(dr.Sqr()))*dr;////////////////////////////////////////////////////////////
					}
					else
						if (dr.Sqr() > pow(l->at(jj).eq_st_sq*e0, 1))//проверка на состоянии пружины
						{//случай разрыва связи
							double xCurrBond = (as->at(kk).r.x + as->at(numSecPart).r.x)*0.5;
							double yCurrBond = (as->at(kk).r.y + as->at(numSecPart).r.y)*0.5;
							bool f = true;

							/*int numgridX = int((xCurrBond + periodic_border) / bond);//номер ячейки по х куда попала частица
							int numgridY = int((yCurrBond + periodic_border) / bond);
							int arrNumCellX[3];//массив номеров ячеек проверяем по х
							int arrNumCellY[3];
							int X_num = 3;//количество ячеек проверяем по х
							int Y_num = 3;//количество ячеек проверяем по н
							if (numgridX == numCellsX - 1) {
								//arrNumCellX = new int [2];
								arrNumCellX[0] = numgridX - 1;
								arrNumCellX[1] = numgridX;
								X_num = 2;
							}
							else if (numgridX == 0) {
								//arrNumCellX = new int [2];
								arrNumCellX[0] = 0;
								arrNumCellX[1] = 1;
								X_num = 2;
							}
							else {
								//arrNumCellX = new int [3];
								arrNumCellX[0] = numgridX - 1;
								arrNumCellX[1] = numgridX;
								arrNumCellX[2] = numgridX + 1;
							}
							if (numgridY == numCellsY - 1) {
								//arrNumCellY = new int [2];
								arrNumCellY[0] = numgridY - 1;
								arrNumCellY[1] = numgridY;
								Y_num = 2;
							}
							else if (numgridY == 0) {
								//arrNumCellY = new int [2];
								arrNumCellY[0] = 0;
								arrNumCellY[1] = 1;
								Y_num = 2;
							}
							else {
								//arrNumCellY = new int [3];
								arrNumCellY[0] = numgridY - 1;
								arrNumCellY[1] = numgridY;
								arrNumCellY[2] = numgridY + 1;
							}
							arrNumCellX[0] = numgridX - 1;
							arrNumCellX[1] = numgridX;
							arrNumCellX[2] = numgridX + 1;
							arrNumCellY[0] = numgridY - 1;
							arrNumCellY[1] = numgridY;
							arrNumCellY[2] = numgridY + 1;
							for (int i = 0; i < X_num;i++)
							{
								for (int j = 0; j < Y_num;j++)
								{
									int ii = arrNumCellX[i];
									int jj = arrNumCellY[j];
									if (gridPressBond[ii][jj])
									{
										f = true;
										break;
									}
								}
								if (f)
								{
									break;
								}
							}
							if (f)
							{
								for (int i = 0; i < X_num;i++)
								{
									for (int j = 0; j < Y_num;j++)
									{
										int ii = arrNumCellX[i];
										int jj = arrNumCellY[j];
										gridPressBond[ii][jj] = true;
									}
								}
							}*/
							/*int NN = gridPressBond[ii][jj].numberPart;
							for (int k = 0;k<NN;k++)
							{
							int nomer = gridPressBond[ii][jj].arrNumPart[k];
							int b_i = l->at(nomer).i;
							int b_j = l->at(nomer).j;
							double x = (as->at(b_i).r.x + as->at(b_j).r.x)*0.5;
							double y = (as->at(b_i).r.y + as->at(b_j).r.y)*0.5;
							double dist = (x - xCurrBond)*(x - xCurrBond) + (y - yCurrBond)*(y - yCurrBond);
							if (dist < 1)
							{
							f = true;
							gridPressBond[ii][jj].Add(jj);
							//
							break;
							}
							}*/
							/*if (f)
							{
							break;
							}
							}
							if (f)
							{
							break;
							}
							}*/
							//bool f = false;
							if ((true))
							{
								l->at(jj).pressure = true;
								int ii1 = l->at(jj).i;
								int jj1 = l->at(jj).j;
								as->at(ii1).numSpringBond--;
								as->at(jj1).numSpringBond--;
								as->at(ii1).pressure = true;
								as->at(jj1).pressure = true;
								PressBond->push_back(jj);
								numPressBond++;
								num_errase++;
								//cout << num_errase << "\n";
							}
							else
							{
								l->at(jj).free = true;
								int ii1 = l->at(jj).i;
								int jj1 = l->at(jj).j;
								as->at(ii1).numSpringBond--;
								as->at(jj1).numSpringBond--;
								//num_errase++;
								//cout << num_errase << "\n";
							}
							//l->erase(l->begin()+ii);
							//system("pause");
						}
						else if (!(l->at(jj).free))
						{
							absForce = length*k_bond;//it->koef;//модуль силы
							if (ff != 0) {
								//system("pause");
							}
							as->at(kk).f -= (absForce / sqrt(dr.Sqr()))*dr;
							as->at(numSecPart).f += (absForce / sqrt(dr.Sqr()))*dr;
						}
					l->at(jj).cheek = true;
				}
				else {
					l->at(jj).cheek = false;
				}
			}

		}
		else//Удаляю частицы которые только с давлением
		{
			for (int it_b = 0;it_b <as->at(kk).arrNumBond.size(); it_b++)//прохожу по связям которые имеет данная частица и удаляю их
			{
				l->at(it_b).pulled = true;//c этой связью больше не работаем
			}
			as->at(kk).pressure = true;
		}
	}
	for (int i = 0; i<as->size();i++)//интегрируем по каждой частице скороксти
	{
		if (!(as->at(i).im))
		{
			as->at(i).f -= as->at(i).v*visc_k;
		}
	}
	for (int i = 0; i<as->size();i++)//интегрируем по каждой частице скороксти
	{
		if (!(as->at(i).im))
		{
			as->at(i).v += as->at(i).f*dt;
		}
	}
	for (int i = 0; i<as->size();i++)//интегрируем перемещения каждой частице
	{
		
		
		if ((!(as->at(i).im))&&(!(as->at(i).pressure)))
		{
			as->at(i).r += as->at(i).v*dt;
			
		}
		as->at(i).f = VECT2D_0;
	}
	stepp++;
	//    std::cout<<"step finished \n";
	delete numBondNodes;
	for (int ipr3 = 0;ipr3 < 2;ipr3++)
	{
		delete[] prW[ipr3];
	}
	delete[]prW;
}
void AproximationPress(vector<double>* Press, double deltaX,int numNodes)
{
	int numCurNodes = 0;
	int numCurNodes2 = 0;
	int i = 0;
	int i2 = 0;
	bool NotZeroPress = true;
	bool NotZeroPress2 = true;
	while ((NotZeroPress) && (i < Press->size()))
	{
		if (Press->at(i)> 1e-5)
		{
			
			numCurNodes = i;

		}
		else {
			NotZeroPress = false;
		}
		i++;
	}
	while ((NotZeroPress2) && (i2 < Press->size()))
	{
		if (Press->at(i2+numNodes)> 1e-5)
		{

			numCurNodes2 = i2;

		}
		else {
			NotZeroPress2 = false;
		}
		i2++;
	}
	/*double A[3][3];
	double B[3];
	for (int j1 = 0;j1 < 3;j1++)
		for (int j2 = 0; j2 < 3;j2++)
		{
			A[j1][j2] = 0;
			B[j2] = 0;
		}
	for (int it = 0; it < numCurNodes;it++)
	{
		double x_i = deltaX*it;
		A[0][0] += pow(x_i, 4);
		A[0][1] += pow(x_i, 3);
		A[0][2] += pow(x_i, 2);
		A[1][1] += pow(x_i, 2);
		A[1][2] += x_i;
		B[0] += pow(x_i, 2)*Press->at(it);
		B[1] += x_i*Press->at(it);
		B[2] += Press->at(it);
	}
	A[1][0] = A[0][1];
	A[2][0] = A[0][2];
	A[2][1] = A[1][2];
	A[2][2] = numPressBond;
	double** matrix = new double*[3];
	for (int ii = 0;ii < 3;ii++)
	{
		matrix[ii] = new double[3];
	}
	for (int j1 = 0;j1 < 3;j1++)
		for (int j2 = 0; j2 < 3;j2++)
		{
			matrix[j1][j2] = A[j1][j2];
		}
	double delta = GetDeterminant(matrix);
	matrix[0][0] = B[0];
	matrix[1][0] = B[1];
	matrix[2][0] = B[2];
	double delta0 = GetDeterminant(matrix);
	for (int j1 = 0;j1 < 3;j1++)
		for (int j2 = 0; j2 < 3;j2++)
		{
			matrix[j1][j2] = A[j1][j2];
		}
	matrix[0][1] = B[0];
	matrix[1][1] = B[1];
	matrix[2][1] = B[2];
	double delta1 = GetDeterminant(matrix);
	for (int j1 = 0;j1 < 3;j1++)
		for (int j2 = 0; j2 < 3;j2++)
		{
			matrix[j1][j2] = A[j1][j2];
		}
	matrix[0][2] = B[0];
	matrix[1][2] = B[1];
	matrix[2][2] = B[2];
	double delta2 = GetDeterminant(matrix);
	if (abs(delta) < 1e-5) { system("pause"); }
	
	double a = delta0/ delta;
	double b = delta1 / delta;
	double c = delta2 / delta;
	for (int i2 = 0;i2 < 3;i2++)
	{
		delete[]matrix[i2];
	}
	delete[] matrix;*/
	double b = numCurNodes*deltaX;
	double a = Press->at(0) / sqrt(b);
	double b2 = numCurNodes2*deltaX;
	double a2 = Press->at(numNodes) / sqrt(b2);
	for (int it = 0;it < numCurNodes;it++)
	{
		double x_i = deltaX*it;
		Press->at(it) = a*pow((b - x_i),0.5);
		Press->at(it+numNodes) = a2*pow((b2 - x_i), 0.5);
	}
	for (int it = 0;it < numCurNodes2;it++)
	{
		double x_i = deltaX*it;
		Press->at(it + numNodes) = a2*pow((b2 - x_i), 0.5);
	}
	cout << a << "  " << a2 << "  " << b << "  " << b2 << "  ";
}
double GetDeterminant(double ** matrix)
{
	double res = matrix[0][0] * matrix[1][1] * matrix[2][2] + matrix[0][1] * matrix[1][2] * matrix[2][0] + matrix[1][0] * matrix[2][1] * matrix[0][2] -
		matrix[2][0] * matrix[1][1] * matrix[0][2] - matrix[1][0] * matrix[0][1] * matrix[2][2] - matrix[2][1] * matrix[1][2] * matrix[0][0];
	return res;
}
void Clear_all(vector<Atom2D>*as)
{
	for (int j = 0;j<as->size();j++)
	{

	}
}
void SavePressBondDisplacement(vector<Atom2D>* as, vector<connection>* l, vector<int>* PressBond)
{
	ofstream savePress;
	string name = "BondPress.txt";
	savePress.open(name, std::ofstream::out | std::ofstream::app);
	savePress << "step  " << stepp << "\n";
	for (int i1 = 0;i1 < PressBond->size();i1++)
	{
		int num = PressBond->at(i1);
		int a1 = l->at(num).i;
		int a2 = l->at(num).j;
		savePress << as->at(a1).r.x << " " << as->at(a1).r.y << "\n";
		savePress << as->at(a2).r.x << " " << as->at(a2).r.y << "\n";
	}
}
void extention_system(vector<Atom2D>*as,double e,int ax,double &w,double &h)
{
	switch (ax)
	{
	case axis_x:{
		for (int ii = 0; ii<as->size();ii++)
		{
			if(1){
				as->at(ii).r.x *= (1 + e);
			}
			
		}
		w *= (1 + e);

		break;
	}
	case axis_y:{
		for (int ii = 0; ii<as->size();ii++)
		{
			if (1) {
				as->at(ii).r.y *= (1 + e);
			}
				
		}
		h *= (1 + e);
		break;
	}
	case volume_ext:
		for (int ii = 0; ii<as->size();ii++)
		{
			if (!(as->at(ii).im))
			{
				/*as->at(ii).r *=(1+e);
				as->at(ii).r_obr_x *=(1+e);
				as->at(ii).r_obr_y *=(1+e);
				as->at(ii).r_obr_xy *=(1+e);
				if (as->at(ii).obr_x)
				{
					as->at(as->at(ii).number_x).r = as->at(ii).r+as->at(ii).r_obr_x;
				}
				if (as->at(ii).obr_y)
				{
					as->at(as->at(ii).number_y).r = as->at(ii).r+as->at(ii).r_obr_y;
				}
				if (as->at(ii).obr_xy)
				{
					as->at(as->at(ii).number_xy).r = as->at(ii).r+as->at(ii).r_obr_xy;
				}*/
			}
		}
		break;
	default:
		break;
	}


}
void dublicate_config(vector<Atom2D>*Initial,vector<Atom2D>*Final)
{
	Final->clear();
	for (int ii = 0;ii<Initial->size();ii++)
	{
		double x = Initial->at(ii).r.x;
		double y = Initial->at(ii).r.y;
		Final->push_back(*(new Atom2D(x, y, false)));
		Final->at(ii) = Initial->at(ii);
	}
}
void Save_StressField(vector<Atom2D>* P, vector<connection>* B)
{
	string name = "Sigma " + std::to_string(stepp) + ".txt";
	ofstream fileStress;
	fileStress.open(name);
	fileStress << P->size()<<"\n";
	double V_a = 1;
	double sigma[4] = { 0,0,0,0 };//xx xy yx yy
	for (int i = 0;i < 4;i++) {
		sigma[i] = 0;
	}
	double force = 0;
	for (int iA = 0;iA < P->size();iA++) {
		for (int iC = 0;iC < P->at(iA).arrNumBond.size();iC++) {
			double numBond = P->at(iA).arrNumBond.at(iC);
			int jA = B->at(numBond).j;
			if (jA == iA) {
				jA = B->at(numBond).j;
			}
			Vect2D dr = P->at(jA).r - P->at(iA).r;
			if (B->at(numBond).pressure) {
				force = press;
			}
			else if ((B->at(numBond).free)) {
				force = 0;
			}
			else {
				force = k_bond*(LenConnection(P->at(jA), P->at(iA))- sqrt(B->at(numBond).eq_st_sq));
			}
			Vect2D f = force*dr / dr.Abs();
			if (dr.Abs() < 1e-9) {
				cout << "i " << B->at(numBond).i << " j " << B->at(numBond).j << "\n";
				cout << "i " << iA << " j " << jA << "\n";
				cout << "BAD ABS " << dr.Abs() << "\n";
				system("pause");
			}
			sigma[0] += f.x*dr.x;
			sigma[1] += f.x*dr.y;
			sigma[2] += f.y*dr.x;
			sigma[3] += f.y*dr.y;
		}
		fileStress << P->at(iA).r.x << " " << P->at(iA).r.y << " " <<
			sigma[0] << " " << sigma[1] << " " << sigma[2] << " " << sigma[3] << "\n";
		for (int i = 0;i < 4;i++) {
			sigma[i] = 0;
		}
	}
	fileStress.close();
}
void new_free_fracture(vector<Atom2D>*P,vector<connection>*B, double x, double y, double alf_rad, double len)
{
	double x1 = x-0.5*len*cos(alf_rad);
	double y1 = y-0.5*len*sin(alf_rad);
	double x2 = x+0.5*len*cos(alf_rad);
	double y2 = y+0.5*len*sin(alf_rad);

	for (int p = B->size()-1;p>0;p--)//
	{
		int i1 = B->at(p).i;
		int j1 = B->at(p).j;
		double x_i = P->at(i1).r.x;
		double y_i = P->at(i1).r.y;
		double x_j = P->at(j1).r.x;
		double y_j = P->at(j1).r.y;
		if (intersection(x1,y1,x2,y2,x_i,y_i,x_j,y_j))
		{
			B->at(p).free = true;
			int ii = B->at(p).i;
			int jj = B->at(p).j;
			P->at(ii).numSpringBond--;
			P->at(jj).numSpringBond--;
		}
	}
}
bool intersection(double A_x,double A_y,double B_x, double B_y,double C_x, double C_y,double D_x, double D_y)
{
	bool f = false;
	Vect2D AB = Vect2D(float(B_x-A_x),float(B_y-A_y));
	Vect2D CD = Vect2D(float(D_x-C_x),float(D_y-C_y));
	Vect2D AC = Vect2D(float(C_x-A_x),float(C_y-A_y));
	Vect2D AD = Vect2D(float(D_x-A_x),float(D_y-A_y));
	Vect2D CA = -1*AC;
	Vect2D CB = Vect2D(float(B_x-C_x),float(B_y-C_y));
	if ((((AB%AC)*(AB%AD))<0)&&((CD%CA)*(CD%CB)<0))
	{
		f = true;
	}
	return f;

}
void add_rand_free_fractures(vector<Atom2D>*P,vector<connection>*B,int n,double l,double w, double h)
{
	uniform_real_distribution<double> unif(0,1);
    mt19937 mt(time(NULL));
	for (int i = 0;i<n;i++)
	{
		double x = unif(mt) * w;
		double y = unif(mt) * h;
		double alf = unif(mt) * pi-pi*0.5;
		new_free_fracture(P,B,x,y,alf,l);
	}
}
void output_pressured_area(vector<Atom2D>* P, vector<connection> *B, int j, int number_iter, double w, double h, double r) {
	vector<Atom2D> Particles_out;
	Particles_out.clear();
	for (int ii = 0;ii < P->size();ii++)
	{
		if ((P->at(ii).im))
		{
			Atom2D At_pr = P->at(ii);
			Particles_out.push_back(At_pr);
		}
	}
	cout << Particles_out.size()<<"\n";
	for (int it = 0;it < B->size();it++)
	{
		if ((B->at(it).pressure)|| (B->at(it).free))
		{
			int ai = B->at(it).i;
			int aj = B->at(it).j;
			Atom2D At_pr = P->at(ai);
			if (P->at(ai).Spring)
			{
				
				Particles_out.push_back(At_pr);
			}
			
			if (P->at(aj).Spring)
			{
				At_pr = P->at(aj);
				Particles_out.push_back(At_pr);
			}
		}
	}

	//
	//for (int k = 0; k < kol; k++) {

	//	//cout << particles_nums[k] << "\n";
	//	
	//	Atom2D At_pr = P->at(particles_nums[k]);
	//	
	//	Particles_out.push_back(At_pr);
	//}
	//delete[]particles_nums;
	char file_name[20];	sprintf(file_name, "a%06u%06u.a3r", j + 1, number_iter);
	Save_A3R(file_name, &Particles_out, Particles_out.size(), float(r) / 2, w, h);
	Particles_out.clear();
	//particles_nums.clear();
}
double FindFractureSquare(vector<Atom2D>* P, vector<connection> *B, vector<int> *PressBond, double bond_len)
{
	std::vector<int>::iterator it;
	double Square = 0;
	double length = numPressBond * bond_len;
	double width = 0;
	double dr = 0;
	/*for (it = PressBond->begin();it != PressBond->end();it++)
	{
		int i = B->at(*it).i;
		int j = B->at(*it).j;
		if ((P->at(i).numSpringBond>0) && (P->at(j).numSpringBond > 0))
		{
			dr = sqrt(pow((P->at(i).r.x - P->at(j).r.x), 2) + pow((P->at(i).r.y - P->at(j).r.y), 2)) / numPressBond;
			width += dr;
		}
		
	}*/
	return length*length;
}
double LenConnection(Atom2D A, Atom2D B)
{
	double len = sqrt(pow((A.r.x - B.r.x), 2) + pow((A.r.y - B.r.y), 2));
	return len;
}

double LenBetweenConnection(Atom2D A1, Atom2D A2, Atom2D B1, Atom2D B2)
{
	double xA = A1.r.x + A2.r.x;
	double xB = B1.r.x + B2.r.x;
	double yA = A1.r.y + A2.r.y;
	double yB = B1.r.y + B2.r.y;
	return sqrt(pow(xA-xB,2)+pow(yA-yB,2));
}

void CreateAvtoModelSolution(vector<Atom2D>* P, vector<connection>* B, vector<int>* PressBond, double LenFracture, double * w_avt,  int numNodes,double h, double w)
{
	double dy = LenFracture / (numNodes - 1);
	int* masCr = new int[2*PressBond->size()];
	int kol = 0;
	for (int i = 0;i < PressBond->size();i++)
	{
		/*int numBond = PressBond->at(i);
		int numNode = 0;
		int crI = B->at(numBond).i;
		int crJ = B->at(numBond).j;
		for (int jj = 0; jj < kol;jj++)
		{
			if (masCr[jj] == crI)
			{
				crI = -1;
			}
			if (masCr[jj] == crJ)
			{
				crJ = -1;
			}
		}
		if (crI)
		{
			masCr[kol] = crI;
			kol++;
		}
		if (crJ)
		{
			masCr[kol] = crJ;
			kol++;
		}*/
		/*if ((crI>0)&&(abs(P->at(crI).r.y - h*0.5) < LenFracture*0.5))
		{
			numNode = int((abs(P->at(crI).r.y - h*0.5) + dy*0.5) / dy);
			if (P->at(crI).r.x >w*0.5)
			{
				P->at(crI).r.x += w_avt[numNode];
			}
			else
			{
				P->at(crI).r.x -= w_avt[numNode];
			}
		}
		else if ((crJ>0)&&(abs(P->at(crJ).r.y - h*0.5) < LenFracture*0.5))
		{
			numNode = int((abs(P->at(crJ).r.y - h*0.5) + dy*0.5) / dy);
			if (P->at(crJ).r.x > w*0.5)
			{
				P->at(crJ).r.x += w_avt[numNode];
			}
			else
			{
				P->at(crJ).r.x -= w_avt[numNode];
			}
		}*/
		
	}
	delete[]masCr;
}
void StepManyFrackPress(vector<Atom2D>* as, vector<connection>* l, bool ** gridPressBond, int numCellsX, int numCellsY, double w, double h, double bond, double periodic_border, vector<int>* PressBond,vector<int>*FrackPress, double * masPress, int numFrack)
{
	std::vector<connection>::iterator it;
	double length;
	double ff;
	Vect2D dr;
	Vect2D ww = Vect2D(float(w), 0);
	Vect2D hh = Vect2D(0, float(h));
	Vect2D wh = Vect2D(float(w), float(h));
	for (int kk = 0; kk < as->size();kk++)//проходим по частицам
	{
		for (int j1 = 0;j1 < as->at(kk).arrNumBond.size();j1++)//проходим по номерам связей имеющихся у частицы
		{
			int jj = as->at(kk).arrNumBond.at(j1);//номер связи
			if (!(l->at(jj).cheek))
			{
				int numSecPart = l->at(jj).j;
				if (numSecPart == kk)
				{
					numSecPart = l->at(jj).i;
				}
				dr = as->at(kk).r - as->at(numSecPart).r;
				length = sqrt(dr.Sqr()) - sqrt(l->at(jj).eq_st_sq); //растяжение-сжатие пружины
				if (l->at(jj).pressure)
				{
					bool fl = true;
					int numFrack = -1;
					int it = 0;
					while (fl) {
						if (jj == FrackPress->at(it)) {
							fl = false;
							numFrack = FrackPress->at(it + 1);
						}
						else {
							it += 2;
						}
					}
					ff = -masPress[numFrack];
					as->at(kk).f -= (float(ff) / sqrt(dr.Sqr()))*dr;
					as->at(numSecPart).f += (float(ff) / sqrt(dr.Sqr()))*dr;
				}
				else
					if (dr.Sqr() > pow(l->at(jj).eq_st_sq*e0, 1))//проверка на состоянии пружины
					{//случай разрыва связи
						double xCurrBond = (as->at(kk).r.x + as->at(numSecPart).r.x)*0.5;
						double yCurrBond = (as->at(kk).r.y + as->at(numSecPart).r.y)*0.5;
						bool f = false;

						int numgridX = int((xCurrBond + periodic_border) / bond);//номер ячейки по х куда попала частица
						int numgridY = int((yCurrBond + periodic_border) / bond);
						int arrNumCellX[3];//массив номеров ячеек проверяем по х
						int arrNumCellY[3];
						int X_num = 3;//количество ячеек проверяем по х
						int Y_num = 3;//количество ячеек проверяем по н
						if (numgridX == numCellsX - 1) {
							//arrNumCellX = new int [2];
							arrNumCellX[0] = numgridX - 1;
							arrNumCellX[1] = numgridX;
							X_num = 2;
						}
						else if (numgridX == 0) {
							//arrNumCellX = new int [2];
							arrNumCellX[0] = 0;
							arrNumCellX[1] = 1;
							X_num = 2;
						}
						else {
							//arrNumCellX = new int [3];
							arrNumCellX[0] = numgridX - 1;
							arrNumCellX[1] = numgridX;
							arrNumCellX[2] = numgridX + 1;
						}
						if (numgridY == numCellsY - 1) {
							//arrNumCellY = new int [2];
							arrNumCellY[0] = numgridY - 1;
							arrNumCellY[1] = numgridY;
							Y_num = 2;
						}
						else if (numgridY == 0) {
							//arrNumCellY = new int [2];
							arrNumCellY[0] = 0;
							arrNumCellY[1] = 1;
							Y_num = 2;
						}
						else {
							//arrNumCellY = new int [3];
							arrNumCellY[0] = numgridY - 1;
							arrNumCellY[1] = numgridY;
							arrNumCellY[2] = numgridY + 1;
						}
						arrNumCellX[0] = numgridX - 1;
						arrNumCellX[1] = numgridX;
						arrNumCellX[2] = numgridX + 1;
						arrNumCellY[0] = numgridY - 1;
						arrNumCellY[1] = numgridY;
						arrNumCellY[2] = numgridY + 1;
						if (true)//((f))
						{
							l->at(jj).pressure = true;
							int ii1 = l->at(jj).i;
							int jj1 = l->at(jj).j;
							bool fl = false;
							/*Определяю какая трещина ближайшая*/
							for (int ii = 0;ii < FrackPress->size();ii+=2) {
								int numBond = FrackPress->at(ii);
								int numFrack = FrackPress->at(ii+1);
								int at1 = l->at(numBond).i;
								int at2 = l->at(numBond).j;
								double dist = LenBetweenConnection(as->at(ii1), as->at(jj1), as->at(at1), as->at(at2));
								if (dist < 5 * bond) {//Находится достаточно близок трещина
									FrackPress->push_back(jj);
									FrackPress->push_back(numFrack);
									//cout << numFrack << "\n";
									break;
								}
								
							}

							
							as->at(ii1).numSpringBond--;
							as->at(ii1).pressure = true;
							as->at(jj1).numSpringBond--;
							as->at(jj1).pressure = true;
							PressBond->push_back(jj);
							numPressBond++;
							num_errase++;
							//cout << num_errase << "\n";
						}
						else
						{
							l->at(jj).free = true;
							int ii1 = l->at(jj).i;
							int jj1 = l->at(jj).j;
							as->at(ii1).numSpringBond--;
							as->at(jj1).numSpringBond--;

						}

					}
					else if (!(l->at(jj).free))
					{
						ff = length*k_bond;//it->koef;//модуль силы
						as->at(kk).f -= (float(ff) / sqrt(dr.Sqr()))*dr;
						as->at(numSecPart).f += (float(ff) / sqrt(dr.Sqr()))*dr;
					}
				l->at(jj).cheek = true;
			}
			else {
				l->at(jj).cheek = false;
			}
		}
	}
	if (typeBC == 1) {
		RecalculateVelosityPositionFixed(as, w, h);
	}
	else {
		RecalculateVelosityPositionPeriodic(as, w, h);
	}
	stepp++;
}
void StepRungeKutta(vector<double>* w, vector<double>* p, vector<double>* v, int numElements, double deltaX, double x_c0, double x_c1, double h, double* v_cur)
{
	double E_ung = 2 * sqrt(3)*k_bond / 3;
	double time_s = pow(60, gamma);
	double w_n = pow((time_s*E_ung / visc_liq), (1 / (gamma + 2)));
	double** K = new double*[4];//Массив коэффициентов для метода Рунге-кутта 4 порядка
	for (int i1 = 0;i1 < 4;i1++)
	{
		
			K[i1] = new double[numElements];
		
	}
	double ** W_y = new double*[2];//временный массив для решения системы дифуров
	for (int i1 = 0;i1 < 2;i1++)
	{
		W_y[i1] = new double[numElements];
	}
	if (stepp > 101)
	{
		for (int pp = 0;pp < 2;pp++)
		{
			for (int j = 0;j < numElements-1;j++)
			{
				if (j == 0)
				{
					K[0][0] = h*(3 * q_0 - (4 * w->at(pp*numElements+1) * v->at(pp*numElements + 1) - w->at(pp*numElements + 2) * v->at(pp*numElements + 2))) / 2 / deltaX;
					//cout << "q  " << q << "deltaX " << deltaX << " h " << h << " " << K[0][0] << " " << K[0][0] << "\n";
				}
				else if (j == 1)
				{
					K[0][j] = h*(-w->at(pp*numElements +j+ 1) * v->at(pp*numElements +j+ 1) + q_0) / 2 / deltaX;
				}
				else
				{
					K[0][j] = h*(-w->at(pp*numElements + j + 1) * v->at(pp*numElements + j + 1) + w->at(pp*numElements + j - 1) * v->at(pp*numElements + j - 1)) / 2 / deltaX;
				}
			}
			cout << K[0][numElements - 1] << "  " << K[0][numElements - 1] << "\n";
			for (int j = 0;j < numElements-1;j++)
			{
				if (j == 0)
				{
					K[1][0] = h*(3 * q_0 - (4 * (w->at(pp*numElements+1) + 0.5*h*K[0][1]) * v->at(pp*numElements+1) - (w->at(pp*numElements + 2) + 0.5*h*K[0][2]) * v->at(pp*numElements + 2))) / 2 / deltaX;
				}
				else if (j == 1)
				{
					K[1][j] = h*(-(w->at(pp*numElements +j+ 1) + K[0][j + 1] * 0.5*h) * v->at(pp*numElements +j+ 1) + (q_0)) / 2 / deltaX;
				}
				else
				{
					K[1][j] = h*(-(w->at(pp*numElements + j + 1) + K[0][j + 1] * 0.5*h) * v->at(pp*numElements + j + 1) + (w->at(pp*numElements + j - 1) + 0.5*h*K[0][j - 1]) * v->at(pp*numElements + j - 1)) / 2 / deltaX;
				}
			}
			for (int j = 0;j < numElements-1;j++)
			{
				if (j == 0)
				{
					K[2][0] = h*(3 * q_0 - (4 * (w->at(pp*numElements + 1) + 0.5*h*K[1][1]) * v->at(pp*numElements  + 1) - (w->at(pp*numElements + 2) + 0.5*h*K[1][2]) * v->at(pp*numElements + 2))) / 2 / deltaX;
				}
				else if (j == 1)
				{
					K[2][j] = h*(-(w->at(pp*numElements+j + 1) + K[1][j + 1] * 0.5*h) * v->at(pp*numElements + j + 1) + q_0) / 2 / deltaX;
				}
				else
				{
					K[2][j] = h*(-(w->at(pp*numElements + j + 1) + K[1][j + 1] * 0.5*h) * v->at(pp*numElements + j + 1) + (w->at(pp*numElements + j - 1) + 0.5*h*K[1][j - 1]) * v->at(pp*numElements + j - 1)) / 2 / deltaX;
				}
			}
			for (int j = 0;j < numElements-1;j++)
			{
				if (j == 0)
				{
					K[3][0] = h*(3 * q_0 - (4 * (w->at(pp*numElements  + 1) + h*K[2][1]) * v->at(pp*numElements  + 1) - (w->at(pp*numElements + 2) + h*K[2][2]) * v->at(pp*numElements + 2))) / 2 / deltaX;
				}
				else if (j == 1)
				{
					K[3][j] = h*(-(w->at(pp*numElements + j + 1) + h* K[2][j + 1]) * v->at(pp*numElements + j + 1) + q_0) / 2 / deltaX;
				}
				else
				{
					K[3][j] = h*(-(w->at(pp*numElements + j + 1) + h* K[2][j + 1]) * v->at(pp*numElements + j + 1) + (w->at(pp*numElements + j - 1) + h*K[2][j - 1]) * v->at(pp*numElements + j - 1)) / 2 / deltaX;
				}
			}

			int curr0 = int((abs(x_c0)) / deltaX);
			int curr1 = int((abs(x_c1)) / deltaX);
			if (pp == 1) {
				curr0 = curr1;
			}
			for (int ii = 0;ii < 2;ii++)//определяем новые раскрытия в узлах
			{
				for (int j = 0;j < numElements - 1;j++)
				{
					if (j == 0)
					{
						W_y[pp][0] = w->at(pp*numElements + 0) + h*(K[0][0] + 2 * K[1][0] + 2 * K[2][0] + K[3][0]);
						w->at(pp*numElements + 0) = 0.1*W_y[pp][0];
					}
					else if ((j <= curr0))
					{
						W_y[pp][j] = w->at(pp*numElements + j) + 1e-1*h*(K[0][j] + 2 * K[1][j] + 2 * K[2][j] + K[3][j]);
						w->at(pp*numElements + j) = W_y[pp][j];
					}
					else {
						W_y[pp][j] = 0;//w[0][j] + h*(K[0][0][j] + 2 * K[0][1][j] + 2 * K[0][2][j] + K[0][3][j]);
						w->at(pp*numElements + j) = W_y[pp][j];
					}
				}
			}
		}
	/*for (int ii = 0;ii < 2;ii++)//определяем новые раскрытия в узлах
	{
		//w[0][0] = 0.5*(W_y[0][0] + W_y[1][0]);
		//w[1][0] = w[0][0];
		for (int j = 0;j < numElements - 1;j++)
		{
			//W_y[0][j] = w[0][j] + (K[0][0][j] + 2 * K[0][1][j] + 2 * K[0][2][j] + K[0][3][j]);
			//w[0][j] = W_y[0][j];
			if (j == 0)
			{
				W_y[1][0] = w[1][0] + h*(K[1][0][0] + 2 * K[1][1][0] + 2 * K[1][2][0] + K[1][3][0]);
				w[1][0] = W_y[1][0];
			}
			else if ((j <= curr1))
			{
				W_y[1][j] = w[1][j] + 1e-1*h*(K[1][0][j] + 2 * K[1][1][j] + 2 * K[1][2][j] + K[1][3][j]);
				w[1][j] = W_y[1][j];
			}
			else
			{
				W_y[1][j] = 0;//w[1][j] + h*(K[1][0][j] + 2 * K[1][1][j] + 2 * K[1][2][j] + K[1][3][j]);
				w[1][j] = (W_y[1][j]);
			}
		}
	}*/
	}
	/*ofstream fileFirst;
	string namefile = "K_and_width.txt";
	fileFirst.open(namefile);
	fileFirst << "deltaX = " << deltaX << "  numNodes" << numElements <<" h="<<h<<"  q="<<q<< "\n";
	fileFirst << "K0\n";
	for (int i3 = 0;i3 < numElements;i3++)
	{
		fileFirst << K[0][0][i3] << "  " << K[1][0][i3] << "\n";
	}
	fileFirst << "K1\n";
	for (int i3 = 0;i3 < numElements;i3++)
	{
		fileFirst << K[0][1][i3] << "  " << K[1][1][i3] << "\n";
	}
	fileFirst << "K2\n";
	for (int i3 = 0;i3 < numElements;i3++)
	{
		fileFirst << K[0][2][i3] << "  " << K[1][2][i3] << "\n";
	}
	fileFirst << "K3\n";
	for (int i3 = 0;i3 < numElements;i3++)
	{
		fileFirst << K[0][3][i3] << "  " << K[1][3][i3] << "\n";
	}
	fileFirst << "width\n";
	for (int i3 = 0;i3 < numElements;i3++)
	{
		fileFirst << w[0][i3] << "  " << w[1][i3] << "\n";
	}
	fileFirst.close();*/
	double p0 = pow((abs(p->at(0) - p->at(1)) / deltaX), 1/gamma);
	double p1 = pow((abs(p->at(numElements + 0) - p->at(numElements + 1)) / deltaX), 1/gamma);
	v->at(0) = (p0)*w_n*1e-2;//;*E_ung;
	v->at(numElements) = v->at(0);
	int curNode0 = int((abs(x_c0))/ deltaX);
	int curNode1 = int((abs(x_c1)) / deltaX);
	// = deltaX*(numCurentPoint - 1) + 0.001;
	double alfa = 2 / (gamma + 2);
	double B = 0.25*alfa / tan(pi*(1 - alfa));
	double A_u = pow(((1 - alfa)*B), (-0.5*alfa));
	if (curNode0>0)
	{
		for (int j = 1;j < curNode0-1;j++)
		{
			double mnoj = 10*pow(w_n, (gamma + 1) / gamma)*pow(w_n*E_ung, 1 / gamma);
 			double df = (pow((x_c0 - deltaX*j), (alfa*(gamma + 1))) / (alfa*(gamma + 1) - 1))*((1 / (pow((x_c0 - deltaX*(j + 1)), (alfa*(gamma + 1) - 1))) - (1 / (pow((x_c0 - deltaX*(j - 1)), (alfa*(gamma + 1) - 1))))));
			v->at(j) = mnoj*pow((pow(w->at(j), gamma + 1)*abs((p->at(j - 1) - p->at(j + 1))) / df), 1 / gamma);
		}
		double mm = 1*pow(w_n, 1 / (1 - alfa));
		double d_x = x_c0 - deltaX*(curNode0);
		if (abs(d_x) < 1e-5) {
			d_x = x_c0 - deltaX*(curNode0 - 1);
		}
		v_cur[0] = mm* abs(pow(A_u, (-1 / (1 - alfa)))*pow(abs(w->at(curNode0)) / pow((d_x), alfa), (1 / (1 - alfa))));//скорость распространения конца трещины
		if (v_cur[0] < 1e-5)
		{
			cout << w->at(curNode0) << " " << x_c0 << " " << deltaX*(curNode0 - 1) << "\n";
			//system("pause");
		}
		cout << w->at(curNode0 - 1) << "  " << w->at(curNode0) << " \n ";
		double a_v = -gamma / (gamma + 1) / (gamma + 4);
		v->at(curNode0-1) = 1e-1*v_cur[0]*(1 + a_v*(1 - deltaX*(curNode0 - 1) / x_c0));
	}
	else
	{
		v_cur[0] = 1e-2;//1e10*pow(A_u, (-1 / (1 - alfa)))*pow(w[0][0] / pow((x_c0 - deltaX*(0)), alfa), (1 / (1 - alfa)));//скорость распространения конца трещины
		cout << w[0][0] << "  " << x_c0 <<"  "<< v_cur[0] <<"  \n";
		//system("pause");
																											   //v[0][0] = 0;
	}
	if(curNode1>0)
	{
		for (int j = 1;j < curNode1 - 1;j++)
		{
			double mnoj = 10*pow(w_n, (gamma + 1) / gamma)*pow(w_n*E_ung, 1 / gamma);
			double df = (pow((x_c1 - deltaX*j), (alfa*(gamma + 1))) / (alfa*(gamma + 1) - 1))*((1 / (pow((x_c1 - deltaX*(j + 1)), (alfa*(gamma + 1) - 1))) - (1 / (pow((x_c1 - deltaX*(j - 1)), (alfa*(gamma + 1) - 1))))));
			v->at(numElements+j) = mnoj*pow((pow(w->at(numElements + j), gamma + 1)*abs(p->at(numElements + j - 1) - p->at(numElements + j + 1)) / df), 1 / gamma);
		}
		double mm = 1*pow(w_n, 1 / (1 - alfa));
		double d_x = x_c1 - deltaX*(curNode1);
		if (abs(d_x) < 1e-5) {
			d_x = x_c1 - deltaX*(curNode1 - 1);
		}
		v_cur[1] = mm*pow(A_u, (-1 / (1 - alfa)))*pow(w->at(numElements + curNode1)/ pow((d_x), alfa), (1 / (1 - alfa)));//скорость распространения конца трещины
		
		if (v_cur[1] < 1e-5)
		{
			cout << w->at(numElements + curNode1) << " " << x_c0 << " " << deltaX*(curNode1 - 1) << "\n";
			//system("pause");
		}
		double a_v = -gamma / (gamma + 1) / (gamma + 4);
		cout << w->at(numElements + curNode1 - 1) << "  " << w->at(numElements + curNode1) <<"  "<< v_cur[1]<< " \n ";
		v->at(numElements + curNode1-1) = 1e-1*v_cur[1]*(1 + a_v*(1 - deltaX*(curNode1) / x_c0));
	}
	else
	{
		v_cur[1] = 1e-2;//1e10*pow(A_u, (-1 / (1 - alfa)))*pow(w[1][0] / pow((x_c1 - deltaX*(0)), alfa), (1 / (1 - alfa)));;// pow(A_u, (-1 / 1 - alfa))*pow(w[1][curNode1] / pow((x_c1 - deltaX*(curNode1)), alfa), (1 / (1 - alfa)));//скорость распространения конца трещины
		cout << w->at(numElements + 0) << "  " << x_c1 << "  " << v_cur[1] << "  \n";
		//system("pause");																									//v[1][0] = 0;
	}
	if ((abs(curNode0) > 0) || (abs(curNode1) > 0))
	{
		
		//cout << v[0][curNode0] << "   " << v[0][curNode1] << "\n";
		//system("pause");
	}
	/*for (int i = 0;i < numElements;i++)
	{
		p->at(i) *= (w_n*E_ung);
		w->at(i) *= w_n;
	}*/
	ofstream velos;
	ofstream widt;
	ofstream Press;
	
	string nameVel = "Velosity.txt";
	string nameWit = "Widht.txt";
	string namePress = "Press.txt";
	
	Press.open(namePress, std::ofstream::out | std::ofstream::app);
	velos.open(nameVel,std::ofstream::out | std::ofstream::app);
	widt.open(nameWit, std::ofstream::out | std::ofstream::app);
	/*Press<<"Step " << stepp + 1 << "\n";
	velos << "Step " << stepp + 1 << "\n";
	widt << "Step " << stepp + 1 << "\n";
	*/
	for (int i0 = 0;i0 < numElements;i0++)
	{
		
		velos << v->at(i0) << "  " << v->at(numElements + i0)<<"  ";// << "\n";
		widt << w->at(i0) << "  " << w->at(numElements + i0) << "  ";// << "\n";
		Press << p->at(i0) << "   " << p->at(numElements + i0) << "  ";// << "\n";
	}
	velos << "  \n";//v->at(i0) << " " << v->at(numElements + i0);// << "\n";
	widt << " \n";//<< w->at(i0) << " " << w->at(numElements + i0);// << "\n";
	Press << " \n";//<< p->at(i0) << " " << p->at(numElements + i0);// << "\n";
	velos.close();
	widt.close();
	

	for (int j = 0;j < 4;j++)
	{
		delete[] K[j];
	}
	delete[] K;
	

	for (int i = 0;i < 2;i++)
	{
		delete[] W_y[i];
	}
	delete[]W_y;
}
void RecalculateVelosityPositionFixed(vector<Atom2D>* as,double w, double h)
{
	std::vector<connection>::iterator it;
	Vect2D dr;
	Vect2D ww = Vect2D(w, 0);
	Vect2D hh = Vect2D(0, h);
	Vect2D wh = Vect2D(w, h);
	for (int i = 0; i<as->size();i++)//интегрируем по каждой частице скороксти
	{
		if (1)//((!(as->at(i).im))||((as->at(i).Spring)))
		{
			as->at(i).f -= as->at(i).v*visc_k;
		}
	}
	for (int i = 0; i<as->size();i++)//интегрируем по каждой частице скороксти
	{
		if (1)//((!(as->at(i).im)) || ((as->at(i).Spring)))
		{
			as->at(i).v += as->at(i).f*dt;
		}
	}
	for (int i = 0; i<as->size();i++)//интегрируем перемещения каждой частице
	{
		
		if ((!(as->at(i).im)))//if ((!(as->at(i).im)) || ((as->at(i).Spring)))
		{
			as->at(i).r += as->at(i).v*dt;
			double ul = abs(0 - as->at(i).r.x), ur = w - as->at(i).r.x, ut = as->at(i).r.y, ub = h - as->at(i).r.y;
			if ((ul<w*part_fixed) || (ut<h*part_fixed) || (ur<w*part_fixed) || (ub<h*part_fixed))//if ((ul<10*part_fixed) || (ut<10*part_fixed) || (ur<10*part_fixed) || (ub<10*part_fixed))
			//if ((ul<w*part_fixed) || (ur<w*part_fixed))
			//if ((ut<h*part_fixed) || (ub<h*part_fixed))
			//if ((ub<h*part_fixed))
			{
				as->at(i).im = true;
			}
			
		}
		else {
			/*double ub = h - as->at(i).r.y;//моделирую ниднюю подложку(образец стоит на столе)
			if (ub < w*part_fixed) {
					as->at(i).r.x += as->at(i).v.x*dt;
			}
			else {
				as->at(i).r += as->at(i).v*dt;
			}*/
			
		}
		
		as->at(i).f = VECT2D_0;
	}
	
}
void RecalculateVelosityPositionPeriodic(vector<Atom2D>* as,double w, double h)
{
	std::vector<connection>::iterator it;
	Vect2D dr;
	Vect2D ww = Vect2D(w, 0);
	Vect2D hh = Vect2D(0, h);
	Vect2D wh = Vect2D(w, h);
	for (int i = 0; i<as->size();i++)//интегрируем по каждой частице скороксти
	{
		//if (!(as->at(i).im))
		{
			as->at(i).f -= as->at(i).v*visc_k;
		}
	}
	for (int i = 0; i<as->size();i++)//интегрируем по каждой частице скороксти
	{
		if ((!(as->at(i).im)) || ((as->at(i).Spring)))
		{
			as->at(i).v += as->at(i).f*dt;
		}
	}
	for (int i = 0; i<as->size();i++)//интегрируем по каждой частице скороксти
	{
		if ((!(as->at(i).im)))
		{
			as->at(i).r += as->at(i).v*dt;

			if (((as->at(i).obr_x) || (as->at(i).obr_y) || (as->at(i).obr_xy)) && (as->at(i).r.x < 0 || as->at(i).r.x > w || as->at(i).r.y < 0 || as->at(i).r.y > h))//проверяем стала ли частица мнимой
			{
				if (as->at(i).obr_x)//пересчет координат
				{
					as->at(as->at(i).number_x).r = as->at(i).r + as->at(i).r_obr_x*ww;
				}
				if (as->at(i).obr_y)
				{
					as->at(as->at(i).number_y).r = as->at(i).r + as->at(i).r_obr_y*hh;
				}
				if (as->at(i).obr_xy)
				{
					int px = 0;
					int py = 0;
					if (as->at(i).r_obr_xy == 22)
					{
						px = 1;
						py = 1;
					}
					if (as->at(i).r_obr_xy == 12)
					{
						px = -1;
						py = 1;
					}
					if (as->at(i).r_obr_xy == 21)
					{
						px = 1;
						py = -1;
					}
					if (as->at(i).r_obr_xy == 11)
					{
						px = -1;
						py = -1;
					}
					Vect2D pr = wh;
					pr.x *= px;
					pr.y *= py;
					as->at(as->at(i).number_xy).r = as->at(i).r + pr;
				}
				as->at(i).im = true;
				if ((as->at(i).obr_x) && (!(as->at(as->at(i).number_x).r.x < 0 || as->at(as->at(i).number_x).r.x > w || as->at(as->at(i).number_x).r.y < 0 || as->at(as->at(i).number_x).r.y > h)))//проверяем стал ли образ по оси х ральной частицей
				{
					as->at(as->at(i).number_x).im = false;//меняем образ х который теперь стал реальной с рассматриываемой частицей, которая стала теперь мнимой
					as->at(as->at(i).number_x).obr_x = true;
					as->at(as->at(i).number_x).number_x = i;
					if ((as->at(i).r - as->at(as->at(i).number_x).r).x > 0)
					{
						as->at(as->at(i).number_x).r_obr_x = 1;
					}
					else {
						as->at(as->at(i).number_x).r_obr_x = -1;
					}
					//as->at(as->at(i).number_x).r_obr_x = as->at(i).r - as->at(as->at(i).number_x).r;
					if (as->at(i).obr_y)//у рассматриваемой частицы 3 образа, оставшиеся меняем местами
					{
						//if ((abs(as->at(i).number_y)>2000)||(abs(as->at(i).number_xy)>2000)){system("pause");}
						as->at(as->at(i).number_x).number_xy = as->at(i).number_y;
						as->at(as->at(i).number_x).number_y = as->at(i).number_xy;
						as->at(as->at(i).number_x).obr_xy = true;
						as->at(as->at(i).number_x).obr_y = true;
						int px, py;
						if ((as->at(as->at(i).number_y).r - as->at(as->at(i).number_x).r).x > 0)
						{
							px = 2;
						}
						else {
							px = 1;
						}
						if ((as->at(as->at(i).number_y).r - as->at(as->at(i).number_x).r).y > 0)
						{
							py = 2;
						}
						else {
							py = 1;
						}
						as->at(as->at(i).number_x).r_obr_xy = px * 10 + py;//as->at(as->at(i).number_y).r - as->at(as->at(i).number_x).r;
						if (py == 2)
						{
							as->at(as->at(i).number_x).r_obr_y = 1;
						}
						else
						{
							as->at(as->at(i).number_x).r_obr_y = -1;
						}
						//as->at(as->at(i).number_x).r_obr_y = as->at(as->at(i).number_xy).r - as->at(as->at(i).number_x).r;
					}
					as->at(i).obr_x = false;
					as->at(i).obr_y = false;
					as->at(i).obr_xy = false;

				}
				else if ((as->at(i).obr_y) && (!(as->at(as->at(i).number_y).r.x < 0 || as->at(as->at(i).number_y).r.x > w || as->at(as->at(i).number_y).r.y < 0 || as->at(as->at(i).number_y).r.y > h)))//проверяем стал ли образ по оси y ральной частицей
				{
					as->at(as->at(i).number_y).im = false;//меняем образ y который теперь стал реальной с рассматриываемой частицей, которая стала теперь мнимой
					as->at(as->at(i).number_y).obr_y = true;
					as->at(as->at(i).number_y).number_y = i;
					if ((as->at(i).r - as->at(as->at(i).number_y).r).y > 0)
					{
						as->at(as->at(i).number_y).r_obr_y = 1;
					}
					else {
						as->at(as->at(i).number_y).r_obr_y = -1;
					}
					//as->at(as->at(i).number_y).r_obr_y = as->at(i).r - as->at(as->at(i).number_y).r;
					if (as->at(i).obr_x)//у рассматриваемой частицы 3 образа, оставшиеся меняем местами
					{
						//if ((abs(as->at(i).number_x)>2000)||(abs(as->at(i).number_xy)>2000)){system("pause");}
						as->at(as->at(i).number_y).number_xy = as->at(i).number_x;
						as->at(as->at(i).number_y).number_x = as->at(i).number_xy;
						as->at(as->at(i).number_y).obr_xy = true;
						as->at(as->at(i).number_y).obr_x = true;
						int px, py;
						if ((as->at(as->at(i).number_x).r - as->at(as->at(i).number_y).r).x > 0)
						{
							px = 2;
						}
						else {
							px = 1;
						}
						if ((as->at(as->at(i).number_x).r - as->at(as->at(i).number_y).r).y > 0)
						{
							py = 2;
						}
						else {
							py = 1;
						}
						as->at(as->at(i).number_y).r_obr_xy = px * 10 + py;//as->at(as->at(i).number_x).r - as->at(as->at(i).number_y).r;
						if (px == 2)
						{
							as->at(as->at(i).number_y).r_obr_x = 1;
						}
						else {
							as->at(as->at(i).number_y).r_obr_x = -1;
						}
						//as->at(as->at(i).number_y).r_obr_x = as->at(as->at(i).number_xy).r - as->at(as->at(i).number_y).r;
					}
					as->at(i).obr_x = false;
					as->at(i).obr_y = false;
					as->at(i).obr_xy = false;
				}
				else if (as->at(i).obr_xy)//реальной стал образ по ху
				{
					as->at(as->at(i).number_xy).im = false;//меняем образ х который теперь стал реальной с рассматриываемой частицей, которая стала теперь мнимой
					as->at(as->at(i).number_xy).obr_xy = true;
					as->at(as->at(i).number_xy).number_xy = i;
					if ((as->at(i).r - as->at(as->at(i).number_xy).r).x > 0)
					{
						as->at(as->at(i).number_xy).r_obr_x = 1;
					}
					else {
						as->at(as->at(i).number_xy).r_obr_x = -1;
					}
					//as->at(as->at(i).number_xy).r_obr_x = as->at(i).r - as->at(as->at(i).number_xy).r;
					if (as->at(i).obr_x)//у рассматриваемой частицы 3 образа, оставшиеся меняем местами
					{
						//if ((abs(as->at(i).number_xy)>2000)||(abs(as->at(i).number_y)>2000)){system("pause");}
						as->at(as->at(i).number_xy).number_y = as->at(i).number_x;
						as->at(as->at(i).number_xy).number_x = as->at(i).number_y;
						as->at(as->at(i).number_xy).obr_x = true;
						as->at(as->at(i).number_xy).obr_y = true;
						if ((as->at(as->at(i).number_x).r - as->at(as->at(i).number_xy).r).y > 0)
						{
							as->at(as->at(i).number_xy).r_obr_y = 1;
						}
						else {
							as->at(as->at(i).number_xy).r_obr_y = -1;
						}
						//as->at(as->at(i).number_xy).r_obr_y = as->at(as->at(i).number_x).r - as->at(as->at(i).number_xy).r;
						if ((as->at(as->at(i).number_y).r - as->at(as->at(i).number_xy).r).x > 0)
						{
							as->at(as->at(i).number_xy).r_obr_x = 1;
						}
						else {
							as->at(as->at(i).number_xy).r_obr_x = -1;
						}
						//as->at(as->at(i).number_xy).r_obr_x = as->at(as->at(i).number_y).r - as->at(as->at(i).number_xy).r;
					}
					as->at(i).obr_x = false;
					as->at(i).obr_y = false;
					as->at(i).obr_xy = false;
				}


			}

			else //частица не пересекала границы
			{
				if (as->at(i).obr_x)
				{
					as->at(as->at(i).number_x).r = as->at(i).r + as->at(i).r_obr_x*ww;
				}
				if (as->at(i).obr_y)
				{
					as->at(as->at(i).number_y).r = as->at(i).r + as->at(i).r_obr_y*hh;
				}
				if (as->at(i).obr_xy)
				{
					Vect2D pr = wh;
					int px = 0;
					int py = 0;
					if ((as->at(i).r_obr_xy) == 22) {
						px = 1;
						py = 1;
					}
					if ((as->at(i).r_obr_xy) == 12) {
						px = -1;
						py = 1;
					}
					if ((as->at(i).r_obr_xy) == 21) {
						px = 1;
						py = -1;
					}
					if ((as->at(i).r_obr_xy) == 11) {
						px = -1;
						py = -1;
					}
					pr.x *= px;
					pr.y *= py;
					as->at(as->at(i).number_xy).r = as->at(i).r + pr;
				}
			}

		}
	}
	for (int i = 0; i<as->size();i++)//интегрируем перемещения каждой частице
	{
		as->at(i).r += as->at(i).v*dt;
		as->at(i).f = VECT2D_0;
	}
}
void StepWithSimpleModle(vector<Atom2D>* as, vector<connection>* l, vector<int>* PressBond,double w,double h, int numNodes, int numNodesDisplacement, double lenBetweenNodes, vector<double>* Press, vector<double>* width, vector<double>* q,vector<double>*u, vector<double>*pr_u, vector<int>*num_u)
{
	pr_u->clear();
	num_u->clear();
	for (int iii = 0;iii < 2 * numNodesDisplacement;iii++)
	{
		pr_u->push_back(double(0));
		num_u->push_back(int(0));
	}
	ofstream particl;
	string nammm = "Particles.txt";
	particl.open(nammm);
	for (int i = 0;i < PressBond->size();i++)//Определение раскрытия
	{
		int num = PressBond->at(i);
		int left_at = l->at(num).i;
		int right_at = l->at(num).j;
		if (as->at(left_at).r.x>as->at(right_at).r.x)
		{
			int pr = left_at;
			left_at = right_at;
			right_at = pr;
		}
		double x = 0.5*(as->at(left_at).r.x + as->at(right_at).r.x);
		double y = 0.5*(as->at(left_at).r.y + as->at(right_at).r.y);
		particl << x - w / 2 << "  " << y - h / 2 << "\n";
		int ff = 0;
		if (y - h*0.5 < 0)
		{
			ff = 1;
		}
		else
		{
			ff = 0;
		}
	
		
		int nomerNode = int(abs(y - h*0.5) / lenBetweenNodes);
		if (nomerNode >= numNodesDisplacement)
		{
			cout << "Vihod za predeli" << "\n";
			cout << y << "  " << as->at(left_at).r.y << " " << as->at(right_at).r.y << "  " << left_at << "  " << right_at << "\n";
			system("pause");
			//nomerNode--;
		}
		else{
		Vect2D drr = (as->at(left_at).r - as->at(right_at).r);
		double cos_alff = abs(drr.y / drr.x);
		double pp1 = LenConnection(as->at(left_at), as->at(right_at));
		double pp2 = sqrt(l->at(num).eq_st_sq);
		double ppppa = pp1 - pp2;
		double ur = (LenConnection(as->at(left_at), as->at(right_at)) - sqrt(l->at(num).eq_st_sq))*cos_alff;
		if ((stepp == 5) && (nomerNode == 0))
		{ 
			cout << ur << "\n";
			if (ur > 10)
			{
				cout << "left " << left_at << "right " << right_at << "\n";
				int iiiiiiii12 = 0;
			}
		}
		pr_u->at(ff*numNodesDisplacement + nomerNode) += ur;
		num_u->at(ff*numNodesDisplacement + nomerNode) = num_u->at(ff*numNodesDisplacement + nomerNode) +1;
		}
	}
	particl.close();
	if ((q->at(0)<1e-7)){
		q->at(0) += 1e-9;
		q->at(numNodes) = q->at(0);
		cout << q->at(0) << "\n";
	}
	for (int it = 0;it < numNodesDisplacement;it++)
	{
		double pr = 0;
		if (num_u->at(it) != 0)
		{
			pr = pr_u->at(it) / num_u->at(it);
			if (pr < 1e-11)
			{
				u->at(it) = 0;
			}
			else {
				u->at(it) = pr;
			}
		}
		else {
			pr = pr_u->at(it);
			if (pr < 1e-11)
			{
				u->at(it) = 0;
			}
			else {
				u->at(it) = pr;
			}
			
		}
		if (num_u->at(it + numNodesDisplacement) != 0)
		{
			pr = pr_u->at(it + numNodesDisplacement) / num_u->at(it + numNodesDisplacement);
			if (pr < 1e-11)
			{
				u->at(it + numNodesDisplacement) = 0;
			}
			else {
				u->at(it + numNodesDisplacement) = pr;
			} 
		}
		else {
			pr = pr_u->at(it + numNodesDisplacement);
			if (pr < 1e-11)
			{
				u->at(it + numNodesDisplacement) = 0;
			}
			else {
				u->at(it + numNodesDisplacement) = pr;
			}
		}
	}
	pr_u->clear();
	num_u->clear();
	for (int it = 1;it < numNodes - 1;it++)//Обновления потоков
	{
		if (u->at(it) > 1e-10)
		{
			q->at(it) = koef_q*pow((width->at(it) - width->at(it - 1)), 3)*(width->at(it) - width->at(it - 1) - (u->at(it) - u->at(it - 1)));
			q->at(it + numNodes) = koef_q*pow((width->at(it + numNodesDisplacement) - width->at(it + numNodesDisplacement - 1)), 3)*(width->at(it + numNodesDisplacement) - width->at(it + numNodesDisplacement - 1) - (u->at(it + numNodesDisplacement) - u->at(it + numNodesDisplacement - 1)));

		}
	}
	for (int it = 0;it < numNodesDisplacement;it++)//Обновления jобъемов жидкости в узлах
	{
		double pr1 = -(dt_liq / lenBetweenNodes)*(q->at(it+1)-q->at(it));
		width->at(it) += pr1*dt;
		double pr2 = -(dt_liq / lenBetweenNodes)*(q->at(it + numNodes + 1) - q->at(it + numNodes));
		width->at(it + numNodesDisplacement) += pr2*dt;
	}
	for (int it = 0;it < numNodesDisplacement;it++)//Обновления давлений в узлах
	{
		double tt = width->at(it) - u->at(it);
		if (width->at(it)>0)
		{
			Press->at(it) = koef_pp*(tt);
		}
		else
		{
			Press->at(it) = 0;
		}
		tt = width->at(it + numNodesDisplacement) - u->at(it + numNodesDisplacement);
		if (width->at(it + numNodesDisplacement)>0)
		{
			Press->at(it + numNodesDisplacement) = koef_pp*(tt);
		}
		else
		{
			Press->at(it + numNodesDisplacement) = 0;
		}
	}
	ofstream flow;
	ofstream widt;
	ofstream Press_t;
	ofstream u_t;
	string nameVel = "Flow.txt";
	string nameWit = "Widht.txt";
	string namePress = "Press.txt";
	string nameU = "U_.txt";
	Press_t.open(namePress, std::ofstream::out | std::ofstream::app);
	flow.open(nameVel, std::ofstream::out | std::ofstream::app);
	widt.open(nameWit, std::ofstream::out | std::ofstream::app);
	u_t.open(nameU, std::ofstream::out | std::ofstream::app);
	for (int st = 0;st < 2*numNodes;st++)
	{
		flow << q->at(st)<<"  ";
	}
	flow << "\n";
	
	for (int st = 0;st < 2*numNodesDisplacement;st++)
	{
		u_t << u->at(st) << "  ";
		Press_t << Press->at(st) << "  ";
		widt << width->at(st) << "  ";
	}
	u_t << "\n";
	Press_t << "\n";
	widt << "\n";
	Press_t.close();
	u_t.close();
	widt.close();
	flow.close();
	for (int ii = 0;ii < as->size();ii++)
	{
		Vect2D pppp = Vect2D(0, 0);
		as->at(ii).f = pppp;
	}
	for (int i = 0;i < PressBond->size();i++)//Определение давления
	{
		int num = PressBond->at(i);
		int left_at = l->at(num).i;
		int right_at = l->at(num).j;
		if (as->at(left_at).r.x>as->at(right_at).r.x)
		{
			int pr = left_at;
			left_at = right_at;
			right_at = pr;
		}
		double x = 0.5*(as->at(left_at).r.x + as->at(right_at).r.x);
		double y = 0.5*(as->at(left_at).r.y + as->at(right_at).r.y);
		int ff = 0;
		if (y - h*0.5 < 0)
		{
			ff = 1;
		}
		int nomerNode = int(abs(y - h*0.5) / lenBetweenNodes);
		float len_e = (as->at(right_at).r - as->at(left_at).r).Abs();
		Vect2D e = (as->at(right_at).r - as->at(left_at).r)/len_e;
		/*ЗАДАЕМ ДАВЛЕНИЯ ПОСТОЯННЫМ*/
		as->at(left_at).f -= e*press;//e*Press->at(ff*numNodesDisplacement + nomerNode);
		as->at(right_at).f += e*press;//e*Press->at(ff*numNodesDisplacement + nomerNode);
	}
	
	float absForce = 0;
	Vect2D dr;
	double length = 0;
	for (int kk = 0; kk < as->size();kk++)//проходим по частицам
	{
		if (as->at(kk).numSpringBond>0)
		{

			for (int j1 = 0;j1 < as->at(kk).arrNumBond.size();j1++)//проходим по номерам связей имеющихся у частицы
			{
				int jj = as->at(kk).arrNumBond.at(j1);//номер связи
				if (!(l->at(jj).cheek))
				{
					int numSecPart = l->at(jj).j;
					if (numSecPart == kk)
					{
						numSecPart = l->at(jj).i;
					}
					dr = as->at(kk).r - as->at(numSecPart).r;

					length = sqrt(dr.Sqr()) - sqrt(l->at(jj).eq_st_sq); //растяжение-сжатие пружины
					if (length != 0) {
						//		system("pause");
					}
					if (l->at(jj).pressure)
					{
						/*ff = -press;
						if (ff != 0) {
							//system("pause");
						}
						as->at(kk).f -= (ff / sqrt(dr.Sqr()))*dr;
						as->at(numSecPart).f += (ff / sqrt(dr.Sqr()))*dr;*/
					}
					else
						if (dr.Sqr() > pow(l->at(jj).eq_st_sq*e0, 1))//проверка на состоянии пружины
						{//случай разрыва связи
							int i1 = l->at(jj).i;
							int j1 = l->at(jj).j;
							int fi = 1;
							int fj = 1;
							double x1 = as->at(i1).r.x;
							double x2 = as->at(j1).r.x;
							double y1 = as->at(i1).r.y;
							double y2 = as->at(j1).r.y;
							double x_tr1 = w / 2;//- lenInitialFracture*0.5*e*cos(pi / 2);
							double x_tr2 = w / 2; //+ lenInitialFracture*0.5*e*cos(pi / 2);
							double y_tr1 = 0;//h / 2 + lenInitialFracture*0.5*e*sin(pi / 2);
							double y_tr2 = h;// -lenInitialFracture*0.5*e*sin(pi / 2);
							if (intersection(x1, y1, x2, y2, x_tr1, y_tr1, x_tr2, y_tr2))
							{
								bool f = true;
								if ((f))
								{
									l->at(jj).pressure = true;
									int ii1 = l->at(jj).i;
									int jj1 = l->at(jj).j;
									as->at(ii1).numSpringBond--;
									as->at(ii1).pressure = true;
									as->at(jj1).numSpringBond--;
									as->at(jj1).pressure = true;
									PressBond->push_back(jj);
									numPressBond++;
									num_errase++;
									//cout << num_errase << "\n";
								}
								else
								{
									l->at(jj).free = true;
									int ii1 = l->at(jj).i;
									int jj1 = l->at(jj).j;
									as->at(ii1).numSpringBond--;
									as->at(jj1).numSpringBond--;
									//num_errase++;
									//cout << num_errase << "\n";
								}
								//l->erase(l->begin()+ii);
								//system("pause");
							}

						}
						else if (!(l->at(jj).free))
						{
							absForce = float(length*k_bond);//it->koef;//модуль силы
							/*if (ff != 0) {
								//system("pause");
							}*/
							as->at(kk).f -= (absForce / sqrt(dr.Sqr()))*dr;
							as->at(numSecPart).f += (absForce / sqrt(dr.Sqr()))*dr;
						}
					l->at(jj).cheek = true;
				}
				else {
					l->at(jj).cheek = false;
				}
			}

		}
		
	}
	/*ofstream rrr;
	string nn = "force.txt";
	rrr.open(nn, std::ofstream::out | std::ofstream::app);
	if ((stepp > 100)&&(stepp%10 ==0))
	{
		rrr << "Step " << stepp << "\n";
		for (int iii = 0;iii < as->size();iii++)
		{
			rrr << as->at(iii).f.x << "  " << as->at(iii).f.y << "\n";
		}
	}
	rrr.close();*/
	RecalculateVelosityPositionFixed(as, w, h);
	stepp++;
	
}
void StepWithLiquidPipe(vector<Atom2D>* as, vector<connection>* l, vector<int>* PressBond, Cell**gridBonds, double w, double h, double periodic_border, double bond_len, int numCellsX, int numCellsY)// int numNodes, int numNodesDisplacement, double lenBetweenNodes, vector<double>* Press, vector<double>* widht, vector<double>* q, vector<double>* u, vector<double>* pr_u, vector<int>* num_u)
{
	ofstream dpp;//Файлы для сохранения значений давлений и раскрытий
	char file_press[20];
	sprintf(file_press, "dp%08u.xyz", stepp + 1);
	ofstream dqq;
	char file_qqq[20];
	sprintf(file_qqq, "dq%08u.xyz", stepp + 1);
	
	if (stepp % save_step == 0)
	{
		dqq.open(file_qqq);
		dpp.open(file_press);
		//duu.open(file_uuu);
		dqq << numPressBond << "\n" << "comment\n";
		dpp << numPressBond << "\n" << "comment\n";
	}
	/*пересчет потоков*/
	for (int it_l = 0;it_l < l->size();it_l++)
	{
			if (l->at(it_l).origin)//Добавка от источника
			{
				l->at(it_l).Q += q_0*dt_liq;
			}
			int numNodePipe = 0;
			int ai2 = l->at(it_l).i;
			int aj2= l->at(it_l).j;
			double u_it = LenConnection(as->at(ai2),as->at(aj2))-sqrt(l->at(it_l).eq_st_sq);

			if ((l->at(it_l).pressure)&&(u_it>0))
			{
				for (int i = 0;i < l->at(it_l).numLiquidBond.size();i++)
				{
					numNodePipe = l->at(it_l).numLiquidBond.at(i);
					double dop_koef = 0;
					int ai1 = l->at(numNodePipe).i;
					int aj1= l->at(numNodePipe).j;
					double u_numNodePipe = LenConnection(as->at(ai1),as->at(aj1))-sqrt(l->at(numNodePipe).eq_st_sq);
					
					if ((l->at(numNodePipe).pressure)&&(u_numNodePipe>0))
					{
						//cout<<u_it<<" "<<u_numNodePipe<<"\n";
						double distt = 0.5*(u_it+u_numNodePipe);
						dop_koef = distt;
					}
					double dq = -dop_koef*koef_dq * (l->at(it_l).pressInBond - l->at(numNodePipe).pressInBond);
					l->at(it_l).Q += dq*dt_liq;
					if (l->at(it_l).Q<0)
					{
						l->at(it_l).Q = 0;//std::cout<<dq<<"\n";
					}
					
				}
				if (l->at(it_l).Q < 0)
				{
					//l->at(it_l).Q = 0;
				}
			}
	}
	int numberLiquid = 1;
	/*пересчет объемов*/
	for (int it_l = 0;it_l < l->size();it_l++)//Сохранение значений объемов для последующей итерации
	{

		int numNodePipe = 0;
		
		if ((l->at(it_l).pressure))// && (!(l->at(it_l).pulled)))
		{
			//l->at(it_l).p_Q = l->at(it_l).Q;
			
		}
		if ((l->at(it_l).pressure) && (stepp % save_step == 0))
		{

			int ai = l->at(it_l).i;
			int aj = l->at(it_l).j;
			dqq << numberLiquid << " " << l->at(it_l).Q << " " << (as->at(ai).r.x + as->at(aj).r.x)*0.5 << " " << (as->at(ai).r.y + as->at(aj).r.y)*0.5 << "\n";
			numberLiquid++;

		}
	}
	vector<int>numChekingCells;
	numChekingCells.clear();
	//#pragma omp parallel
	//{
	//#pragma omp for
	numberLiquid = 1;
	for (int it = 0;it < as->size();it++)
	{
		if(true)//((as->at(it).Spring))//Не явлеется ли эта частица только с давленчискими связями
		{
			for (int j1 = 0;j1 < as->at(it).arrNumBond.size();j1++)//проходим по номерам связей имеющихся у частицы
			{
				int jj = as->at(it).arrNumBond.at(j1);//номер связи
				/*ДОБАВКА К СИЛАМ*/
				if (true)//(!(l->at(jj).pulled))//Не связывает ли связь частицу с только с далвением
				{
					int ai = l->at(jj).i;
					int aj = l->at(jj).j;
					//Vect2D pr = as->at(ai).r - as->at(aj).r;
					double distt = LenConnection(as->at(ai), as->at(aj)) - sqrt(l->at(jj).eq_st_sq);//pr.Abs() - sqrt(l->at(jj).eq_st_sq);

					double delta_P = (l->at(jj).Q-distt)*koef_pp;//(l->at(jj).Q - distt)*koef_pp;//жидкость - длина связи на коэф.
					l->at(jj).pressInBond = delta_P;
					if(abs(delta_P) > 1e6)
					{
						cout << "press\n";
						//system("pause");
					}
					if ((l->at(jj).pressure) && (stepp%save_step == 0))
					{

						double ur = LenConnection(as->at(ai), as->at(aj)) - sqrt(l->at(jj).eq_st_sq);

						dpp << numberLiquid << " " << l->at(jj).pressInBond << " " << (as->at(ai).r.x+ as->at(aj).r.x)*0.5 << " " << (as->at(ai).r.y + as->at(aj).r.y)*0.5 <<" "<<ur<< "\n";// << " " << as->at(aj).r.x << " " << as->at(aj).r.y << " " << ur << "\n";
						numberLiquid++;
					}
					if (!(l->at(jj).cheek))
					{
						int numSecPart = l->at(jj).j;
						if (numSecPart == it)
						{
							numSecPart = l->at(jj).i;
						}
						Vect2D dr = as->at(it).r - as->at(numSecPart).r;

						double length = sqrt(dr.Sqr()) - sqrt(l->at(jj).eq_st_sq); //растяжение-сжатие пружины

						if (l->at(jj).pressure)
						{
							as->at(it).f += float(delta_P)*dr / dr.Abs();
							as->at(numSecPart).f -= float(delta_P)*dr / dr.Abs();
						}
						else
							if (dr.Sqr() > pow(l->at(jj).eq_st_sq*e0, 1))//проверка на состоянии пружины
							{//случай разрыва связи
								int i1 = l->at(jj).i;
								int j1 = l->at(jj).j;
								int fi = 1;
								int fj = 1;
								bool f = true;
								if ((f))
								{
									l->at(jj).pressure = true;
									int ii1 = l->at(jj).i;
									int jj1 = l->at(jj).j;
									as->at(ii1).numSpringBond--;
									as->at(ii1).pressure = true;
									as->at(jj1).numSpringBond--;
									as->at(jj1).pressure = true;
									/*Проверяю не стали ли частицы только с давлением*/
									if (as->at(ii1).numSpringBond == 0)
									{
										as->at(ii1).Spring = false;
										l->at(jj).pulled = true;
									}
									if (as->at(jj1).numSpringBond == 0)
									{
										as->at(jj1).Spring = false;
										l->at(jj).pulled = true;
									}
									if (!(as->at(ii1).Spring) && (!(as->at(jj1).Spring)))
									{
										
									}
									PressBond->push_back(jj);
									double x = 0.5*(as->at(ii1).r.x + as->at(jj1).r.x);
									double y = 0.5*(as->at(ii1).r.y + as->at(jj1).r.y);
									int numgridX = int((x + periodic_border) / bond_len);//номер ячейки по х куда попала частица
									int numgridY = int((y + periodic_border) / bond_len);
									numChekingCells.push_back(numgridX);
									numChekingCells.push_back(numgridY);
									numPressBond++;
									num_errase++;
									//cout << num_errase << "\n";
								}
								else
								{
									l->at(jj).free = true;
									int ii1 = l->at(jj).i;
									int jj1 = l->at(jj).j;
									as->at(ii1).numSpringBond--;
									as->at(jj1).numSpringBond--;
									//num_errase++;
									//cout << num_errase << "\n";
								}
								//l->erase(l->begin()+ii);
								//system("pause");


							}
							else if (!(l->at(jj).free))
							{
								float absForce = float(length*k_bond);//it->koef;//модуль силы
														 /*if (ff != 0) {
														 //system("pause");
														 }*/
								as->at(it).f -= (absForce / sqrt(dr.Sqr()))*dr;

								as->at(numSecPart).f += (absForce / sqrt(dr.Sqr()))*dr;
							}
						l->at(jj).cheek = true;
					}
					else {
					l->at(jj).cheek = false;
				}
				}
			}
		}

	}
	//}
	if (nodes)
	{
		//cout <<"num " <<numChekingCells.size()<<"\n";
		if (numChekingCells.size() > 0)
		{
			for (int i1 = 0;i1 < numChekingCells.size() - 1;i1 = i1 + 2)
			{
				//system("pause");
				int xgrid = numChekingCells.at(i1);
				int ygrid = numChekingCells.at(i1 + 1);
				ConnectLiquidNodes(as, l, gridBonds, xgrid, ygrid, numCellsX, numCellsY, bond_len);
			}
		}
	}
	dpp.close();
	RecalculateVelosityPositionFixed(as, w, h);
	stepp++;
}
void CreateLiquidSourse(vector<Atom2D>*as,vector<connection> *l,double x,double y)
{
	bool presec = true;
	for (int it_liq = 0;it_liq < l->size();it_liq++)
	{
		if ((l->at(it_liq).pressure)&&(presec))
		{
			int at1 = l->at(it_liq).i;
			int at2 = l->at(it_liq).j;
			double x11 = (as->at(at1).r.x + as->at(at2).r.x) / 2;
			double y11 = (as->at(at1).r.y + as->at(at2).r.y) / 2;
			if ((abs(x11 - x) < 5 * 1e-2) && (abs(y11 - y) < 5 * 1e-2))
			{
				l->at(it_liq).origin= true;
				presec = false;

			}
		}
	}
}
void CreateFracture(vector<Atom2D>*as,vector<connection> *l,vector<int> *PressBond, Cell **gridBonds,double x, double y,
	double alfa,double lenght,double periodic_border,double bond_len,int numCellsX,int numCellsY)
{
	vector<int>numChekingCells;
	bool** CheckGrid= new bool*[numCellsX];
	for (int ii1 = 0;ii1 < numCellsX;ii1++)
	{
		CheckGrid[ii1] = new bool[numCellsY];
		for (int jj1 = 0;jj1 < numCellsY;jj1++)
		{
			CheckGrid[ii1][jj1] = false;
		}
	}
	for (int p = 0;p<l->size();p++)//удаляю связи в области
	{
		int i1 = l->at(p).i;
		int j1 = l->at(p).j;
		int fi=1;
		int fj=1;
		double x1 = as->at(i1).r.x;
		double x2 = as->at(j1).r.x;
		double y1 = as->at(i1).r.y;
		double y2 = as->at(j1).r.y;
		double x_tr1 = x- 0.5*lenght*cos(alfa);
		double x_tr2 = x+ lenght*0.5*cos(alfa);
		double y_tr1 = y+ 0.5*lenght*sin(alfa);
		double y_tr2 = y- 0.5*lenght*sin(alfa);
		if (intersection(x1,y1,x2,y2,x_tr1,y_tr1,x_tr2,y_tr2))
		{
			l->at(p).pressure = true;
			int ii = l->at(p).i;
			int jj = l->at(p).j;
			as->at(ii).numSpringBond--;
			as->at(jj).numSpringBond--;
			as->at(ii).pressure = true;
			as->at(jj).pressure = true;
			
			numPressBond++;
			PressBond->push_back(p);
			double x;
			double y;
			x = 0.5*(x1 + x2);
			y = 0.5*(y1 + y2);
			int numgridX = int((x + periodic_border) / bond_len);//номер ячейки по х куда попала частица
			int numgridY = int((y + periodic_border) / bond_len);
			if (!CheckGrid[numgridX][numgridY])
			{
				numChekingCells.push_back(numgridX);
				numChekingCells.push_back(numgridY);
				CheckGrid[numgridX][numgridY] = true;
			}
			

			//gridPressBond[numgridX][numgridY] = true;
			//kol_bond++;
			
		}
		
	}
	for (int ii1 = 0;ii1 < numCellsX;ii1++)
	{
		delete[] CheckGrid[ii1];
	}
	delete[] CheckGrid;
	if (nodes)
	{
		for (int i1 = 0;i1 < numChekingCells.size() - 1;i1 = i1 + 2)
		{
			int xgrid = numChekingCells.at(i1);
			int ygrid = numChekingCells.at(i1 + 1);
			cout << xgrid << " " << ygrid << "\n";
			ConnectLiquidNodes(as, l, gridBonds, xgrid, ygrid, numCellsX, numCellsY, bond_len);
		}
	}
}
void ConnectLiquidNodes(vector<Atom2D>*P, vector<connection>* B, Cell** gridBonds, int numgridX, int numgridY, int numCellsX, int numCellsY,double bond_len)
{
	int N_bonds = gridBonds[numgridX][numgridY].arrNumPart.size();
	for (int it = 0;it < N_bonds;it++)//Проходим по связям в ячейке.
	{
		int it_bond = gridBonds[numgridX][numgridY].arrNumPart.at(it);
		int ai = B->at(it_bond).i;
		int  aj = B->at(it_bond).j;
		double x = (P->at(ai).r.x + P->at(aj).r.x) / 2;
		double y = (P->at(ai).r.y + P->at(aj).r.y) / 2;
		int* arrNumCellX;//массив номеров ячеек проверяем по х
		int* arrNumCellY;
		int X_num = 3;//количество ячеек проверяем по х
		int Y_num = 3;//количество ячеек проверяем по н
		if (numgridX == numCellsX - 1) {
			arrNumCellX = new int[2];
			arrNumCellX[0] = numgridX - 1;
			arrNumCellX[1] = numgridX;
			X_num = 2;
		}
		else if (numgridX == 0) {
			arrNumCellX = new int[2];
			arrNumCellX[0] = 0;
			arrNumCellX[1] = 1;
			X_num = 2;
		}
		else {
			arrNumCellX = new int[3];
			arrNumCellX[0] = numgridX - 1;
			arrNumCellX[1] = numgridX;
			arrNumCellX[2] = numgridX + 1;
		}
		if (numgridY == numCellsY - 1) {
			arrNumCellY = new int[2];
			arrNumCellY[0] = numgridY - 1;
			arrNumCellY[1] = numgridY;
			Y_num = 2;
		}
		else if (numgridY == 0) {
			arrNumCellY = new int[2];
			arrNumCellY[0] = 0;
			arrNumCellY[1] = 1;
			Y_num = 2;
		}
		else {
			arrNumCellY = new int[3];
			arrNumCellY[0] = numgridY - 1;
			arrNumCellY[1] = numgridY;
			arrNumCellY[2] = numgridY + 1;
		}
		for (int iiCon1 = 0;iiCon1 < X_num;iiCon1++)//По ячейкам Х
		{
			for (int iiCon2 = 0;iiCon2 < Y_num;iiCon2++)//По ячейкам У
			{
				int ii = arrNumCellX[iiCon1];
				int jj = arrNumCellY[iiCon2];
				int NN = gridBonds[ii][jj].numberPart;
				for (int k = 0;k < NN;k++)
				{
					int nomer = gridBonds[ii][jj].arrNumPart.at(k);
					int ai1 = ai;//итый атом 1 связи
					int ai2 = B->at(nomer).i;//итый атом 2 связи
					int aj1 = aj;
					int aj2 = B->at(nomer).j;
					auto res = std::find(B->at(it_bond).numLiquidBond.begin(),B->at(it_bond).numLiquidBond.end(),nomer);
					if((it_bond==nomer)||((res!=B->at(it_bond).numLiquidBond.end())))
					{
						break;
					}
					double x1 = double(0.5*(P->at(ai1).r.x + P->at(aj1).r.x));
					double x2 = double(0.5*(P->at(ai2).r.x + P->at(aj2).r.x));
					double y1 = double(0.5*(P->at(ai1).r.y + P->at(aj1).r.y));
					double y2 = double(0.5*(P->at(ai2).r.y + P->at(aj2).r.y));
					double dist = pow((x1 - x2), 2) + pow((y1 - y2), 2);
					if (dist < pow(part_bond_liquid* bond_len, 2))
					{
						B->at(nomer).numLiquidBond.push_back(it_bond);
						B->at(it_bond).numLiquidBond.push_back(nomer);

					}
				}
			}
		}
		delete[]arrNumCellX;
		delete[]arrNumCellY;

	}
	/*Прошли по всем связям в ячейки, значит они уже все соединены. Можно удалить из памяти*/
	//gridBonds[numgridX][numgridY].arrNumPart.clear();
	gridBonds[numgridX][numgridY].Delete();
}

void SaveConfiguration(vector<Atom2D>* P, vector<connection>* B)//Функция сохранения грунта
{
	ofstream Part_list;//Файлы для сохранения частиц
	char file_part_list[20];
	string name = "Particles_list" + std::to_string(stepp)+".txt";
	//sprintf(file_part_list,name);
	Part_list.open(name);
	Part_list << P->size() << "\n";
	for (int i = 0;i < P->size();i++)
	{
		if ((P->at(i).pressure)) {
			Part_list << P->at(i).r.x << " " << P->at(i).r.y << " " << 0 << "\n";
		}
	}
	Part_list.close();
/*	ofstream Bond_list;
	char file_bond_list[20];
	sprintf(file_bond_list, "Bonds_list.txt");
	Bond_list.open(file_bond_list);
	Bond_list << B->size() << "\n";
	for (int j = 0;j < B->size();j++)
	{
		Bond_list << B->at(j).i << " " << B->at(j).j <<" "<< B->at(j).eq_st_sq<<"\n";
	}
	Bond_list.close();*/
}

void LoadConfiguration(vector<Atom2D>* P, vector<connection>* B, Cell** gridBonds,int numCellsX, int numCellsY, double bond_len,double periodic_border)
{
	
	FILE *file;
	file = fopen("Particles_list.txt", "r");
	/*if ((file = fopen("Particles_list.txt", "r")) == NULL) {
		printf("Ошибка открытия файла.\n");
		system("pause");
	}*/
	int NN;
	fscanf(file, "%i", &(NN));
	cout << NN << endl;
	float x = 0;
	float y = 0;
	int fl = 0;
	bool im = false;
	for (int i = 0;i < NN;i++)
	{
		fscanf(file, "%g", &(x));
		fscanf(file, "%g", &(y));
		fscanf(file, "%u", &(fl));
		if (fl == 1) { im = true; }
		else { im = false; }
		Atom2D pr(x, y, fl);
		P->push_back(pr);
	}
	fclose(file);
	file = fopen("Bonds_list.txt", "r");
	if ((file = fopen("Bonds_list.txt", "r")) == NULL) {
		printf("Ошибка открытия файла.\n");
		system("pause");
	}
	//string pr;
	//char str[80];
	//int NN = 0;
	fscanf(file, "%i", &(NN));
	int i;
	int j;
	float eq;
	cout << NN << endl;
	for (int ii = 0;ii < NN;ii++)
	{
		fscanf(file, "%u", &(i));
		fscanf(file, "%u", &(j));
		fscanf(file, "%g", &(eq));
		connection p(i, j, eq);
		p.pulled = true;
		B->push_back(p);
		P->at(i).arrNumBond.push_back(B->size()-1);
		P->at(j).arrNumBond.push_back(B->size()-1);
		P->at(i).numSpringBond++;
		P->at(j).numSpringBond++;
	}
	for (int ii = 0;ii<B->size();ii++)//Заполняю сетку
	{
		int a2i = B->at(ii).i;
		int  a2j = B->at(ii).j;
		double xx = (P->at(a2i).r.x + P->at(a2j).r.x) / 2;
		double yy = (P->at(a2i).r.y + P->at(a2j).r.y) / 2;
		int numgridX = int((xx + periodic_border) / bond_len);//номер ячейки по х куда попала частица
		int numgridY = int((yy + periodic_border) / bond_len);
		gridBonds[numgridX][numgridY].Add(ii);

	}
	

}


void SaveLiquidConnection(vector<Atom2D>*P,vector<connection>*B)
{
	ofstream LiquidConn;//Файлы для сохранения частиц
	char file_part_list[20];
	sprintf(file_part_list, "Liquid_conn.txt");
	LiquidConn.open(file_part_list);
	//LiquidConn << P->size() << "\n";
	for (int i = 0;i < B->size();i++)
	{
		if (B->at(i).pressure)
		{
			LiquidConn<<i<<"  ";
			for (int it = 0;it<B->at(i).numLiquidBond.size();it++)
			{
				LiquidConn << B->at(i).numLiquidBond.at(it) << " " ;
		
			}
			LiquidConn<<"\n";
		}
	}
	LiquidConn.close();
}

void FindPressInFractures(double V_liq, vector<int>* FrackPress, double * masPress, int numFracks,double bond)
{
	int* numBond = new int[numFracks];
	for (int ii = 0;ii < numFracks;ii++) {
		numBond[ii] = 0;
	}
	for (int ii = 0;ii < FrackPress->size();ii+=2) {//считаю количество связей в каждой трещине
		int cur_fr = FrackPress->at(ii + 1);
		numBond[cur_fr]++;
	}
	for (int ii = 0;ii < numFracks;ii++) {
		masPress[ii] = V_liq / pow(numBond[ii]*bond,2);
	}
	delete[] numBond;
}
