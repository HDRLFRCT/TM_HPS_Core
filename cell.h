#include<vector>

class Cell
{
public:
	Cell(int N);
	Cell();
	int numberPart;
	int maxPart;
	//int* arrNumPart;
	std::vector<int>arrNumPart;
	void Add(int numPart);
	void Delete();
	Cell& operator=(const Cell &at2d );
	~Cell();
};
