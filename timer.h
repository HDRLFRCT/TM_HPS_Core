//-----------------------------------------------------------------------------------
//
//	Timer class
//
//	Vladimir Bratov
//
//	20.06.2001
//
//-----------------------------------------------------------------------------------

#ifndef ___Timer_h___
#define ___Timer_h___

#include <time.h>

//-----------------------------------------------------------------------------------

class Timer
{
public:

	Timer()									{ Start(); }

	void Start()							{ start = clock(); }
	
	double GetTime()
	{
		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC; // in Vc++<6.0 clocks_per_sec constant is called CLK_TCK
  		return duration;
	}

	void Wait( clock_t wait )
	{
		clock_t goal;
		goal = wait + clock();
		while( goal > clock() );
	}

private:
	
	clock_t start, finish;
	double  duration;
};

//-----------------------------------------------------------------------------------

#endif //___Timer_h___
