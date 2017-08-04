#include <stdio.h>  
#include <stdlib.h>  

#include <windows.h>   
#include <stdio.h>   
#include <mmsystem.h>   
#include <omp.h>
#include <vector>
#include <algorithm>

#include<numeric>
#include"./mainpart/part2.h"
#include"./draw/draw.h"

#include"./basicfunction/basicfunction.h"
#include"./fourier/fourier.h"

class cfreq :public complexfunction
{
public:
	cpoly x; //系统输入部分
	cpoly y;//系统输出部分
	cfreq() {};
	cfreq(cpoly x, cpoly y) { this->x = x; this->y = y; };
	virtual complex operator ()(double xx)
	{
		return value(y, complex(0, xx)) / value(x, complex(0, xx));
	}
};

template<class T>
class mylog
{
public:
	T input;
	mylog() {};
	mylog(T in) { input = in; };
	double operator() (double x)
	{
		return 20*log10(input.cabs(x/6000));
	}
};

int main()
{
	cpoly but = Butterworth(5);
	cpoly xx;
	xx.push_back(1);
	cfreq butfilt = cfreq( but,xx);
	mylog<cfreq> l = mylog<cfreq>(butfilt);
	showc(l, 0, 12000, -50, 0);

	getchar();
}