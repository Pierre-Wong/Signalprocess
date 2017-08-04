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



class freq :public complexfunction
{
public:
	poly x; //系统输入部分
	poly y;//系统输出部分
	freq() {};
	freq(poly x, poly y) { this->x = x; this->y = y; };
	virtual complex operator ()(double xx)
	{
		return value(y, complex(0, xx)) / value(x, complex(0, xx));
	}
};


int main()
{
	poly x, y;
	x.push_back(2);
	x.push_back(3);
	x.push_back(1);
	y.push_back(3);
	y.push_back(2);
	freq zz = freq(x, y);
	showc(zz.cabs, 0, 10, 0, 2);


	getchar();
}