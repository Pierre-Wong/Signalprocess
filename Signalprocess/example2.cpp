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



//对应第三章，利用matlab进行时域分析
int main()
{
	
	poly aa, bb;
	aa.push_back(2);
	aa.push_back(3);
	aa.push_back(1);

	bb.push_back(0);
	bb.push_back(0.5);
	
	//vecfunc im = impulse(aa, bb);
	
	//vecfunc2show showfun = vecfunc2show(im);
	//showc(showfun.re, 0, 10, -10, 10);
	vevecom im = dzeroinput(aa, bb);

	/*
	poly aa, bb;
	aa.push_back(4);
	aa.push_back(4);
	aa.push_back(1);

	bb.push_back(1);
	bb.push_back(2);

	vevecom im = zeroinput(aa, bb);
	*/


	getchar();
}