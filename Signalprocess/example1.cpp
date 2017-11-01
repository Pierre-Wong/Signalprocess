
#include <stdio.h>  
#include <stdlib.h>  
#include"./mainpart/part2.h"
#include"./draw/draw.h"

#include"./basicfunction/basicfunction.h"
#include"./fourier/fourier.h"


using namespace std;

//对应第二章，信号的matlab表示

//信号Ae^(ax)
double expfun(double x)
{
	double A = 1;
	double a = 1;
	return A*exp(a*x);
}

double expseq(int x)
{
	double a = 2;
	return pow(a, x);
}



int main()
{
	
	//showc(expfun, -10, 10, 0, 100);
	showd(expseq, -10, 10, 0, 1000);
	getchar();
}