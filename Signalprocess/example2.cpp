#include <stdio.h>  
#include <stdlib.h>  


#include"./mainpart/part2.h"
#include"./draw/draw.h"

#include"./basicfunction/basicfunction.h"
#include"./fourier/fourier.h"



//��Ӧ�����£�����matlab����ʱ�����
int main()
{
	/*
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
	*/
	
	poly aa, bb;
	aa.push_back(8);
	aa.push_back(6);
	aa.push_back(1);

	bb.push_back(1);
	bb.push_back(2);


	cout <<"aa="<< aa<<endl;
	vecfunc im = zeroinput(aa, bb);
	
	
	cout <<"zeroinput="<< im<<endl;
	vecfunc2show showfun = vecfunc2show(im);
	showc(showfun.re, 0, 3, -10, 10);
	
	
	getchar();
}