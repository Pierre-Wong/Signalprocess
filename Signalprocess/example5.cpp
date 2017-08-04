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

//拉普拉斯变换
int main()
{
	cpoly nu, de;
	nu.push_back(2);
	nu.push_back(1);
	de.push_back(3);
	de.push_back(4);
	de.push_back(1);
	fracfunction lap = fracfunction(nu, de);
	vevecom dec= decomp(lap);
	vecfunc la = vevec2vecf(dec);

	getchar();


}