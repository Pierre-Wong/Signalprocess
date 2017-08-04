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


double cyclerect(double x)
{
	if (x<0.5 && x>-0.5)
		return 1;
	else return 0;
}

double gibbs(double x)
{
	return fser(cyclerect, x, -1, 1).real;
}

int main()
{
	decom seq;

	for (int i = 0; i < 32; i++)
	{
		if (i < 5)
			seq.push_back(complex(1));
		else if (i < 28)
			seq.push_back(complex(0));
		else if (i < 32)
			seq.push_back(complex(1));

	}
	decom iff = fft(seq, 5);
	designal dis = designal(iff,1,32,0);
	//showd(dis.re,-31,31,-10,10);
	showc(gibbs, -10, 10, -2, 2);
	
	

	getchar();

	


	



}
