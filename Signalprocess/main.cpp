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
//#include "msp.h"
using namespace std;

//#define pi 3.14159265358979323846





double fre(int x)
{
	return fouierser(juxing, x, -2, 2).real;
}


double db(double x)
{
	return log(DTFT(blackmanwin, x*pie, -100, 100).cabs()) / log(10) * 20;
}



double df(int x)
{
	return DFT(rec, -10,10, x).real;
}


double re(double x)
{
	//return fouier(juxing,x).real;
	//return fser(juxing, x, -2, 2).real;
	return fouier(juxing, x).real;
	//return DTFT(sinc, x, -1000, 1000).real;
	//return sinc(x);
	//return fouier(sinc, x).real;
	//return DTFT(kaiserwin, x, -100, 100).cabs();
}




int main()
{
	//showc(showva, -10, 10, 0, 12);
	//showd(df, -10, 10, -30, 30);
	//showc(db, 0, 1, -1000, 0);
	//showc(db, 0, 1, -150, 0);

	
	poly aa, bb;
	aa.push_back(1);
	aa.push_back(3);
	aa.push_back(3);
	aa.push_back(1);
	bb.push_back(0);
	bb.push_back(1);
	bb.push_back(6);
	vevecom im = zeroinput(aa, bb);
	
	showc(re, -5, 5, -30, 30);
	/*showfunction a;
	printf("%lf", a(1));
	getchar();
	*/
/*	printf("lww\n");
	
	// 定义一个每个颜色8位(bit)的640x400的彩色图像   
	CImg<unsigned char> img(640, 400, 1,3);
	img._is_shared = 1;
	//将像素值设为0（黑色）   
	img.fill(255);
	// 定义一个紫色   
	const unsigned char purple[] = { 255, 0, 255 };
	// 在坐标(100, 100)处画一个紫色的“Hello world”
	Point a;
	a.x = 100;
	a.y = 100;
	putText(img,"hello world", a, 1, 1, Scalar(255,0,255));
	//img.draw_text(100, 100, "Hello World", purple,0,1,10);
	// 在一个标题为“My first CImg code”的窗口中显示这幅图像   
	
	//const unsigned char color[] = { 0, 0, 0 };
	img.draw_line(100, 100, 200, 200, purple,3);
	img.draw_circle(200, 200, 10, purple);
	img.display("My first CImg code");
	
	//showc(db, 0, 1, -120, 0);
	//showd(df,0,10,-10,10);
	*/
//	vector<unsigned int> dd =ibj (3);

	/*

	decom x ;
	x.push_back(complex(1));
	x.push_back(complex(2));
	x.push_back(complex(3));
	x.push_back(complex(4));
	x.push_back(complex(5));
	x.push_back(complex(6));
	x.push_back(complex(7));
	x.push_back(complex(8));
	x.push_back(complex(1));
	x.push_back(complex(2));
	x.push_back(complex(3));
	x.push_back(complex(4));
	x.push_back(complex(5));
	x.push_back(complex(6));
	x.push_back(complex(7));
	x.push_back(complex(8));
	
	
	decom y = fft(x, 4);
	decom z = DFT(x);
	decom zz = ifft(y, 4);
	
	*/
	

	/*
	CImg<double> image("lena.bmp"), visu(512, 512, 1, 1, 0);
	
	double *p1 = image.data(0, 0, 0, 0), *p2 = image.data(0, 0, 0, 1), *p3 = image.data(0, 0, 0, 2), *p4 = visu.data(0, 0, 0, 0);
	for (unsigned int i = 0; i < visu._width*visu._height; i++) {
		*(p4++) = (*(p1++) + *(p2++) + *(p3++)) / 3;
	//	*(p4++) = 0;
	}
	
	vector<vector<complex> >ve = fft(visu);
	CImg<double> im = ifft(ve);
	bool d = 1;
	*/
	// 定义一个每个颜色8位(bit)的640x400的彩色图像   
/*	CImg<unsigned char> img(640, 400, 1, 3);
	img._is_shared = 1;
	//将像素值设为0（黑色）   
	img.fill(255);
	// 定义一个紫色   
	const unsigned char purple[] = { 255, 0, 255 };
	// 在坐标(100, 100)处画一个紫色的“Hello world”
/*	Point a;
	a.x = 100;
	a.y = 100;
	putText(img, "hello world", a, 1, 1, Scalar(255, 0, 255));
	//img.draw_text(100, 100, "Hello World", purple,0,1,10);
	// 在一个标题为“My first CImg code”的窗口中显示这幅图像   

	//const unsigned char color[] = { 0, 0, 0 };
	img.draw_line(100, 100, 200, 200, purple, 3);*/

//	getchar();    
/*	//showc(re, -2, 2, -20, 20);
	CImg<unsigned char> img(600, 600, 1, 3,0);
	img._is_shared = 1;
	img.fill(255);
	CImg<unsigned char>img2=img.resize(300, 300);
	//img.resize_halfXY();
	img2.display("a");
	*/



/*	record();
	double T1=0;
	vector<double> energy = fra(buffer, BUFFER_SIZE);
	vector<double> zeros = ze(buffer, BUFFER_SIZE);
	for (int i = 0; i < energy.size(); i++)
		T1 += energy[i] / energy.size();
	
	for (start = 0; start < energy.size(); start++)
	{
		if (energy[start]>T1)
			break;
	}
	if (start>0)
		T1 = energy[0];
	for (start = 0; start < energy.size(); start++)
	{
		if (energy[start]>2*T1)
			break;
	}
	for (endd = energy.size()-1; endd>=0; endd--)
	{
		if (energy[endd]>2 * T1)
			break;
	}
	

	for (int i = start; i <= 0; i--)
	{
		if (zeros[i] < 10)
		{
			
		
				start = i;
				break;
		}
	}

	for (int i = endd; i <zeros.size(); i++)
	{
		if (zeros[i] > 10)
		{
			
				endd = i;
				break;
		}
	}
	//start = 0; endd = energy.size();
	buffer2 = new unsigned char[(endd - start) * 80];
	int count = 0;
	
	for (int i = start * 80; i < endd * 80; i++)
	{
		buffer2[count] = buffer[i];
		count++;
	}

//	showd(rsf, 0, 319, -1, 1);
//	sound(buffer, 0, 15999, 0, 255);
	lpc(buffer2, (endd - start) * 80);
	play();
	
	


	
	
	*/

	
	return 0;
}

