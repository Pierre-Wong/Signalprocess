#pragma once
#include "../complex/mycomplex.h"
#include "../Def/Definitions.h"
#include "../draw/draw.h"


complex fser(double(*ptr)(double), double t, double Tstart, double Tend);

complex fouierser(double(*ptr)(double), int k, double Tstart, double Tend);


complex fouier(double(*ptr)(double), double omega,                     //傅里叶变换函数
	double wide = 100);                                             //这是一个反常积分，为了运用计算机计算，需要限定范围，这里默认为-100到+100范围内。

complex DTFT(double(*ptr)(int), double omega, int st, int ed);             //离散时间傅里叶变换


complex DFT(double(*ptr)(int), int st, int ed, int k);            //离散傅里叶变换


decom DFT(decom x);             //离散傅里叶变换

vector<unsigned int> ibj(int number); //n阶倒序输出



decom fft(decom x, unsigned int n);



decom ifft(decom x, unsigned int n);


vector<vector<complex> > fft(CImg<double> image);


CImg<double> veve2im(vector<vector<complex> > veve);


CImg<double> ifft(vector<vector<complex> > ff);
