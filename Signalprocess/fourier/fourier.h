#pragma once
#include "../complex/mycomplex.h"
#include "../Def/Definitions.h"
#include "../draw/draw.h"


complex fser(double(*ptr)(double), double t, double Tstart, double Tend);

complex fouierser(double(*ptr)(double), int k, double Tstart, double Tend);


complex fouier(double(*ptr)(double), double omega,                     //����Ҷ�任����
	double wide = 100);                                             //����һ���������֣�Ϊ�����ü�������㣬��Ҫ�޶���Χ������Ĭ��Ϊ-100��+100��Χ�ڡ�

complex DTFT(double(*ptr)(int), double omega, int st, int ed);             //��ɢʱ�丵��Ҷ�任


complex DFT(double(*ptr)(int), int st, int ed, int k);            //��ɢ����Ҷ�任


decom DFT(decom x);             //��ɢ����Ҷ�任

vector<unsigned int> ibj(int number); //n�׵������



decom fft(decom x, unsigned int n);



decom ifft(decom x, unsigned int n);


vector<vector<complex> > fft(CImg<double> image);


CImg<double> veve2im(vector<vector<complex> > veve);


CImg<double> ifft(vector<vector<complex> > ff);
