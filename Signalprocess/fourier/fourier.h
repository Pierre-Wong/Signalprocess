#pragma once
#include "../complex/mycomplex.h"
#include "../Def/Definitions.h"
#include "../draw/draw.h"

template<class T>
complex fser(T ptr , double t, double Tstart, double Tend);

template<class T>
complex fouierser(T ptr, int k, double Tstart, double Tend);

template<class T>
complex fouier(T ptr , double omega,                     //����Ҷ�任����
	double wide = 100);                                             //����һ���������֣�Ϊ�����ü�������㣬��Ҫ�޶���Χ������Ĭ��Ϊ-100��+100��Χ�ڡ�

template<class T>
complex DTFT( T ptr, double omega, int st, int ed);             //��ɢʱ�丵��Ҷ�任

template<class T>
complex DFT(T ptr, int st, int ed, int k);            //��ɢ����Ҷ�任





decom DFT(decom x);             //��ɢ����Ҷ�任

vector<unsigned int> ibj(int number); //n�׵������



decom fft(decom x, unsigned int n);



decom ifft(decom x, unsigned int n);


vector<vector<complex> > fft(CImg<double> image);


CImg<double> veve2im(vector<vector<complex> > veve);


CImg<double> ifft(vector<vector<complex> > ff);

template<class T>
complex DFT(T ptr , int st, int ed, int k)             //��ɢ����Ҷ�任
{
	complex dt = complex(0);
	for (int x = st; x <ed; x++)
	{
		dt = ptr(x) * exp(complex(0, -1 * 2 * pie*x*k / (ed - st))) + dt;
	}
	return dt;
}

template<class T>
complex DTFT(T ptr , double omega, int st, int ed)             //��ɢʱ�丵��Ҷ�任
{
	complex dt = complex(0);
	for (int x = st; x <= ed; x++)
	{
		dt = ptr(x) * exp(complex(0, -1 * omega*x)) + dt;
	}
	return dt;
}


template<class T>
complex fouier( T ptr, double omega,                     //����Ҷ�任����
	double wide)                                             //����һ���������֣�Ϊ�����ü�������㣬��Ҫ�޶���Χ������Ĭ��Ϊ-100��+100��Χ�ڡ�
{
	complex cof = complex(0);
	double ree = 0;
	double imm = 0;
	for (int x = -1 * (int)(wide / 0.01); x < (int)(wide / 0.01); x++)
	{
		ree = ptr(x*0.01)*exp(complex(0, -1 * 2 * pie* x*0.01*omega)).real*0.01 + ree;
		imm = ptr(x*0.01)*exp(complex(0, -1 * 2 * pie* x*0.01*omega)).image*0.01 + imm;
		//cof = ptr(x*0.01)*exp(complex(0, -1 * 6.28* x*0.01*omega))*0.01 + cof;
	}
	cof.real = ree;
	cof.image = imm;
	return cof;
}

template<class T>
complex fouierser(T ptr , int k, double Tstart, double Tend)
{
	complex a;
	complex ret = complex(0);
	a = complex(0);
	for (double x = Tstart; x <= Tend; x = x + 0.01)
		a = ptr(x)*exp(complex(0, -1 * k * 2 * pie / (Tend - Tstart)*x))*0.01 + a;
	return a;

}



#define num 30
template<class T>
complex fser(T ptr, double t, double Tstart, double Tend)
{

	complex a[2 * num + 1];
	complex ret = complex(0);
	for (int k = -1 * num; k <= num; k++)
	{
		a[k + num] = complex(0);
		for (double x = Tstart; x <= Tend; x = x + 0.01)
			a[k + num] = ptr(x)*exp(complex(0, -1 * k * 2 * pie / (Tend - Tstart)*x))*0.01 + a[k + num];
		ret = ret + a[k + num] * exp(complex(0, k*6.28 / (Tend - Tstart)*t)) / (Tend - Tstart);
	}
	return ret;
}