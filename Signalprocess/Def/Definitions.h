#pragma once

#include<vector>
#include <math.h>
#include <algorithm>
#include"../complex/mycomplex.h"
#include<Windows.h>


extern double epison;
using namespace std;

class functionobject
{
public:

	virtual double operator() (double x) = 0;

};

class showfunction :public functionobject
{
public:
	double objectfunction(double x) { return x; };
	virtual double operator()(double x) {
		return objectfunction(x);
	}
};



typedef vector<double>  poly;   //����ʽ������������ʽ�洢���±����x��ָ����
typedef vector<complex> cpoly;

complex value(cpoly f, complex x);

cpoly operator +(cpoly a, cpoly b);           //����ʽ�����������

template <class T>cpoly operator *(T a, cpoly b);

template <class T> cpoly operator *(cpoly b, T a);


cpoly operator -(cpoly a, cpoly b);




cpoly operator *(cpoly a, cpoly b);



typedef vector<double>  lfvector;
typedef vector <vector<double> > matrix; //��˫vector�������������������ʱ��������һ����ά������ʹ�ã�
										 //������Ҫ���ö�����ϤC++ STL��vector��ʹ�ú������������ķ�����

typedef vector<vector<complex> >cmatrix; //�������󣬱���ֻ���������ϵͳ���и�����ֵ��ʱ��ʹ�ø�������

typedef vector < vector<complex> > vevecom;


typedef struct comvec
{
	lfvector real;
	lfvector image;
}comvec;

class expfunction        //ָ��������
{
public:
	complex expo;    //e��ָ��
	expfunction() {};
	expfunction(complex t) { expo = t; };
	expfunction(double t) { expo = complex(t); };
	complex val(complex t) { return exp(expo*t); }
};
class fracfunction     //��ʽ������
{
public:
	cpoly numer;     //���Ӷ���ʽ
	cpoly deno;      //��ĸ����ʽ
	fracfunction() { numer = cpoly(), deno = cpoly(1,complex(1)); }
	fracfunction(cpoly a, cpoly b) { numer = a, deno = b; };
	complex val(complex t) {
		if (numer.size() == 0)
			return 0;
		return value(numer, t) / value(deno, t);
	}
};

class func    //ָ�������ͷ�ʽ�����˻��ĺ����࣬����������˹��Z�任
{
public:
	expfunction ex;
	fracfunction fra;
	func() { };
	func(expfunction e, fracfunction f) {
		ex = e, fra = f;
	};
	func(expfunction e) { ex = e, fra = fracfunction(cpoly(1, 1), cpoly(1, 1)); };
	func(fracfunction f) { ex = expfunction(complex(0)); fra = f; }
	complex val(complex t) { return ((complex)(ex.val(t)) * (fra.val(t))); }
};

typedef vector<func> vecfunc; //ָ�������ͷ�ʽ�����˻��Ķ���ʽ



func operator *(func x, func y);

func operator *(expfunction x, fracfunction y);

func operator *(func x, fracfunction y);

func operator *(func x, expfunction y);



func operator /(func x, func y);

expfunction operator / (expfunction x, expfunction y);

fracfunction operator / (fracfunction x, fracfunction y);

expfunction operator * (expfunction x, expfunction y);

fracfunction operator * (fracfunction x, fracfunction y);

fracfunction operator * (fracfunction x, complex y);

fracfunction operator + (fracfunction x, fracfunction y);

vecfunc operator +(vecfunc x, vecfunc y);

vecfunc operator +(func x, func y);

poly operator +(poly a, poly b);             //����ʽ�����������

poly operator -(poly a, poly b);

template <class T> poly operator *(T a, poly b);


template <class T> poly operator *(poly b, T a);



poly operator *(poly a, poly b);


typedef vector<vevecom> ve3com;


#ifndef pie

#define pie 3.1415926535



#endif
typedef vector<complex> decom;


template <class T>cpoly operator *(T a, cpoly b)
{
	if (b.empty())
		return b;
	cpoly product;
	for (int i = 0; i < b.size(); i++)
	{
		product.push_back(a*b[i]);
	}
	return product;
}
template <class T> cpoly operator *(cpoly b, T a)
{
	if (b.empty())
		return b;
	cpoly product;
	for (int i = 0; i < b.size(); i++)
	{
		product.push_back(a*b[i]);
	}
	return product;
}


template <class T> poly operator *(T a, poly b)
{
	poly product;
	for (int i = 0; i < b.size(); i++)
	{
		product.push_back(a*b[i]);
	}
	return product;
}
template <class T> poly operator *(poly b, T a)
{
	poly product;
	for (int i = 0; i < b.size(); i++)
	{
		product.push_back(a*b[i]);
	}
	return product;
}

template<class T,class D>
class realpart
{
public:
	T*com;
	double operator() (D x)
	{
		return (*com)(x).real;
	}

};

template<class T,class D>
class abspart
{
public:
	T* com;
	double operator() (D x)
	{
		return (*com)(x).cabs();
	}

};

class vecfunc2show
{
public:
	vecfunc f;
	
	realpart<vecfunc2show,double> re;
	abspart<vecfunc2show,double> cabs;
	vecfunc2show() { re.com = this; cabs.com = this; };
	vecfunc2show(vecfunc f) { this->f = f; re.com = this; cabs.com = this;};
	complex operator ()(double x)
	{
		complex sum = complex(0);
		for (int i = 0; i < f.size(); i++)
		{
			sum = f[i].val(x)+sum;
		}
		return sum;
	}
	
};

class complexfunction
{
public:
	realpart<complexfunction, double> re;
	abspart<complexfunction, double> cabs;
	complexfunction() { re.com = this; cabs.com = this; };
	virtual complex operator() (double x) = 0;
};

//��ɢ�ź�
class designal
{
public:
	decom data;
	bool cycle;//�Ƿ������ں���
	int cytime;//����
	int initial;//�����ں���ʱ��vec[0]��Ӧ����ɢʱ��
	realpart<designal,int> re;
	abspart<designal,int> cabs;
	designal() { cycle = 0; cytime = 0; initial = 0; re.com = this; cabs.com = this;};
	//designal(decom data, bool cycle) { this->data = data; this->cycle = cycle; re.com = this; cabs.com = this;};
	designal(decom data, bool cycle, int cytime, int initial) { this->data = data; this->cycle = cycle; this->cytime = cytime; this->initial = initial; re.com = this; cabs.com = this; };

	complex operator ()(int x)
	{
		if (!cycle)
		{
			int y = x - initial;
			if (y < 0 || y >= data.size())
				return complex(0);
			else
				return data[y];
		}
		else
		{
			while (x <= 0)
			{
				x += cytime;
			}
			return(data[x%cytime]);
		}
	}


};



class powfunction        //ָ��������
{
public:
	complex r;    //e��ָ��
	powfunction() {};
	powfunction(complex t) { r = t; };
	powfunction(double t) { r = complex(t); };
	complex val(int t) { return power(r,t); }
};


class dfunc    //��ɢʱ������ϵͳ�Ľ�
{
public:
	powfunction ex;
	fracfunction fra;
	dfunc() { };
	dfunc(powfunction e, fracfunction f) {
		ex = e, fra = f;
	};
	dfunc(powfunction e) { ex = e, fra = fracfunction(cpoly(1, 1), cpoly(1, 1)); };
	dfunc(fracfunction f) { ex = powfunction(1); fra = f; }
	complex val(int t) { return ((complex)(ex.val(t)) * (fra.val(t))); }
};

typedef vector<dfunc> dvecfunc; //ָ�������ͷ�ʽ�����˻��Ķ���ʽ



dfunc operator *(dfunc x, dfunc y);

dfunc operator *(powfunction x, fracfunction y);

dfunc operator *(dfunc x, fracfunction y);

dfunc operator *(dfunc x, powfunction y);



dfunc operator /(dfunc x, dfunc y);

powfunction operator / (powfunction x, powfunction y);

fracfunction operator / (fracfunction x, fracfunction y);

powfunction operator * (powfunction x, powfunction y);

fracfunction operator * (fracfunction x, fracfunction y);

fracfunction operator * (fracfunction x, complex y);

fracfunction operator + (fracfunction x, fracfunction y);

dvecfunc operator +(dvecfunc x, dvecfunc y);

dvecfunc operator +(dfunc x, dfunc y);