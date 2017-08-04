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



typedef vector<double>  poly;   //多项式以向量容器方式存储，下标代表x的指数。
typedef vector<complex> cpoly;

complex value(cpoly f, complex x);

cpoly operator +(cpoly a, cpoly b);           //多项式的运算符重载

template <class T>cpoly operator *(T a, cpoly b);

template <class T> cpoly operator *(cpoly b, T a);


cpoly operator -(cpoly a, cpoly b);




cpoly operator *(cpoly a, cpoly b);



typedef vector<double>  lfvector;
typedef vector <vector<double> > matrix; //用双vector来构造矩阵容器，访问时把它当做一个二维数组来使用，
										 //这里主要是让读者熟悉C++ STL中vector的使用和用其替代数组的方法。

typedef vector<vector<complex> >cmatrix; //复数矩阵，本书只在求解线性系统（有复特征值）时才使用复数矩阵。

typedef vector < vector<complex> > vevecom;


typedef struct comvec
{
	lfvector real;
	lfvector image;
}comvec;

class expfunction        //指数函数类
{
public:
	complex expo;    //e上指数
	expfunction() {};
	expfunction(complex t) { expo = t; };
	expfunction(double t) { expo = complex(t); };
	complex val(complex t) { return exp(expo*t); }
};
class fracfunction     //分式函数类
{
public:
	cpoly numer;     //分子多项式
	cpoly deno;      //分母多项式
	fracfunction() { numer = cpoly(), deno = cpoly(1,complex(1)); }
	fracfunction(cpoly a, cpoly b) { numer = a, deno = b; };
	complex val(complex t) {
		if (numer.size() == 0)
			return 0;
		return value(numer, t) / value(deno, t);
	}
};

class func    //指数函数和分式函数乘积的函数类，用于拉普拉斯和Z变换
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

typedef vector<func> vecfunc; //指数函数和分式函数乘积的多项式



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

poly operator +(poly a, poly b);             //多项式的运算符重载

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

//离散信号
class designal
{
public:
	decom data;
	bool cycle;//是否是周期函数
	int cytime;//周期
	int initial;//非周期函数时，vec[0]对应的离散时间
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



class powfunction        //指数函数类
{
public:
	complex r;    //e上指数
	powfunction() {};
	powfunction(complex t) { r = t; };
	powfunction(double t) { r = complex(t); };
	complex val(int t) { return power(r,t); }
};


class dfunc    //离散时间线性系统的解
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

typedef vector<dfunc> dvecfunc; //指数函数和分式函数乘积的多项式



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