#pragma once

#include<vector>
#include <math.h>
#include <algorithm>
#include"../complex/mycomplex.h"
#include<Windows.h>


extern double epison;
using namespace std;

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
	fracfunction() { numer = cpoly(), deno = cpoly(); }
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

expfunction operator * (expfunction x, expfunction y);

fracfunction operator * (fracfunction x, fracfunction y);

fracfunction operator * (fracfunction x, complex y);

func operator *(func x, func y);

func operator *(expfunction x, fracfunction y);

func operator *(func x, fracfunction y);

func operator *(func x, expfunction y);

expfunction operator / (expfunction x, expfunction y);

fracfunction operator / (fracfunction x, fracfunction y);

func operator /(func x, func y);

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

//H文件内容

#endif
typedef vector<complex> decom;