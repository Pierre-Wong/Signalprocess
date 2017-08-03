#pragma once

#include<vector>
#include <math.h>
#include <algorithm>
#include"../complex/mycomplex.h"
#include<Windows.h>


extern double epison;
using namespace std;

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
	fracfunction() { numer = cpoly(), deno = cpoly(); }
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

poly operator +(poly a, poly b);             //����ʽ�����������

poly operator -(poly a, poly b);

template <class T> poly operator *(T a, poly b);


template <class T> poly operator *(poly b, T a);



poly operator *(poly a, poly b);


typedef vector<vevecom> ve3com;


#ifndef pie

#define pie 3.1415926535

//H�ļ�����

#endif
typedef vector<complex> decom;