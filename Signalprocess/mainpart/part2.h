#pragma once
//#include "stdafx.h"
#include"../Def/Definitions.h"



//多项式的导数
poly derivative(poly f);

//多项式的导数
cpoly derivative(cpoly f);

//导数
func derivative(expfunction  f);

func derivative(fracfunction f);

vecfunc derivative(func f);


vecfunc derivative(vecfunc f);

//m次导数
vecfunc derivative(vecfunc f, unsigned int m);


vecfunc derivative(func f, unsigned int m);





/*
劈因子法求解多项式
基本原理是实数域R上每个多项式要么可以分解为一次和二次多项式的乘积
证明：由刘维尔定理，复数域C上每个非常数全纯函数（从而包括多项式函数）都不可能在整个复球面（包括无穷远点处）全纯。
这就是说，一个函数1/f(z)，在z=无穷时不发散，那么在复平面上必有一点z0使其发散，从而f(z)=0，
从而能推导出代数基本定理： 复数域上每个多项式都有解。
设R上有一多项式f，在C上，它有解a，那么因为f是实系数多项式，如果a不是实数，那么a的厄米共轭也是它的一个解，从而f=(x-a)(x-a_bar) g= (x+2re(a) x+|a|^2)-g
得证。
劈因子法就是找到一个二次多项式，不断逼近使得余式足够小，从而求解这个二次多项式即可得到高次多项式的解。
*/
comvec PolynomialEquation(poly A);






//求矩阵行列式的递归函数，按行（列）展开方式计算，简单易懂，注意，这种方法绝对不适用于大规模矩阵。
double detm(matrix a1);


//重载detm函数
complex detm(cmatrix a1);


//保证多项式p的最高此项非0
void preparepoly(poly &p);





/*
结式:
在高等代数中，m次多项式f=a_0X^m+a……+a_m  和n次多项式g=b_0X^n+……+b_n  的如下定义的行列式称为结式：
|a_0 …… …… a_m         |
|  ……                    |
|    ……                  |   n行
|      ……                |
|       a_0 …… …… a_m  |
|b_0 …… …… b_n         |
|  ……                    |
|    ……                  |   m行
|      ……                |
|       b_0 …… …… b_n  |
如果f和g有公共因子的话，这个行列式为0。
本书中此函数主要用于计算多项式和其导数的结式，如果为0，则证明此多项式有重根。
*/
double eliminant(poly f, poly g);


/*对复数重载*/

complex eliminant(cpoly f, cpoly g);









/*利用结式求最高重根重数，不断求导再求结式即可*/

unsigned int mulroot(poly a);



unsigned int mulroot(cpoly a);






/*
多项式赋值函数，更快的方法是秦九昭算法，这里直接赋值是为了让poly类的性质更清楚。
*/
complex value(poly f, complex x);

complex value(vecfunc f, complex x);

/*分式函数的幂*/
fracfunction power(fracfunction f, unsigned int n);

/*分式函数变量替换的函数，即用分式函数s替代fun中的变量，主要用于双线性不变法*/
fracfunction replace(fracfunction fun, fracfunction s);






/*给出高次方程的解和根重数，返回值为vecom类型vecom a=vector<vector<complex> >
a[i]表示i重特征根对应的comple组成的vector, a[0]无意义。
*/

vevecom solve(poly a);


vevecom solve(cpoly a);






/*
克拉默法则
矩阵方程Ax=C利用克拉默法则求解，这里为了方便运用公式，matrix的行和列与结式的实现是相反的
这里matrix是vector<列>，而每个列是一个vector<double>。
*/
//注意，这里是低维矩阵才使用克拉默法则，一般情况下则应该采用雅各比迭代等近似方法。
lfvector cramer(matrix a, lfvector c);

vector<complex> cramer(cmatrix a, vector<complex> c);


complex postivepower(complex x, int i);


/*阶乘函数*/
unsigned gamma(int i);

unsigned comb(int m, int n);

/*零输入响应的求解*/
vevecom zeroinput(poly linearsystem, lfvector initial);


/*冲击响应求解*/
vevecom impulse(poly y, poly x);


cpoly ve2cpoly(vevecom sing);




/*实多项式约分函数，为了应用劈因式法，尽管fracfuction定义为cpoly，但是这里要求多项式系数都是实的*/
ve3com reduce(fracfunction f);






/*实多项式的部分因子分解*/
vevecom decomp(fracfunction f);



/*正规化 Butterworth多项式*/
cpoly Butterworth(int n);

/* sinh 双曲正弦函数*/
double sinh(double x);

/* cosh 双曲余弦函数*/
double cosh(double x);
/*反双曲正弦函数*/
double arsinh(double x);

/*反双曲余弦函数*/
double arcosh(double x);



/*Chebyshev 多项式*/
cpoly Chebyshev(unsigned int n);

/*求值用函数*/
complex Chebyshev(unsigned int n, double fre);





/*

int main()
{
	
	cpoly a, b, c;
//	a.push_back(complex(0));
//	a.push_back(complex(0));
	a.push_back(complex(1));
	
	b.push_back(complex(0));
	b.push_back(complex(0));
//	b.push_back(complex(0.5));
//	b.push_back(complex(-1.5));
	b.push_back(complex(1));
	//b.push_back(complex(3));
	//b.push_back(complex(1));
	//b.push_back(complex(5));
	//b.push_back(complex(1));
	//b.push_back(complex(1));
	//fracfunction ss = reduce(fracfunction(a, b));
//	fracfunction xx = replace(fracfunction(a, b), fracfunction(a,b));
//	vevecom re = decomp(fracfunction(a, b));
	
	cpoly but = Chebyshev(10);

	/*
	double epsion = 0.15262;
	double phi = arsinh(1 / epsion) / 5;
	printf("%lf", sinh(phi));
	*/
	//vevecom ze = zeroinput(aa, bb);
/*
	poly aa, bb;
	aa.push_back(1);
	aa.push_back(3);
	aa.push_back(3);
	aa.push_back(1);
	bb.push_back(0);
	bb.push_back(1);
	bb.push_back(6);
	vevecom im = zeroinput(aa, bb);
*/	
/*	vecfunc s = derivative(func(expfunction(1),fracfunction(a, b)));
	for (int i = 0; i < s.size(); i++)
	{
		for (int j = 0; j < s[i].fra.deno.size(); j++)
			printf("%lf", s[i].fra.deno[j].real);
	}
	*/


	//c = derivative(c);
	//printf("%d", mulroot(a));
	//vevecom e = solve(a);
	//for (int i = 0; i < e[2].size(); i++)
	//	printf("%lf ", e[2][i].real);
	//vector<complex> y = impulse(a, b);
	//for (int i = 0; i < y.size(); i++)
	//	printf("%lf", y[i].real);
/*
	poly cc;
	cc.push_back(0);
	cc.push_back(3);
	cc.push_back(7);
	cc.push_back(5);
	cc.push_back(0);
	preparepoly(cc);
	*/
//	comvec dd=PolynomialEquation(cc);
	//for (int i = 0; i < dd.real.size(); i++)
	//	printf("%lf+%lf\n", dd.real[i],dd.image[i]);
	//for (int i = 0; i < c.size(); i++)
	/*matrix d = matrix(3, lfvector(3));
	d[0][0] = 1;
	d[0][1] = 1;
	d[1][0] = 0;
	d[1][1] = 2;
	d[2][2] = 1;
	c = cramer(d, a);*/
//	sort(a.begin(), a.end());
//	comvec D;
//	D = PolynomialEquation(a);
//	for (int i = 0; i < D.real.size();i++)
//		printf("%lf  \n", D.real[i]);
//	for (int i = 0; i < a.size();i++)
//	printf("%lf", a[i]);
	//matrix m=matrix(10, lfvector(10));
	/*
	getchar();
	return 0;
}
*/