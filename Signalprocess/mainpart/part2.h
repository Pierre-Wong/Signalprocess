#pragma once
//#include "stdafx.h"
#include"../Def/Definitions.h"



//����ʽ�ĵ���
poly derivative(poly f);

//����ʽ�ĵ���
cpoly derivative(cpoly f);

//����
func derivative(expfunction  f);

func derivative(fracfunction f);

vecfunc derivative(func f);


vecfunc derivative(vecfunc f);

//m�ε���
vecfunc derivative(vecfunc f, unsigned int m);


vecfunc derivative(func f, unsigned int m);





/*
�����ӷ�������ʽ
����ԭ����ʵ����R��ÿ������ʽҪô���Էֽ�Ϊһ�κͶ��ζ���ʽ�ĳ˻�
֤��������ά������������C��ÿ���ǳ���ȫ���������Ӷ���������ʽ�������������������������棨��������Զ�㴦��ȫ����
�����˵��һ������1/f(z)����z=����ʱ����ɢ����ô�ڸ�ƽ���ϱ���һ��z0ʹ�䷢ɢ���Ӷ�f(z)=0��
�Ӷ����Ƶ��������������� ��������ÿ������ʽ���н⡣
��R����һ����ʽf����C�ϣ����н�a����ô��Ϊf��ʵϵ������ʽ�����a����ʵ������ôa�Ķ��׹���Ҳ������һ���⣬�Ӷ�f=(x-a)(x-a_bar) g= (x+2re(a) x+|a|^2)-g
��֤��
�����ӷ������ҵ�һ�����ζ���ʽ�����ϱƽ�ʹ����ʽ�㹻С���Ӷ����������ζ���ʽ���ɵõ��ߴζ���ʽ�Ľ⡣
*/
comvec PolynomialEquation(poly A);






//���������ʽ�ĵݹ麯�������У��У�չ����ʽ���㣬���׶���ע�⣬���ַ������Բ������ڴ��ģ����
double detm(matrix a1);


//����detm����
complex detm(cmatrix a1);


//��֤����ʽp����ߴ����0
void preparepoly(poly &p);





/*
��ʽ:
�ڸߵȴ����У�m�ζ���ʽf=a_0X^m+a����+a_m  ��n�ζ���ʽg=b_0X^n+����+b_n  �����¶��������ʽ��Ϊ��ʽ��
|a_0 ���� ���� a_m         |
|  ����                    |
|    ����                  |   n��
|      ����                |
|       a_0 ���� ���� a_m  |
|b_0 ���� ���� b_n         |
|  ����                    |
|    ����                  |   m��
|      ����                |
|       b_0 ���� ���� b_n  |
���f��g�й������ӵĻ����������ʽΪ0��
�����д˺�����Ҫ���ڼ������ʽ���䵼���Ľ�ʽ�����Ϊ0����֤���˶���ʽ���ظ���
*/
double eliminant(poly f, poly g);


/*�Ը�������*/

complex eliminant(cpoly f, cpoly g);









/*���ý�ʽ������ظ������������������ʽ����*/

unsigned int mulroot(poly a);



unsigned int mulroot(cpoly a);






/*
����ʽ��ֵ����������ķ������ؾ����㷨������ֱ�Ӹ�ֵ��Ϊ����poly������ʸ������
*/
complex value(poly f, complex x);

complex value(vecfunc f, complex x);

/*��ʽ��������*/
fracfunction power(fracfunction f, unsigned int n);

/*��ʽ���������滻�ĺ��������÷�ʽ����s���fun�еı�������Ҫ����˫���Բ��䷨*/
fracfunction replace(fracfunction fun, fracfunction s);






/*�����ߴη��̵Ľ�͸�����������ֵΪvecom����vecom a=vector<vector<complex> >
a[i]��ʾi����������Ӧ��comple��ɵ�vector, a[0]�����塣
*/

vevecom solve(poly a);


vevecom solve(cpoly a);






/*
����Ĭ����
���󷽳�Ax=C���ÿ���Ĭ������⣬����Ϊ�˷������ù�ʽ��matrix���к������ʽ��ʵ�����෴��
����matrix��vector<��>����ÿ������һ��vector<double>��
*/
//ע�⣬�����ǵ�ά�����ʹ�ÿ���Ĭ����һ���������Ӧ�ò����Ÿ��ȵ����Ƚ��Ʒ�����
lfvector cramer(matrix a, lfvector c);

vector<complex> cramer(cmatrix a, vector<complex> c);


complex postivepower(complex x, int i);


/*�׳˺���*/
unsigned gamma(int i);

unsigned comb(int m, int n);

/*��������Ӧ�����*/
vevecom zeroinput(poly linearsystem, lfvector initial);


/*�����Ӧ���*/
vevecom impulse(poly y, poly x);


cpoly ve2cpoly(vevecom sing);




/*ʵ����ʽԼ�ֺ�����Ϊ��Ӧ������ʽ��������fracfuction����Ϊcpoly����������Ҫ�����ʽϵ������ʵ��*/
ve3com reduce(fracfunction f);






/*ʵ����ʽ�Ĳ������ӷֽ�*/
vevecom decomp(fracfunction f);



/*���滯 Butterworth����ʽ*/
cpoly Butterworth(int n);

/* sinh ˫�����Һ���*/
double sinh(double x);

/* cosh ˫�����Һ���*/
double cosh(double x);
/*��˫�����Һ���*/
double arsinh(double x);

/*��˫�����Һ���*/
double arcosh(double x);



/*Chebyshev ����ʽ*/
cpoly Chebyshev(unsigned int n);

/*��ֵ�ú���*/
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