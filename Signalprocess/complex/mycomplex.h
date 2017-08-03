#pragma once
#include<math.h>
class complex{
public:
	double real;
	double image;
	double cabs();
	double phase();
	double square();
	complex(){ real = 0, image = 0; }
	complex(double x){ real = x, image = 0; }
	complex(double x, double y){
		real = x, image = y;}
};
complex con(complex a);   //复数的厄米共轭运算函数



complex operator+(complex  a, complex b);//重载+运算符

complex operator-(complex a, complex b); //重载-运算符
complex operator*(complex  a, complex b); //重载*运算符
complex operator*(complex a, double b);
complex operator*(double a, complex b);


complex operator/(complex a, double b); //重载/运算符

complex operator/(double a, complex b);

complex operator/(complex  a, complex b);

bool operator == (complex  a, complex b);


complex exp(complex x);

complex power(complex x, unsigned int i);
