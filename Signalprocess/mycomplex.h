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
complex con(complex a);   //�����Ķ��׹������㺯��



complex operator+(complex  a, complex b);//����+�����

complex operator-(complex a, complex b); //����-�����
complex operator*(complex  a, complex b); //����*�����
complex operator*(complex a, double b);
complex operator*(double a, complex b);


complex operator/(complex a, double b); //����/�����

complex operator/(double a, complex b);

complex operator/(complex  a, complex b);

bool operator == (complex  a, complex b);


complex exp(complex x);

complex power(complex x, unsigned int i);
