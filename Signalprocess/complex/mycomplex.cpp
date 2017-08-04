#include"mycomplex.h"
complex con(complex a)    //复数的厄米共轭运算函数
{
	complex cc;
	cc.real = a.real;
	cc.image = -1 * a.image;
	return cc;
}
double complex::cabs()
{
	return sqrt((real)*(real)+(image)*(image));
}
double complex::phase()
{
	return atan(image / real);
}
double complex::square()
{
	return real*real + image*image;
}
complex operator+(complex  a, complex b)//重载+运算符
{
	return complex(a.real + b.real, a.image + b.image);
}
complex operator-(complex a, complex b)//重载-运算符
{
	return complex(a.real - b.real, a.image - b.image);
}
complex operator*(complex  a, complex b)//重载*运算符
{
	return complex(a.real*b.real - a.image*b.image, a.real*b.image + a.image*b.real);
}
complex operator*(complex a, double b)
{
	return complex(a.real*b, a.image*b);
}
complex operator*(double a, complex b)
{
	return complex(b.real*a, b.image*a);
}

complex operator/(complex a, double b)//重载/运算符
{
	return complex(a.real / b, a.image / b);
}
complex operator/(double a, complex b)
{
	return (a)*(con(b)) / b.square();
}
complex operator/(complex  a, complex b)
{
	return a*((double)1 / b);
}
bool operator == (complex  a, complex b)
{
	return (a.real == b.real) && (a.image == b.image);
}

complex exp(complex x)
{
	complex euler;
	euler.real = exp(x.real)*cos(x.image);
	euler.image = exp(x.real)*sin(x.image);
	return euler;
}

complex pluspower(complex x, unsigned int i)
{
	if (i == 0)
		return 1;
	else
		return x*pluspower(x, i - 1);
}

complex power(complex x, int i)
{
	if (i >= 0)
		return pluspower(x, i);
	else
		return 1 / pluspower(x, -i);
}

