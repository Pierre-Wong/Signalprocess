#include<vector>
#include <math.h>
#include <algorithm>
#include"../complex/mycomplex.h"
#include<Windows.h>
#include "Definitions.h"

double epison=0.00001;
using namespace std;





complex value(cpoly f, complex x)
{
	if (f.empty())
		return 0;
	complex sum = complex(0);
	for (unsigned int i = 0; i < f.size(); i++)
		sum = f[i] * power(x, i) + sum;
	return sum;
}

cpoly operator +(cpoly a, cpoly b)             //多项式的运算符重载
{
	if (a.empty() || b.empty())
		return a.empty() ? b : a;
	cpoly sum;
	for (int i = 0; i < a.size() || i < b.size(); i++)
	{
		sum.push_back((i < a.size() ? a[i] : complex(0)) + (i < b.size() ? b[i] : complex(0)));//在C++中，? :是唯一的一个三目运算符
																							   //a?b:c 的意义为，如果a为真，则执行b，否则执行c。
	}
	return sum;
}


cpoly operator -(cpoly a, cpoly b)
{
	if (a.empty())
	{
		cpoly c = b*-1.0;
		return c;
	}
	if (b.empty())
		return a;
	cpoly diff;
	for (int i = 0; i < a.size() || i < b.size(); i++)
	{
		diff.push_back((i < a.size() ? a[i] : complex(0)) - (i < b.size() ? b[i] : complex(0)));
	}
	return diff;
}




cpoly operator *(cpoly a, cpoly b)
{
	if (a.empty() || b.empty())
		return cpoly();
	cpoly c, d;

	for (int i = 0; i < a.size(); i++)
	{
		d = a[i] * b;
		for (int j = 0; j < i; j++)
			d.insert(d.begin(), 0);
		c = c + d;

	}
	return c;
}








expfunction operator * (expfunction x, expfunction y)
{
	return(expfunction(x.expo + y.expo));
}

fracfunction operator * (fracfunction x, fracfunction y)
{
	if (x.numer.empty() || y.numer.empty())
		return fracfunction();
	return(fracfunction(x.numer*y.numer, x.deno*y.deno));
}

fracfunction operator * (fracfunction x, complex y)
{
	if (x.numer.empty() || y == 0)
		return fracfunction();
	return(fracfunction(x.numer*y, x.deno));
}

func operator *(func x, func y)
{
	return(func(x.ex*y.ex, x.fra*y.fra));
}
func operator *(expfunction x, fracfunction y)
{
	return(func(x, y));
}
func operator *(func x, fracfunction y)
{
	return(func(x.ex, y*x.fra));
}
func operator *(func x, expfunction y)
{
	return(func(x.ex*y, x.fra));
}
expfunction operator / (expfunction x, expfunction y)
{
	return(expfunction(x.expo - y.expo));
}
fracfunction operator / (fracfunction x, fracfunction y)
{
	return(fracfunction(x.numer*y.deno, x.deno * y.numer));

}
func operator /(func x, func y)
{
	return(func(x.ex / y.ex, x.fra / y.fra));
}
fracfunction operator + (fracfunction x, fracfunction y)
{
	if (x.numer.empty() || y.numer.empty())
		return x.numer.empty() ? y : x;
	return(fracfunction(x.numer*y.deno + y.numer*x.deno, x.deno * y.deno));
}

vecfunc operator +(vecfunc x, vecfunc y)
{
	vecfunc z = x;
	for (int j = 0; j < y.size(); j++)
	{
		bool find = 0;
		for (int i = 0; i < x.size(); i++)
		{
			if (x[i].ex.expo == y[j].ex.expo)
				z[i].fra = z[i].fra + y[j].fra;
			find = 1;
			break;
		}
		if (!find)
			z.push_back(y[j]);
	}
	return z;
}
vecfunc operator +(func x, func y)
{
	if (!x.fra.numer.empty() && !y.fra.numer.empty())
		return (vecfunc(1, x) + vecfunc(1, y));
	else if (x.fra.numer.empty())
		return vecfunc(1, y);
	else return vecfunc(1, x);
}


poly operator +(poly a, poly b)             //多项式的运算符重载
{
	poly sum;
	for (int i = 0; i < a.size() || i < b.size(); i++)
	{
		sum.push_back((i < a.size() ? a[i] : 0) + (i < b.size() ? b[i] : 0));//在C++中，? :是唯一的一个三目运算符
																			 //a?b:c 的意义为，如果a为真，则执行b，否则执行c。
	}
	return sum;
}
poly operator -(poly a, poly b)
{
	poly diff;
	for (int i = 0; i < a.size() || i < b.size(); i++)
	{
		diff.push_back((i < a.size() ? a[i] : 0) - (i < b.size() ? b[i] : 0));
	}
	return diff;
}




poly operator *(poly a, poly b)
{
	poly c, d;
	for (int i = 0; i < a.size(); i++)
	{
		d = a[i] * b;
		for (int j = 0; j < i; j++)
			d.insert(d.begin(), 0);
		c = c + d;

	}
	return c;
}


powfunction operator * (powfunction x, powfunction y)
{
	return(powfunction(x.r*y.r));
}

dfunc operator *(powfunction x, fracfunction y)
{
	return(dfunc(x, y));
}
dfunc operator *(dfunc x, powfunction y)
{
	return(dfunc(x.ex*y, x.fra));
}

powfunction operator / (powfunction x, powfunction y)
{
	return(powfunction(x.r /y.r));
}


dfunc operator *(dfunc x, dfunc y)
{
	return(dfunc(x.ex*y.ex, x.fra*y.fra));
}
/*
dfunc operator *(powfunction x, fracfunction y)
{
	return(dfunc(x, y));
}
*/
dfunc operator *(dfunc x, fracfunction y)
{
	return(dfunc(x.ex, y*x.fra));
}
/*
dfunc operator *(dfunc x, powfunction y)
{
	return(dfunc(x.ex*y, x.fra));
}
*/

dfunc operator /(dfunc x, dfunc y)
{
	return(dfunc(x.ex / y.ex, x.fra / y.fra));
}


dvecfunc operator +(dvecfunc x, dvecfunc y)
{
	dvecfunc z = x;
	for (int j = 0; j < y.size(); j++)
	{
		bool find = 0;
		for (int i = 0; i < x.size(); i++)
		{
			if (x[i].ex.r == y[j].ex.r)
				z[i].fra = z[i].fra + y[j].fra;
			find = 1;
			break;
		}
		if (!find)
			z.push_back(y[j]);
	}
	return z;
}
dvecfunc operator +(dfunc x, dfunc y)
{
	if (!x.fra.numer.empty() && !y.fra.numer.empty())
		return (dvecfunc(1, x) + dvecfunc(1, y));
	else if (x.fra.numer.empty())
		return dvecfunc(1, y);
	else return dvecfunc(1, x);
}

std::ostream & operator<<(std::ostream &out, poly &obj)
{
	bool zero = true;
	for (int i = 0; i < obj.size(); i++)
	{
		if (obj[i] != 0)
		{
			if (!zero)
				out << "+";
			out << "(" << obj[i] << ")";
			if (i)
			{
				out << "x^" << i;
			}
			zero = false;
		}
	}
	return out;
}

std::ostream & operator<<(std::ostream &out, cpoly &obj)
{
	bool zero = true;
	for (int i = 0; i < obj.size(); i++)
	{
		if (!(obj[i] ==complex(0)))
		{
			if (!zero)
				out << "+";
			out << "("<<obj[i]<<")";
			if (i)
			{
				out << "x^" << i;
			}
			zero = false;
		}
	}
	return out;
}

std::ostream & operator<<(std::ostream &out, expfunction &obj)
{
	if (!(obj.expo == complex(0)))
	{
		out << "e^" << obj.expo<<"x";
	}
	else
		out << 1;
	return out;
}

std::ostream & operator<<(std::ostream &out, fracfunction &obj)
{
	if ((obj.deno.size() == 1)&&(obj.deno[0] == complex(1)))
	{
		
		out << obj.numer;
	}
	else
	{
		out << "(" << obj.numer << ")/(" << obj.deno << ")";
	}
	return out;
}

std::ostream & operator<<(std::ostream &out, func &obj)
{
	out << obj.fra << " * " << obj.ex;
	return out;
}
std::ostream & operator<<(std::ostream &out, vecfunc &obj)
{
	bool zero = true;
	for (int i = 0; i < obj.size(); i++)
	{
		if (!zero)
		{
			out << " + ";
		}
		out << obj[i];
		zero = false;
	}
	return out;
}