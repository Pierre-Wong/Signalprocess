#include"basicfunction.h"
#include<cmath>
#include<math.h>
#ifndef pie

#define pie 3.1415926535

#endif // !pie
double show(double x)
{
	complex c = complex(0, 1);
	return(exp(x*c).real);
}
complex frac(complex s)
{
	return  (complex(1) + s*s) / (complex(1) - s*s);
}
double sfra(double x)
{
	return frac(complex(0, x)).cabs();
}
double gauss(double x)
{
	return 10 * exp(-x*x);
}

double f(double x)
{
	return cos(x);
}



double juxing(double x)
{
	if (abs(x) < 1)
		return 10;
	else return 0;
}

double delta(int i)
{
	if (i == 2)
		return 1;
	else return 0;
}
double sinc(double x)               //sinc函数，注意x=0时要使用洛必达法则直接算出结果
{
	if (abs(x) < 0.001)
		return 1;
	else return sin(x) / x;
}
double sinc(int x)
{
	if (x == 0)
		return 1;
	else return sin(x) / x;
}

#define N 51 
#define omegac 0.5
double rectwindow(int n)
{
	if (n >= 0 && n < N)
	{
		if ((n - (N - 1)*0.5) != 0)
			return omegac*sin(omegac*pie*(n - (N - 1)*0.5)) / (omegac*pie*(n - (N - 1)*0.5));
		else return omegac;
	}
	else return 0;

}
double bartlettwin(int n)
{
	if (n >= 0 && n <= (N - 1)*0.5)
		return rectwindow(n) * 2.0 * n / (N - 1);
	else return rectwindow(n)*(2 - 2 * n / (N - 1));
}
double hanningwin(int n)
{
	return((0.5 - 0.5*cos((2 * pie*n) / (N - 1)))*rectwindow(n));
}
double hammingwin(int n)
{
	return (0.54 - 0.46*cos(2 * pie*n / (N - 1)))*rectwindow(n);
}
double blackmanwin(int n)
{
	return (0.42 - 0.5*cos(2 * pie*n / (N - 1)) + 0.08*cos(4 * pie*n / (N - 1)))*rectwindow(n);
}

unsigned int fact(unsigned int k)
{
	if (k == 0)
	{
		return 1;
	}
	else return k*fact(k - 1);
}

#define baseln 20 //贝塞尔函数取20项，读者可任意修改
double basel(double x)
{
	double sum = 0;
	for (int i = 1; i < baseln; i++)
	{
		sum = pow((pow((x / 2), i) / fact(i)), 2) + sum;
	}
	return sum + 1;
}

#define alpha 7.8
double I0 = basel(alpha);
double kaiserwin(int n)
{
	if (n<N&&n>0)
	{
		return rectwindow(n)*basel(alpha*(sqrt(1 - (2.0 * n / (N - 1) - 1)*(2.0 * n / (N - 1) - 1)))) / I0;
	}
	else return 0;
}





double rec(int x)
{
	if (abs(x) <= 5)
		return 1;
	else return 0;
}


double drect(int x)
{
	if (abs(x) <= 5)
		return 10;
	else return 0;
}
double gauss2(double x)
{
	return gauss(3 - 2 * x);
}
double crec(double x)
{
	if (abs(x) < 2)
		return 1;
	else return 0;
}

double cconv(double(*f)(double), double(*g)(double), double t, double flow, double fhigh, double glow, double ghigh, double accur)
//只在有界区间上计算，t是输入参数，f的下界是flow，上界是fhigh，g同理，accur代表黎曼积分的最小矩形区间
{
	double min, max; //积分区间
	double sum = 0; //返回值
	max = ghigh - t<fhigh ? ghigh - t : fhigh;
	min = glow - t>flow ? glow - t : flow;
	if (max>min)
	{
		for (double i = min; i < max; i = i + accur)
		{
			sum += f(i)*g(t - i)*accur;    //黎曼积分公式
		}
	}
	return sum;
}

double showva(double x)
{
	return cconv(crec, crec, x, -2, 2, -2, 2);
}

double lshi(int r, int b, double x)
{
	return 1 - pow((1 - pow(x, r)), b);
}
double lsh(double x)
{
	return lshi(5, 50, x);
}



#define min(a,b) (((a)<(b))?(a):(b))


#define max(a,b) (((a)>(b))?(a):(b))


void conv(double u[], double v[], double w[], int m, int n)
{
	int i, j;

	int k = m + n - 1;

	for (int l = 0; l < k; l++)
		w[l] = 0;
	for (i = 0; i < k; i++)
	{
		for (j = max(0, i + 1 - n); j <= min(i, m - 1); j++)
		{
			w[i] += u[j] * v[i - j];
		}
	}
}

void filter(int ord, double *a, double *b, int np, double *x, double *y)
{
	int i, j;
	y[0] = b[0] * x[0];
	for (i = 1; i<ord + 1; i++)
	{
		y[i] = 0.0;
		for (j = 0; j < i + 1; j++)
		{
			if (i - j >= 0)
			{
				y[i] = y[i] + b[j] * x[i - j];
			}
		}
		for (j = 0; j < i; j++)
		{
			if (i - j - 1 >= 0)
			{
				y[i] = y[i] - a[j + 1] * y[i - j - 1];
			}
		}

	}
	/* end of initial part */
	for (i = ord + 1; i<np + 1; i++)
	{
		y[i] = 0.0;
		for (j = 0; j<ord + 1; j++)
		{
			if (i - j >= 0)
			{
				y[i] = y[i] + b[j] * x[i - j];
			}
		}
		for (j = 0; j<ord; j++)
		{
			if (i - j - 1 >= 0)
			{
				y[i] = y[i] - a[j + 1] * y[i - j - 1];
			}
		}
	}
} /* end of filter */