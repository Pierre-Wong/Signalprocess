

#include"part2.h"


//多项式的导数
poly derivative(poly f)
{
	if (f.size() <= 1)
		return poly();
	poly diff;
	for (int i = 1; i < f.size(); i++)
	{
		diff.push_back(f[i] * i);
	}
	return diff;
}
//多项式的导数
cpoly derivative(cpoly f)
{
	if (f.size() <= 1)
		return cpoly();
	cpoly diff;
	for (int i = 1; i < f.size(); i++)
	{
		diff.push_back(f[i] * i);
	}
	return diff;
}

//导数
func derivative(expfunction  f)
{
	//expfunction diff;
	if (f.expo.cabs() != 0)
		return func(f, fracfunction(cpoly(1, f.expo), cpoly(1, 1)));
	else return func();
}
func derivative(fracfunction f)
{
	fracfunction fr;
	fr.numer = derivative(f.numer)*f.deno - derivative(f.deno)*f.numer;
	fr.deno = f.deno*f.deno;
	return fr;
}
vecfunc derivative(func f)
{
	if (f.fra.numer.empty())
		return vecfunc();
	func a = (derivative(f.fra) * f.ex);
	func b = (derivative(f.ex)* f.fra);
	return a + b;
}

vecfunc derivative(vecfunc f)
{
	vecfunc sum;
	for (int i = 0; i < f.size(); i++)
		sum = sum + derivative(f[i]);
	return sum;
}
//m次导数
vecfunc derivative(vecfunc f, unsigned int m)
{
	if (m == 0)
		return f;
	if (m == 1)
	{
		return derivative(f);
	}
	else return derivative(derivative(f, m - 1));
}

vecfunc derivative(func f, unsigned int m)
{
	if (m == 0)
		return vecfunc(1, f);
	if (m == 1)
	{
		return derivative(f);
	}
	else return derivative(derivative(f, m - 1));
}




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
comvec PolynomialEquation(poly A)
{
	int r = 0;

	int n = A.size() - 1;//方程的最高次数

	if (n < 1) return comvec();
	poly a, x;
	comvec X;
	a = poly(A);
	x = poly(2 * n);
	int m = (n + (n % 2)) / 2;//需要计算的次数
	long double u = 0.0;//U(x)的一次项
	long double v = 0.0;//U(X)的常数项
	long double du = 0.0;
	long double dv = 0.0;
	int count = 0;

	poly b, cu, cv;
	b = poly(n + 1);//Q(x)的系数
	cu = poly(n + 1);//偏导数u
	cv = poly(n + 1);//偏导数v


	long double r0, r1;
	long double pr0u, pr0v;
	long double pr1u, pr1v;

	long double det, p1, p2;
	for (int k = 0; k < m; k++)
	{
		u = 4.0;
		v = 3.0;//这里固定了取值，对于一般计算满足要求
		if (n >= 3)//三次及三次以上的方程
		{
			count = 0;
			b[n] = b[n - 1] = 0.0;
			cu[n] = cu[n - 1] = cu[n - 2] = 0.0;
			cv[n] = cv[n - 1] = cv[n - 2] = 0.0;
			do
			{
				for (int i = n - 2; i >= 0; i--)
				{
					b[i] = a[i + 2] - u * b[i + 1] - v * b[i + 2];
				}
				r1 = a[1] - u * b[0] - v * b[1];
				r0 = a[0] - v * b[0];
				for (int i = n - 3; i >= 0; i--)
				{
					cu[i] = -b[i + 1] - u * cu[i + 1] - v * cu[i + 2];
					cv[i] = -u * cv[i + 1] - b[i + 2] - v * cv[i + 2];
				}

				pr1u = -b[0] - u * cu[0] - v * cu[1];
				pr1v = -u * cv[0] - b[1] - v * cv[1];
				pr0u = -v * cu[0];
				pr0v = -b[0] - v * cv[0];

				while (pr1v * pr0u - pr1u * pr0v == 0)
				{
					u = u + 0.001;
					v = v + 0.001;
					pr1u = -b[0] - u * cu[0] - v * cu[1];
					pr1v = -u * cv[0] - b[1] - v * cv[1];
					pr0u = -v * cu[0];
					pr0v = -b[0] - v * cv[0];
					count = 0;
				}


				du = (r1 * pr0v - r0 * pr1v) / (pr1v * pr0u - pr1u * pr0v);
				dv = (r0 * pr1u - r1 * pr0u) / (pr1v * pr0u - pr1u * pr0v);

				u += du;
				v += dv;

				count++;
				if (count > 10000) break;//这里可以适当调整,u、v取值适当的话，收敛很快
			} while (fabs(du) > 1.0e-10 || fabs(dv) > 1.0e-10);
			//确定下一步计算的数据
			n -= 2;
			a = poly(n + 1);
			for (int i = n; i >= 0; i--) a[i] = b[i];
		}
		else
		{
			if (n == 2)//剩下一个二次式
			{
				if (fabs(a[2]) > 0.0)
				{
					u = a[1] / a[2];
					v = a[0] / a[2];
				}
				else
				{
					MessageBox(0, "1", "提示", 0);//	_Debug_message//不应该出现这种情况
				}
			}
			else
				if (n == 1)//剩下一个一次式
				{
					r += 11;
					if (fabs(a[1]) > 0.0)
					{
						x[4 * k] = -a[0] / a[1];
					}
					else
					{
						MessageBox(0, "2", "提示", 0);//不应该出现这种情况
					}
					break;
				}
		}

		//下面求解一元二次方程
		det = u * u - 4.0 * v;
		p1 = -0.5 * u;
		p2 = 0.5 * sqrt(fabs(det));
		if (det >= 0.0)
		{
			r += 22;
			x[4 * k] = p1 + (p1 > 0.0 ? 1.0 : p1 < 0.0 ? -1.0 : 0.0) * p2;
			if (fabs(x[4 * k]) > 0.0)
			{
				x[4 * k + 2] = v / x[4 * k];
			}
			else
			{
				x[4 * k + 2] = p1 - (p1 > 0.0 ? 1.0 : p1 < 0.0 ? -1.0 : 0.0) * p2;
			}
		}
		else
		{
			r += 20;
			x[4 * k] = p1;
			x[4 * k + 1] = p2;
			x[4 * k + 2] = p1;
			x[4 * k + 3] = -p2;
		}
	}

	//X = poly(x.size());
	for (int i = 0; i < x.size(); i = i + 2)
	{
		X.real.push_back(x[i]);
		X.image.push_back(x[i + 1]);
	}
	return X;

}






//求矩阵行列式的递归函数，按行（列）展开方式计算，简单易懂，注意，这种方法绝对不适用于大规模矩阵。
double detm(matrix a1)
{
	if (a1.empty()) return 0;
	if (a1.size() != a1[0].size())
		return 0;
	int n1 = a1.size();
	int i, j, c;            //c为b的行
	matrix  b = matrix(n1 - 1, lfvector(n1 - 1));        //按行（列）展开的余子式    
	int p = 0, q = 0;
	double sum = 0;
	if (n1 == 1) return a1[0][0];
	for (i = 0; i<n1; i++)
	{
		for (c = 0; c<n1 - 1; c++)
		{
			if (c<i)  p = 0;
			else  p = 1;
			for (j = 0; j<n1 - 1; j++)
			{
				b[c][j] = a1[c + p][j + 1];
			}

		}
		q = i % 2 ? -1 : 1;
		sum = sum + a1[i][0] * q*detm(b);
	}

	return sum;

}

//重载detm函数
complex detm(cmatrix a1)
{
	if (a1.empty()) return 0;
	if (a1.size() != a1[0].size())
		return 0;
	int n1 = a1.size();
	int i, j, c;            //c为b的行
	cmatrix  b = cmatrix(n1 - 1, vector<complex>(n1 - 1, complex()));        //按行（列）展开的余子式    
	int p = 0, q = 0;
	complex sum = complex();
	if (n1 == 1) return a1[0][0];
	for (i = 0; i<n1; i++)
	{
		for (c = 0; c<n1 - 1; c++)
		{
			if (c<i)  p = 0;
			else  p = 1;
			for (j = 0; j<n1 - 1; j++)
			{
				b[c][j] = a1[c + p][j + 1];
			}

		}
		q = i % 2 ? -1 : 1;
		sum = sum + a1[i][0] * q*detm(b);
	}

	return sum;

}


//保证多项式p的最高次项系数非0
void preparepoly(poly &p)
{


	vector<double>::reverse_iterator del;
	for (del = p.rbegin(); del != p.rend();)
	{

		if ((*del) <= epison)
		{

			del = vector<double>::reverse_iterator(p.erase((++del).base()));

		}
		else break;
	}
}





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
double eliminant(poly f, poly g)
{
	preparepoly(f);
	preparepoly(g);
	int nl = f.size() - 1 + g.size() - 1;
	int m = f.size() - 1;            //m
	int n = g.size() - 1;         //n
	matrix e;
	//matrix e = matrix(nl, lfvector(nl));
	for (int i = 0; i < n; i++)
	{

		f.resize(nl);             //补零
		e.push_back(f);
		f.insert(f.begin(), 0);
		f.resize(nl);
	}
	for (int i = 0; i < m; i++)
	{
		g.resize(nl);
		e.push_back(g);
		g.insert(g.begin(), 0);
		g.resize(nl);
	}
	return detm(e);
}


/*对复数重载*/

complex eliminant(cpoly f, cpoly g)
{
	int nl = f.size() - 1 + g.size() - 1;
	int m = f.size() - 1;            //m
	int n = g.size() - 1;         //n
	cmatrix e;
	//matrix e = matrix(nl, lfvector(nl));
	for (int i = 0; i < n; i++)
	{

		f.resize(nl);             //补零
		e.push_back(f);
		f.insert(f.begin(), complex(0));
		f.resize(nl);
	}
	for (int i = 0; i < m; i++)
	{
		g.resize(nl);
		e.push_back(g);
		g.insert(g.begin(), complex(0));
		g.resize(nl);
	}
	return detm(e);
}









/*利用结式求最高重根重数，不断求导再求结式即可*/

unsigned int mulroot(poly a)
{
	unsigned count = 1;
	if (a.size() - 1 == 1)
		return 1;
	poly diff = a;
	for (int i = 0; i < a.size() - 1; i++)
	{
		diff = derivative(diff);
		if (eliminant(a, diff) != 0)
			return count;
		count++;
	}
	return 0;
}




unsigned int mulroot(cpoly a)
{
	unsigned count = 1;
	if (a.size() - 1 == 1)
		return 1;
	cpoly diff = a;
	for (int i = 0; i < a.size() - 1; i++)
	{
		diff = derivative(diff);
		if (eliminant(a, diff).cabs() != 0)
			return count;
		count++;
	}
	return 0;
}






/*
多项式赋值函数，更快的方法是秦九昭算法，这里直接赋值是为了让poly类的性质更清楚。
*/
complex value(poly f, complex x)
{
	if (f.empty())
		return 0;
	complex sum = complex(0);
	for (unsigned int i = 0; i < f.size(); i++)
		sum = f[i] * power(x, i) + sum;
	return sum;
}

complex value(vecfunc f, complex x)
{
	if (f.empty())
		return 0;
	complex sum = complex(0);
	for (unsigned int i = 0; i < f.size(); i++)
		sum = f[i].val(x) + sum;
	return sum;
}

/*分式函数的幂*/
fracfunction power(fracfunction f, unsigned int n)
{
	if (n == 0)
	{
		fracfunction o;
		o.numer = cpoly(1, complex(1));
		o.deno = cpoly(1, complex(1));
		return o;
	}
	if (n == 1)
	{
		return f;
	}
	return f*power(f, n - 1);
}
/*分式函数变量替换的函数，即用分式函数s替代fun中的变量，主要用于双线性不变法*/
fracfunction replace(fracfunction fun, fracfunction s)
{
	fracfunction nu, de;
	nu.numer = cpoly(1, complex(0));
	nu.deno = cpoly(1, complex(1));
	de.numer = cpoly(1, complex(0));
	de.deno = cpoly(1, complex(1));
	for (int i = 0; i < fun.numer.size(); i++)
		nu = nu + power(s, i)*fun.numer[i];
	for (int i = 0; i < fun.deno.size(); i++)
		de = de + power(s, i)*fun.deno[i];
	return nu / de;
}





/*给出高次方程的解和根重数，返回值为vecom类型vecom a=vector<vector<complex> >
a[i]表示i重特征根对应的comple组成的vector, a[0]无意义。
*/

vevecom solve(poly a)
{
	if (a.empty())
		return vevecom();
	unsigned int max = mulroot(a);       //最高重数
	vevecom com = vevecom(1, vector<complex>());                        //最后的返回值
	comvec so = PolynomialEquation(a);           //所有解都保存在这里
	poly der = a;
	for (int i = 1; i <= max; i++)
	{
		der = derivative(der);
		comvec temp;
		for (int j = 0; j < so.real.size(); j++)
		{
			complex comval = complex(so.real[j], so.image[j]);  //第j个解的复数表达式
			if (value(der, comval).cabs() > 0.01)
			{
				temp.real.push_back(comval.real);
				temp.image.push_back(comval.image);
				lfvector::iterator it = so.real.begin() + j;
				so.real.erase(it);
				it = so.image.begin() + j;
				j--;
				so.image.erase(it);
			}
		}
		sort(temp.real.begin(), temp.real.end());
		sort(temp.image.begin(), temp.image.end());
		vector<complex> temproot;
		for (int k = 0; k < temp.real.size() / i; k++)
		{
			complex pushcom = complex(0);
			for (int m = 0; m < i; m++)
				pushcom = complex(temp.real[k*i + m], temp.image[k*i + m]) + pushcom;
			temproot.push_back(pushcom / (double)i);
			//	temproot.push_back(complex(temp.real[k*i ], temp.image[k*i]));
		}
		com.push_back(temproot);


	}
	return com;
}



vevecom solve(cpoly a)
{
	if (a.empty())
		return vevecom();
	unsigned int max = mulroot(a);       //最高重数
	vevecom com = vevecom(1, vector<complex>());                        //最后的返回值
	poly areal;
	for (int i = 0; i < a.size(); i++)
		areal.push_back(a[i].real);
	comvec so = PolynomialEquation(areal);           //所有解都保存在这里
	cpoly der = a;
	for (int i = 1; i <= max; i++)
	{
		der = derivative(der);
		comvec temp;
		for (int j = 0; j < so.real.size(); j++)
		{
			complex comval = complex(so.real[j], so.image[j]);  //第j个解的复数表达式
			if (value(der, comval).cabs() > 0.01)
			{
				temp.real.push_back(comval.real);
				temp.image.push_back(comval.image);
				lfvector::iterator it = so.real.begin() + j;
				so.real.erase(it);
				it = so.image.begin() + j;
				j--;
				so.image.erase(it);
			}
		}
		sort(temp.real.begin(), temp.real.end());
		sort(temp.image.begin(), temp.image.end());
		vector<complex> temproot;
		for (int k = 0; k < temp.real.size() / i; k++)
		{
			complex pushcom = complex(0);
			for (int m = 0; m < i; m++)
				pushcom = complex(temp.real[k*i + m], temp.image[k*i + m]) + pushcom;
			temproot.push_back(pushcom / (double)i);
			//	temproot.push_back(complex(temp.real[k*i ], temp.image[k*i]));
		}
		com.push_back(temproot);


	}
	return com;
}






/*
克拉默法则
矩阵方程Ax=C利用克拉默法则求解，这里为了方便运用公式，matrix的行和列与结式的实现是相反的
这里matrix是vector<列>，而每个列是一个vector<double>。
*/
//注意，这里是低维矩阵才使用克拉默法则，一般情况下则应该采用雅各比迭代等近似方法。
lfvector cramer(matrix a, lfvector c)
{
	double det = detm(a);
	if (det == 0)
		return lfvector();   //vector的默认构造函数，返回空vector。
	lfvector x;
	//x.reserve(a.size());
	for (int i = 0; i < a.size(); i++)
	{
		matrix b = a;
		b[i] = c;
		x.push_back(detm(b) / det);
	}
	return x;
}
vector<complex> cramer(cmatrix a, vector<complex> c)
{
	complex det = detm(a);
	if (det.cabs() == 0)
		return vector<complex>();   //vector的默认构造函数，返回空vector。
	vector<complex> x;
	//x.reserve(a.size());
	for (int i = 0; i < a.size(); i++)
	{
		cmatrix b = a;
		b[i] = c;
		x.push_back(detm(b) / det);
	}
	return x;
}


complex postivepower(complex x, int i)
{
	if (i < 0)
		return complex(0);
	else
		return power(x, i);
}


/*阶乘函数*/
unsigned gamma(int i)
{
	if (i < 0)
		return 0; //error
	if (i == 0)
		return 1;
	else return i*gamma(i - 1);

}
unsigned comb(int m, int n)
{
	if (m < 0 || n < 0 || m<n)
		return 0;
	return gamma(m) / gamma(n) / gamma(m - n);
}

/*零输入响应的求解*/
vecfunc zeroinput(poly linearsystem, lfvector initial)
{
	if (linearsystem.size() - 1 != initial.size())
		return vecfunc();          //初始状态和系统阶数要匹配
		
	
								   //cmatrix ma = cmatrix(initial.size(), vector<complex>(initial.size(), complex(0)));
	cmatrix coef;                  //转换成线性方程组后的系数存储在coef矩阵里
	vevecom s = solve(linearsystem);    //通过solve函数求解特征值，solve函数利用劈因子法求解特征值并且能够求出每个特征值的阶数。
	int count = 0;
	for (int i = 1; i < s.size(); i++)       //i循环代表遍历所有阶数的重根，i重根存储在s[i]中
	{
		for (int j = 0; j < s[i].size(); j++)  //遍历每个根
		{
			for (int m = 0; m <= i - 1; m++)
			{
				vector<complex> temp;
				for (int k = 0; k < initial.size(); k++)
				{
					func f;
					f.ex = expfunction(s[i][j]);
					f.fra = fracfunction();
					f.fra.deno = cpoly(1, complex(1));
					f.fra.numer = cpoly(m + 1, complex(0));
					f.fra.numer[m] = complex(1);                               //f=exp(lambda*t)*t^m 
					temp.push_back(value(derivative(f, k), complex(0)));      //对f求k次导数再求值得到系数
																			  //temp.push_back(comb(k, m) *postivepower(s[i][j], k - m));
				}
				coef.push_back(temp);   //根据根的重数带入系数
			}


		}
	}



	//带入初始状态，用克拉默法则求解
	vector<complex> cinitial;
	for (int i = 0; i < initial.size(); i++)
		cinitial.push_back(complex(initial[i]));  //初始状态向量
	vector<complex> solv = cramer(coef, cinitial);   //cramer函数，给定矩阵求解未知数
	vevecom zeros;
	unsigned int count1 = 0;
	unsigned int count2 = 0;
	for (int i = 1; i < s.size(); i++)
	{
		for (int j = 0; j < s[i].size(); j++)
		{
			zeros.push_back(vector<complex>());
			zeros[count1].push_back(s[i][j]);      //首先把特征值放在最前面
			for (int k = 0; k < i; k++)
			{

				zeros[count1].push_back(solv[count2]); //然后放入解的系数
				count2++;
			}
			count1++;
		}
	}
	
	
	vecfunc zerrepo;
	for (int i = 0; i < zeros.size(); i++)
	{
		func f;
		f.ex = zeros[i][0];
		for (int j = 1; j < zeros[i].size(); j++)
		{
			f.fra.numer.push_back(zeros[i][j]);
		}
		zerrepo.push_back(f);
	}
	
	return zerrepo;

}
/*冲击响应求解*/
vecfunc impulse(poly y, poly x)
{
	int n = y.size() - 1, m = x.size() - 1;
	if (y.size() <= x.size())
		return vecfunc();
	lfvector E = x;
	E.resize(n);
	//	lfvector initial=lfvector(y.size()-1);
	//	for (int i = 0; n - 1 - i >= 0 && i <= m;i++)
	//		initial[n - 1 - i] = x[]
	matrix coef;
	lfvector col = y;
	lfvector::iterator it = col.begin();
	col.erase(it);
	for (int i = 0; i < n; i++)
		coef.push_back(col);
	for (int i = 0; i < n; i++)
	{
		for (int k = 0; k < n - 1 - i; k++)
		{
			lfvector::iterator it = coef[i].begin();
			coef[i].erase(it);
			//coef[i].insert(coef[i].begin(), 0);
			coef[i].resize(n);
		}
	}
	lfvector initial = cramer(coef, E);
	reverse(initial.begin(), initial.end());
	vecfunc result = zeroinput(y, initial);
	return result;
}



cpoly ve2cpoly(vevecom sing)
{
	cpoly ss = cpoly(1, complex(1));
	for (int j = 1; j < sing.size(); j++)
	{
		for (int jj = 0; jj < sing[j].size(); jj++)
		{
			cpoly linear;

			linear.push_back(-1 * sing[j][jj]);
			linear.push_back(complex(1));
			ss = ss*linear;
		}
	}
	return ss;
}




/*实多项式约分函数，为了应用劈因式法，尽管fracfuction定义为cpoly，但是这里要求多项式系数都是实的*/
ve3com reduce(fracfunction f)
{
	//if (eliminant(f.numer, f.deno).cabs() != 0)
	//	return f;
	vevecom soln = solve(f.numer);
	vevecom sode = solve(f.deno);
	//	complex coef = f.numer[f.numer.size() - 1] / f.deno[f.deno.size() - 1]; //最高项系数
	for (int i = 1; i < soln.size(); i++)
	{
		for (int ii = 0; ii < soln[i].size(); ii++)
		{
			for (int j = 1; j < sode.size(); j++)
			{
				for (int jj = 0; jj < sode[j].size(); jj++)
				{
					if ((soln[i][ii] - sode[j][jj]).cabs() < 0.01)
					{
						if (i == j)
						{
							vector<complex>::iterator it = soln[i].begin() + ii;
							soln[i].erase(it);
							//	i--;
							it = sode[j].begin() + jj;
							sode[j].erase(it);
							goto out; //用goto跳出两层循环
						}
						if (i < j)
						{
							vector<complex>::iterator it = soln[i].begin() + ii;
							soln[i].erase(it);

							sode[j - i].push_back(sode[j][jj]);
							//	i--;
							it = sode[j].begin() + jj;
							sode[j].erase(it);
							goto out; //用goto跳出两层循环

						}
						if (i > j)
						{
							soln[i - j].push_back(sode[j][jj]);
							vector<complex>::iterator it = soln[i].begin() + ii;
							soln[i].erase(it);
							//	i--;
							it = sode[j].begin() + jj;
							sode[j].erase(it);
							goto out; //用goto跳出两层循环
						}

					}
				}
			}

		out:   _asm {};
		}
	}
	ve3com p(2, vevecom());
	p[0] = soln;
	p[1] = sode;
	return p;
	/*	cpoly f1 = cpoly(1, complex(1));
	for (int i = 1; i < soln.size(); i++)
	{
	for (int ii = 0; ii < soln[i].size(); ii++)
	{
	cpoly linear;
	linear.push_back(-1*soln[i][ii]);
	linear.push_back(complex(1));
	for (int k = 0; k < i; k++)
	f1 = f1*linear;
	}
	}
	cpoly f2 = cpoly(1, complex(1));
	for (int i = 1; i < sode.size(); i++)
	{
	for (int ii = 0; ii < sode[i].size(); ii++)
	{
	cpoly linear;
	linear.push_back(sode[i][ii]);
	linear.push_back(complex(1));
	for (int k = 0; k < i; k++)
	f2 = f2*linear;
	}
	}*/
	//return fracfunction(f1*coef, f2);
}



vecfunc vevec2vecf(vevecom zeros)
{
	vecfunc zerrepo;
	for (int i = 0; i < zeros.size(); i++)
	{
		func f;
		f.ex = zeros[i][0];
		for (int j = 1; j < zeros[i].size(); j++)
		{
			f.fra.numer.push_back(zeros[i][j]);
		}
		zerrepo.push_back(f);
	}

	return zerrepo;
}



/*实多项式的部分因子分解*/
vevecom decomp(fracfunction f)
{
	vevecom sing;// = solve(f.fra.deno);
	bool flag = 0;
	complex coef = f.numer[f.numer.size() - 1] / f.deno[f.deno.size() - 1]; //最高项系数
	cpoly num;
	ve3com fr = reduce(f);
	if (eliminant(f.numer, f.deno).cabs() == 0)  //只求解不可约多项式
	{
		flag = 1;

		//sing = reduce(f.fra)[1];
		//num = ve2cpoly(reduce(f.fra)[0])*coef;
		sing = fr[1];
		num = ve2cpoly(fr[0])*coef;
		//delete fr;
	}
	else
	{
		sing = solve(f.deno);
		num = f.numer *(1 / f.deno[f.deno.size() - 1]);
	}

	vevecom laplace;
	unsigned int count = 0;
	for (int i = 1; i < sing.size(); i++)
	{
		for (int ii = 0; ii < sing[i].size(); ii++)      //先遍历每个极点
		{
			cpoly ss = cpoly(1, complex(1));

			for (int j = 1; j < sing.size(); j++)
			{
				for (int jj = 0; jj < sing[j].size(); jj++)
				{
					if (!(sing[j][jj] == sing[i][ii]))
					{
						for (int m = 0; m < j; m++)
						{
							cpoly linear;
							linear.push_back(-1 * sing[j][jj]);
							linear.push_back(complex(1));
							ss = ss*linear;
						}
					}


				}
			}

			laplace.push_back(vector<complex>());
			laplace[count].push_back(sing[i][ii]);    //返回值第一个存储特征值
			for (int k = 1; k <= i; k++)
			{
				fracfunction ff = fracfunction(num, ss);
				func fc = func(ff);
				vecfunc mder = derivative(fc, k - 1);
				laplace[count].push_back(value(mder, sing[i][ii]) / gamma(k - 1));
			}

			count++;
		}

	}
	return laplace;


	//vecfunc lapfunc = vevec2vecf(laplace);
	//return lapfunc;


}



/*正规化 Butterworth多项式*/
/*
cpoly Butterworth(int n)
{
	cpoly mul(1, complex(1));
	if (n % 2 == 0)  //偶数情况
	{
		cpoly two;
		two.push_back(complex(1));
		two.push_back(complex(0));
		two.push_back(complex(1));
		for (int k = 1; k <= n / 2; k++)
		{
			two[1] = -2 * cos((2 * k + n - 1)*pie / 2 / n);
			mul = mul*two;
		}
		return mul;
	}
	else {
		cpoly one;
		one.push_back(complex(1));
		one.push_back(complex(1));
		return one*Butterworth(n - 1);
	}
}
*/
cpoly Butterworth(int n)
{
	cpoly mul;
	mul.push_back(complex(1));
	for (int k = 0; k < n; k++)
	{
		cpoly p; //一次多项式
		p.push_back(complex(0));
		p.push_back(complex(1));
		double z = 0.5 + (2 * k + 1)*1.0 / (2 * n);
		p[0] = -1*exp(complex(0, pie)*z) ;
		mul = mul*p;
	}
	return mul;
}





/* sinh 双曲正弦函数*/
double sinh(double x)
{
	return (exp(x) - exp(-1 * x)) / 2;
}

/* cosh 双曲余弦函数*/
double cosh(double x)
{
	return (exp(x) + exp(-1 * x)) / 2;
}
/*反双曲正弦函数*/
double arsinh(double x)
{
	return log(x + sqrt(x*x + 1));
}

/*反双曲余弦函数*/
double arcosh(double x)
{
	return log(x + sqrt(x*x - 1));
}


/*Chebyshev 多项式*/
cpoly Chebyshev(unsigned int n)
{

	if (n == 0)
	{
		return cpoly(1, complex(1));
	}
	if (n == 1)
	{
		cpoly one;
		one.push_back(complex(0));
		one.push_back(complex(1));
		return one;
	}
	cpoly tw;
	tw.push_back(complex(0));
	tw.push_back(complex(2));
	return (tw*Chebyshev(n - 1) - Chebyshev(n - 2));
}
/*求值用函数*/
complex Chebyshev(unsigned int n, double fre)
{
	if (abs(fre) <= 1)
		return cos(n*acos(fre));
	else
		return cosh(n*arcosh(fre));
}

//离散时间线性系统零输入响应
vevecom dzeroinput(poly linearsystem, lfvector initial)
{
	if (linearsystem.size() - 1 != initial.size())
		return vevecom();          //初始状态和系统阶数要匹配


								   //cmatrix ma = cmatrix(initial.size(), vector<complex>(initial.size(), complex(0)));
	cmatrix coef;                  //转换成线性方程组后的系数存储在coef矩阵里
	vevecom s = solve(linearsystem);    //通过solve函数求解特征值，solve函数利用劈因子法求解特征值并且能够求出每个特征值的阶数。
	int count = 0;
	for (int i = 1; i < s.size(); i++)       //i循环代表遍历所有阶数的重根，i重根存储在s[i]中
	{
		for (int j = 0; j < s[i].size(); j++)  //遍历每个根
		{
			for (int m = 0; m <= i - 1; m++)
			{
				vector<complex> temp;
				for (int k = 0; k < initial.size(); k++)
				{
					dfunc f;
					f.ex = powfunction(s[i][j]);
					f.fra = fracfunction();
					f.fra.deno = cpoly(1, complex(1));
					f.fra.numer = cpoly(m + 1, complex(0));
					f.fra.numer[m] = complex(1);                               //f=exp(lambda*t)*t^m 
					temp.push_back(f.val(-k-1));      //对f求k次导数再求值得到系数
																			  //temp.push_back(comb(k, m) *postivepower(s[i][j], k - m));
				}
				coef.push_back(temp);   //根据根的重数带入系数
			}


		}
	}



	//带入初始状态，用克拉默法则求解
	vector<complex> cinitial;
	for (int i = 0; i < initial.size(); i++)
		cinitial.push_back(complex(initial[i]));  //初始状态向量
	vector<complex> solv = cramer(coef, cinitial);   //cramer函数，给定矩阵求解未知数
	vevecom zeros;
	unsigned int count1 = 0;
	unsigned int count2 = 0;
	for (int i = 1; i < s.size(); i++)
	{
		for (int j = 0; j < s[i].size(); j++)
		{
			zeros.push_back(vector<complex>());
			zeros[count1].push_back(s[i][j]);      //首先把特征值放在最前面
			for (int k = 0; k < i; k++)
			{

				zeros[count1].push_back(solv[count2]); //然后放入解的系数
				count2++;
			}
			count1++;
		}
	}


	
	return zeros;

}