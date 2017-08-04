

#include"part2.h"


//����ʽ�ĵ���
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
//����ʽ�ĵ���
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

//����
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
//m�ε���
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
�����ӷ�������ʽ
����ԭ����ʵ����R��ÿ������ʽҪô���Էֽ�Ϊһ�κͶ��ζ���ʽ�ĳ˻�
֤��������ά������������C��ÿ���ǳ���ȫ���������Ӷ���������ʽ�������������������������棨��������Զ�㴦��ȫ����
�����˵��һ������1/f(z)����z=����ʱ����ɢ����ô�ڸ�ƽ���ϱ���һ��z0ʹ�䷢ɢ���Ӷ�f(z)=0��
�Ӷ����Ƶ��������������� ��������ÿ������ʽ���н⡣
��R����һ����ʽf����C�ϣ����н�a����ô��Ϊf��ʵϵ������ʽ�����a����ʵ������ôa�Ķ��׹���Ҳ������һ���⣬�Ӷ�f=(x-a)(x-a_bar) g= (x+2re(a) x+|a|^2)-g
��֤��
�����ӷ������ҵ�һ�����ζ���ʽ�����ϱƽ�ʹ����ʽ�㹻С���Ӷ����������ζ���ʽ���ɵõ��ߴζ���ʽ�Ľ⡣
*/
comvec PolynomialEquation(poly A)
{
	int r = 0;

	int n = A.size() - 1;//���̵���ߴ���

	if (n < 1) return comvec();
	poly a, x;
	comvec X;
	a = poly(A);
	x = poly(2 * n);
	int m = (n + (n % 2)) / 2;//��Ҫ����Ĵ���
	long double u = 0.0;//U(x)��һ����
	long double v = 0.0;//U(X)�ĳ�����
	long double du = 0.0;
	long double dv = 0.0;
	int count = 0;

	poly b, cu, cv;
	b = poly(n + 1);//Q(x)��ϵ��
	cu = poly(n + 1);//ƫ����u
	cv = poly(n + 1);//ƫ����v


	long double r0, r1;
	long double pr0u, pr0v;
	long double pr1u, pr1v;

	long double det, p1, p2;
	for (int k = 0; k < m; k++)
	{
		u = 4.0;
		v = 3.0;//����̶���ȡֵ������һ���������Ҫ��
		if (n >= 3)//���μ��������ϵķ���
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
				if (count > 10000) break;//��������ʵ�����,u��vȡֵ�ʵ��Ļ��������ܿ�
			} while (fabs(du) > 1.0e-10 || fabs(dv) > 1.0e-10);
			//ȷ����һ�����������
			n -= 2;
			a = poly(n + 1);
			for (int i = n; i >= 0; i--) a[i] = b[i];
		}
		else
		{
			if (n == 2)//ʣ��һ������ʽ
			{
				if (fabs(a[2]) > 0.0)
				{
					u = a[1] / a[2];
					v = a[0] / a[2];
				}
				else
				{
					MessageBox(0, "1", "��ʾ", 0);//	_Debug_message//��Ӧ�ó����������
				}
			}
			else
				if (n == 1)//ʣ��һ��һ��ʽ
				{
					r += 11;
					if (fabs(a[1]) > 0.0)
					{
						x[4 * k] = -a[0] / a[1];
					}
					else
					{
						MessageBox(0, "2", "��ʾ", 0);//��Ӧ�ó����������
					}
					break;
				}
		}

		//�������һԪ���η���
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






//���������ʽ�ĵݹ麯�������У��У�չ����ʽ���㣬���׶���ע�⣬���ַ������Բ������ڴ��ģ����
double detm(matrix a1)
{
	if (a1.empty()) return 0;
	if (a1.size() != a1[0].size())
		return 0;
	int n1 = a1.size();
	int i, j, c;            //cΪb����
	matrix  b = matrix(n1 - 1, lfvector(n1 - 1));        //���У��У�չ��������ʽ    
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

//����detm����
complex detm(cmatrix a1)
{
	if (a1.empty()) return 0;
	if (a1.size() != a1[0].size())
		return 0;
	int n1 = a1.size();
	int i, j, c;            //cΪb����
	cmatrix  b = cmatrix(n1 - 1, vector<complex>(n1 - 1, complex()));        //���У��У�չ��������ʽ    
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


//��֤����ʽp����ߴ���ϵ����0
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

		f.resize(nl);             //����
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


/*�Ը�������*/

complex eliminant(cpoly f, cpoly g)
{
	int nl = f.size() - 1 + g.size() - 1;
	int m = f.size() - 1;            //m
	int n = g.size() - 1;         //n
	cmatrix e;
	//matrix e = matrix(nl, lfvector(nl));
	for (int i = 0; i < n; i++)
	{

		f.resize(nl);             //����
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









/*���ý�ʽ������ظ������������������ʽ����*/

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
����ʽ��ֵ����������ķ������ؾ����㷨������ֱ�Ӹ�ֵ��Ϊ����poly������ʸ������
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

/*��ʽ��������*/
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
/*��ʽ���������滻�ĺ��������÷�ʽ����s���fun�еı�������Ҫ����˫���Բ��䷨*/
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





/*�����ߴη��̵Ľ�͸�����������ֵΪvecom����vecom a=vector<vector<complex> >
a[i]��ʾi����������Ӧ��comple��ɵ�vector, a[0]�����塣
*/

vevecom solve(poly a)
{
	if (a.empty())
		return vevecom();
	unsigned int max = mulroot(a);       //�������
	vevecom com = vevecom(1, vector<complex>());                        //���ķ���ֵ
	comvec so = PolynomialEquation(a);           //���нⶼ����������
	poly der = a;
	for (int i = 1; i <= max; i++)
	{
		der = derivative(der);
		comvec temp;
		for (int j = 0; j < so.real.size(); j++)
		{
			complex comval = complex(so.real[j], so.image[j]);  //��j����ĸ������ʽ
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
	unsigned int max = mulroot(a);       //�������
	vevecom com = vevecom(1, vector<complex>());                        //���ķ���ֵ
	poly areal;
	for (int i = 0; i < a.size(); i++)
		areal.push_back(a[i].real);
	comvec so = PolynomialEquation(areal);           //���нⶼ����������
	cpoly der = a;
	for (int i = 1; i <= max; i++)
	{
		der = derivative(der);
		comvec temp;
		for (int j = 0; j < so.real.size(); j++)
		{
			complex comval = complex(so.real[j], so.image[j]);  //��j����ĸ������ʽ
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
����Ĭ����
���󷽳�Ax=C���ÿ���Ĭ������⣬����Ϊ�˷������ù�ʽ��matrix���к������ʽ��ʵ�����෴��
����matrix��vector<��>����ÿ������һ��vector<double>��
*/
//ע�⣬�����ǵ�ά�����ʹ�ÿ���Ĭ����һ���������Ӧ�ò����Ÿ��ȵ����Ƚ��Ʒ�����
lfvector cramer(matrix a, lfvector c)
{
	double det = detm(a);
	if (det == 0)
		return lfvector();   //vector��Ĭ�Ϲ��캯�������ؿ�vector��
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
		return vector<complex>();   //vector��Ĭ�Ϲ��캯�������ؿ�vector��
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


/*�׳˺���*/
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

/*��������Ӧ�����*/
vecfunc zeroinput(poly linearsystem, lfvector initial)
{
	if (linearsystem.size() - 1 != initial.size())
		return vecfunc();          //��ʼ״̬��ϵͳ����Ҫƥ��
		
	
								   //cmatrix ma = cmatrix(initial.size(), vector<complex>(initial.size(), complex(0)));
	cmatrix coef;                  //ת�������Է�������ϵ���洢��coef������
	vevecom s = solve(linearsystem);    //ͨ��solve�����������ֵ��solve�������������ӷ��������ֵ�����ܹ����ÿ������ֵ�Ľ�����
	int count = 0;
	for (int i = 1; i < s.size(); i++)       //iѭ������������н������ظ���i�ظ��洢��s[i]��
	{
		for (int j = 0; j < s[i].size(); j++)  //����ÿ����
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
					temp.push_back(value(derivative(f, k), complex(0)));      //��f��k�ε�������ֵ�õ�ϵ��
																			  //temp.push_back(comb(k, m) *postivepower(s[i][j], k - m));
				}
				coef.push_back(temp);   //���ݸ�����������ϵ��
			}


		}
	}



	//�����ʼ״̬���ÿ���Ĭ�������
	vector<complex> cinitial;
	for (int i = 0; i < initial.size(); i++)
		cinitial.push_back(complex(initial[i]));  //��ʼ״̬����
	vector<complex> solv = cramer(coef, cinitial);   //cramer�����������������δ֪��
	vevecom zeros;
	unsigned int count1 = 0;
	unsigned int count2 = 0;
	for (int i = 1; i < s.size(); i++)
	{
		for (int j = 0; j < s[i].size(); j++)
		{
			zeros.push_back(vector<complex>());
			zeros[count1].push_back(s[i][j]);      //���Ȱ�����ֵ������ǰ��
			for (int k = 0; k < i; k++)
			{

				zeros[count1].push_back(solv[count2]); //Ȼ�������ϵ��
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
/*�����Ӧ���*/
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




/*ʵ����ʽԼ�ֺ�����Ϊ��Ӧ������ʽ��������fracfuction����Ϊcpoly����������Ҫ�����ʽϵ������ʵ��*/
ve3com reduce(fracfunction f)
{
	//if (eliminant(f.numer, f.deno).cabs() != 0)
	//	return f;
	vevecom soln = solve(f.numer);
	vevecom sode = solve(f.deno);
	//	complex coef = f.numer[f.numer.size() - 1] / f.deno[f.deno.size() - 1]; //�����ϵ��
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
							goto out; //��goto��������ѭ��
						}
						if (i < j)
						{
							vector<complex>::iterator it = soln[i].begin() + ii;
							soln[i].erase(it);

							sode[j - i].push_back(sode[j][jj]);
							//	i--;
							it = sode[j].begin() + jj;
							sode[j].erase(it);
							goto out; //��goto��������ѭ��

						}
						if (i > j)
						{
							soln[i - j].push_back(sode[j][jj]);
							vector<complex>::iterator it = soln[i].begin() + ii;
							soln[i].erase(it);
							//	i--;
							it = sode[j].begin() + jj;
							sode[j].erase(it);
							goto out; //��goto��������ѭ��
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



/*ʵ����ʽ�Ĳ������ӷֽ�*/
vevecom decomp(fracfunction f)
{
	vevecom sing;// = solve(f.fra.deno);
	bool flag = 0;
	complex coef = f.numer[f.numer.size() - 1] / f.deno[f.deno.size() - 1]; //�����ϵ��
	cpoly num;
	ve3com fr = reduce(f);
	if (eliminant(f.numer, f.deno).cabs() == 0)  //ֻ��ⲻ��Լ����ʽ
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
		for (int ii = 0; ii < sing[i].size(); ii++)      //�ȱ���ÿ������
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
			laplace[count].push_back(sing[i][ii]);    //����ֵ��һ���洢����ֵ
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



/*���滯 Butterworth����ʽ*/
/*
cpoly Butterworth(int n)
{
	cpoly mul(1, complex(1));
	if (n % 2 == 0)  //ż�����
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
		cpoly p; //һ�ζ���ʽ
		p.push_back(complex(0));
		p.push_back(complex(1));
		double z = 0.5 + (2 * k + 1)*1.0 / (2 * n);
		p[0] = -1*exp(complex(0, pie)*z) ;
		mul = mul*p;
	}
	return mul;
}





/* sinh ˫�����Һ���*/
double sinh(double x)
{
	return (exp(x) - exp(-1 * x)) / 2;
}

/* cosh ˫�����Һ���*/
double cosh(double x)
{
	return (exp(x) + exp(-1 * x)) / 2;
}
/*��˫�����Һ���*/
double arsinh(double x)
{
	return log(x + sqrt(x*x + 1));
}

/*��˫�����Һ���*/
double arcosh(double x)
{
	return log(x + sqrt(x*x - 1));
}


/*Chebyshev ����ʽ*/
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
/*��ֵ�ú���*/
complex Chebyshev(unsigned int n, double fre)
{
	if (abs(fre) <= 1)
		return cos(n*acos(fre));
	else
		return cosh(n*arcosh(fre));
}

//��ɢʱ������ϵͳ��������Ӧ
vevecom dzeroinput(poly linearsystem, lfvector initial)
{
	if (linearsystem.size() - 1 != initial.size())
		return vevecom();          //��ʼ״̬��ϵͳ����Ҫƥ��


								   //cmatrix ma = cmatrix(initial.size(), vector<complex>(initial.size(), complex(0)));
	cmatrix coef;                  //ת�������Է�������ϵ���洢��coef������
	vevecom s = solve(linearsystem);    //ͨ��solve�����������ֵ��solve�������������ӷ��������ֵ�����ܹ����ÿ������ֵ�Ľ�����
	int count = 0;
	for (int i = 1; i < s.size(); i++)       //iѭ������������н������ظ���i�ظ��洢��s[i]��
	{
		for (int j = 0; j < s[i].size(); j++)  //����ÿ����
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
					temp.push_back(f.val(-k-1));      //��f��k�ε�������ֵ�õ�ϵ��
																			  //temp.push_back(comb(k, m) *postivepower(s[i][j], k - m));
				}
				coef.push_back(temp);   //���ݸ�����������ϵ��
			}


		}
	}



	//�����ʼ״̬���ÿ���Ĭ�������
	vector<complex> cinitial;
	for (int i = 0; i < initial.size(); i++)
		cinitial.push_back(complex(initial[i]));  //��ʼ״̬����
	vector<complex> solv = cramer(coef, cinitial);   //cramer�����������������δ֪��
	vevecom zeros;
	unsigned int count1 = 0;
	unsigned int count2 = 0;
	for (int i = 1; i < s.size(); i++)
	{
		for (int j = 0; j < s[i].size(); j++)
		{
			zeros.push_back(vector<complex>());
			zeros[count1].push_back(s[i][j]);      //���Ȱ�����ֵ������ǰ��
			for (int k = 0; k < i; k++)
			{

				zeros[count1].push_back(solv[count2]); //Ȼ�������ϵ��
				count2++;
			}
			count1++;
		}
	}


	
	return zeros;

}