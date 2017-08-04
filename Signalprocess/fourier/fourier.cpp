#include"fourier.h"









decom DFT(decom x)             //离散傅里叶变换
{
	decom sum;
	for (int j = 0; j < x.size(); j++)
	{
		complex dt = complex(0);
		for (int i = 0; i < x.size(); i++)
		{
			dt = x[i] * exp(complex(0, -1 * 2 * pie*i*j / x.size())) + dt;
		}
		sum.push_back(dt);

	}
	return sum;
}
vector<unsigned int> ibj(int number) //n阶倒序输出
{
	int m = 1 << number;
	vector<unsigned int> index;
	for (int i = 0; i < m; i++)
	{
		int sum = 0, count = 0;
		for (int k = number - 1; k >= 0; k--)
		{
			sum += ((i >> k) % 2)*(1 << count);
			count++;
		}
		index.push_back(sum);
	}
	return index;

}


decom fft(decom x, unsigned int n)
{
	x.resize(1 << n);
	decom reindex;
	vector<unsigned int> index = ibj(n);

	for (int i = 0; i < x.size(); i++)
	{
		reindex.push_back(x[index[i]]);
	}

	/*
	for (int i = 0; i < n; i++)   //n次蝶型算法
	{
	decom w;
	for (int k = 0; k < x.size(); k++)
	{
	w.push_back(exp(complex(0, -1 *k* 2 * pie /x.size())));
	}
	vector<bool> flag = vector<bool>(x.size(), 0);//用于标记蝶形运算的上下节点
	int count = 0;
	for (int j = 0; j < x.size(); j++)
	{


	if (flag[(j)% x.size()]==0)
	{
	complex temp = reindex[j];
	reindex[j] = temp + reindex[(j + (1 << i)) % x.size()] * w[(j << (n - i - 1)) % x.size()];
	reindex[(j + (1 << i)) % x.size()] = temp - reindex[j + (1 << i) % x.size()] * w[(j << (n - i - 1)) % x.size()];

	flag[(j + (1 << i))%x.size()] = 1;
	}
	}
	}
	*/
	for (int l = 0; l < n; l++)
	{
		int d = 2 << l;
		complex u = complex(1);
		complex wn = exp(complex(0, -1 * 2 * pie / d));
		for (int j = 0; j < d / 2; j++)
		{
			for (int k = j; k < x.size(); k = k + d)
			{
				int kp = k + d / 2;
				complex t = reindex[kp] * u;
				reindex[kp] = reindex[k] - t;
				reindex[k] = reindex[k] + t;
			}
			u = u*wn;
		}

	}


	return reindex;
}


decom ifft(decom x, unsigned int n)
{
	for (int i = 0; i < x.size(); i++)
	{
		x[i] = con(x[i]) / x.size();
	}
	decom y = fft(x, n);
	for (int i = 0; i < y.size(); i++)
	{
		y[i] = con(y[i]);
	}
	return y;
}

vector<vector<complex> > fft(CImg<double> image)
{
	vector<vector<complex> > fimg;
	vector<vector<complex> > fimg2;
	int n = ceil(log2(image._width));
	for (int i = 0; i < image._height; i++)
	{
		vector<complex> col;

		for (int j = 0; j < image._width; j++)
			col.push_back(complex(image.atXY(j, i)));

		vector<complex> ff = fft(col, n);
		fimg.push_back(ff);

	}

	for (int i = 0; i < image._width; i++)
	{
		vector<complex> col;
		for (int j = 0; j < image._height; j++)
			col.push_back(fimg[j][i]);
		vector<complex> ff = fft(col, n);
		fimg2.push_back(ff);
	}


	return fimg2;


}


CImg<double> veve2im(vector<vector<complex> > veve)
{
	CImg<double>image(veve.size(), veve[0].size(), 1, 1, 0);
	image._is_shared = 1;
	image.fill(255);
	double max = 255;
	for (int i = 0; i < veve[0].size(); i++)
	{
		for (int j = 0; j < veve.size(); j++)
		{
			//	image.atXY(j, i) = veve[j][i].cabs();
			double m = veve[j][i].cabs();
			if (m>max)
				max = m;


		}
	}

	for (int i = 0; i < veve[0].size(); i++)
	{
		for (int k = 0; k < veve.size(); k++)
			image.atXY(k, i) = veve[k][i].cabs() / max * 255;

	}
	/*double *p1 = image.data(0, 0, 0, 0);
	for (int i = 0; i < veve[0].size(); i++)
	{
	for (int j = 0; j < veve.size(); j++)
	{
	double c = veve[j][i].cabs() / max * 255;

	*(p1++) = c;


	}
	}*/
	return image;
}

CImg<double> ifft(vector<vector<complex> > ff)
{
	vector<vector<complex> > ma;
	CImg<double> image(ff.size(), ff[0].size(), 1, 1, 0);
	unsigned int n = log2(ff[0].size());
	for (int i = 0; i < ff.size(); i++)
	{
		vector<complex>temp = ifft(ff[i], n);
		ma.push_back(temp);
	}

	for (int i = 0; i < ff[0].size(); i++)
	{
		vector<complex> col;

		for (int j = 0; j < ff.size(); j++)
			col.push_back(ma[j][i]);

		vector<complex> ff = ifft(col, n);
		for (int k = 0; k < ff.size(); k++)
			image.atXY(k, i) = ff[k].cabs();

	}
	return image;

}