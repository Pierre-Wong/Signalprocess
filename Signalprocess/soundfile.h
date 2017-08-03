#pragma once

using namespace std;
#include <windows.h>  
#include "draw.h"
#include <vector>

#pragma comment(lib, "winmm.lib")  

#define BUFFER_SIZE (8000*8*1/8*2)    // 录制声音长度  
#define FRAGMENT_SIZE 1024              // 缓存区大小  
#define FRAGMENT_NUM 4                  // 缓存区个数  

unsigned char buffer[BUFFER_SIZE] = { 0 };
static int buf_count = 0;



// 录音、播音回调函数声明 
void CALLBACK waveInProc(HWAVEIN hwi,
	UINT uMsg,
	DWORD_PTR dwInstance,
	DWORD_PTR dwParam1,
	DWORD_PTR dwParam2);
void CALLBACK waveOutProc(HWAVEOUT hwo,
	UINT uMsg,
	DWORD_PTR dwInstance,
	DWORD_PTR dwParam1,
	DWORD_PTR dwParam2);


int sound(unsigned char *ptr,//要显示的函数
	int xlow,//x坐标初始位置
	int xhigh,//x坐标结束位置
	int ylow,//y坐标初始位置
	int yhigh,//y坐标结束位置
	int dotted = 1 //是否画框虚线,默认为1，代表是,0则不画。
)//连续函数作图函数，五个参数分别为，要显示的函数指针，x初始坐标，x结束坐标，y初始坐标，y结束坐标。
{
	char wndname[] = "信号与系统";
	char str[10];

	int  width = 600, height = 600;


	CImg<unsigned char> image2((xhigh - xlow), height, 1, 3, 0);
	image2._is_shared = 1;
	image2.fill(255);
	//imshow(wndname, image);




	Point pt1, pt2, pt3, pt4;
	pt1.x = 40;
	pt1.y = 600 - 40;          //左下端点  注意opencv当中左上角的坐标是(0,0)自左向右x值增加，自上到下y值增加。
	pt2.x = 600 - 40;
	pt2.y = 600 - 40;		 //右下端点
	pt3.x = 40;
	pt3.y = 40;                  //左上端点
	pt4.x = 600 - 40;
	pt4.y = 40;                  //右下端点


								 //根据函数画出函数图像
	for (int i = 0; i < 520.0 * (xhigh - xlow) / 600 - 2; i = i + 2)
	{
		Point pta, ptb, ptc;
		pta.x = 40 * (xhigh - xlow) / 600 + i;
		ptb.x = 40 * (xhigh - xlow) / 600 + i + 1;
		pta.y = pt1.y - ptr[i + xlow] * 480.0 / (yhigh - ylow) + ylow * 480.0 / (yhigh - ylow);
		ptb.y = pt1.y - ptr[(i + 1) + xlow] * 480.0 / (yhigh - ylow) + ylow * 480.0 / (yhigh - ylow);
		ptc.x = 40 * (xhigh - xlow) / 600 + i + 2;
		ptc.y = pt1.y - ptr[(i + 2) + xlow] * 480.0 / (yhigh - ylow) + ylow * 480.0 / (yhigh - ylow);
		line(image2, pta, ptb, Scalar(255, 0, 0));
		line(image2, ptb, ptc, Scalar(255, 0, 0));
		line(image2, Point(pta.x - 1, pta.y), Point(ptb.x - 1, ptb.y), Scalar(255, 0, 0));
		line(image2, Point(ptb.x - 1, ptb.y), Point(ptc.x - 1, ptc.y), Scalar(255, 0, 0));
	}
	CImg<unsigned char>image = image2.get_resize(600, 600, -100, -100, 2, -1, 1);
	image._is_shared = 1;


	//根据xlow和xhigh的值的正负设置Y轴位置
	Point orial, orial0;
	int texflag = 1;
	orial.y = height - 40;
	orial0.y = 40;
	if (xlow >= 0)
	{
		orial.x = pt1.x;
		orial0.x = orial.x;
		line(image, orial, orial0, Scalar(0, 0, 0));
	}
	else if (xhigh <= 0)
	{
		orial.x = pt1.x + 480;
		orial0.x = orial.x;
		line(image, orial, orial0, Scalar(0, 0, 0));
		texflag = 0;
	}
	else
	{
		orial.x = (double)(0 - xlow) / (xhigh - xlow) * 480 + pt1.x;
		orial0.x = orial.x;
		line(image, orial, orial0, Scalar(0, 0, 0));
	}

	//根据ylow和yhigh的值的正负设置X轴位置
	Point orialy, orialy0;
	orialy.x = 40;
	orialy0.x = 560;
	if (ylow >= 0)
	{
		orialy.y = pt1.y;
		orialy0.y = orialy.y;
		line(image, orialy, orialy0, Scalar(0, 0, 0));
	}
	else if (yhigh <= 0)
	{
		orialy.y = pt3.y;
		orialy0.y = orial.y;
		line(image, orialy, orialy0, Scalar(0, 0, 0));
	}
	else
	{
		orialy.y = pt1.y - (double)(0 - ylow) / (yhigh - ylow) * 480;
		orialy0.y = orialy.y;
		line(image, orialy, orialy0, Scalar(0, 0, 0));
	}


	int count = 0;
	for (int i = pt1.x; i < pt2.x; i = i + 40)       //画横轴坐标循环
	{
		if (dotted&&i != pt1.x)                      //纵虚线
		{
			for (int j = pt1.y; j - 10 >= pt3.y; j = j - 10)
			{
				Point pta, ptb;
				pta.x = i;
				ptb.x = i;
				pta.y = j;
				ptb.y = j - 5;
				line(image, pta, ptb, Scalar(0, 0, 0), 1);
			}
		}
		Point tex, axis1, axis2;
		tex.x = i - 10;
		tex.y = orialy.y + 20;
		axis1.x = i;
		axis1.y = orialy.y;
		axis2.x = i;
		axis2.y = orialy.y - 5;
		circle(image, axis1, 2, Scalar(0, 0, 0));
		line(image, axis1, axis2, Scalar(0, 0, 0));
		double a = (double)(xhigh - xlow) / 12.0 *count + xlow;
		if (abs(xhigh - xlow)>1)
			sprintf(str, "%.1lf", a);
		else
			sprintf(str, "%.3lf", a);
		putText(image, str, tex, 1,
			0.8, Scalar(0, 0, 0));
		count++;
	}



	count = 0;
	for (int i = pt1.y; i > pt3.y; i = i - 40)//画纵轴坐标循环
	{
		if (dotted&&i != pt1.y)                   //横虚线
		{
			for (int j = pt1.x; j + 10 <= pt2.x; j = j + 10)
			{
				Point pta, ptb;
				pta.y = i;
				ptb.y = i;
				pta.x = j;
				ptb.x = j + 5;
				line(image, pta, ptb, Scalar(0, 0, 0), 1);
			}
		}
		Point tex;
		tex.y = i - 6;
		if (texflag)
			tex.x = orial.x - 30;
		else tex.x = orial.x;

		Point  axis1, axis2;
		axis1.y = i;
		axis1.x = orial.x;
		axis2.y = i;
		axis2.x = orial.x + 5;
		line(image, axis1, axis2, Scalar(0, 0, 0));
		circle(image, axis1, 2, Scalar(0, 0, 0));
		double a = (double)(yhigh - ylow) / 12.0 *count + ylow;
		if (abs(yhigh - ylow)>1)
			sprintf(str, "%.1lf", a);
		else
			sprintf(str, "%.3lf", a);
		putText(image, str, tex, 1,
			0.8, Scalar(0, 0, 0));
		count++;
	}





	image.display(wndname);

	return 0;
}


int nReturn = 0;

void record()
{
	// Device   
	nReturn = waveInGetNumDevs();
	printf("输入设备数目：%d\n", nReturn);
	for (int i = 0; i<nReturn; i++)
	{
		WAVEINCAPS wic;
		waveInGetDevCaps(i, &wic, sizeof(WAVEINCAPS));
		printf("#%d\t设备名：%s\n", i, wic.szPname);
	}

	// open   
	HWAVEIN hWaveIn;
	WAVEFORMATEX wavform;
	wavform.wFormatTag = WAVE_FORMAT_PCM;
	wavform.nChannels = 1;
	wavform.nSamplesPerSec = 8000;
	wavform.nAvgBytesPerSec = 8000 * 8 * 1 / 8;
	wavform.nBlockAlign = 1;
	wavform.wBitsPerSample = 8;
	wavform.cbSize = 0;

	waveInOpen(&hWaveIn, WAVE_MAPPER, &wavform, (DWORD_PTR)waveInProc, 0, CALLBACK_FUNCTION);

	WAVEINCAPS wic;
	waveInGetDevCaps((UINT_PTR)hWaveIn, &wic, sizeof(WAVEINCAPS));
	printf("打开的输入设备：%s\n", wic.szPname);

	// prepare buffer   
	static WAVEHDR wh[FRAGMENT_NUM];
	for (int i = 0; i<FRAGMENT_NUM; i++)
	{
		wh[i].lpData = new char[FRAGMENT_SIZE];
		wh[i].dwBufferLength = FRAGMENT_SIZE;
		wh[i].dwBytesRecorded = 0;
		wh[i].dwUser = NULL;
		wh[i].dwFlags = 0;
		wh[i].dwLoops = 1;
		wh[i].lpNext = NULL;
		wh[i].reserved = 0;

		waveInPrepareHeader(hWaveIn, &wh[i], sizeof(WAVEHDR));
		waveInAddBuffer(hWaveIn, &wh[i], sizeof(WAVEHDR));
	}

	// record   
	printf("Start to Record...\n");

	buf_count = 0;
	waveInStart(hWaveIn);

	while (buf_count < BUFFER_SIZE)
	{
		Sleep(1);
	}

	printf("Record Over!\n\n");

	// clean   
	waveInStop(hWaveIn);
	waveInReset(hWaveIn);
	for (int i = 0; i<FRAGMENT_NUM; i++)
	{
		waveInUnprepareHeader(hWaveIn, &wh[i], sizeof(WAVEHDR));
		delete wh[i].lpData;
	}
	waveInClose(hWaveIn);
	//unsigned char buffer2[8000];
	//for (int j = 0; j < 8000; j++)
	//{
	//	buffer2[j] = buffer[2 * j];
	//}

}
unsigned char* buffer2;
unsigned int start;
unsigned int endd;

void play()
{

	HWAVEIN hWaveIn = 0;
	WAVEFORMATEX wavform;
	wavform.wFormatTag = WAVE_FORMAT_PCM;
	wavform.nChannels = 1;
	wavform.nSamplesPerSec = 8000;
	wavform.nAvgBytesPerSec = 8000 * 8 * 1 / 8;
	wavform.nBlockAlign = 1;
	wavform.wBitsPerSample = 8;
	wavform.cbSize = 0;
	WAVEINCAPS wic;
	waveInGetDevCaps((UINT_PTR)hWaveIn, &wic, sizeof(WAVEINCAPS));



	//    放音
	nReturn = waveOutGetNumDevs();
	printf("\n输出设备数目：%d\n", nReturn);
	for (int i = 0; i<nReturn; i++)
	{
		WAVEOUTCAPS woc;
		waveOutGetDevCaps(i, &woc, sizeof(WAVEOUTCAPS));
		printf("#%d\t设备名：%s\n", i, wic.szPname);
	}

	// open   
	HWAVEOUT hWaveOut;
	waveOutOpen(&hWaveOut, WAVE_MAPPER, &wavform, (DWORD_PTR)waveOutProc, 0, CALLBACK_FUNCTION);

	WAVEOUTCAPS woc;
	waveOutGetDevCaps((UINT_PTR)hWaveOut, &woc, sizeof(WAVEOUTCAPS));
	printf("打开的输出设备：%s\n", wic.szPname);

	// prepare buffer   
	WAVEHDR wavhdr;
	wavhdr.lpData = (LPSTR)buffer2;
	wavhdr.dwBufferLength = (endd - start) * 80;
	wavhdr.dwFlags = 0;
	wavhdr.dwLoops = 0;

	waveOutPrepareHeader(hWaveOut, &wavhdr, sizeof(WAVEHDR));

	// play   
	printf("Start to Play...\n");

	buf_count = 0;
	waveOutWrite(hWaveOut, &wavhdr, sizeof(WAVEHDR));
	while (buf_count < (endd - start) * 80)
	{
		Sleep(1);
	}

	// clean   
	waveOutReset(hWaveOut);
	waveOutUnprepareHeader(hWaveOut, &wavhdr, sizeof(WAVEHDR));
	waveOutClose(hWaveOut);

	printf("Play Over!\n\n");

	//sound(buffer, 0, 6000, 0, 255);

}

vector<double> fra(unsigned char *buffera, unsigned int n)
{
	vector<double> su;
	for (int i = 0; i + 80 <= n; i = i + 80)     //10ms
	{
		double sum = 0;
		for (int j = 0; j < 80; j++)
			sum += ((double)buffera[i + j] - 126)*((double)buffera[i + j] - 126) / 80;
		su.push_back(sum);
	}
	return su;
}


vector<double> ze(unsigned char *buffera, unsigned int n)
{
	vector<double> su;
	for (int i = 0; i + 80 <= n; i = i + 80)     //10ms
	{
		double sum = 0;
		for (int j = 1; j < 80; j++)
			sum += abs((double)(buffera[i + j]> 126) - (double)(buffera[i + j - 1] >126));
		su.push_back(sum);
	}
	return su;
}




vector < double> xiaobo()
{
	vector<double> xiao;
	double max = 0;
	for (int i = 0; i < (endd - start) * 80; i++)
	{
		if (abs(buffer2[i] - 126)>max)
			max = abs(buffer2[i] - 126);
	}
	for (int i = 0; i < (endd - start) * 80; i++)
	{
		double x = buffer2[i] - 126;
		if (x>0.7*max)
			x = x - 0.7*max;
		else if (x < -0.7*max)
			x = x + 0.7*max;
		else x = 0;
		xiao.push_back(x);
	}
	return xiao;
}
vector<double> Rself(vector<double> b)
{
	vector<double> A;
	int ii = 0;
	while (b[ii] == 0)
	{
		ii++;
	}
	for (int k = 0; k < 320; k++)
	{
		double sum = 0;
		for (int m = ii; m < 320 + ii; m++)
		{
			if (m + k - ii >= 320)
				sum += 0;
			else
				sum += b[m] * b[m + k];
		}
		A.push_back(sum);
	}
	double m = A[0];
	for (int i = 0; i < A.size(); i++)
	{

		A[i] = A[i] / m;
	}
	return A;
}

vector<double> xx;

void setxx(vector<double> A)
{
	xx = A;
}

double showfun(int x)
{
	if (x >= 0 && x < xx.size())
		return xx[x];
	else return 0;
}
double rsf(int x)
{
	vector<double> b = Rself(xiaobo());
	if (x >= 0 && x < b.size())
		return b[x];
	else return 0;
}


double uniform(double a, double b, long int *seed)
{
	double t;
	*seed = 2045 * (*seed) + 1;
	*seed = *seed - (*seed / 1048576) * 1048576;
	t = *seed / 1048576.0;
	t = a + (b - a)*t;
	return t;
}

double boxmuller(double mean, double var)
{
	double u1 = (double)rand() / (RAND_MAX + 1);
	double u2 = (double)rand() / (RAND_MAX + 1);
	double z = sqrt(-2 * log(u2))*cos(2 * pie*u1);
	return mean + (z*sqrt(var));
}
void whitenoise(unsigned M, double *gauss, double SigmaNeed = 1, long d = 0)
{
	for (int i = 0; i < M; i++)
		gauss[i] = boxmuller(0, SigmaNeed);
}

/*void whitenoise(unsigned M, double *gauss, double SigmaNeed = 1,long d = 0)
{
//利用中心极限定理产生高斯白噪声
//#define M 1000
//#define N_perpoint 50
//#define MeanNeed 1
//#define SigmaNeed sqrt(2)
//该程序产生M个均值为MeanNeed，方差为SigmaNeed平方的高斯随机数
//每产生M个高斯点中的一个点需要N_perpoint个均匀分布的随机数。N_perpoint越大越精确
int N_perpoint = 50, MeanNeed = 0;
int         n, i;
double       y[50], mean = 0, sigma = 0;
double a, b;
long int s;

a = 0.0;
b = 1.0;
s = 13576 + d;
for (i = 0; i<M; i++)
{
gauss[i] = 0;
s = s + 1463;//修改每次的种子，使产生不同的变量
for (n = 0; n<N_perpoint; n++)
{
y[n] = uniform(a, b, &s);//产生均匀分布的随机变量
gauss[i] = gauss[i] + (double)sqrt((double)12 / N_perpoint)*y[n];
}
gauss[i] = gauss[i] - (double)sqrt((double)12 / N_perpoint)*(N_perpoint / 2);
gauss[i] = (double)(MeanNeed + sqrt(SigmaNeed)*gauss[i]);
mean = mean + gauss[i] / M;
}
for (i = 0; i<M; i++)
sigma = sigma + (gauss[i] - mean)*(gauss[i] - mean) / M;

}
*/


#define lpclen  10 //10 lpc

typedef vector < vector<double> >  fmartrix;
void lpc(unsigned char* buff, unsigned int length)
{
	//double pre[2] = { 1,-0.95 };


	double *one = new double[11];
	one[0] = 1;
	for (int l = 1; l < 11; l++)
		one[l] = 0;
	unsigned int Time = 80;
	for (int i = 0; i + 160 < length; i = i + 160)
	{
		double *temp2 = new double[160];
		for (int j = 0; j < 160; j++)
		{
			temp2[j] = (int)buff[i + j] - 126;
		}
		//		conv(pre, temp, temp2, 2, 160);

		vector<double> R;
		for (int jj = 0; jj < lpclen + 1; jj++)
		{
			double sum = 0;
			int total = 160 - jj;
			for (int j = jj; j < 160; j++)
			{
				sum += temp2[j] * temp2[j - jj] / total;
			}
			R.push_back(sum);
		}
		double E0 = R[0];
		vector<double> Ep(lpclen, 0);
		vector<double> k(lpclen, 0);;
		fmartrix a(lpclen, vector<double>(lpclen, 0));

		/*	for (int x = 0; x < lpclen; x++)
		{
		a.push_back(vector<float>());
		for (int y = 0; y < lpclen; y++)
		a[x].push_back(0);
		}*/


		double* lpc = new double[lpclen + 1];
		vector<double> E_V(lpclen + 1, 0);
		vector<double> a_yuce(lpclen + 1, 0);
		vector<double> gam(lpclen + 1, 0);
		vector<double>D(lpclen, 0);
		lpc[0] = 1;
		lpc[1] = -1 * R[1] / R[0];
		for (int l = 2; l < lpclen + 1; l++)
			lpc[l] = 0;

		for (int k = 0; k < lpclen - 1; k++)
		{
			E_V[k + 1] = R[0];
			for (int j = 1; j <= k + 1; j++)
				E_V[k + 1] += lpc[j] * R[j];
			for (int j = 0; j <= k + 1; j++)
				D[k + 1] += lpc[j] * R[k + 2 - j];

			gam[k + 2] = D[k + 1] / E_V[k + 1];

			E_V[k + 2] = E_V[k + 1] * (1 - gam[k + 2] * gam[k + 2]);
			a_yuce[0] = 1;
			for (int q = 0; q <= k; q++)
				a_yuce[q + 1] = lpc[q + 1] - gam[k + 2] * lpc[k - q + 2 - 1];

			a_yuce[k + 2] = -gam[k + 2];
			for (int j = 0; j < lpclen + 1; j++)
				lpc[j] = a_yuce[j];
		}
		double EE = E_V[lpclen];
		if (EE < 0)
			EE = 0;
		vector<double> xiao;
		/*	double max = 0;
		for (int i = 0; i < 160; i++)
		{
		if (abs(temp2[i])>max)
		max = abs(temp2[i]);
		}
		for (int i = 0; i <160; i++)
		{
		double x = temp2[i];
		if (x>0.5*max)
		x = x - 0.5*max;
		else if (x < -0.5*max)
		x = x + 0.5*max;
		else x = 0;
		xiao.push_back(x);
		}
		*/
		for (int i = 0; i <160; i++)
		{
			double x = temp2[i];
			xiao.push_back(x);
		}
		double *A = new double[160];

		for (int k = 0; k < 160; k++)
		{
			double sum = 0;
			for (int m = 0; m < 160; m++)
			{
				if (m + k >= 160)
					sum += 0;
				else
					sum += xiao[m] * xiao[m + k];
			}
			A[k] = sum;
		}


		double m = 0.25;
		/*for (int i = 0; i < A.size(); i++)
		{

		A[i] = A[i] / m;
		}*/
		Time = 0;
		for (int jj = 24; jj < 160; jj++)
		{
			if (A[jj] / A[0]>m)
			{
				m = A[jj] / A[0];
				Time = jj;
			}
		}
		delete[]A;

		//while (Time > 50 && A[Time / 2] > 0.8*A[Time])
		//	Time = Time / 2;

		//	setxx(A);
		//	showd(showfun, 0, 159, -1, 1);

		/*	vector<double> amdf;
		for (int i = 0; i<160; i++)
		{
		double sum = 0;
		for (int j = 0; j < 160; j++)
		{
		if (j + i>=160)
		sum +=abs(temp2[j]);
		else
		sum += abs(temp2[j] - temp2[j + i]);
		}
		amdf.push_back(sum);
		}
		setxx(amdf);
		showd(showfun, 0, 159, 0, 1000);
		*/
		/*	int count = 0;
		for (int i = 1; i<160; i++)
		{
		count += abs((int)(temp2[i]> 0) - (int)(temp2[i- 1] > 0));
		}*/
		//printf("过零率为 %d", count);
		if (Time == 0)        //清音
		{
			double *noise = new double[160];
			double* buffd = new double[160];
			whitenoise(160, buffd, EE*EE);
			for (int l = 0; l < 160; l++)
			{
				buffd[l] = 0;
			}
			//	conv(noise, lpc, buffd, 150, 11);
			//filter(10, lpc, one, 160, noise, buffd);
			for (int kk = 0; kk < 160; kk++)
			{
				buffer2[i + kk] = (int)buffd[kk] + 126;
			}

			delete[]buffd, noise;
		}
		else
		{
			double *priord = new double[160];
			double *buffd = new double[160];
			bool sign = 0;
			for (int kk = 0; kk < 160; kk++)
			{
				if (kk%Time == 0)
				{
					priord[kk] = sign ? -1 : 1 * EE * 10 / ceil(160.0 / Time);//*sqrt(EE / ceil(160.0 / Time) );
					sign = ~sign;
				}
				else priord[kk] = 0;
			}
			//conv(priord, lpc, buffd, 150, 11);
			for (int l = 0; l < 160; l++)
			{
				buffd[l] = 0;
			}

			filter(10, lpc, one, 160, priord, buffd);
			for (int kk = 0; kk < 160; kk++)
			{
				buffer2[i + kk] = (int)buffd[kk] + 126;
			}
			//	delete[]buffd,priord;
		}
		delete[]temp2;
	}
	delete[] one;
}

// 录音回调函数   
void CALLBACK waveInProc(HWAVEIN hwi,
	UINT uMsg,
	DWORD_PTR dwInstance,
	DWORD_PTR dwParam1,
	DWORD_PTR dwParam2)
{
	LPWAVEHDR pwh = (LPWAVEHDR)dwParam1;

	if ((WIM_DATA == uMsg) && (buf_count<BUFFER_SIZE))
	{
		int temp = BUFFER_SIZE - buf_count;
		temp = (temp>pwh->dwBytesRecorded) ? pwh->dwBytesRecorded : temp;
		memcpy(buffer + buf_count, pwh->lpData, temp);
		buf_count += temp;

		waveInAddBuffer(hwi, pwh, sizeof(WAVEHDR));
	}
}

// 放音回调函数   
void CALLBACK waveOutProc(HWAVEOUT hwo,
	UINT uMsg,
	DWORD_PTR dwInstance,
	DWORD_PTR dwParam1,
	DWORD_PTR dwParam2)
{
	if (WOM_DONE == uMsg)
	{
		//buf_count = BUFFER_SIZE;
		buf_count = (endd - start) * 80;
	}
}