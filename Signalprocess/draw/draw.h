#pragma once

#ifndef CImg.h
#include "../CImg/CImg.h"
#endif // !CImg.h



using namespace cimg_library;



struct Point                                   //opencv2��Point��
{
	int x;
	int y;
	Point() {}
	Point(int a, int b) { x = a, y = b; }
};

struct Scalar                               //opencv2��Scalar�࣬ע��Scalar����ɫ˳��ΪBGR�洢����CImg��RGB�෴��
{
	unsigned char blue;
	unsigned char green;
	unsigned char red;
	Scalar(unsigned char x, unsigned char y, unsigned char z) { blue = x, green = y, red = z; }
};
void line(CImg<unsigned char> image, Point a, Point b, Scalar color, int opacity = 1);    //������������Ϊopencv��ʽ��Ϊ�˴��븴�ã���д����ʹ���߽ӿ�ͳһ��

void circle(CImg<unsigned char>image, Point a, int radius, Scalar color, int flag = 1);

void putText(CImg<unsigned char>image, const char str[], Point a, int thickness, float frontheight, Scalar color);



class functionobject
{
public:

	virtual double operator() (double x) = 0;

};

class showfunction :public functionobject
{
public:
	double objectfunction(double x) { return x; };
	virtual double operator()(double x) {
		return objectfunction(x);
	}
};

int showc(double(*ptr)(double),//Ҫ��ʾ�ĺ���
	int xlow,//x�����ʼλ��
	int xhigh,//x�������λ��
	int ylow,//y�����ʼλ��
	int yhigh,//y�������λ��
	int dotted = 1 //�Ƿ񻭿�����,Ĭ��Ϊ1��������,0�򲻻���
);//����������ͼ���������������ֱ�Ϊ��Ҫ��ʾ�ĺ���ָ�룬x��ʼ���꣬x�������꣬y��ʼ���꣬y�������꣬�Ƿ����߿�



int showd(double(*ptr)(int),//Ҫ��ʾ����ɢ�ź�
	int xlow,//x�����ʼλ��
	int xhigh,//x�������λ��
	int ylow,//y�����ʼλ��
	int yhigh,//y�������λ��
	int dotted = 0     //��ɢ�ź�Ĭ�ϲ�������
);//��ɢ�ź���ͼ���������������ֱ�Ϊ��Ҫ��ʾ���źŵĺ���ָ�룬x��ʼ���꣬x�������꣬y��ʼ���꣬y�������꣬�Ƿ�����
