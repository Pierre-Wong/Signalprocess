#pragma once

#ifndef CImg.h
#include "../CImg/CImg.h"
#endif // !CImg.h



using namespace cimg_library;



struct Point                                   //opencv2的Point类
{
	int x;
	int y;
	Point() {}
	Point(int a, int b) { x = a, y = b; }
};

struct Scalar                               //opencv2的Scalar类，注意Scalar的颜色顺序为BGR存储，和CImg的RGB相反。
{
	unsigned char blue;
	unsigned char green;
	unsigned char red;
	Scalar(unsigned char x, unsigned char y, unsigned char z) { blue = x, green = y, red = z; }
};
void line(CImg<unsigned char> image, Point a, Point b, Scalar color, int opacity = 1);    //以下三个函数为opencv格式，为了代码复用，编写函数使两者接口统一。

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

int showc(double(*ptr)(double),//要显示的函数
	int xlow,//x坐标初始位置
	int xhigh,//x坐标结束位置
	int ylow,//y坐标初始位置
	int yhigh,//y坐标结束位置
	int dotted = 1 //是否画框虚线,默认为1，代表是,0则不画。
);//连续函数作图函数，六个参数分别为，要显示的函数指针，x初始坐标，x结束坐标，y初始坐标，y结束坐标，是否画虚线框。



int showd(double(*ptr)(int),//要显示的离散信号
	int xlow,//x坐标初始位置
	int xhigh,//x坐标结束位置
	int ylow,//y坐标初始位置
	int yhigh,//y坐标结束位置
	int dotted = 0     //离散信号默认不画虚线
);//离散信号作图函数，六个参数分别为，要显示的信号的函数指针，x初始坐标，x结束坐标，y初始坐标，y结束坐标，是否画虚线
