#include"draw.h"



void line(CImg<unsigned char> image, Point a, Point b, Scalar color, int opacity)    //以下三个函数为opencv格式，为了代码复用，编写函数使两者接口统一。
{
	const unsigned char colour[3] = { color.red, color.green, color.blue };
	image.draw_line(a.x, a.y, b.x, b.y, colour, (float)opacity);
}
void circle(CImg<unsigned char>image, Point a, int radius, Scalar color, int flag )
{
	const unsigned char colour[3] = { color.red, color.green, color.blue };
	if (flag)                            //flag不为0时为实心圆，否则为空心圆
		image.draw_circle((const int)a.x, (const int)a.y, radius, colour);
	else
		image.draw_circle((const int)a.x, (const int)a.y, radius, colour, 1.0, 1); //加入最后一个参数后draw_circle函数重载为画空心圆函数，具体见CImg原代码。
}
void putText(CImg<unsigned char>image, const char str[], Point a, int thickness, float frontheight, Scalar color)
{
	const unsigned char colour[3] = { color.red, color.green, color.blue };
	image.draw_text((const int)a.x, (const int)a.y, str, colour);
	//	image.display("a");
}


 

