#include"draw.h"



void line(CImg<unsigned char> image, Point a, Point b, Scalar color, int opacity)    //������������Ϊopencv��ʽ��Ϊ�˴��븴�ã���д����ʹ���߽ӿ�ͳһ��
{
	const unsigned char colour[3] = { color.red, color.green, color.blue };
	image.draw_line(a.x, a.y, b.x, b.y, colour, (float)opacity);
}
void circle(CImg<unsigned char>image, Point a, int radius, Scalar color, int flag )
{
	const unsigned char colour[3] = { color.red, color.green, color.blue };
	if (flag)                            //flag��Ϊ0ʱΪʵ��Բ������Ϊ����Բ
		image.draw_circle((const int)a.x, (const int)a.y, radius, colour);
	else
		image.draw_circle((const int)a.x, (const int)a.y, radius, colour, 1.0, 1); //�������һ��������draw_circle��������Ϊ������Բ�����������CImgԭ���롣
}
void putText(CImg<unsigned char>image, const char str[], Point a, int thickness, float frontheight, Scalar color)
{
	const unsigned char colour[3] = { color.red, color.green, color.blue };
	image.draw_text((const int)a.x, (const int)a.y, str, colour);
	//	image.display("a");
}


 

