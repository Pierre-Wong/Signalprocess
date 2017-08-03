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



int showc(double(*ptr)(double),//Ҫ��ʾ�ĺ���
	int xlow,//x�����ʼλ��
	int xhigh,//x�������λ��
	int ylow,//y�����ʼλ��
	int yhigh,//y�������λ��
	int dotted  //�Ƿ񻭿�����,Ĭ��Ϊ1��������,0�򲻻���
)//����������ͼ���������������ֱ�Ϊ��Ҫ��ʾ�ĺ���ָ�룬x��ʼ���꣬x�������꣬y��ʼ���꣬y�������꣬�Ƿ����߿�
{
	char wndname[] = "�ź���ϵͳ";
	char str[10];

	int  width = 600, height = 600;


	//Mat image = Mat::zeros(height, width, CV_8UC3);
	CImg<unsigned char> image(height, width, 1, 3, 0);
	image._is_shared = 1;
	image.fill(255);
	//	imshow(wndname, image);


	Point pt1, pt2, pt3, pt4;
	pt1.x = 40;
	pt1.y = 600 - 40;          //���¶˵�  ע��opencv�������Ͻǵ�������(0,0)��������xֵ���ӣ����ϵ���yֵ���ӡ�
	pt2.x = 600 - 40;
	pt2.y = 600 - 40;		 //���¶˵�
	pt3.x = 40;
	pt3.y = 80;                  //���϶˵�
	pt4.x = 600 - 40;
	pt4.y = 80;                  //���¶˵�


								 //����xlow��xhigh��ֵ����������Y��λ��
	Point orial, orial0;
	Scalar black = Scalar(0, 0, 0);
	int texflag = 1;
	orial.y = height - 40;
	orial0.y = 40;
	if (xlow >= 0)
	{
		orial.x = pt1.x;
		orial0.x = orial.x;
		line(image, orial, orial0, black);
	}
	else if (xhigh <= 0)
	{
		orial.x = pt1.x + 480;
		orial0.x = orial.x;
		line(image, orial, orial0, black);
		texflag = 0;
	}
	else
	{
		orial.x = (double)(0 - xlow) / (xhigh - xlow) * 480 + pt1.x;
		orial0.x = orial.x;
		line(image, orial, orial0, black);
	}

	//����ylow��yhigh��ֵ����������X��λ��
	Point orialy, orialy0;
	orialy.x = 40;
	orialy0.x = 560;
	if (ylow >= 0)
	{
		orialy.y = pt1.y;
		orialy0.y = orialy.y;
		line(image, orialy, orialy0, black);
	}
	else if (yhigh <= 0)
	{
		orialy.y = pt3.y;
		orialy0.y = orialy.y;
		line(image, orialy, orialy0, black);
	}
	else
	{
		orialy.y = pt1.y - (double)(0 - ylow) / (yhigh - ylow) * 480;
		orialy0.y = orialy.y;
		line(image, orialy, orialy0, black);
	}


	int count = 0;
	for (int i = pt1.x; i < pt2.x; i = i + 40)       //����������ѭ��
	{
		if (dotted&&i != pt1.x)                      //������
		{
			for (int j = pt1.y; j - 10 >= pt3.y; j = j - 10)
			{
				Point pta, ptb;
				pta.x = i;
				ptb.x = i;
				pta.y = j;
				ptb.y = j - 5;
				line(image, pta, ptb, black, 1);
			}
		}
		Point tex, axis1, axis2;
		tex.x = i - 10;
		tex.y = orialy.y + 20;
		axis1.x = i;
		axis1.y = orialy.y;
		axis2.x = i;
		axis2.y = orialy.y - 5;
		circle(image, axis1, 2, black);
		line(image, axis1, axis2, black);
		double a = (double)(xhigh - xlow) / 12.0 *count + xlow;
		if (abs(xhigh - xlow)>1)
			sprintf(str, "%.1lf", a);
		else
			sprintf(str, "%.3lf", a);
		putText(image, str, tex, 1,
			0.8, black);
		count++;
	}



	count = 0;
	for (int i = pt1.y; i > pt3.y; i = i - 40)//����������ѭ��
	{
		if (dotted&&i != pt1.y)                   //������
		{
			for (int j = pt1.x; j + 10 <= pt2.x; j = j + 10)
			{
				Point pta, ptb;
				pta.y = i;
				ptb.y = i;
				pta.x = j;
				ptb.x = j + 5;
				line(image, pta, ptb, black, 1);
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
		line(image, axis1, axis2, black);
		circle(image, axis1, 2, black);
		double a = (double)(yhigh - ylow) / 12.0 *count + ylow;
		if (abs(yhigh - ylow)>1)
			sprintf(str, "%.1lf", a);
		else
			sprintf(str, "%.3lf", a);
		putText(image, str, tex, 1,
			0.8, black);
		count++;
	}



	Scalar red = Scalar(0, 0, 255);
	//���ݺ�����������ͼ��
	for (int i = 0; i < 480; i = i + 2)
	{
		Point pta, ptb, ptc;
		pta.x = 40 + i;
		ptb.x = 40 + i + 1;
		pta.y = pt1.y - (*ptr)(i*(xhigh - xlow) / 480.0 + xlow) * 480.0 / (yhigh - ylow) + ylow * 480.0 / (yhigh - ylow);
		ptb.y = pt1.y - (*ptr)((i + 1)*(xhigh - xlow) / 480.0 + xlow) * 480.0 / (yhigh - ylow) + ylow * 480.0 / (yhigh - ylow);
		ptc.x = 40 + i + 2;
		ptc.y = pt1.y - ((*ptr)((i + 2)*(xhigh - xlow) / 480.0 + xlow)) * 480.0 / (yhigh - ylow) + ylow * 480.0 / (yhigh - ylow);
		line(image, pta, ptb, red, 1);
		line(image, ptb, ptc, red, 1);
	}
	image.display(wndname);

	return 0;
}


int showd(double(*ptr)(int),//Ҫ��ʾ����ɢ�ź�
	int xlow,//x�����ʼλ��
	int xhigh,//x�������λ��
	int ylow,//y�����ʼλ��
	int yhigh,//y�������λ��
	int dotted     //��ɢ�ź�Ĭ�ϲ�������
)//��ɢ�ź���ͼ���������������ֱ�Ϊ��Ҫ��ʾ���źŵĺ���ָ�룬x��ʼ���꣬x�������꣬y��ʼ���꣬y�������꣬�Ƿ�����
{
	char wndname[] = "�ź���ϵͳ";
	char str[10];

	int  width = 600, height = 600;


	CImg<unsigned char> image(height, width, 1, 3, 0);
	image._is_shared = 1;
	image.fill(255);



	Point pt1, pt2, pt3, pt4;
	pt1.x = 40;
	pt1.y = 600 - 40;          //���¶˵�  ע��opencv�������Ͻǵ�������(0,0)��������xֵ���ӣ����ϵ���yֵ���ӡ�
	pt2.x = 600 - 40;
	pt2.y = 600 - 40;		 //���¶˵�
	pt3.x = 40;
	pt3.y = 40;                  //���϶˵�
	pt4.x = 600 - 40;
	pt4.y = 40;                  //���¶˵�



								 //resize(image2, image, Size(600, 600), 0, 0, INTER_AREA);


								 //����xlow��xhigh��ֵ����������Y��λ��
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
		orial.x = (double)(-xlow) / (xhigh - xlow) * 520 + pt1.x;
		orial0.x = orial.x;
		line(image, orial, orial0, Scalar(0, 0, 0));
	}

	//����ylow��yhigh��ֵ����������X��λ��
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

	//���ݺ�����������ͼ��
	for (int i = xlow; i <= xhigh; i++)
	{
		Point pta, ptb, tex;
		pta.x = 40 + 520.0 / (xhigh - xlow)* (i - xlow);
		pta.y = pt1.y - ptr(i) * 480.0 / (yhigh - ylow) + ylow * 480.0 / (yhigh - ylow);
		ptb.x = pta.x;
		ptb.y = orialy.y;

		tex.x = ptb.x - 5;
		tex.y = orialy.y + 10;
		sprintf(str, "%d", i);
		putText(image, str, tex, 1,
			0.8, Scalar(0, 0, 0));
		circle(image, pta, 3, Scalar(0, 0, 0));

		line(image, pta, ptb, Scalar(255, 0, 0), 1);
	}



	int count = 0;
	for (int i = pt1.x; i < pt2.x; i = i + 40)       //����������ѭ��
	{
		if (dotted&&i != pt1.x)                      //������
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

	}



	count = 0;
	for (int i = pt1.y; i > pt3.y; i = i - 40)//����������ѭ��
	{
		if (dotted&&i != pt1.y)                   //������
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