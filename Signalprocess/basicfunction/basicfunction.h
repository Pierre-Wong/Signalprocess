#pragma once

#include "../complex/mycomplex.h"



#include <windows.h>   
#include <stdio.h>   
#include <mmsystem.h>   
#include <vector>
#include <algorithm>

#include<numeric>


double show(double x);

complex frac(complex s);

double sfra(double x);

double gauss(double x);


double f(double x);




double juxing(double x);


double delta(int i);

double sinc(double x);            //sinc函数，注意x=0时要使用洛必达法则直接算出结果

double sinc(int x);


double rectwindow(int n);

double bartlettwin(int n);

double hanningwin(int n);

double hammingwin(int n);

double blackmanwin(int n);

unsigned int fact(unsigned int k);

double basel(double x);

double kaiserwin(int n);


double rec(int x);


double drect(int x);

double gauss2(double x);

double crec(double x);

double cconv(double(*f)(double), double(*g)(double), double t, double flow, double fhigh, double glow, double ghigh, double accur = 0.01);
//只在有界区间上计算，t是输入参数，f的下界是flow，上界是fhigh，g同理，accur代表黎曼积分的最小矩形区间


double showva(double x);

double lshi(int r, int b, double x);

double lsh(double x);

void conv(double u[], double v[], double w[], int m, int n);

void filter(int ord, double *a, double *b, int np, double *x, double *y);