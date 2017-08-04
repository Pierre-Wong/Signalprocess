#include "stdio.h"





void f(int x)
{
	return;
}


template<class T> void g(T a)
{
	a(1);
	return;
}





int main()
{

	g(f);
	getchar();


}