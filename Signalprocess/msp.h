#include<math.h>
#define abs_error 1.e-10

#ifndef _MSP_H_
#define _MSP_H_


float randnu(long *iseed);
void meavar(float u[], int *n, float *pum, float *puv);
void mrandom(float u[], int *n, long *piseed, int ITYPE, float p);


/*-------------------------------------------------------------------*/
#endif 