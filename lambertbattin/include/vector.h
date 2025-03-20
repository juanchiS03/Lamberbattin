#ifndef _VECTOR_
#define _VECTOR_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double *v_create(int n);
double v_show(double *v, int n);
double v_norm(double *v, int n);
double v_dot(double *v, double *w, int n);
double *v_cross(double *v, double *w, int n);

#endif