#include "..\include\vector.h"

double *v_create(int n){
	if(n<=0)
		exit(EXIT_FAILURE);
	double *v = (double*)calloc(n, sizeof(double));
	if(v==NULL)
		exit(EXIT_FAILURE);
	return v;

}

double v_show(double *v, int n){
	if(v==NULL || n<=0)
		exit(EXIT_FAILURE);
	for(int i=0; i<n; i++)
		printf("%lf ",v[i]);
	printf("\n");
}


double v_norm(double *v, int n){
	if(v==NULL || n<=0)
		exit(EXIT_FAILURE);
	double r = 0;
	
	for(int i=0; i<n; i++)
		r += pow(v[i],2);
	return sqrt(r);
	printf("\n");

}

double v_dot(double *v, double *w, int n){
	if(v==NULL || n<=0)
		exit(EXIT_FAILURE);
	double r = 0;
	
	for(int i=0; i<n; i++)
		r += v[i] * w[i];
	return r;
	printf("\n");

}

double *v_cross(double *v, double *w, int n){
	if(v==NULL || n<=0)
		exit(EXIT_FAILURE);
	 double *r = v_create(3);
    if (r == NULL) {
        fprintf(stderr, "Error al asignar memoria.\n");
        exit(EXIT_FAILURE);
    }

    r[0] = v[1] * w[2] - v[2] * w[1];
    r[1] = v[2] * w[0] - v[0] * w[2];
    r[2] = v[0] * w[1] - v[1] * w[0];

    return r;
	printf("\n");

}
