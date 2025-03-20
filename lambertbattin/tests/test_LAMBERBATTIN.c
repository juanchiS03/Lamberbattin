#include "..\include\vector.h"
#include "..\include\seebatt.h"
#include "..\include\seebattk.h"

void main(){
	int n=3;
	double *v = v_create(n);
	v[0] = 1;
	v[1] = 2;
	v[2] = 3;
	v_show(v,n);
	
	double *w = v_create(n);
	w[0] = 4;
	w[1] = 5;
	w[2] = 6;
	v_show(w,n);
	
	printf("%lf ", v_norm(v,n));
	printf("\n");
	printf("%lf ", v_dot(v,w,n));
	printf("\n");

	v_show(v_cross(v,w,n),n);
	printf("\n");
	printf("%lf", seebatt(3));
	printf("\n");
	printf("%lf", seebattk(3));
	
	getchar();





	double r1[] = {20.0e6, 20.0e6, 0};  // [m]
	double r2[] = {-20.0e6, 10.0e6, 0}; // [m]
	double tof = 1.0 * 86400;
	
	double *v1 = v_create(n);
	double *v2 = v_create(n);
	LAMBERTBATTIN(r1, r2, "retro", tof, v1, v2);
	v_show(v1,n);
	v_show(v2,n);
}