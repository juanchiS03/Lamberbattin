#include "..\include\vector.h"
#include <math.h>

int tests_run = 0;
double eps = 1e-10;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

int v_equals(double *v1, int n1, double *v2, int n2) {
	int equal = 1;
	
	if (n1 == n2) {
		for(int i = 0; i < n1; i++)
			if(fabs(v1[i] - v2[i]) > eps) {
				equal = 0;
				break;
			}
	}
	else
		equal = 0;
	
	return equal;
}

int v_create_01() {
	int n = 5;
    double v[] = {0,0,0,0,0};
    
    _assert(v_equals(v_create(n), n, v, n));
    
    return 0;
}
int v_norm_01() {
    double v[] = {3, 4};
    _assert(fabs(v_norm(v, 2) - 5.0) < eps);
    return 0;
}

int v_dot_01() {
    double v[] = {1, 2, 3};
    double w[] = {4, 5, 6};
    _assert(fabs(v_dot(v, w, 3) - 32.0) < eps);
    return 0;
}

int v_cross_01() {
    double v[] = {1, 2, 3};
    double w[] = {4, 5, 6};
    double expected[] = {-3, 6, -3};
    _assert(v_equals(v_cross(v, w, 3), 3, expected, 3));
    return 0;
}

int v_seebatt_01() {
    _assert(fabs(seebatt(0.5) - 0.238) < 0.01);
    return 0;
}

int v_seebattk_01() {
    _assert(fabs(seebattk(0.5) - 0.182) < 0.01);
    return 0;
}

int all_tests() {
    _verify(v_create_01);
    _verify(v_norm_01);
    _verify(v_dot_01);
    _verify(v_cross_01);
    _verify(v_seebatt_01);
    _verify(v_seebattk_01);
    return 0;
}

int main() {
    int result = all_tests();
    if (result == 0)
        printf("PASSED\n");
    printf("Tests run: %d\n", tests_run);
    getchar();
}