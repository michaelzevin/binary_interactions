#include <gsl/gsl_sf_ellint.h>
void c_wrapper_gsl_sf_ellint_rj_(double* R, double* x, double* y, double* z, double* p, double* err){
        *R = gsl_sf_ellint_RJ(*x, *y, *z, *p, *err);
}
