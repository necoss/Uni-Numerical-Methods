#include "trapezoid.h"

double function(double x) {
    return 0.1 * pow(x, 2) / log10(x);
}
double trapezoid(double a, double b, double eps) {
    int n = 1;
    double h = (b - a) / n; // step
    double prev_step = (function(a) + function(b)) / 2; // initialize previous step
    double step = prev_step;
    double e = eps + 1;
    int max_n = 1000000; // limit on number of segments
    while (e > eps && n < max_n) {
        n *= 2;
        h /= 2;
        prev_step = step;
        step = 0;
        for (int i = 1; i <= n - 1; i++) {
            double x = a + i * h;
            step += function(x);
        }
        step += (function(a) + function(b)) / 2;
        step *= h;
        e = fabs(step - prev_step) / 3.0;
    }
    return step;
}