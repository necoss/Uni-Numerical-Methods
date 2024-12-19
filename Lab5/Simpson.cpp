#include "simpson.h"

double function_1(double x) {
    return 0.1 * pow(x, 2) / log10(x);
}

namespace simpson {
    double function_simpson(double x, double y) {
        return 4 - x * x - y * y;
    }

    double one(double a, double b, double eps) {
        int n = 2; //num of segments, should be even (%2 = 0)
        double h = (b - a) / n;
        double prev_s =
            (function_1(a) + 4 * function_1((a + b) / 2) + function_1(b)) * h / 3; // initialize previous simpson
        double s = prev_s;
        double e = eps + 1;
        int max_n = 1000000;
        while (e > eps && n < max_n) {
            n *= 2;
            h /= 2;
            prev_s = s;
            s = 0;
            for (int i = 1; i <= n - 1; i += 2) {
                double x = a + i * h;
                s += 4 * function_1(x);
            }
            for (int i = 2; i <= n - 2; i += 2) {
                double x = a + i * h;
                s += 2 * function_1(x);
            }
            s += function_1(a) + function_1(b);
            s *= h / 3;
            e = fabs(s - prev_s) / 15.0;
        }
        return s;
    }

    double twin(double a, double b, double c, double d, int n, int m) {
        double h = (b - a) / n;
        double k = (d - c) / m;
        double sum = 0.0;

        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= m; j++) {
                double x = a + i * h;
                double y = c + j * k;
                double z = function_simpson(x, y);

                if (i == 0 || i == n) {
                    z /= 2;
                }
                if (j == 0 || j == m) {
                    z /= 2;
                }
                sum += z;
            }
        }

        return sum * h * k;
    }
}