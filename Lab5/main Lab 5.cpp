#include "trapezoid.h"
#include "simpson.h"

int main() {
    double a = 2.0;
    double b = 3.104;
    double eps = pow(10, -6);
    double a2 = -1.0, b2 = 1.0;
    double c = -1.0, d = 1.0;
    int n = 1000, m = 1000;

    cout << "The value of the integral of the function f(x) on the interval [" << a << ", " << b << "] by the trapezoidal method with accuracy " << eps << " is " << trapezoid(a, b, eps) << endl;
    cout << "The value of the integral of the function f(x) on the interval [" << a << ", " << b << "] using the Simpson method with accuracy " << eps << " is " << simpson::one(a, b, eps) << endl;
    cout << "The value of the double integral of the function f(x) on the interval [" << a2 << ", " << b2 << "] and [" << c << ", " << d << "]using the Simpson method with num of partitions: n = " << n << " m = " << m << " is " << simpson::twin(a2, b2, c, d, n, m) << endl;
    return 0;
}