#ifndef LABA_5_SIMPSON_H
#define LABA_5_SIMPSON_H

#endif //LABA_5_SIMPSON_H
#include <iostream>
#include <cmath>
using namespace std;
double function_1(double x);
double function_simpson(double x, double y);
namespace simpson {
    double one(double a, double b, double eps);
    double twin(double a, double b, double c, double d, int n, int m);
}