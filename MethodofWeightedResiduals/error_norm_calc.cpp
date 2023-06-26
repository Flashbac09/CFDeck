/*
Jia [2023.04]

Calculate specific error produced by several Weighted Residual Method.
L[u(x)]=u''(x)+u;f(x)=-x;0<x<1,x(0)=0,x(1)=0.
In this program you can't set the dimension of galerkin method, it's 2 fixed.

Check [equation settings] & [numerical settings] before compiling & running the program.

*/
#define _USE_MATH_DEFINES // to utilize cmath lib
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <fstream>
using namespace std;
double u0[10000], u1[10000], u2[10000], u3[10000];
double dx = 0.001;
double e1, e2, e3, e4, e5; // 配置法、子区域法、最小二乘法、矩法、Galerkin法.
int grid_x = 1000;
double x;
inline double theoretical_result(double x)
{
    return ((1.0 / sin(1.0)) * sin(x) - x);
}
int main()
{
    for (int i = 0; i < grid_x; ++i)
    {
        x = ((double)i) / grid_x;
        e1 += pow((x * (1 - x) * (0.173076923 + 0.194711539 * x) - theoretical_result(x)), 2);
        e2 += pow((x * (1 - x) * (0.187620890 + 0.170212766 * x) - theoretical_result(x)), 2);
        e3 += pow((x * (1 - x) * (0.187541897 + 0.169470661 * x) - theoretical_result(x)), 2);
        e4 += pow((x * (1 - x) * (0.18798151 + 0.169491525 * x) - theoretical_result(x)), 2);
        e5 += pow((x * (1 - x) * (0.192411924 + 0.170731707 * x) - theoretical_result(x)), 2);
    }
    e1 *= dx;
    e2 *= dx;
    e3 *= dx;
    e4 *= dx;
    e5 *= dx;
    e1 = sqrt(e1);
    e2 = sqrt(e2);
    e3 = sqrt(e3);
    e4 = sqrt(e4);
    e5 = sqrt(e5);
    cout << "比较各个误差范数：\n"
         << "配置法:" << e1 << endl
         << "子区域法:" << e2 << endl
         << "最小二乘法:" << e3 << endl
         << "矩法:" << e4 << endl
         << "Galerkin法:" << e5 << endl;
}