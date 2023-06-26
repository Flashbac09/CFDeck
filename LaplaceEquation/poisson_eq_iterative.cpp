/*

[Jia/Flashbac09] [2023.04]

Numercial Solution of 2-dimensional
poisson equation
with rectangular area & boundary conditions.

methods:
    1. Fourier Transform [haven't done yet]
    2. Iterative method, Gauss-Seidel & SOR.

If you just want to solve laplace equation, you can delete those parts related to func(),
which contains data on RHS of the poisson equation.

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

// -----numerical settings----- //
const int T_LIMIT = 100000;
const double EPS_0 = 0.0001;
const double L = 0.0, R = 4.0, D = 0.0, U = 3.0; // support rectangle region only
const int grid_x = 400, grid_y = 300;
double dx = 0.01, dy = 0.01;
double u[grid_x+10][grid_y+10];
// double f[grid_x][grid_y];
double beta;
// -----numerical settings----- //

// -----equation settings----- //
double func(int i, int j) // RHS f(x,y)
{
    return 0;
}
void poisson_boundary_value()
{
    for (int i = 0; i < grid_x; ++i)
    {
        u[i][0] = 100;     // d
        u[i][grid_y] = 20; // u
    }
    for (int i = 0; i < grid_y; ++i)
    {
        u[0][i] = 20;      // l
        u[grid_x][i] = 20; // r
    }
}
// -----equation settings----- //

// -----generated variables----- //

// -----generated variables----- //

// -----manipulation----- //
void func_init()
{
    for (int i = 0; i <= grid_x; ++i)
    {
        for (int j = 0; j <= grid_y; ++j)
        {
            // f[i][j] = func(i, j);
        }
    }
}
void output_u()
{
    ofstream fout;
    string ss="x"+to_string(grid_x)+"y"+to_string(grid_y)+".txt";
    fout.open(ss);
    for(int j=0;j<=grid_y;++j)
    {
        for(int i=0;i<=grid_x;++i)
        {
            fout<<u[i][j]<<"   ";
        }
        fout<<endl;
    }
    fout.close();
}
// -----manipulation----- //

int main()
{

    // 2. Iterative, Gauss Seidel,SOR
    beta = dx / dy;
    // for poisson:
    //   func_init();  //
    poisson_boundary_value();
    int t = 0;
    double eps = 0xffffff, temp = 0;
    while (t <= T_LIMIT && eps > EPS_0)
    {
        eps=0;
        for (int i = 1; i < grid_x; ++i)
        {
            for (int j = 1; j < grid_y; ++j)
            {
                // for laplace:
                temp = u[i][j];
                u[i][j] = (u[i + 1][j] + u[i - 1][j] + beta * beta * (u[i][j + 1] + u[i][j - 1])) / (2 * (1 + beta * beta));
                if (fabs(u[i][j] - temp) > eps)
                    eps = fabs(u[i][j] - temp);
                // for poisson:
                // u[i][j]=(u[i+1][j]+u[i-1][j]+beta*beta*(u[i][j+1]+u[i][j-1])+dx*dx*f[i][j])/(2*(1+beta*beta));
            }
        }
        t+=1;
        if (t % 100 == 0)
            cout << setprecision(2) << ((double)t * 100.0) / T_LIMIT << "%, eps=" << eps << endl;
    }
    output_u();
}