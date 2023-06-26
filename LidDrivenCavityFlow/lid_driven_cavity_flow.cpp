/*
Flashbac09 [2023.05]

Lid-driven Cavity Flow
Case:
    Lid Speed:u(x)=16*x^2*(1-x^2).
    Region:[0,1]*[0,1].
Note: Check [equation settings] & [numerical settings] before compiling & running the program.


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
// -----equation settings----- //
double nju = 0.01;
// -----equation settings----- //

// -----numerical settings----- //
int S_LIM = 100000; // threshold of cycle
const double h = 0.01;
const int N = 100;
double dt = 0.001;
// -----numerical settings----- //

// -----physics variables----- //
double psi[N + 1][N + 1];
double omega[N + 1][N + 1];
double old_omega[N + 1][N + 1];
double t[N-1][N-1]; // for fourier
double f[N-1][N-1]; // in poisson, -omega
double u[N-1][N-1]; // in poisson, psi
double temp[N-1][N-1];
// -----physics variables----- //

void region_update()
{
    double psi_dx, psi_dy, ome_dx, ome_dy, ome_laplace;
    for (int i = 0; i <= N; ++i)
    {
        for (int j = 0; j <= N; ++j)
        {
            old_omega[i][j] = omega[i][j];
        }
    }
    for (int i = 1; i < N; ++i)
    {
        for (int j = 1; j < N; ++j)
        {
            psi_dx = (psi[i + 1][j] - psi[i - 1][j]) / (2.0 * h);
            psi_dy = (psi[i][j + 1] - psi[i][j - 1]) / (2.0 * h);
            ome_dx = (old_omega[i + 1][j] - old_omega[i - 1][j]) / (2.0 * h);
            ome_dy = (old_omega[i][j + 1] - old_omega[i][j - 1]) / (2.0 * h);
            ome_laplace = (old_omega[i + 1][j] + old_omega[i - 1][j] + old_omega[i][j + 1] + old_omega[i][j - 1] - 4.0 * old_omega[i][j]) / (h * h);
            omega[i][j] = old_omega[i][j] + dt * (nju * ome_laplace + psi_dx * ome_dy - psi_dy * ome_dx);
        }
    }
}

void fourier_set_mat()
{
    for (int i = 0; i < N - 1; ++i)
    {
        for (int j = 0; j < N - 1; ++j)
        {
            t[i][j] = sin((i + 1) * (j + 1) * (M_PI / (double)N));
        }
    }
}
void fourier_on_f()
{
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = 0; j < N - 1; j++)
        {
            f[i][j] = -omega[i + 1][j + 1];
        }
    }
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = 0; j < N - 1; j++)
        {
            temp[i][j] = 0;
            for (int k = 0; k < N - 1; k++)
            {
                temp[i][j] += t[i][k] * f[k][j];
            }
        }
    }
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = 0; j < N - 1; j++)
        {
            f[i][j] = 0;
            for (int k = 0; k < N - 1; k++)
            {
                f[i][j] += temp[i][k] * t[k][j];
            }
        }
    }
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = 0; j < N - 1; j++)
        {
            f[i][j] = f[i][j] * 4.0 / pow((double)N, 2);
        }
    }
}
void fourier_calc_u()
{
    double lambda = 0, ilambda = 0, jlambda = 0;
    for (int i = 0; i < N - 1; i++)
    {
        ilambda = (-4) / pow(h, 2.0) * pow(sin((i + 1) * M_PI / (2 * (double)N)), 2.0);
        for (int j = 0; j < N - 1; j++)
        {
            jlambda = (-4) / pow(h, 2.0) * pow(sin((j + 1) * M_PI / (2 * (double)N)), 2.0);
            lambda = ilambda + jlambda;
            u[i][j] = f[i][j] / lambda;
        }
    }
}
void refourier_on_u()
{
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = 0; j < N - 1; j++)
        {
            temp[i][j] = 0;
            for (int k = 0; k < N - 1; k++)
            {
                temp[i][j] += t[i][k] * u[k][j];
            }
        }
    }
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = 0; j < N - 1; j++)
        {
            u[i][j] = 0;
            for (int k = 0; k < N - 1; k++)
            {
                u[i][j] += temp[i][k] * t[k][j];
            }
        }
    }
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = 0; j < N - 1; j++)
        {
            psi[i + 1][j + 1] = u[i][j];
        }
    }
}
void fourier_solve_poisson()
{
    fourier_on_f();
    fourier_calc_u();
    refourier_on_u();
}
void thom_bound()
{
    double u;
    double x;
    for (int i = 1; i < N; i++)
    {
        omega[i][0] = -2 * (psi[i][1] - psi[i][0]) / (h * h);
    }
    for (int j = 1; j < N; j++)
    {
        omega[0][j] = -2 * (psi[1][j] - psi[0][j]) / (h * h);
        omega[N][j] = -2 * (psi[N - 1][j] - psi[N][j]) / (h * h);
    }
    for (int i = 1; i < N; i++)
    {
        x = (double)i * h;
        u = 16.0 * x * x * (1.0 - x) * (1.0 - x);
        omega[i][N] = -2 * (psi[i][N - 1] - psi[i][N] + u * h) / (h * h);
    }
    return;
}
void output()
{
    ofstream fout("omega_nju_" + to_string(nju) + "h_" + to_string(h) + "S_" + to_string(S_LIM) + ".txt");
    for (int j = 0; j <= N; j++)
    {
        for (int i = 0; i <= N; i++)
        {
            fout << omega[i][j] << " ";
        }
        fout << endl;
    }
    fout.close();
    fout.open("psi_nju_" + to_string(nju) + "h_" + to_string(h) + "S_" + to_string(S_LIM) + ".txt");
    for (int j = 0; j <= N; j++)
    {
        for (int i = 0; i <= N; i++)
        {
            fout << psi[i][j] << " ";
        }
        fout << endl;
    }
    fout.close();
}
int main()
{
    fourier_set_mat();
    for (int i = 0; i < S_LIM; ++i)
    {
        region_update();
        fourier_solve_poisson();
        thom_bound();
    }
    output();
}