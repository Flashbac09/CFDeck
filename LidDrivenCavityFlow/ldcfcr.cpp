/*
[Flashbac09] [2023.06]

Lid-driven Cavity Flow in Complex Region
Case:
    Lid Speed:u(x)=16*x^2*(1-x^2).
    Region:[0,1]*[0,1].

    Complex Region:
        1. 0.50*0.50 tile
        2. 0.15*0.25 tile in middle
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
double nju = 0.0001;
// -----equation settings----- //

// -----optimized variables----- //
double w;
// -----optimized variables----- //

// -----numerical settings----- //
int SOR_CYCLE_LIM = 300; // SOR cycle
int SOR_OPT_LIM = 500;   // To optimize SOR factor
int S_LIM = 100000;       // threshold of cycle
const double h = 0.01;
const int N = 100; // grid
double dt = 0.001; // time steps
// -----numerical settings----- //

// -----physics variables----- //
double psi[N + 1][N + 1];
double old_psi[N + 1][N + 1];
double omega[N + 1][N + 1];
double old_omega[N + 1][N + 1];
int ZoneCate[N + 1][N + 1];
// -----physics variables----- //
void preset_SOR_bound()
{
    for (int i = 0; i < N; ++i)
    {
        psi[i][0] = 100; // d
        psi[i][N] = 20;  // u
    }
    for (int i = 0; i < N; ++i)
    {
        psi[0][i] = 20; // l
        psi[N][i] = 20; // r
    }
}
void SOR_opt()
{
    cout << "Obtain the best relaxation factor w of SOR Iteration Method." << endl;
    ofstream fout;
    fout.open("SOR-w-test.txt");
    fout << "To obtain the best relaxation factor of Successive Over Relaxation Iteration Method" << endl;
    fout << "[ w - iteration rounds - final eps ]" << endl;
    double w_test, EPS_0_SOR = 1e-4;
    int cycle_min = SOR_OPT_LIM;
    for (w_test = 0.9; w_test <= 2; w_test += 0.01)
    {
        int cycle = 0;
        double eps = 0xfffffff, temp = 0;
        memset(psi, 0, sizeof(psi));
        preset_SOR_bound();
        while (cycle < SOR_OPT_LIM && eps > EPS_0_SOR)
        {
            eps = 0;
            for (int i = 0; i <= N; i++)
            {
                for (int j = 0; j <= N; j++)
                {

                    // simplified formula:
                    // beta=dx/dy=dh/dh=1
                    temp = psi[i][j];
                    psi[i][j] += 0.25 * w_test * (psi[i + 1][j] + psi[i - 1][j] + psi[i][j + 1] + psi[i][j - 1] + h * h * omega[i][j]) - w_test * psi[i][j];
                    if (fabs(psi[i][j] - temp) > eps)
                        eps = fabs(psi[i][j] - temp);
                }
            }
            cycle += 1;
        }
        fout << w_test << "  " << cycle << "  " << eps << endl;
        // cout << w_test << "  " << cycle << "  " << eps << endl;
        if (cycle < cycle_min)
        {
            w = w_test;
            cycle_min = cycle;
        }
    }
    fout << "The best factor is " << w << "." << endl;
    cout << "The best factor is " << w << "." << endl;
    memset(psi, 0, sizeof(psi));
    memset(omega, 0, sizeof(omega));
    cout << "Run vorticity-flow function iteration." << endl;
}
void Complex_Region_A()
{
    for (int i = 1; i <= 50; i++)
    {
        for (int j = 51; j < 100; j++)
        {
            ZoneCate[i][j] = 1;
        }
    }
    for (int i = 51; i < 100; i++)
    {
        for (int j = 1; j < 100; j++)
        {
            ZoneCate[i][j] = 1;
        }
    }
    for (int i = 0; i < 50; i++)
    {
        for (int j = 0; j < 50; j++)
        {
            ZoneCate[i][j] = -1;
        }
    }
    return;
}
void Complex_Region_B()
{
    for (int i = 1; i < 25; i++)
    {
        for (int j = 1; j < 100; j++)
        {
            ZoneCate[i][j] = 1;
        }
    }
    for (int i = 25; i < 51; i++)
    {
        for (int j = 1; j < 60; j++)
        {
            ZoneCate[i][j] = 1;
        }
        for (int j = 76; j < 100; j++)
        {
            ZoneCate[i][j] = 1;
        }
    }
    for (int i = 51; i < 100; i++)
    {
        for (int j = 1; j < 100; j++)
        {
            ZoneCate[i][j] = 1;
        }
    }
    for (int i = 26; i < 50; i++)
    {
        for (int j = 61; j < 75; j++)
        {
            ZoneCate[i][j] = -1;
        }
    }
    return;
}
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

void thom_bound()
{
    // lid:
    double u = 0, x = 0;
    for (int i = 1; i < N; i++)
    {
        x = (double)i * h;
        u = 16.0 * x * x * (1.0 - x) * (1.0 - x);
        omega[i][N] = -2 * (psi[i][N - 1] - psi[i][N] + u * h) / (h * h);
    }
    // for y<1:
    for (int i = 0; i <= N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (ZoneCate[i][j] == 0)
            {
                if (i == 0 && j != 0)
                {
                    omega[i][j] = 2.0 * (psi[i][j] - psi[i + 1][j]) / (h * h);
                }
                else if (i == 100 && j != 0)
                {
                    omega[i][j] = 2.0 * (psi[i][j] - psi[i - 1][j]) / (h * h);
                }
                else if (j == 0)
                {
                    omega[i][j] = 2.0 * (psi[i][j] - psi[i][j + 1]) / (h * h);
                }
                else
                {
                    int I = i;
                    int J = j + 1;
                    if (ZoneCate[i + 1][j] == 1)
                    {
                        I = i + 1;
                        J = j;
                    }
                    else if (ZoneCate[i][j - 1] == 1)
                    {
                        I = i;
                        J = j - 1;
                    }
                    else if (ZoneCate[i - 1][j] == 1)
                    {
                        I = i - 1;
                        J = j;
                    }
                    omega[i][j] = 2.0 * (psi[i][j] - psi[I][J]) / (h * h);
                }
            }
        }
    }
}
void output()
{
    ofstream fout("omega_nju_" + to_string(nju) + "S_" + to_string(S_LIM) + ".txt");
    for (int j = 0; j <= N; j++)
    {
        for (int i = 0; i <= N; i++)
        {
            fout << omega[i][j] << " ";
        }
        fout << endl;
    }
    fout.close();
    fout.open("psi_nju_" + to_string(nju) + "S_" + to_string(S_LIM) + ".txt");
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

void iter_solve_poisson()
{
    double EPS = 1e-8;
    int cycle = 0;
    while (cycle < SOR_CYCLE_LIM)
    {
        for (int i = 0; i <= N; i++)
        {
            for (int j = 0; j <= N; j++)
            {
                old_psi[i][j] = psi[i][j];
            }
        }
        for (int i = 0; i < 101; i++)
        {
            for (int j = 0; j < 101; j++)
            {
                if (ZoneCate[i][j] == 1)
                {
                    psi[i][j] += w / 4.0 * (psi[i + 1][j] + psi[i - 1][j] + psi[i][j + 1] + psi[i][j - 1] + h * h * omega[i][j]) - w * psi[i][j];
                }
            }
        }
        double eps_m = 0, eps = 0;
        for (int i = 0; i <= N; i++)
        {
            for (int j = 0; j <= N; j++)
            {
                if (ZoneCate[i][j] == 1)
                {
                    eps = psi[i][j] - old_psi[i][j];
                    if (eps < 0)
                        eps = 0 - eps;
                    if (eps > eps_m)
                        eps_m = eps;
                }
            }
        }
        if (eps_m < EPS)
            break;
        cycle += 1;
    }
}
int main()
{
    // 1. Optimize SOR factor w
    // 2. Set Complex Region and Boundary
    // 3. Iteration of Vorticity-Stream
    // 4. Output
    // 5. (Python)Draw

    SOR_opt();
    Complex_Region_A();
    //Complex_Region_B();
    for (int i = 0; i < S_LIM; ++i)
    {
        region_update();
        iter_solve_poisson();
        thom_bound();
        if (i % (S_LIM / 10) == 0)
            cout << i / (S_LIM / 100) << "%" << endl;
    }
    output();
}

/* Fourier method to solve Possion equation
   I didn't make it to utilize it in complex region.

//double t[N - 1][N - 1]; // for fourier
//double f[N - 1][N - 1]; // in poisson, -omega
//double u[N - 1][N - 1]; // in poisson, psi
//double temp[N - 1][N - 1];


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

*/