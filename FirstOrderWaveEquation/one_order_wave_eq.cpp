/*
[Flashbac09] [2023.04]

Numerical Solution of one-order wave equation

Preset Functions:
    TheoreticalSolution:theoretical_sol
    BoundaryConditionSet:bound_cond
    SolutionRegionSet&:Global Variables

SchemeSet:
    1. Forward Time Central Space(FTCS) scheme
    2. Lax Scheme
    3. Lax-Wendroff Scheme
    4. First Order Upwind Scheme
    5. Second Order Upwind Scheme (a>0)
    6. Upwind Warming-Beam Scheme (a>0)

Analysis:
    1. test stable c limit
    2. get order
    3. watch num/phase cut.
Theoretical Solution: u(x,t)=sin(2*PI*(x-t))
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cmath>

using namespace std;
const double PI = 3.14159265358979;

// settings
double dt = 0.001, dx = 0.001;
const double L = 0.0, R = 3.0;
double TIME_ZONE = 40;
const double a = 1;
double c = a * dt / dx;
//

int SLIM = (R - L) / dx;
int TLIM = TIME_ZONE / dt;
double u[80010][3010]; // u(t,x)
double p[20];
double x_temp;
int n;
/*
StackOverFlow:
stack 1M --128000 double
overall 2G --2.4e8 double --> 10^5 * 10^4 / 4*10^4 * 4*10^4
0-249000000 normal
249000000-269000000 silent
269000000-inf error
c=a*dt/dx,a=1;
*/
inline void FTCS_region_update()
{
    for (int j = 1; j < SLIM - 1; ++j) // in region
    {
        u[n + 1][j] = u[n][j] + (u[n][j - 1] - u[n][j + 1]) * (dt) / (2 * dx);
    }
}
inline void FTCS_bound_update()
{
    // periodic boundary
    u[n + 1][0] = u[n][0] + (u[n][SLIM - 1] - u[n][1]) * (dt) / (2 * dx);
    u[n + 1][SLIM - 1] = u[n][SLIM - 1] + (u[n][SLIM - 2] - u[n][0]) * (dt) / (2 * dx);
}
inline void Lax_region_update()
{
    for (int j = 1; j < SLIM - 1; ++j)
    {
        u[n + 1][j] = 0.5 * (1 - dt / dx) * u[n][j + 1] + 0.5 * (1 + dt / dx) * u[n][j - 1];
    }
}
inline void Lax_bound_update()
{
    u[n + 1][0] = 0.5 * (1 - dt / dx) * u[n][1] + 0.5 * (1 + dt / dx) * u[n][SLIM - 1];
    u[n + 1][SLIM - 1] = 0.5 * (1 - dt / dx) * u[n][0] + 0.5 * (1 + dt / dx) * u[n][SLIM - 2];
}
inline void Lax_Wendroff_region_update()
{
    for (int j = 1; j < SLIM - 1; ++j)
    {
        u[n + 1][j] = u[n][j] + 0.5 * dt * dt * (u[n][j + 1] - 2 * u[n][j] + u[n][j - 1]) / (dx * dx) - dt * (u[n][j + 1] - u[n][j - 1]) / (2 * dx);
    }
}
inline void Lax_Wendroff_bound_update()
{
    u[n + 1][0] = u[n][0] + 0.5 * dt * dt * (u[n][1] - 2 * u[n][0] + u[n][SLIM - 1]) / (dx * dx) - dt * (u[n][SLIM - 1]);
    u[n + 1][SLIM - 1] = u[n][SLIM - 1] + 0.5 * dt * dt * (u[n][0] - 2 * u[n][SLIM - 1] + u[n][SLIM - 2]) / (dx * dx) - dt * (u[n][0] - u[n][SLIM - 2]) / (2 * dx);
}
inline void first_order_upwind_region_update()
{
    for (int j = 1; j < SLIM; ++j)
    {
        u[n + 1][j] = u[n][j] + a * (u[n][j - 1] - u[n][j]) * dt / dx;
    }
}
inline void first_order_upwind_bound_update()
{
    u[n + 1][0] = u[n][0] + a * (u[n][SLIM - 1] - u[n][0]) * dt / dx;
}
inline void second_order_upwind_region_update()
{
    for (int j = 2; j < SLIM; ++j)
    {
        u[n + 1][j] = u[n][j] - 0.5 * dt * (3 * u[n][j] - 4 * u[n][j - 1] + u[n][j - 2]) / dx;
    }
}
inline void second_order_upwind_bound_update()
{
    u[n + 1][0] = u[n][0] - 0.5 * dt * (3 * u[n][0] - 4 * u[n][SLIM - 1] + u[n][SLIM - 2]) / dx;
    u[n + 1][1] = u[n][1] - 0.5 * dt * (3 * u[n][1] - 4 * u[n][0] + u[n][SLIM - 1]) / dx;
}
inline void WB_region_update()
{
    for (int j = 2; j < SLIM; ++j)
    {
        u[n + 1][j] = u[n][j] - c * (u[n][j] - u[n][j - 1]) - 0.5 * c * (1 - c) * (u[n][j] - 2 * u[n][j - 1] + u[n][j - 2]);
    }
}
inline void WB_bound_update()
{
    u[n + 1][0] = u[n][0] - c * (u[n][0] - u[n][SLIM - 1]) - 0.5 * c * (1 - c) * (u[n][0] - 2 * u[n][SLIM - 1] + u[n][SLIM - 2]);
    u[n + 1][1] = u[n][1] - c * (u[n][1] - u[n][0]) - 0.5 * c * (1 - c) * (u[n][1] - 2 * u[n][0] + u[n][SLIM - 1]);
}
inline void leap_frog_region_update()
{
    for (int j = 1; j < SLIM - 1; ++j)
    {
        u[n + 1][j] = u[n - 1][j] - c * (u[n][j + 1] - u[n][j - 1]);
    }
}
inline void leap_frog_bound_update()
{
    u[n + 1][0] = u[n - 1][0] - c * (u[n][1] - u[n][SLIM - 1]);
    u[n + 1][SLIM - 1] = u[n - 1][SLIM - 1] - c * (u[n][0] - u[n][SLIM - 2]);
}
inline double bound_cond(double x)
{
    return sin(2 * PI * x);
}
double theoretical_sol(double x, double t)
{
    return sin(2 * PI * (x - t));
}
void bound_init()
{
    for (int i = 0; i < SLIM; ++i)
    {
        x_temp = i * (R - L) / SLIM;
        u[0][i] = bound_cond(x_temp);
    }
}
void OutputTime(int time)
{
    ofstream fout;
    fout.open("output_t_fixed.txt");

    for (int j = 0; j < SLIM; ++j)
    {
        fout << j * (R - L) / (SLIM) << "   " << u[time][j] << endl;
    }

    fout.close();
}
void OutputSpace(int x)
{
    ofstream fout;
    fout.open("output_x_fixed.txt");

    for (int n = 0; n < TLIM; ++n)
    {
        fout << n * TIME_ZONE / (TLIM) << "   " << u[n][x] << endl;
    }

    fout.close();
}
void Iterative(string method)
{
    if (method.compare("FTCS") == 0 || method.compare("1") == 0)
    {
        for (n = 0; n <= TLIM; ++n)
        {
            FTCS_region_update();
            FTCS_bound_update();
        }
    }
    else if (method.compare("Lax") == 0 || method.compare("2") == 0)
    {
        for (n = 0; n <= TLIM; ++n)
        {
            Lax_region_update();
            Lax_bound_update();
        }
    }
    else if (method.compare("Lax-Wendroff") == 0 || method.compare("3") == 0)
    {
        for (n = 0; n <= TLIM; ++n)
        {
            Lax_Wendroff_region_update();
            Lax_Wendroff_bound_update();
        }
    }
    else if (method.compare("1-upwind") == 0 || method.compare("4") == 0)
    {
        for (n = 0; n <= TLIM; ++n)
        {
            first_order_upwind_region_update();
            first_order_upwind_bound_update();
        }
    }
    else if (method.compare("2-upwind") == 0 || method.compare("5") == 0)
    {
        for (n = 0; n <= TLIM; ++n)
        {
            second_order_upwind_region_update();
            second_order_upwind_bound_update();
        }
    }
    else if (method.compare("Warming-Beam") == 0 || method.compare("6") == 0)
    {
        for (n = 0; n <= TLIM; ++n)
        {
            WB_region_update();
            WB_bound_update();
        }
    }
    else if (method.compare("Leap-Frog") == 0 || method.compare("7") == 0)
    {
        for (n = 0; n <= TLIM; ++n)
        {
            leap_frog_region_update();
            leap_frog_bound_update();
        }
    }
    else
        cout << "scheme not set." << endl;
}
void cal_p(string method)
{
    dt = 0.0001;
    dx = 0.001;
    c = a * dt / dx;
    TIME_ZONE = 0.01;
    SLIM = (R - L) / dx;
    TLIM = TIME_ZONE / dt;
    for (int i = 0; i < 5; ++i)
    {
        double p1 = 0, p2 = 0;
        bound_init();
        Iterative(method);
        for (int j = 0; j < SLIM; ++j)
        {
            p1 += pow((u[TLIM][j] - theoretical_sol(j * (R - L) / SLIM, TIME_ZONE)), 2) * dx;
            // if(i==0)cout<<j<<"   "<<p1<<endl;
        }
        p1 = log(sqrt(p1));

        dt *= 2;
        dx *= 2;
        SLIM = (R - L) / dx;
        TLIM = TIME_ZONE / dt;
        bound_init();
        Iterative(method);
        for (int j = 0; j < SLIM; ++j)
        {
            p2 += pow((u[TLIM][j] - theoretical_sol(j * (R - L) / SLIM, TIME_ZONE)), 2) * dx;
            // if(i==0)cout<<j<<"   "<<p2<<endl;
        }
        p2 = log(sqrt(p2));
        // cout << (p2 - p1) / log(2) << endl;
        p[2 * i] = p1 / log(2);
        p[2 * i + 1] = p2 / (log(2));
    }
    ofstream fout;
    fout.open("Error-Order-log.txt");
    for (int i = 0; i < 5; i++)
    {
        fout << p[2 * i] << "   " << p[2 * i + 1] << "   " << p[2 * i + 1] - p[2 * i] << endl;
    }
    fout.close();
}
int main()
{
    bound_init();
    Iterative("Warming-Beam");
    OutputTime(80000);   // 1/2000 t
    OutputSpace(1500); // 1/1000 x
    cal_p("Warming-Beam");
    return 0;
}