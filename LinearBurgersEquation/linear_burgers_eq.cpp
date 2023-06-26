/*
[Flashbac09] 2023/04/12 - 2023/04/15

Numercial Solution of 1-dimensional
convection-diffusion equation(linear Burgers equation)
with certain initial value
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
// -----equation settings----- //
const double a = 1;      // convection term,wave velocity
const double nju = 0.01; // diffusion term, 10^-6 for water at 300k or so.
inline double Burgers_LW_init_value(double x)
{
    return sin(2 * M_PI * x);
}
// -----equation settings----- //

// -----numerical settings----- //
const double L = 0.0, R = 3.0, T = 1;
int grid_x = 3000;
int grid_t = 100000;
// -----numerical settings----- //

// -----generated variables----- //
double u1[1000000], u2[1000000], u_s[1000000];
double u_test_std[10010], u_test_rgh[10010], u_test_key1[1000000], u_test_key2[1000000];
double dt = T / grid_t;
double dx = (R - L) / grid_x;
double c = a * dt / dx;
double d = (nju * dt) / (dx * dx);
double precision_order = 0, key_pre_ord = 0;
int grid_x_test = 300, grid_t_test = 100000, grid_x_ori = grid_x, grid_t_ori = grid_t;
int grid_x_key = 2499, grid_t_key = 41650;
// -----generated variables----- //

// -----internals----- //
int t_, x_;
int dataflow_flag;
// -----internals----- //

// -----manipulation----- //
inline void init_set()
{
    for (x_ = 0; x_ < grid_x; ++x_)
    {
        u1[x_] = Burgers_LW_init_value((double)(R - L) * x_ / grid_x);
    }
}

inline void Burgers_LW_region()
{
    if (dataflow_flag)
        for (x_ = 1; x_ < grid_x - 1; ++x_)
        {
            u2[x_] = u1[x_] - 0.5 * c * (u1[x_ + 1] - u1[x_ - 1]) + (d + 0.5 * c * c) * (u1[x_ + 1] - 2 * u1[x_] + u1[x_ - 1]);
        }
    else
    {
        for (x_ = 1; x_ < grid_x - 1; ++x_)

        {
            u1[x_] = u2[x_] - 0.5 * c * (u2[x_ + 1] - u2[x_ - 1]) + (d + 0.5 * c * c) * (u2[x_ + 1] - 2 * u2[x_] + u2[x_ - 1]);
        }
    }
}
inline void Burgers_LW_boundary()
{
    if (dataflow_flag)
    {
        u2[0] = u1[0] - 0.5 * c * (u1[1] - u1[grid_x - 1]) + (d + 0.5 * c * c) * (u1[1] - 2 * u1[0] + u1[grid_x - 1]);
        u2[grid_x - 1] = u1[grid_x - 1] - 0.5 * c * (u1[0] - u1[grid_x - 2]) + (d + 0.5 * c * c) * (u1[0] - 2 * u1[grid_x - 1] + u1[grid_x - 2]);
    }
    else
    {
        u1[0] = u2[0] - 0.5 * c * (u2[1] - u2[grid_x - 1]) + (d + 0.5 * c * c) * (u2[1] - 2 * u2[0] + u2[grid_x - 1]);
        u1[grid_x - 1] = u2[grid_x - 1] - 0.5 * c * (u2[0] - u2[grid_x - 2]) + (d + 0.5 * c * c) * (u2[0] - 2 * u2[grid_x - 1] + u2[grid_x - 2]);
    }
}
void stability_check(string methods)
{
    if (methods == "Burgers_LW")
    {
        string Burgers_LW_condition = "(2*d+c^2) <= 1";
        double c = a * dt / dx;
        double d = nju * dt / (dx * dx);
        cout << "Settings:" << endl;
        cout << "Convection term: a=" << a << "     Diffusion term: nju=" << nju << endl;
        cout << "dt=" << dt << "     "
             << "dx=" << dx << endl;
        cout << "Checking Stability:     " << Burgers_LW_condition << endl;
        cout << "c=" << c << "     "
             << "d=" << d << endl;
        cout << "(2*d+c^2) = " << 2 * d + c * c << endl;
        if ((6 * d + c * c) == 1)
            cout << "Significantly Stable." << endl;
        else if (2 * d + c * c < 1 && 2 * d + c * c > 0)
            cout << "Stable." << endl;
        else if (2 * d + c * c == 1)
            cout << "Marginally Stable." << endl;
        else
            cout << "Unstable." << endl;
    }
}
void save_high_accuracy()
{
    if (!dataflow_flag)
    {
        memcpy(u_s, u1, sizeof(u1));
        cout << "Result Saved." << endl;
    }
    else
    {
        memcpy(u_s, u2, sizeof(u2));
        cout << "Result Saved." << endl;
    }
}
void compare_std()
{
    grid_x = grid_x_test;
    grid_t = grid_t_test;
    dt = T / grid_t;
    dx = (R - L) / grid_x;
    c = a * dt / dx;
    d = (nju * dt) / (dx * dx);
    memset(u1, 0, sizeof(u1));
    memset(u2, 0, sizeof(u2));
    init_set();
    for (t_ = 1; t_ <= grid_t; ++t_) // high-accuracy solution
    {
        dataflow_flag = t_ % 2; // [u1 - u2 circulation]
        Burgers_LW_region();
        Burgers_LW_boundary();
        if (t_ % 100000 == 0)
            cout << t_ << "/" << grid_t << endl;
    }
    for (x_ = 0; x_ < grid_x; ++x_)
    {
        u_test_std[x_] = u1[x_];
    }
}
void compare_rougher()
{
    grid_x = grid_x_test / 2;
    grid_t = grid_t_test / 2;
    dt = T / grid_t;
    dx = (R - L) / grid_x;
    c = a * dt / dx;
    d = (nju * dt) / (dx * dx);
    memset(u1, 0, sizeof(u1));
    memset(u2, 0, sizeof(u2));
    init_set();
    for (t_ = 1; t_ <= grid_t; ++t_) // high-accuracy solution
    {
        dataflow_flag = t_ % 2; // [u1 - u2 circulation]
        Burgers_LW_region();
        Burgers_LW_boundary();
        if (t_ % 100000 == 0)
            cout << t_ << "/" << grid_t << endl;
    }
    for (x_ = 0; x_ < grid_x; ++x_)
    {
        u_test_rgh[x_] = u1[x_];
    }
}
void calculate_order()
{
    grid_x = grid_x_test;
    grid_t = grid_t_test;
    dt = T / grid_t;
    dx = (R - L) / grid_x;
    double norm1 = 0, norm2 = 0;
    for (x_ = 0; x_ < grid_x; ++x_)
    {
        norm1 += (u_test_std[x_] - u_s[x_ * grid_x_ori / grid_x]) * (u_test_std[x_] - u_s[x_ * grid_x_ori / grid_x]);
    }
    norm1 *= dx;
    norm1 = sqrt(norm1);
    grid_x = grid_x_test / 2;
    grid_t = grid_t_test / 2;
    dt = T / grid_t;
    dx = (R - L) / grid_x;
    for (x_ = 0; x_ < grid_x; ++x_)
    {
        norm2 += (u_test_rgh[x_] - u_s[x_ * grid_x_ori / grid_x]) * (u_test_rgh[x_] - u_s[x_ * grid_x_ori / grid_x]);
    }
    norm2 *= dx;
    norm2 = sqrt(norm2);
    precision_order = (1.0 / log(2.0)) * (log(norm2) - log(norm1));
}
void verify_order_of_accuracy()
{
    cout << "Verify Order of Precision:" << endl;
    compare_std();
    compare_rougher();
    calculate_order();
    cout << "Order of Precision:" << precision_order << endl;
}
void three_order_term_elimination()
{
    double norm_key1 = 0, norm_key2 = 0;
    // 1/2
    grid_x = grid_x_key;
    grid_t = grid_t_key;
    dt = T / grid_t;
    dx = (R - L) / grid_x;
    c = a * dt / dx;
    d = (nju * dt) / (dx * dx);
    memset(u1, 0, sizeof(u1));
    memset(u2, 0, sizeof(u2));
    init_set();
    for (t_ = 1; t_ <= grid_t; ++t_) // high-accuracy solution
    {
        dataflow_flag = t_ % 2; // [u1 - u2 circulation]
        Burgers_LW_region();
        Burgers_LW_boundary();
        if (t_ % 100000 == 0)
            cout << t_ << "/" << grid_t << endl;
    }
    for (x_ = 0; x_ < grid_x; ++x_)
    {
        u_test_key1[x_] = u1[x_];
    }
    for (x_ = 0; x_ < grid_x; ++x_)
    {
        norm_key1 += (u_test_key1[x_] - u_s[x_ * grid_x_ori / grid_x]) * (u_test_key1[x_] - u_s[x_ * grid_x_ori / grid_x]);
    }
    norm_key1 *= dx;
    norm_key1 = sqrt(norm_key1);
    // 2/2
    grid_x = grid_x_key / 7;
    grid_t = grid_t_key / 7;
    dt = T / grid_t;
    dx = (R - L) / grid_x;
    c = a * dt / dx;
    d = (nju * dt) / (dx * dx);
    memset(u1, 0, sizeof(u1));
    memset(u2, 0, sizeof(u2));
    init_set();
    for (t_ = 1; t_ <= grid_t; ++t_) // high-accuracy solution
    {
        dataflow_flag = t_ % 2; // [u1 - u2 circulation]
        Burgers_LW_region();
        Burgers_LW_boundary();
        if (t_ % 100000 == 0)
            cout << t_ << "/" << grid_t << endl;
    }
    for (x_ = 0; x_ < grid_x; ++x_)
    {
        u_test_key2[x_] = u1[x_];
    }
    for (x_ = 0; x_ < grid_x; ++x_)
    {
        norm_key2 += (u_test_key2[x_] - u_s[x_ * grid_x_ori / grid_x]) * (u_test_key2[x_] - u_s[x_ * grid_x_ori / grid_x]);
    }
    norm_key2 *= dx;
    norm_key2 = sqrt(norm_key2);
    key_pre_ord = (1.0 / log(7.0)) * (log(norm_key2) - log(norm_key1));
    if(grid_x_ori%grid_x==0)
    cout << "Significantly Higher Precision Order:"<<key_pre_ord<<endl;
}
void Output_high_accuracy()
{
    ofstream fout;
    string ss = "a" + to_string(a) + "nju" + to_string(nju) + "x" + to_string(grid_x) + "t" + to_string(grid_t) + ".txt";
    fout.open(ss);
    for (x_ = 0; x_ < grid_x; ++x_)
    {
        fout << ((R - L) * x_) / ((double)grid_x) << "   " << u_s[x_] << endl;
    }

    fout.close();
}
// -----manipulation----- //
int main()
{
    stability_check("Burgers_LW");
    init_set();
    for (t_ = 1; t_ <= grid_t; ++t_) // high-accuracy solution
    {
        dataflow_flag = t_ % 2; // [u1 - u2 circulation]
        Burgers_LW_region();
        Burgers_LW_boundary();
        if (t_ % 100000 == 0)
            cout << t_ << "/" << grid_t << endl;
    }
    save_high_accuracy();
    Output_high_accuracy();
    verify_order_of_accuracy();
    three_order_term_elimination();
    return 0;
}