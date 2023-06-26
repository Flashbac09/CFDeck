/*
[Flashbac09]

// 2023/03/19 differential_ftest
// This program will export txt files into folder "test_float" as it runs.
// Sorry to mention that the program isn't responsible for creating the folder if it does not exist.
*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <math.h>

using namespace std;

// These variables can be adjusted on demand.
int n_init = 8;
// int LIMIT=0xFFFFFFFF;
int LIMIT = 200000;
float L = -3.2, R = 3.2; // interval

float x = 0, xf = 0, xb = 0, xff = 0; // variable
float r1[100], r2[100], r3[100], r4[100], r5[100], r6[100], r7[100];
float std1, std2, std3, std4;

inline void theoretical_set()
{
    std1 = 2 * x + 1;                                               // f1求一阶导,f1_d1
    std2 = 0;                                                       // f2_d2
    std3 = 4.6 * cos(4.6 * x) - 0.16 * sin(0.16 * x);               // f3_d1
    std4 = -4.6 * 4.6 * sin(4.6 * x) - 0.16 * 0.16 * cos(0.16 * x); // f3_d2
}
inline float f3(float x)
{
    return sin(4.6 * x) + cos(0.16 * x);
}
void FileOut(int cycle)
{
    ofstream fout;
    // 1
    fout.open("test_float/r1.txt");
    int n = n_init;
    for (int i = 1; i <= cycle; ++i)
    {
        fout << log10(n) << "  " << r1[i] << endl;
        n *= 2;
    }
    fout.close();
    // 2
    n = n_init;
    fout.open("test_float/r2.txt");
    for (int i = 1; i <= cycle; ++i)
    {
        fout << log10(n) << "  " << r2[i] << endl;
        n *= 2;
    }
    fout.close();
    // 3
    n = n_init;
    fout.open("test_float/r3.txt");
    for (int i = 1; i <= cycle; ++i)
    {
        fout << log10(n) << "  " << r3[i] << endl;
        n *= 2;
    }
    fout.close();
    // 4
    n = n_init;
    fout.open("test_float/r4.txt");
    for (int i = 1; i <= cycle; ++i)
    {
        fout << log10(n) << "  " << r4[i] << endl;
        n *= 2;
    }
    fout.close();
    // 5
    n = n_init;
    fout.open("test_float/r5.txt");
    for (int i = 1; i <= cycle; ++i)
    {
        fout << log10(n) << "  " << r5[i] << endl;
        n *= 2;
    }
    fout.close();
    // 6
    n = n_init;
    fout.open("test_float/r6.txt");
    for (int i = 1; i <= cycle; ++i)
    {
        fout << log10(n) << "  " << r6[i] << endl;
        n *= 2;
    }
    fout.close();
    // 7
    n = n_init;
    fout.open("test_float/r7.txt");
    for (int i = 1; i <= cycle; ++i)
    {
        fout << log10(n) << "  " << r7[i] << endl;
        n *= 2;
    }
    fout.close();
}
int main()
{
    unsigned int n = n_init, cycle = 0;
    cout.setf(ios::fixed);
    cout << setprecision(40);
    while (n < LIMIT)
    {
        /*
        1. 多项式，一阶中心。f1(x)=x^2+x,这保证了中心差分求1阶导数无截断误差。
        2. 多项式，二阶中心。f2(x)=x+10,这保证了中心差分求2阶导数无截断误差。
        3. 三角，一阶向前。f3(x)=sin(4.7*x)+cos(0.16*x),有截断误差和舍入误差，下同。
        4. 三角，一阶中心。
        5. 三角，一阶三点向前
        6. 三角，二阶向前
        7. 三角，二阶中心
        */
        cycle += 1;
        float h = (R - L) / n;
        for (int i = 0; i <= n; ++i)
        {
            x = L + i * h;
            xf = x + h;
            xb = x - h;
            xff = x + 2 * h;
            theoretical_set();
            r1[cycle] += fabs((((xf * xf + xf) - (xb * xb + xb)) / (2 * h)) - std1);
            r2[cycle] += fabs((xf - 2 * x + xb) / (h * h) - std2);
            r3[cycle] += fabs((f3(xf) - f3(x)) / (h)-std3);
            r4[cycle] += fabs((f3(xf) - f3(xb)) / (2 * h) - std3);
            r5[cycle] += fabs((-3 * f3(x) + 4 * f3(xf) - f3(xff)) / (2 * h) - std3);
            r6[cycle] += fabs((f3(x) - 2 * f3(xf) + f3(xff)) / (h * h) - std4);
            r7[cycle] += fabs((f3(xf) - 2 * f3(x) + f3(xb)) / (h * h) - std4);
        }
        r1[cycle] = log10(r1[cycle] / n);
        r2[cycle] = log10(r2[cycle] / n);
        r3[cycle] = log10(r3[cycle] / n);
        r4[cycle] = log10(r4[cycle] / n);
        r5[cycle] = log10(r5[cycle] / n);
        r6[cycle] = log10(r6[cycle] / n);
        r7[cycle] = log10(r7[cycle] / n);
        n *= 2;
    }
    FileOut(cycle);
}