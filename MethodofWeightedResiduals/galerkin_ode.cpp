/*

Jia 2023.04

Solve specific ODE:L[u(x)]=f by galerkin method.
L[u(x)]=u''(x)+u;f(x)=-x;0<x<1,x(0)=0,x(1)=0.
In this program you can set the dimension of galerkin method from 1 to N_LIM.

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
const int N_LIM = 5000;
const int N = 8;
int grid_x = 1000; // for norm calculation
double count_as_zero_threshold=1e-12;
// -----numerical settings----- //

// -----equation settings----- //
double phi[N_LIM][N_LIM], lphi[N_LIM][N_LIM], fphi[N_LIM][N_LIM], lphij_phii[2 * N_LIM + 1], A[N_LIM][N_LIM], x[N_LIM], b[N_LIM];
// -----equation settings----- //

// -----generated settings----- //
 ofstream fout;
// -----generated settings----- //

// -----internals----- //
class Mat
{
public:
    int m = 1, n = 1;                 // 行数和列数
    double mat[N+5][N+5] = {0}; // 矩阵开始的元素

    Mat() {}
    Mat(int mm, int nn)
    {
        m = mm;
        n = nn;
    }
    void copy(double A[][N_LIM]); // 创建矩阵
    void Print();
    void augmat(Mat a, Mat b); // 求矩阵 a 和向量 b 的增广矩阵
    bool solve(Mat a, Mat b);  // 求线性方程组的解
};

void Mat::copy(double A[][N_LIM])
{
    for (int i = 1; i <= N; i++)
    {
        for (int j = 1; j <= N; j++)
        {
            mat[i][j] = A[i][j];
        }
    }
}
void Mat::augmat(Mat a, Mat b) // 求矩阵 a 和向量 b 的增广矩阵
{
    m = a.m;
    n = a.n + 1;
    for (int i = 1; i <= a.m; i++)
    {
        for (int j = 1; j <= a.n; j++)
        {
            mat[i][j] = a.mat[i][j];
        }
        mat[i][n] = b.mat[i][1];
    }
}
void Mat::Print()
{
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            cout << mat[i][j] << "\t";
        }
        cout << endl;
    }
}
bool Mat::solve(Mat a, Mat b) // a 为方阵 ，b 为列向量
// 求线性方程组的解(ax=b ,求 x)，矩阵 a 为方阵并且方程组有唯一解时返回 true
{
    if (a.n != a.m) // 只求解是方阵时的情形
    {
        cout << "系数矩阵不是方阵" << endl;
        return false; // 返回false
    }
    m = a.n;
    n = 1; // 解向量中必定有 a.n（ a.m ）个分量,是 a.n * 1 的列向量
    Mat aa;
    aa.augmat(a, b); // 求增广矩阵

    // 下面代码将增广矩阵化为上三角矩阵，并判断增广矩阵秩是否为 n
    for (int i = 1; i <= aa.m; i++)
    {
        // 寻找第 i 列不为零的元素
        int k;
        for (k = i; k <= aa.m; k++)
        {
            if (fabs(aa.mat[k][i]) > count_as_zero_threshold) // 满足这个条件时，认为这个元素不为0
                break;
        }
        if (k <= aa.m) // 说明第 i 列有不为0的元素
        {
            // 交换第 i 行和第 k 行所有元素
            for (int j = i; j <= aa.n; j++) // 从第 i 个元素交换即可，因为前面的元素都为0
            {                               // 使用aa.mat[0][j]作为中间变量交换元素
                aa.mat[0][j] = aa.mat[i][j];
                aa.mat[i][j] = aa.mat[k][j];
                aa.mat[k][j] = aa.mat[0][j];
            }
            double c; // 倍数
            for (int j = i + 1; j <= aa.m; j++)
            {
                c = -aa.mat[j][i] / aa.mat[i][i];
                for (k = i; k <= aa.n; k++)
                {
                    aa.mat[j][k] += c * aa.mat[i][k]; // 第 i 行 a 倍加到第 j 行
                }
            }
        }
        else // 没有找到则说明系数矩阵秩不为 n ，说明方程组中有效方程的个数小于 n
        {
            cout << "系数矩阵奇异，线性方程组无解或有无数解" << endl;
            return false;
        }
    }
    // 自下而上求解
    for (int i = a.m; i >= 1; i--)
    {
        mat[i][1] = aa.mat[i][aa.n];
        for (int j = a.m; j > i; j--)
        {
            mat[i][1] -= mat[j][1] * aa.mat[i][j];
        }
        mat[i][1] /= aa.mat[i][i];
    }
    return true;
}
// -----internals----- //

// -----manipulation----- //
void given_phi()
{
    for (int i = 1; i <= N; ++i)
    {
        phi[i][i] = 1;
        phi[i][i + 1] = -1;
    }
}
void get_lphi()
{
    for (int i = 1; i <= N; ++i)
    {
        lphi[i][i] = 1;
        lphi[i][i + 1] = -1;
        if ((i - 2) >= 0)
            lphi[i][i - 2] = i * (i - 1);
        if ((i - 1) >= 0)
            lphi[i][i - 1] = (i + 1) * i * (-1);
    }
}
void get_lphij_phii(int i, int j)
{
    for (int m = 0; m <= N + 1; ++m) // control lphij
    {
        if (lphi[j][m] != 0)
            for (int n = 0; n <= N + 1; ++n)
            {
                if (phi[i][n] != 0)
                    lphij_phii[m + n] += lphi[j][m] * phi[i][n];
            }
    }
}
void get_integral(int i, int j)
{
    for (int k = 0; k <= 2 * N + 2; ++k)
    {
        if (lphij_phii[k] != 0)
            A[i][j] += lphij_phii[k] * 1.0 / (k + 1);
    }
}
void print_all_polynomial()
{
    for (int i = 0; i <= 2 * N + 2; ++i)
    {
        cout << lphij_phii[i] << "  ";
    }
    cout << endl;
}
void whole_Aij()
{
    for (int i = 1; i <= N; ++i)
    {
        for (int j = 1; j <= N; ++j)
        {
            memset(lphij_phii, 0, sizeof(lphij_phii));
            get_lphij_phii(i, j);
            // print_all_polynomial();
            get_integral(i, j);
        }
    }
}
void print_equations()
{
    cout << "A:\n";
    cout << setprecision(12);
    for (int i = 1; i <= N; ++i)
    {
        for (int j = 1; j <= N; ++j)
        {
            cout << A[i][j] << "  ";
        }
        cout << endl;
    }
    cout << "b:\n";
    for (int i = 1; i <= N; ++i)
    {
        cout << b[i] << "  ";
    }
    cout << endl;
}
void get_fphi()
{
    for (int i = 1; i <= N; ++i)
    {
        fphi[i][i + 1] = -1;
        fphi[i][i + 2] = 1;
    }
}
void whole_b()
{
    for (int i = 1; i <= N; ++i)
    {
        b[i] += 1.0 / (i + 3) - 1.0 / (i + 2);
    }
}
void solve_Ax_equal_b()
{
    Mat mat_A(N, N), mat_b(N, 1);
    mat_A.copy(A);
    for (int i = 1; i <= N; ++i)
    {
        mat_b.mat[i][1] = b[i];
    }
    Mat mat_x;
    if(mat_x.solve(mat_A, mat_b))
    {
        for(int i=1;i<=N;++i)
        {
            x[i]=mat_x.mat[i][1];
        }
    }
}
void produce_u_n_x()
{
    fout.open("result_for_" + to_string(N) + ".txt");
    fout<<setprecision(12);
    fout << "general:" << endl;
    fout << "x(1-x)(";
    for (int i = 1; i <= N; ++i)
    {
        if (i == 1)
            fout << x[i];
        else if (i == 2)
            fout << x[i] << "*x";
        else
            fout << x[i] << "*x^" << i-1;
        if (i != N)
            fout << "+";
    }
    fout << ")" << endl
         << endl;
    fout << "matlab:" << endl;
    fout << "x=0:0.01:1;" << endl;
    fout << "y=";
    fout << "x.*(1-x).*(";
    for (int i = 1; i <= N; ++i)
    {
        if (i == 1)
            fout << x[i];
        else if (i == 2)
            fout << x[i] << "*x";
        else
            fout << x[i] << "*x.^" << i-1;
        if (i != N)
            fout << "+";
    }

    fout << ");" << endl;
    fout << "plot(x,y);" << endl;
    fout<<"number series:a_1--a_n"<<endl;
    for(int i=1;i<=N;++i)
    {
        fout<<x[i]<<endl;
    }
     cout<<endl<<"results:a_1--a_n"<<endl;
    for(int i=1;i<=N;++i)
    {
        cout<<x[i]<<endl;
    }
}
double u_n_x(double m)
{
    double temp = 0;
    for (int i = 1; i <= N; ++i)
    {
        temp += x[i] * pow(m, i - 1);
    }
    return m * (1 - m) * temp;
}
inline double theoretical_result(double x)
{
    return ((1.0 / sin(1.0)) * sin(x) - x);
}
void error_norm_calc()
{
    double m = 0, e1;
    for (int i = 0; i <= grid_x; ++i)
    {
        m = ((double)i) / grid_x;
        e1 += pow(u_n_x(m) - theoretical_result(m), 2);
        // cout<<e1<<endl;
    }
    e1 /= grid_x;
    e1 = sqrt(e1);
    cout <<endl<< "N=" << N << "     Error:" << e1 << endl;
    fout<<endl<< "N=" << N << "     Error:" << e1 << endl;
    fout.close();
    cout<<"complete result has been written to result_for_"<<N<<".txt"<<endl;
}

// -----manipulation----- //

int main()
{
    given_phi();
    get_lphi();
    whole_Aij();
    whole_b();
    print_equations();
    solve_Ax_equal_b();
    produce_u_n_x();
    error_norm_calc();
}