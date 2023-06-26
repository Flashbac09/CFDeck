/*

[Jia/Flashbac09 2023.05]
Cylinder Flow Calculation by Finite Method 

input: GRID.DAT
output: stream.txt, velocity.txt, pressure.txt

Special Case:
Boundary Point are manually distinguished.
Essential:1-42,51-86,87-105;
Natural:42-50;

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
using namespace std;

//----- numerical settings -----//
const double eps = 1e-12;
int n = 765; // node
int N_node = 765;
int N_unit = 1423;
double X[765], Y[765];
double unit[1423][3]; // unit
double vx[765], vy[765];
double vx_u[1423], vy_u[1423];
double p[765];
double adjacent[765];
//----- numerical settings -----//

//----- variables ----//
const int N = 1000;
double a[N][N];
double Bound[N];
double f[N];
//----- variables ----//

int gauss()
{
    int c, r; // column, row
    for (c = 0, r = 0; c < n; c++)
    {
        int t = r;                  // 从对角线元素开始往下遍历这一列
        for (int i = t; i < n; i++) // 找到该列绝对值最大的元素所在的行号
            if (fabs(a[i][c]) > fabs(a[t][c]))
                t = i;

        if (fabs(a[t][c]) < eps)
            continue; // 最大绝对值是0，那么这一列剩下的元素全是0，不用管这列。这个地方导致后面不能r ++，也意味着底部将会增加一个系数全0的方程

        for (int i = c; i <= n; i++)
            swap(a[t][i], a[r][i]); // 把最大绝对值元素所在的行换到未处理行的最上面(即当前要处理的的第r行)
        for (int i = n; i >= c; i--)
            a[r][i] /= a[r][c];         // 把现在的第r行的数字全部除以一个系数，使得左上角a[r][c]变成1
        for (int i = r + 1; i < n; i++) // 把当前列下的所有数都消成0,要对应两行元素一起变化
            if (fabs(a[i][c]) > eps)    // 已经是0的就不用操作了，省点计算
                for (int j = n; j >= c; j--)
                    a[i][j] -= a[r][j] * a[i][c];

        r++;
    }

    // 上面步骤走完之后，矩阵a[][]扣掉增广的最后一列系数以外，剩下的已经是个上三角阵或者阶梯阵
    if (r < n) // 说明有效的方程个数小于n,那要么无穷解，要么无解
    {
        for (int i = r; i < n; i++)
        {
            if (fabs(a[i][n]) > eps) // a[i][n] = b_i不等于0
                return 2;            // 无解
            return 1;                // 都是0 = 0的方程，无穷解
        }
    }

    // 唯一解，从下往上回代，得到方程的解
    for (int i = n - 1; i >= 0; i--)
        for (int j = i + 1; j < n; j++)
            a[i][n] -= a[i][j] * a[j][n];
    /*
    //output
    for (int i = 0; i < n; i++)
        printf("%.2f\n", a[i][n]);
    */
    return 0;
}
void simple_gauss(int n) // success by default, no warning, no epsilon setting.Ax=f, results stored in f[i]
{
    double c = 0;
    // pt1
    for (int k = 0; k < n - 1; ++k)
    {
        for (int i = k + 1; i < n; ++i)
        {
            c = (-1)*a[i][k] / a[k][k];
            for (int j = k + 1; j < n; ++j)
            {
                a[i][j] += c * a[k][j];
            }
            f[i] += c * f[k];
        }
    }
    // pt2
    for (int i = n - 1; i >= 0; i--)
    {
        for (int j = n - 1; j >= i + 1; j--)
        {
            f[i] -= a[i][j] * f[j];
        }
        f[i] /= a[i][i];
    }
}
double triangle_area(int i, int j, int k)
{
    return fabs(0.5 * ((X[j] - X[i]) * (Y[k] - Y[i]) - (Y[j] - Y[i]) * (X[k] - X[i])));
}
void Input()
{
    ifstream fin("GRID.DAT");
    fin >> N_node >> N_unit;
    for (int i = 0; i < N_node; ++i)
    {
        fin >> X[i] >> Y[i];
    }
    for (int i = 0; i < N_unit; ++i)
    {
        fin >> unit[i][0] >> unit[i][1] >> unit[i][2];
    }
    fin.close();
}
void generate_unit()
{
    double tri_area;
    double b[4], c[4]; // use 1,2,3 as book demostrates
    int tmp_unit[4];   // use 1,2,3
    for (int k = 0; k < N_unit; ++k)
    {
        tmp_unit[1] = unit[k][0] - 1;
        tmp_unit[2] = unit[k][1] - 1;
        tmp_unit[3] = unit[k][2] - 1;
        tri_area = triangle_area(tmp_unit[1], tmp_unit[2], tmp_unit[3]);
        b[1] = Y[tmp_unit[2]] - Y[tmp_unit[3]];
        b[2] = Y[tmp_unit[3]] - Y[tmp_unit[1]];
        b[3] = Y[tmp_unit[1]] - Y[tmp_unit[2]];
        c[1] = X[tmp_unit[3]] - X[tmp_unit[2]];
        c[2] = X[tmp_unit[1]] - X[tmp_unit[3]];
        c[3] = X[tmp_unit[2]] - X[tmp_unit[1]];
        for (int i = 1; i <= 3; ++i)
        {
            adjacent[tmp_unit[i]] += tri_area;
        }

        for (int i = 1; i <= 3; ++i)
        {
            for (int j = 1; j <= 3; ++j)
            {
                a[tmp_unit[i]][tmp_unit[j]] += (b[i] * b[j] + c[i] * c[j]) / (4 * tri_area);
            }
        }
    }
    return;
}
void fix_boundary()
{
    for (int i = 0; i < 41; ++i)
    {
        Bound[i] = 0;
    }
    for (int i = 50; i < 86; ++i)
    {
        Bound[i] = 2.0;
    }
    for (int i = 86; i < 105; ++i)
    {
        Bound[i] = 2.0 - 0.1 * (i - 85);
    }
}
void generate_overall()
{
    for (int i = 0; i < 41; ++i)
    {
        f[i] = Bound[i];

        for (int k = 0; k < N_node; ++k)
        {
            if (k == i)
                a[i][k] = 1.0;
            else
                a[i][k] = 0;
        }
    }
    for (int i = 41; i < 50; ++i)
    {
        for (int j = 0; j < 105; ++j)
        {
            f[i] -= a[i][j] * Bound[j];
        }
        for (int k = 0; k < 41; ++k)
        {
            a[i][k] = 0;
        }
        for (int m = 50; m < 105; ++m)
        {
            a[i][m] = 0;
        }
    }
    for (int i = 50; i < 105; ++i)
    {
        f[i] = Bound[i];
        for (int j = 0; j < N_node; ++j)
        {
            if (i == j)
                a[i][j] = 1;
            else
                a[i][j] = 0;
        }
    }
    for (int i = 105; i < N_node; ++i)
    {
        for (int j = 0; j < 105; ++j)
        {
            f[i] -= a[i][j] * Bound[j];
        }
        for (int j = 0; j < 41; ++j)
        {
            a[i][j] = 0;
        }
        for (int j = 50; j < 105; ++j)
        {
            a[i][j] = 0;
        }
    }
}
void completion()
{
    // pt1:velocity unit/node_ave
    int tmp_unit[4]; // use 1,2,3
    double b[4], c[4];
    double tri_area;
    for (int k = 0; k < N_unit; ++k)
    {
        tmp_unit[1] = unit[k][0] - 1;
        tmp_unit[2] = unit[k][1] - 1;
        tmp_unit[3] = unit[k][2] - 1;
        tri_area = triangle_area(tmp_unit[1], tmp_unit[2], tmp_unit[3]);
        b[1] = 0.5 * (Y[tmp_unit[2]] - Y[tmp_unit[3]]) / tri_area;
        b[2] = 0.5 * (Y[tmp_unit[3]] - Y[tmp_unit[1]]) / tri_area;
        b[3] = 0.5 * (Y[tmp_unit[1]] - Y[tmp_unit[2]]) / tri_area;
        c[1] = 0.5 * (X[tmp_unit[3]] - X[tmp_unit[2]]) / tri_area;
        c[2] = 0.5 * (X[tmp_unit[1]] - X[tmp_unit[3]]) / tri_area;
        c[3] = 0.5 * (X[tmp_unit[2]] - X[tmp_unit[1]]) / tri_area;
        for (int i = 1; i <= 3; ++i)
        {
            vx_u[k] += c[i] * f[tmp_unit[i]];
            vy_u[k] -= b[i] * f[tmp_unit[i]];
        }
        for (int j = 1; j <= 3; ++j)
        {
            vx[tmp_unit[j]] += tri_area * vx_u[k];
            vy[tmp_unit[j]] += tri_area * vy_u[k];
        }
    }
    for (int i = 0; i < N_node; ++i)
    {
        vx[i] /= adjacent[i];
        vy[i] /= adjacent[i];
    }
    // pt2:pressure
    for (int i = 0; i < N_node; ++i)
    {
        p[i] = 0.5 * (1.0 - vx[i] * vx[i] - vy[i] * vy[i]);
    }
}
void Output()
{
    ofstream fout("stream.txt");
    for (int i = 0; i < N_node; ++i)
    {
        fout << X[i] << "  " << Y[i] << "  " << f[i] << endl;
    }
    fout.close();
    fout.open("velocity.txt");
    for (int i = 0; i < N_node; ++i)
    {
        fout << vx[i] << "  " << vy[i] << endl;
    }
    fout.close();
    fout.open("pressure.txt");
    for (int i = 0; i < N_node; ++i)
    {
        fout <<p[i] << endl;
    }
    fout.close();
    fout.open("units.txt");
    for(int i=0;i<N_unit;++i)
    {
        fout<<unit[i][0]<<"  "<<unit[i][1]<<"  "<<unit[i][2]<<endl;
    }
    fout.close();
}
int main()
{
    // node: 105 boundary 660 region 765 total
    // unit: 1423 total

    Input(); // nodes b&r(660*xy+105*xy);unit 1423*(ijk)
    generate_unit();
    fix_boundary();
    generate_overall();
    simple_gauss(N_node);
    completion(); // velocity/pressure
    Output();

    return 0;
}

/*

    cin >> n;

    for (int i = 0; i < n; i++)
        for (int j = 0; j <= n; j++)
            cin >> a[i][j];

    int t = gauss();

    if (t == 0)
        for (int i = 0; i < n; i++)
            printf("%.2f\n", a[i][n]);
    else if (t == 1)
        puts("Infinite group solutions");
    else
        puts("No solution");

*/