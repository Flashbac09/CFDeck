/*
AirfoilMesh
    from specific complex transformation

[Flashbac09] [2023.03]

*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
const int total = 50;
const double h = 1.0 / total;
double L = -0.2, R = 0.8, var = 0;
double list1x[total + 1], list1y[total + 1], list1r[total + 1], list1t[total + 1];
double list2x[total + 1], list2y[total + 1], list2r[total + 1], list2t[total + 1];
double gridwy[total / 2 + 1][2 * (total + 1)], gridwx[total / 2 + 1][2 * (total + 1)];
double f1(double x)
{
    x -= L;
    return 0.6 * (0.2969 * sqrt(x) - 0.126 * x - 0.3516 * x * x + 0.2843 * x * x * x - 0.1015 * x * x * x * x);
}
void output_list()
{
    std::ofstream fout;
    fout.open("line1.txt");
    for (int i = 0; i < total; i++)
    {
        fout << list2x[i] << "   " << list2y[i] << std::endl;
    }

    fout.close();
    return;
}

void output_grid()
{
    std::ofstream fout;
    fout.open("grid1-row.txt");
    for (int i = 0; i <= total / 2; i++)
    {
        for (int j = 0; j <= 2 * total; j++)
        {
            fout << gridwx[i][j] << "   " << gridwy[i][j] << std::endl;
        }
    }
    fout.close();

    fout.open("grid1-col.txt");
    for (int j = 0; j <= 2 * total; j++)
    {
        for (int i = 0; i <= total / 2; i++)
        {
            fout << gridwx[i][j] << "   " << gridwy[i][j] << std::endl;
        }
    }

    fout.close();
    return;
}

void line_transfer_to_w()
{
    for (int i = 0; i < total; i++)
    {
        var = ((double)i / total) + L;
        list1x[i] = var;
        list1y[i] = f1(var);
        list1r[i] = sqrt(list1x[i] * list1x[i] + list1y[i] * list1y[i]);
        list2r[i] = sqrt(list1r[i]);
        list1t[i] = atan2(list1y[i], list1x[i]);
        list2t[i] = list1t[i] * 0.5;
        list2x[i] = list2r[i] * cos(list2t[i]);
        list2y[i] = list2r[i] * sin(list2t[i]);
    }
}

void generate_grid_in_w()
{
    for (int i = 0; i <= total / 2; i++)
    {
        for (int j = 0; j <= 2 * total; j++)
        {
            gridwx[i][j] = ((double)j / total) - 1;
        }
    }
    for (int j = 0; j <= 2 * total; j++)
    {
        for (int i = 0; i <= total / 2; i++)
        {
            gridwy[i][j] = ((double)i / total);
        }
    }
}

void transfer_grid_to_z()
{
    double r = 0, t = 0;
    for (int i = 0; i <= total / 2; i++)
    {
        for (int j = 0; j <= 2 * total; j++)
        {
            r = (gridwx[i][j] * gridwx[i][j] + gridwy[i][j] * gridwy[i][j]);
            t = 2 * atan2(gridwy[i][j], gridwx[i][j]);
            gridwx[i][j] = r * cos(t);
            gridwy[i][j] = r * sin(t);
        }
    }
}
int main()
{
    line_transfer_to_w();
    generate_grid_in_w();
    transfer_grid_to_z();
    output_grid();
    output_list();
    return 0;
}