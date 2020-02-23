//
//  main.cpp
//  HW3-1, explicit finite difference scheme
//
//  Created by JIE QIAN on 2/20/20.
//  Copyright © 2020 JIE QIAN. All rights reserved.
//

//  This code is to implement the explicit finite difference scheme
//  for the transport equation

#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>

using namespace std;

vector<vector<double>> explicit_finite_diff(double x_min, double x_max, double dx,
                                            double t_max, double dt,
                                            double c)
{
    int X = ceil((x_max - x_min) / dx);
    int T = ceil(t_max / dt);
    double p = dt / dx;
    
    vector<vector<double>> mesh;
    vector<double> tmp(X+1);
    mesh.resize(T+1, tmp);
    
    for(int x=0; x<X; x++)
    {
        if( (double) x*dx+x_min < -1)
            mesh[0][x] = 0;
        else if( (double) x*dx+x_min < 0)
            mesh[0][x] = (double) x*dx+x_min + 1.0;
        else
            mesh[0][x] = 1.0;
    }
    
    for(int t=0; t<T; t++)
        for(int x=0; x<X; x++)
        {
            mesh[t+1][x+1] = (1-c*p)*mesh[t][x+1] + c*p*mesh[t][x];
        }
    
    return mesh;
};

int main(int argc, const char * argv[]) {
    ofstream outf1;
    outf1.open("mesh1.csv", ios::out | ios::trunc);
    vector<vector<double>> mesh1 = explicit_finite_diff(-2.0, 3.0, 0.05, 1, 0.01, 1);
    for(int i=0; i<mesh1.size(); i++)
    {
           for(int j=0; j<mesh1[0].size(); j++)
           {
               outf1<<mesh1[i][j]<<",";
           }
        outf1<<endl;
    }
    outf1.close();
    
    ofstream outf2;
    outf2.open("mesh2.csv", ios::out | ios::trunc);
    vector<vector<double>> mesh2 = explicit_finite_diff(-2.0, 3.0, 0.01, 1, 0.01, 1);
    for(int i=0; i<mesh2.size(); i++)
    {
           for(int j=0; j<mesh2[0].size(); j++)
           {
               outf2<<mesh2[i][j]<<",";
           }
        outf2<<endl;
    }
    outf2.close();
    
    ofstream outf3;
    outf3.open("mesh3.csv", ios::out | ios::trunc);
    vector<vector<double>> mesh3 = explicit_finite_diff(-2.0, 3.0, 0.005, 1, 0.01, 1);
    for(int i=0; i<mesh3.size(); i++)
    {
           for(int j=0; j<mesh3[0].size(); j++)
           {
               outf3<<mesh3[i][j]<<",";
           }
        outf3<<endl;
    }
    outf3.close();
    
    return 0;
}
