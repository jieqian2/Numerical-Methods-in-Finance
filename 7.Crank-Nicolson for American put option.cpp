//
//  main.cpp
//  FD American put
//
//  Created by JIE QIAN on 5/2/20.
//  Copyright © 2020 JIE QIAN. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "newmat.h"
#include "newmatap.h"
#include "newmatio.h"

using namespace std;

double pricing_CN_FD(const double& initial_value,
                           const double& r, //risk free rate
                           const double& sigma, //volatility
                           const double& T, //maturity time
                           const double& K,
                           const double& Smax, //P(S>Smax) near 0
                           double ds,
                           double dt)
{
    Tracer tr("name");
    double sigma_sqr = sigma*sigma;
    int M = round(Smax/ds);
    int IN = round(initial_value/ds);
    int N = round(T/dt);

    //construct mesh;
    BandMatrix C(M+1, 1, 1);
    C = 0.0;
    C.element(0, 0) = 1+0.5*dt*(sigma_sqr + r);
    C.element(0, 1) = -0.25*dt*(sigma_sqr + r);
    C.element(M, M-1) = -0.25*dt*((double) sigma_sqr*(M+1)*(M+1)-(double) r*(M+1));
    C.element(M, M) = 1+0.5*dt*((double) sigma_sqr*(M+1)*(M+1)+r);
    
    BandMatrix D(M+1, 1, 1);
    D = 0.0;
    D.element(0, 0) = 1-0.5*dt*(sigma_sqr + r);
    D.element(0, 1) = 0.25*dt*(sigma_sqr + r);
    D.element(M, M-1) = 0.25*dt*((double) sigma_sqr*(M+1)*(M+1)-(double) r*(M+1));
    D.element(M, M) = 1-0.5*dt*((double) sigma_sqr*(M+1)*(M+1)+r);
    
    for(int j=1; j<M; j++)
    {
        double aj = 0.25*dt*((double) sigma_sqr*(j+1)*(j+1)-(double) r*(j+1));
        double bj = -0.5*dt*((double) sigma_sqr*(j+1)*(j+1)+r);
        double cj = 0.25*dt*((double) sigma_sqr*(j+1)*(j+1)+(double) r*(j+1));
        C.element(j, j-1) = -aj;
        C.element(j, j) = 1-bj;
        C.element(j, j+1) = -cj;
        D.element(j, j-1) = aj;
        D.element(j, j) = 1+bj;
        D.element(j, j+1) = cj;
    }
            
    ColumnVector VT(M+1);
    VT = 0.0;
    for(int i=0; i<M; i++){
        VT.element(i) = max(0.0, K - (double) i*ds);
    }
    
    // backward to the initial date
    for(int t=N; t>=0; t--){
        VT = C.i() * D * VT;
        VT.element(0) = 0.0;
        
        for(int i=0; i<M+1; i++)
            VT.element(i) = max(VT.element(i), K - (double) i*ds);
    }
    
    return (VT.element(IN)+VT.element(IN+1))/2;
};


double up_factor, uptick_prob, notick_prob, downtick_prob, risk_free_rate, strike_price;
double initial_stock_price, expiration_time, volatility, R;
int no_of_divisions;

vector<vector<double>> putstore;

double american_put_option(int k, int i, double current_stock_price)
{
    if(putstore[k][no_of_divisions+i] != -1.0)
        return putstore[k][no_of_divisions+i];
    
    if (k == no_of_divisions)
    {
        putstore[k][no_of_divisions+i]=max(0.0, (strike_price - current_stock_price));
        return putstore[k][no_of_divisions+i];
    }
    else
    {
        putstore[k][no_of_divisions+i] = max((strike_price - current_stock_price),
        (uptick_prob*american_put_option(k+1, i+1, current_stock_price*up_factor) + notick_prob*american_put_option(k+1, i, current_stock_price) +
         downtick_prob*american_put_option(k+1, i-1, current_stock_price/up_factor))/R);
        return putstore[k][no_of_divisions+i];
    }
}

int main(int argc, const char * argv[]) {
    initial_stock_price = 50.0;
    no_of_divisions = 1000;
    risk_free_rate = 0.01;
    volatility = 0.4;
    expiration_time = 1.0;
    strike_price = 55.0;
    
    cout<<"American put option pricing"<<endl;
    cout << "Expiration Time (Years) = " << expiration_time << endl;
    cout << "Number of Divisions = " << no_of_divisions << endl;
    cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike Price = " << strike_price << endl;
    cout<<endl;
    cout<<"---As a benchmark (using trinomial model)---"<<endl;
    {
        vector<double> tmp;
        tmp.resize(2*no_of_divisions+2, -1.0);
        putstore.resize(2*no_of_divisions+2, tmp);
        
        up_factor = exp(volatility*sqrt((2*expiration_time)/((float) no_of_divisions)));
        R = exp(risk_free_rate*expiration_time/((float) no_of_divisions));
        uptick_prob = pow((sqrt(R) - 1/sqrt(up_factor))/(sqrt(up_factor)-1/sqrt(up_factor)),2);
        downtick_prob = pow((sqrt(up_factor) - sqrt(R))/(sqrt(up_factor)-1/sqrt(up_factor)),2);
        notick_prob = 1 - uptick_prob - downtick_prob;
        
        double put_price = american_put_option(0, 0, initial_stock_price);
        cout<<"option price = "<<put_price<<endl;
    }
    
    cout<<endl;
    cout<<"---Using C-N finite difference method---"<<endl;
    {
        double Smax = 5.0 * initial_stock_price;
        double dS = 0.01 * initial_stock_price;
        double dt = 1.0/365.0;
        double put_price_cn = pricing_CN_FD(initial_stock_price, risk_free_rate, volatility, expiration_time, strike_price, Smax, dS, dt);
        cout<<"option price = "<<put_price_cn<<endl;
    }
    return 0;
}
