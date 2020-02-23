//
//  main.cpp
//  HW3-2, Crank-Nicolson scheme and down-and-out barrier option
//
//  Created by JIE QIAN on 2/21/20.
//  Copyright © 2020 JIE QIAN. All rights reserved.
//

//European Down-and-out Barrier Option pricing
//with Crank-Nicolon finit difference scheme

#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "newmat.h"
#include "newmatap.h"
#include "newmatio.h"

using namespace std;

double N(const double& z) {
    if (z > 6.0) { return 1.0; }; // this guards against overflow
    if (z < -6.0) { return 0.0; };
    double b1 = 0.31938153;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
    double p = 0.2316419;
    double c2 = 0.3989423;
    double a = fabs(z);
    double t = 1.0/(1.0+a*p);
    double b = c2*exp((-z)*(z/2.0));
    double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
    n = 1.0-b*n;
    if ( z < 0.0 ) n = 1.0 - n;
    return n;
};

//Theorial price for down-and-out put option
double down_and_out_put_option_black_scholes(const double& S,      // spot price
                                             const double& K,      // Strike price
                                             const double& r,      // interest rate
                                             const double& sigma,  // volatility
                                             const double& B,      // barrier
                                             const double& T)
{
    double a = pow((B/S), (-1+(2*r)/(sigma*sigma)));
    double b = pow((B/S), (1+(2*r)/(sigma*sigma)));
    double d1 = (log(S/K) + (r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    double d2 = (log(S/K)+(r-0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    double d3 = (log(S/B)+(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    double d4 = (log(S/B)+(r-0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    double d5 = (log(S/B)-(r-0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    double d6 = (log(S/B)-(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    double d7 = (log((S*K)/(B*B))-(r-0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    double d8 = (log((S*K)/(B*B))-(r+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
    
    double put_price = K*exp(-r*T)*(N(d4)-N(d2)-a*(N(d7)-N(d5)))
                        - S*(N(d3)-N(d1)-b*(N(d8)-N(d6)));
    return put_price;
};

//price for down and out put option with Crank-Nicolson finint difference scheme
double down_and_out_put_option_Crank_Nicolson(const double& S,      // spot price
                                              const double& K,      // Strike price
                                              const double& r,      // interest rate
                                              const double& sigma,  // volatility
                                              const double& B,      // barrier
                                              const double& T,      //time
                                              const double& Smax,
                                              double ds,
                                              double dt)
{
    double sigma_sqr = sigma*sigma;
    int M = round((Smax - B)/ds);
    int N = round(T/dt);
    //ds = (Smax-B)/M;
    //dt = T/N;
    
    BandMatrix C(M+1, 1, 1);
    C = 0.0;
    C.element(0, 0) = 1+0.5*dt*(sigma_sqr + r);
    C.element(0, 1) = -0.25*dt*(sigma_sqr + r);
    C.element(M, M-1) = -0.25*dt*((double) sigma_sqr*(M+1)*(M+1)-(double) r*(M+1));
    C.element(M, M) = 1+0.5*dt*((double) sigma_sqr*(M+1)*(M+1)+r);
    for(int j=1; j<M; j++)
    {
        double aj = 0.25*dt*((double) sigma_sqr*(j+1)*(j+1)-(double) r*(j+1));
        double bj = -0.5*dt*((double) sigma_sqr*(j+1)*(j+1)+r);
        double cj = 0.25*dt*((double) sigma_sqr*(j+1)*(j+1)+(double) r*(j+1));
        C.element(j, j-1) = -aj;
        C.element(j, j) = 1-bj;
        C.element(j, j+1) = -cj;
    }
    
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
        D.element(j, j-1) = aj;
        D.element(j, j) = 1+bj;
        D.element(j, j+1) = cj;
    }
    
    ColumnVector F(M+1);
    F = 0.0;
    for(int m=0; m<M; m++)
    {
        F.element(m) = max(0.0, K - (double) m*ds - B);
    }
    
    for(int t=N; t>=0; t--)
    {
        F = C.i() * D * F;
        F.element(0) = 0.0;
        F.element(M) = 0.0;
    }
    
    //Because in matrix from in Newmat, the index of the first index is 0;
    //Besides, when m=1, S equals B, the first index is then invaild.
    //the whole indexes system is then decrease, so we need to add back 2 indexes
    int price = round((S-B)/ds)+2;
    return F.element(price);
    
};

int main(int argc, const char * argv[]) {
    double exact_price = down_and_out_put_option_black_scholes(50.0, 50.0, 0.1, 0.4, 40.0, (5.0/12.0));
    cout<<"Exact price by Black-Sholes PDE: "<<exact_price<<endl;
    
    double CN_FD_price = down_and_out_put_option_Crank_Nicolson(50.0, 50.0, 0.1, 0.4, 40.0, 5.0/12.0, 100.0, 0.5, 1.0/1200.0);
    cout<<"Price by Crank-Nicolson Scheme:  "<<CN_FD_price<<endl;
    
    return 0;
}
