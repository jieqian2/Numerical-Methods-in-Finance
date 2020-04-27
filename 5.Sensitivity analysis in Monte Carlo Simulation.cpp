//
//  main.cpp
//  Sensitivity Analysis
//
//  Created by JIE QIAN on 4/26/20.
//  Copyright © 2020 JIE QIAN. All rights reserved.
//
// using Pathwise method and Likelihood ratio method to calculate delta and gamma.

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>
#include <algorithm>

using namespace std;

double divide_rate, risk_free_rate, strike_price, initial_stock_price, expiration_time, volatility;
float delta_T, delta_R, delta_SD;

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator (seed);

double normal(){
    normal_distribution<double> distribution (0.0,1.0);
    return distribution(generator);
}

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

double BS_delta_eu_put(const double& S, // spot price
                       const double& K, // Strike (exercise) price,
                       const double& r,  // interest rate
                       const double& sigma,
                       const double& T){
    double time_sqrt = sqrt(T);
    double d1 = (log(S/K)+r*T)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double delta = -N(-d1);
    return delta;
    
}

double BS_gamma_eu_put(const double& S, // spot price
                       const double& K, // Strike (exercise) price,
                       const double& r,  // interest rate
                       const double& sigma,
                       const double& T){
    double time_sqrt = sqrt(T);
    double d1 = (log(S/K)+r*T)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double gamma = exp(-d1*d1/2.0) / (sqrt(2*3.14159265359)*sigma*S*time_sqrt);
    return gamma;
}

// delta = dY/dS0 = dY/dST * dST/dS0 = e^-rT * ST/S0 * 1{K>ST}
double pathwise_delta_eu_put(const double& S, // spot price
                             const double& K, // Strike (exercise) price,
                             const double& r,  // interest rate
                             const double& sigma,
                             const double& T,
                             const int& number_of_simulations){
    double sum = 0;
    for(int i=0; i<number_of_simulations; i++){
        double S_T = S*exp((r-0.5*sigma*sigma)*T + sigma*sqrt(T)*normal());
        if(K > S_T){
            sum += -exp(-r*T) * S_T/S;
        }
    }
    return sum/(double)number_of_simulations;
}

//pathwise method is inapplicable to estimate second derivatives.

double likelihood_ratio_delta_eu_put(const double& S, // spot price
                                     const double& K, // Strike (exercise) price,
                                     const double& r,  // interest rate
                                     const double& sigma,
                                     const double& T,
                                     const int& number_of_simulations){
    double sum = 0;
    for(int i=0; i<number_of_simulations; i++){
        double Z = normal();
        double S_T = S*exp((r-0.5*sigma*sigma)*T + sigma*sqrt(T)*Z);
        //likelihood ratio
        double LR = Z/(S*sigma*sqrt(T));
        double delta = exp(-r*T) * max(K-S_T,0.0) * LR;
        sum += delta;
    }
    return sum/(double)number_of_simulations;
}

double kesi(double x,
            const double& S, // spot price
            const double& r,  // interest rate
            const double& sigma,
            const double& T){
    return (log(x/S) - (r-0.5*sigma*sigma)*T)/(sigma*sqrt(T));
}

double likelihood_ratio_gamma_eu_put(const double& S, // spot price
                                     const double& K, // Strike (exercise) price,
                                     const double& r,  // interest rate
                                     const double& sigma,
                                     const double& T,
                                     const int& number_of_simulations){
    double sum = 0;
    for(int i=0; i<number_of_simulations; i++){
        double S_T = S*exp((r-0.5*sigma*sigma)*T + sigma*sqrt(T)*normal());
        double kesi_ST = kesi(S_T, S, r, sigma, T);
        //likelihood ratio
        double LR = ((kesi_ST*kesi_ST-1)/(S*S*sigma*sigma*T))-(kesi_ST/(S*S*sigma*sqrt(T)));
        double gamma = exp(-r*T) * max(K-S_T,0.0) * LR;
        sum += gamma;
    }
    return sum/(double)number_of_simulations;
    
    
}
int main(int argc, const char * argv[]) {
    vector<double> K = {90, 100, 110};
    cout<<"pathwise method is inapplicable for calculate gamma. Since E[d^2Y/dS0^2] ≠ d^2E[Y]/dS0^2, the pathwise second order derivatives exists with probability 1, but is entirely uninformative."<<endl;
    for(int i=0; i<3; i++){
        cout<<"When K = "<<K[i]<<endl;
        cout<<"BS delta = "<<BS_delta_eu_put(100, K[i], 0.05, 0.3, 0.1)<<endl;
        cout<<"pathwise delta = "<<pathwise_delta_eu_put(100, K[i], 0.05, 0.3, 0.1, 1000000)<<endl;
        cout<<"likelihood delta = "<<likelihood_ratio_delta_eu_put(100, K[i], 0.05, 0.3, 0.1, 1000000)<<endl;
        cout<<endl;
        cout<<"BS gamma = "<<BS_gamma_eu_put(100, K[i], 0.05, 0.3, 0.1)<<endl;
        cout<<"likelihood gamma = "<<likelihood_ratio_gamma_eu_put(100, K[i], 0.05, 0.3, 0.1, 1000000)<<endl;
        cout<<endl;
        cout<<endl;
    }
    return 0;
}
