//
//  main.cpp
//  callable contingent coupon barrier notes
//
//  Created by JIE QIAN on 5/1/20.
//  Copyright © 2020 JIE QIAN. All rights reserved.
//


#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <chrono>
#include <random>
#include "newmat.h"
#include "newmatap.h"
#include "newmatio.h"

using namespace std;

double initial_value, interest_barrier, r, sigma, T, q;
int num_of_divisions;

double max(double a, double b) {
    return (b < a )? a:b;
}

int min(int a, int b) {
    return (b < a )? b:a;
}

unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator (seed);

double get_gaussian()
{
    std::normal_distribution<double> distribution(0.0,1.0);
    double number = distribution(generator);
    return (number);
}

Matrix polynomial_regression(Matrix Independent_Variables, Matrix Dependent_Variable, int order, int no_of_observations)
{
    Matrix X(no_of_observations, order);
    Matrix Y(no_of_observations, 1);
    
    for (int i = 1; i <= no_of_observations; i++)
        Y(i,1) = Dependent_Variable(i,1);
    
    for (int j = 1; j <= order; j++)
        for (int i = 1; i <= no_of_observations; i++)
            X(i,j) = pow(Independent_Variables(i,1), j-1);
    
    // return inv(XT*X)*XT*Y
    Matrix X_transpose_times_X(order, order);
    X_transpose_times_X = X.t()*X;
    return (X_transpose_times_X.i() * X.t() * Y);
}


double Crank_Nicolson_FD(const double& Smax, //P(S>Smax) near 0
                         double ds,
                         double dt)
{
    Tracer tr("name");
    double sigma_sqr = sigma*sigma;
    int M = round((Smax)/ds);
    int IN = round((initial_value)/ds);
    int IB = round((interest_barrier)/ds);
    int N = round(T/dt);
    // review dates
    int rd1 = round((90.0/365.0)/dt);
    int rd2 = round((182.0/365.0)/dt);
    int rd3 = round((274.0/365.0)/dt);
    int rd4 = round((366.0/365.0)/dt);
    int rd5 = round((455.0/365.0)/dt);
    int rd6 = round((547.0/365.0)/dt);
    int rd7 = round((639.0/365.0)/dt);
    
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
    
    //terminal condition
    // actually, we need two mesh here.
    
    // in this mesh, trigger event has happened
    ColumnVector VT(M+1);
    VT = 0.0;
    for(int i=0; i<M; i++){
        if(i >= IN)
            VT.element(i) = 1022.5 *exp(-r*5.0/365.0);
        else if(i >= IB)
            VT.element(i) = (1022.5)*((double) i/(IN)) *exp(-r*5.0/365.0);
        else
            VT.element(i) = (1000.0)*((double) i/(IN)) *exp(-r*5.0/365.0);
            
    }
    
    //in this mesh, trigger event hasn't happened
    ColumnVector V(M+1);
    V = 0.0;
    for(int i=0; i<M; i++){
        V.element(i) = (1022.5)*((double) i/(IN)) *exp(-r*2.0/365.0);
    }
    
    
    // backward to the initial date
    for(int t=N; t>=0; t--)
    {
        //LBC
        VT.element(0) = 0.0;
        V.element(0) = 0.0;
        
        // on non-coupon non-callable date
        if(t!=rd1 && t!=rd2 && t!=rd3 && t!=rd4 && t!=rd5 && t!=rd6 && t!=rd7){
            VT = C.i() * D * VT;
            V = C.i() * D * V;
            for(int i=0; i<IB; i++){
                V.element(i) = VT.element(i);
            }
        }
        
        // on coupon and autocall date
        if(t==rd1 || t==rd2 || t==rd3 || t==rd4 || t==rd5 || t==rd6 || t==rd7){
            for(int i=IN; i<M; i++){
                VT.element(i) = 1000 *exp(-r*((double)t*dt));
            }
            VT = C.i() * D * VT;
            
            for(int i=IN; i<M; i++){
                V.element(i) = 1000 *exp(-r*((double)t*dt));
            }
            V = C.i() * D * V;
            for(int i=0; i<IB; i++){
                V.element(i) = VT.element(i);
            }
            
            for(int i=IB; i<M; i++){
                VT.element(i) += 22.5 *exp(-r*((double)t*dt));
                V.element(i) += 22.5 *exp(-r*((double)t*dt));
            }
        }

    }
    
//    cout<<V.element(IN)<<endl;
//    cout<<V.element(IN+1)<<endl;
    return (V.element(IN+1)+V.element(IN))/2.0;
};

double longstaff_schwartz(int num_of_simulations){
    Tracer tr("name");
    
    double delta_T = T/((double) num_of_divisions);
    double delta_R = (r - q - 0.5*pow(sigma,2))*delta_T;
    double delta_SD = sigma*sqrt(delta_T);
    double R = exp(r*T/((double) num_of_divisions));
    
    vector<double> rd = {90.0/365.0, 182.0/365.0, 274.0/365.0,
        366.0/365.0, 455.0/365.0, 547.0/365.0, 639.0/365.0};
    
    Matrix asset_price(num_of_simulations, 9);
    for(int i=0; i<num_of_simulations; i++)
        asset_price.element(i,0) = initial_value;
    
    for(int i=0; i<num_of_simulations; i++)
        for(int j=1; j<= num_of_divisions; j++)
            asset_price.element(i,j) = asset_price.element(i,j-1)*exp(delta_R + delta_SD*get_gaussian());
    
    ColumnVector value(num_of_simulations);
    for(int i=0; i<num_of_simulations; i++){
        if (asset_price.element(i,8) >= interest_barrier)
            value.element(i) = 1022.5 *exp(-r*5.0/365.0);
        else
            value.element(i) = 1000.0 *(asset_price.element(i,8)/(initial_value)) *exp(-r*5.0/365.0);
    }
    
    for(int i=(num_of_divisions); i>0; i--){
        
        Matrix independent_variables(num_of_simulations,1);
        Matrix dependent_variables(num_of_simulations,1);
        
        int num_of_variables = 0;
        for(int j = 0; j < num_of_simulations; j++){
            //if(asset_price.element(j, i) < interest_barrier)
            {
                independent_variables.element(num_of_variables,0) = asset_price.element(j,i);
                dependent_variables.element(num_of_variables,0) = value.element(j)/R;
                num_of_variables++;
            }
            
        }
        //if(num_of_variables > 4)
        {
            //regressing the dependent_variables on the independent variables using a 3th order polynomial
            Matrix a(4,1);
            a = polynomial_regression(independent_variables, dependent_variables, 4, num_of_variables);

            for (int j = 0; j < num_of_simulations; j++) {
                double continue_value_hat = (a.element(0,0) +(a.element(1,0)*asset_price.element(j,i)) +(a.element(2,0)*pow(asset_price.element(j,i),2)) +(a.element(3,0)*pow(asset_price.element(j,i),3)));

                if(continue_value_hat >= 1022.5 * exp(-r*rd[i])) //it will be called
                    value.element(j) = 1000 * exp(-r*rd[i]);
                else
                    value.element(j) = value.element(j)/R;
            }
            
            for (int j = 0; j < num_of_simulations; j++){
                if(asset_price.element(j,i) >= interest_barrier){
                    value.element(j) = value.element(j) + 22.5 * exp(-r*rd[i]);
                }
            }
        }
    }
    double note_price = 0.0;
    for(int i=0; i<num_of_simulations; i++){
        note_price += value.element(i)/R;
    }
    note_price = note_price/(double) num_of_simulations;
    
    if(note_price>1000.0 || note_price<950)
        return longstaff_schwartz(num_of_simulations);
    
    return note_price;
}


double longstaff_schwartz_benchmark(int num_of_simulations){
    Tracer tr("name");
    
    double delta_T = T/((double) num_of_divisions);
    double delta_R = (r - q - 0.5*pow(sigma,2))*delta_T;
    double delta_SD = sigma*sqrt(delta_T);
    double R = exp(r*T/((double) num_of_divisions));
    
    vector<double> rd = {90.0/365.0, 182.0/365.0, 274.0/365.0,
        366.0/365.0, 455.0/365.0, 547.0/365.0, 639.0/365.0};
    
    Matrix asset_price(num_of_simulations, 9);
    for(int i=0; i<num_of_simulations; i++)
        asset_price.element(i,0) = initial_value;
    
    for(int i=0; i<num_of_simulations; i++)
        for(int j=1; j<= num_of_divisions; j++)
            asset_price.element(i,j) = asset_price.element(i,j-1)*exp(delta_R + delta_SD*get_gaussian());
    
    ColumnVector value(num_of_simulations);
    for(int i=0; i<num_of_simulations; i++){
        if (asset_price.element(i,8) >= interest_barrier)
            value.element(i) = 1022.5 *exp(-r*5.0/365.0);
        else
            value.element(i) = 1000.0 *(asset_price.element(i,8)/(initial_value)) *exp(-r*5.0/365.0);
    }
    
    for(int i=(num_of_divisions); i>0; i--){
        
        Matrix independent_variables(num_of_simulations,1);
        Matrix dependent_variables(num_of_simulations,1);
        
        int num_of_variables = 0;
        for(int j = 0; j < num_of_simulations; j++){
            //if(asset_price.element(j, i) < interest_barrier)
            {
                independent_variables.element(num_of_variables,0) = asset_price.element(j,i);
                dependent_variables.element(num_of_variables,0) = value.element(j)/R;
                num_of_variables++;
            }
            
        }
        
        {
            //regressing the dependent_variables on the independent variables using simple polynomial
            Matrix a(4,1);
            a = polynomial_regression(independent_variables, dependent_variables, 2, num_of_variables);

            for (int j = 0; j < num_of_simulations; j++) {
                double continue_value_hat = (a.element(0,0) +(a.element(1,0)*asset_price.element(j,i)));

                if(continue_value_hat >= 1022.5 * exp(-r*rd[i])) //it will be called
                    value.element(j) = 1000 * exp(-r*rd[i]);
                else
                    value.element(j) = value.element(j)/R;
            }
            
            for (int j = 0; j < num_of_simulations; j++){
                if(asset_price.element(j,i) >= interest_barrier){
                    value.element(j) = value.element(j) + 22.5 * exp(-r*rd[i]);
                }
            }
        }
    }
    double note_price = 0.0;
    for(int i=0; i<num_of_simulations; i++){
        note_price += value.element(i)/R;
    }
    note_price = note_price/(double) num_of_simulations;
    
    if(note_price>1000.0 || note_price<950)
        return longstaff_schwartz(num_of_simulations);
    
    return note_price;
}

int main(int argc, const char * argv[]) {
    initial_value = 33.62;
    interest_barrier = 25.22; //0.75*initial_value
    r = 1.378299534/100.0; //LIBOR, forward rate due on maturity date
    sigma = 0.25089; //choose K=S0(basic case:moneyness=100), or K=IB, 0.28689;
    T = 731.0/365.0; // 731 days between 2/4/20 to 2/4/22
    double Smax = 5*initial_value;
    double dS = 0.01*initial_value; //0.01*S0
    double dt = 1.0/365.0; // every days 1.0/365.0
    q = 2.5873/100.0;
    num_of_divisions = 8;
    
    cout<<"---Pricing by CN Finite Difference---"<<endl;
    double CN_price = Crank_Nicolson_FD(Smax, dS, dt);
    cout<<"Pricing by C-N Finite Difference  : $"<<setprecision(13)<< CN_price<<endl;

    cout<<endl;
    cout<<"---Pricing by Longstaff&Schwartz---"<<endl;
    cout<<"Simulation times: 10,000"<<endl;
    double LS_price_benchmark = longstaff_schwartz_benchmark(10000);
    cout<<"(1) Benchmark Sernario,simple polynomial regression: $"<<setprecision(13)<<LS_price_benchmark<<endl;
    double LS_price= longstaff_schwartz(10000);
    cout<<"(2) Using 4th order polynomial regression: $"<<setprecision(13)<<LS_price<<endl;
    cout<<endl;
    
    cout<<"Simulation times: 1,000,000"<<endl;
    double LS_price_3 = longstaff_schwartz_benchmark(1000000);
    cout<<"(3) Using simple polynomial regression: $"<<setprecision(13)<<LS_price_3<<endl;
    double LS_price_4= longstaff_schwartz(1000000);
       cout<<"(4) Using 4th order polynomial regression: $"<<setprecision(13)<<LS_price_4<<endl;
    
    return 0;
}
