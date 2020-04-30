//
//  main.cpp
//  Pricing Auto-Callable Contingent Interest Notes
//
//  Created by JIE QIAN on 4/21/20.
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

double notes_pricing_CN_FD(const double& initial_value,
                           const double& interest_barrier,
                           const double& r, //risk free rate
                           const double& sigma, //volatility
                           const double& T, //maturity time
                           const double& Smax, //P(S>Smax) near 0
                           const double& Smin, //P(S<Smin) near 0
                           double ds,
                           double dt)
{
    Tracer tr("name");
    double sigma_sqr = sigma*sigma;
    int MIN = round(Smin/ds);
    int M = round((Smax-Smin)/ds);
    int IN = round((initial_value-Smin)/ds);
    int IB = round((interest_barrier-Smin)/ds);
    int N = round(T/dt);
    // review dates
    int rd1 = round((94.0/365.0)/dt);
    int rd2 = round((185.0/365.0)/dt);
    int rd3 = round((276.0/365.0)/dt);
    int rd4 = round((367.0/365.0)/dt);

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
            VT.element(i) = 1019.125 *exp(-r*2.0/365.0);
        else if(i >= IB)
            VT.element(i) = (1019.125)*((double) i/(IN+MIN)) *exp(-r*2.0/365.0);
        else
            VT.element(i) = (1000.00)*((double) i/(IN+MIN)) *exp(-r*2.0/365.0);
            
    }
    
    //in this mesh, trigger event hasn't happened
    ColumnVector V(M+1);
    V = 0.0;
    for(int i=0; i<M; i++){
        V.element(i) = (1019.125)*((double) i/(IN+MIN)) *exp(-r*2.0/365.0);
    }
    
    
    // backward to the initial date
    for(int t=N; t>=0; t--)
    {
        //LBC
        VT.element(0) = 0.0;
        V.element(0) = 0.0;
        
        // on non-coupon non-callable date
        if(t!=rd1 && t!=rd2 && t!=rd3 && t!=rd4){
            VT = C.i() * D * VT;
            V = C.i() * D * V;
            for(int i=0; i<IB; i++){
                V.element(i) = VT.element(i);
            }
        }
        
        // on coupon only date
        if(t==rd1){
            VT = C.i() * D * VT;
            
            for(int i=0; i<IB; i++){
                V.element(i) = VT.element(i);
            }
            V = C.i() * D * V;

            for(int i=IB; i<M; i++){
                VT.element(i) += 19.125 *exp(-r*((double)t*dt));
                V.element(i) += 19.125 *exp(-r*((double)t*dt));
            }
        }
        
        // on coupon and autocall date
        if(t==rd2 || t==rd3 || t==rd4){
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
                VT.element(i) += 19.125 *exp(-r*((double)t*dt));
                V.element(i) += 19.125 *exp(-r*((double)t*dt));
            }
        }
        
//        if(t%10==0){
//            cout<<round(((double) (N-t)/N)*100)<<"%.. ";
//        }
    }
//    cout<<V<<endl;
//    cout<<"in "<<IN<<"th step the value is "<<V.element(IN)<<endl;
    return V.element(IN);
};


int main(int argc, const char * argv[]) {
    
    double initial_value = 1785.00;
    double interest_barrier = 1249.50; //0.7*initial_value
    double r = 0.00443456842; //LIBOR, forward rate due on maturity date
    double sigma = 0.36031; //choose K=S0(basic case:moneyness=100), or K=IB, 0.41321;
    double T = 458.0/365.0; // 458 days between 3/13/20 to 6/14/21
    double Smax = 3*initial_value; //highly unlikely to be touched in 15month
    double Smin = 0.0;
    double dS = 0.01*initial_value; //0.01*S0
    double dt = 1.0/365.0; // every days 1.0/365.0
    
    //benchmark
    double CN_FD_price = notes_pricing_CN_FD(initial_value, interest_barrier, r, sigma, T,
                                             Smax, Smin, dS, dt);
    cout<<"Value: "<<CN_FD_price<<endl;


    //sensitivity analysis
    vector<double> S_step;
    vector<double> t_step;
    for(int i=100; i>=10; i--){
        S_step.push_back( ((double) i/1000.0) * initial_value );
        t_step.push_back( ((double) i/10.0) * (1.0/365.0));
    }

    ofstream outf;
    outf.open("convergence.csv", ios::out | ios::trunc);
    for(int i=0; i<S_step.size(); i++){
        for(int j=0; j<t_step.size(); j++){
            double price = notes_pricing_CN_FD(initial_value, interest_barrier, r, sigma, T, Smax, Smin, S_step[i], t_step[j]);
            outf<<setprecision(13)<<price<<",";
            cout<<S_step[i]<<" "<<t_step[j]<<" "<<setprecision(13)<<price<<endl;
        }
        outf<<endl;
    }
    outf.close();

    //sensitivity to volatility
    // change it here

    vector<double> vololity = {0.41321, 0.40304, 0.39402, 0.38517, 0.37597, 0.36823,0.36031};
    vector<double> moneyness = {70, 75, 80, 85, 90, 95, 100};
    for(int i=0; i<7; i++){
        double CN_FD_price_vo = notes_pricing_CN_FD(initial_value, interest_barrier, r, vololity[i], T, Smax, Smin, dS, dt);
        cout<<"moneyness  "<<moneyness[i]<<", Volitity "<<vololity[i]<<", Value: "<<CN_FD_price_vo<<endl;
    }

    return 0;
}

