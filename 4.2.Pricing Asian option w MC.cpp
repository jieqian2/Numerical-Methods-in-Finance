//
//  main.cpp
//  importance sampling
//
//  Created by JIE QIAN on 4/9/20.
//  Copyright © 2020 JIE QIAN. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>
using namespace std;

int num_of_days;
double divide_rate, risk_free_rate, strike_price, initial_stock_price, expiration_time, volatility;
float delta_T, delta_R, delta_SD;

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator (seed);

double normal(){
    normal_distribution<double> distribution (0.0,1.0);
    return distribution(generator);
}

void initialize(){
    num_of_days = 90;
    divide_rate = 0;
    risk_free_rate = 0.03;
    strike_price = 105;
    initial_stock_price = 99;
    expiration_time = (double) num_of_days/252.0;
    volatility = 0.6;
    
    delta_T = expiration_time/((float) num_of_days);
    delta_R = (risk_free_rate-divide_rate-0.5*pow(volatility,2))*delta_T;
    delta_SD = volatility*sqrt(delta_T);
}

vector<double> plain_monte_carlo(int num_of_simulations){
    vector<double> option_price;
    vector<double> result;
    
      for(int i=0; i<num_of_simulations; i++){
          float current_stock_price = initial_stock_price;
          double asian_sum = 0.0;
          for(int j=0; j<num_of_days; j++){
              current_stock_price = current_stock_price*exp(delta_R + delta_SD*normal());
              asian_sum += current_stock_price;
          }
          double call_price_end = max(0.0, asian_sum /(1.0+ (double)num_of_days) - strike_price);
          option_price.push_back(call_price_end);
      }
    
    double sum = accumulate(begin(option_price), end(option_price), 0.0);
    double mean =  sum / num_of_simulations;
      
    mean = mean * exp(-risk_free_rate*expiration_time);
      
    double accum  = 0.0;
    for(int i=0; i<num_of_simulations; i++){
        accum += (option_price[i]-mean)*(option_price[i]-mean);
    };
    double stdev = sqrt(accum/(num_of_simulations - 1));
    stdev = stdev/sqrt(num_of_simulations);
    
    result.push_back(mean);
    result.push_back(stdev);
    
    return result;
}


int main(int argc, const char * argv[]) {
    initialize();
    
    //Plain vanilla Monte Carlo Simulations
    vector<double> plain_MC = plain_monte_carlo(100000);
    cout<<"-----------------------------------------"<<endl;
    cout<<"Plain vanilla Monte Carlo Simulations"<<endl;
    cout<<"option price estimate: "<<plain_MC[0]<<endl;
    cout<<"standard error:        "<<plain_MC[1]<<endl;
    
    ofstream outf1;
    outf1.open("convergence.csv", ios::out | ios::trunc);
    outf1<<"steps,price"<<endl;
    int repeat=1000;
    while(repeat<=200000){
        vector<double> price;
        price = plain_monte_carlo(repeat);
        outf1<<repeat<<","<<price[0]<<endl;
        repeat += 1000;
    }
    outf1.close();
//
    return 0;
}
