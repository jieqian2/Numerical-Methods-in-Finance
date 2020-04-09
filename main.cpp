//
//  main.cpp
//  3.Variance_Reduction_in_MC_simulation
//
//  Created by JIE QIAN on 4/8/20.
//  Copyright © 2020 Jay QIAN. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>
using namespace std;

int num_of_barriers, num_of_trails;
double divide_rate, risk_free_rate, strike_price, initial_stock_price, expiration_time, volatility, barrier_price;
float delta_T, delta_R, delta_SD;

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator (seed);

double normal(){
    normal_distribution<double> distribution (0.0,1.0);
    return distribution(generator);
}

void initialize(){
    num_of_barriers = 25;
    num_of_trails = 10000;
    divide_rate = 0;
    risk_free_rate = 0.03;
    strike_price = 105;
    initial_stock_price = 99;
    expiration_time = 25.0/252.0;
    volatility = 0.6;
    barrier_price = 90;
    
    delta_T = expiration_time/((float) num_of_barriers);
    delta_R = (risk_free_rate-divide_rate-0.5*pow(volatility,2))*delta_T;
    delta_SD = volatility*sqrt(delta_T);
}

vector<double> plain_monte_carlo(){
    vector<double> call_option_price;
    vector<double> result;
    
      for(int i=0; i<num_of_trails; i++){
          float current_stock_price = initial_stock_price;
          
          for(int j=0; j<num_of_barriers; j++){
              double fi = normal();
              if(current_stock_price*exp(delta_R+delta_SD*fi) > barrier_price){
                  current_stock_price = current_stock_price*exp(delta_R + delta_SD*fi);
              }else{
                  current_stock_price = 0;
              }
          }
          
          double call_price = max(0.0, current_stock_price - strike_price);
          call_option_price.push_back(call_price);
      }
    
    double sum = accumulate(begin(call_option_price), end(call_option_price), 0.0);
    double mean =  sum / num_of_trails;
      
    mean = mean * exp(-risk_free_rate*expiration_time);
      
    double accum  = 0.0;
    double actual_accum  = 0.0;
    for(int i=0; i<num_of_trails; i++){
        accum += (call_option_price[i]-mean)*(call_option_price[i]-mean);
        actual_accum += (call_option_price[i]-4.647650)*(call_option_price[i]-4.647650);
    };
    double stdev = sqrt(accum/(num_of_trails - 1));
    double actual_error = sqrt(accum/(num_of_trails - 1));
    
    result.push_back(mean);
    result.push_back(stdev);
    result.push_back(actual_error);
    
    return result;
}

int main(int argc, const char * argv[]) {
    initialize();
    
    //Plain vanilla Monte Carlo Simulations
    vector<double> plain_MC;
    plain_MC =plain_monte_carlo();
    cout<<"-----------------------------------------"<<endl;
    cout<<"Plain vanilla Monte Carlo Simulations"<<endl;
    cout<<"option price estimate: "<<plain_MC[0]<<endl;
    cout<<"standard error:        "<<plain_MC[1]<<endl;
    cout<<"actual error:          "<<plain_MC[2]<<endl;
    
    
    return 0;
}
