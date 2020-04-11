Autocallable Yield Notes Pricing
====
 
Description
----
 UBS AG Airbag Autocallable Yield Notes (the “Notes”) are unsubordinated, unsecured debt obligations issued by UBS AG (“UBS” or the “issuer”) linked to the common stock of a specific company (the “underlying asset”). The issue price of each Note will be $1,000. UBS will pay you a coupon on each coupon payment date regardless of the performance of the underlying asset unless the Notes were previously subject to an automatic call. If the closing level of the underlying asset is equal to or greater than the initial level on any observation date prior to the final valuation date, UBS will automatically call the Notes (an “automatic call”) and pay you a cash payment per Note equal to your principal amount plus the coupon otherwise due on the applicable coupon payment date following such observation date ( the “call settlement date”), and no further payments will be owed to you under the Notes. If the Notes are not subject to an automatic call and the closing level of the underlying asset on the final valuation date (the “final level”) is equal to or greater than the conversion level, UBS will pay you a cash payment per Note equal to the principal amount. If, however, the Notes are not subject to an automatic call and the final level is less than the conversion level, UBS will deliver to you a number of shares of the underlying asset per Note equal to the quotient of (i) the principal amount divided by (ii) the conversion level (the “share delivery amount”), which is expected to be worth less than your principal amount and, in extreme situations, you could lose all of your initial investment. Any fractional share included in the share delivery amount will be paid in cash at an amount equal to the product of the fractional share and the final level. 
 
Detailed description, see: https://www.sec.gov/Archives/edgar/data/1114446/000091412119003603/ub54526835-424b2.htm

Method
----
Use a binomial tree to value a slightly simplified version of this product for General Motors:

* Face value: $1000

* Payoff at maturity: If the stock price is greater than or equal to $32.78 then you receive the face value plus the final coupon payment. If the stock price is below $32.78 you receive a cash amount equivalent to 1000/32.78 = 30.5064 stocks plus the final coupon payment.

* Autocall feature: If the stock price is greater than or equal to the initial price on any of the observation dates t = 1/4, 1/2, 3/4 then the notes are immediately called for the face value + coupon.

* Coupons: There is a monthly coupon of 8% (this is an annual figure) of the face value, payable at t = 1/12, 2/12 etc. 

* Proportional Dividends: 4.195% (annual figure) quarterly at t = 2/12, 5/12, 8/12, 11/12, so assume 0.25 x 4.195% is paid each quarter.

* Initial stock price = 37.25, T = 1, r = 1.755%, volitatity =  26.125%


```cpp
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <vector>

using namespace std;

double up_factor, down_factor, uptick_prob, downtick_prob, risk_free_rate, strike_price;
double initial_stock_price, expiration_time, dividends_rate, volatility, R;
double period_dividend_rate, face_value;
int no_of_divisions;

vector<vector<double>> callstore;
vector<int> dividend_date;
vector<int> autocall_date;



//change the initial condition here;
void initialize(){
    //user-change parameters
    initial_stock_price = 37.25;
    expiration_time = 1.0;
    risk_free_rate = 0.01755;
    volatility = 0.2615;
    dividends_rate = 0.04195;
    face_value = 1000.0;
    no_of_divisions = 1200;
    
    //auto-initialize all other parameters
    vector<double> tmp;
    tmp.resize(2*no_of_divisions+2, -1.0);
    callstore.resize(2*no_of_divisions+2, tmp);
    
    dividend_date.push_back((int) (no_of_divisions * (2.0/12.0)));
    dividend_date.push_back((int) (no_of_divisions * (5.0/12.0)));
    dividend_date.push_back((int) (no_of_divisions * (8.0/12.0)));
    dividend_date.push_back((int) (no_of_divisions * (11.0/12.0)));
    
    period_dividend_rate = (double)dividends_rate/dividend_date.size();
    
    autocall_date.push_back((int) (no_of_divisions * (1.0/4.0)));
    autocall_date.push_back((int) (no_of_divisions * (2.0/4.0)));
    autocall_date.push_back((int) (no_of_divisions * (3.0/4.0)));
    
    up_factor = exp(volatility*sqrt((expiration_time)/((float) no_of_divisions)));
    down_factor = 1.0/up_factor;
    R = exp(risk_free_rate*expiration_time/((float) no_of_divisions));
    uptick_prob = (R - down_factor)/(up_factor - down_factor);
    downtick_prob = 1.0 - uptick_prob;
}

double autocallable_yield_note(int k, int i, double current_stock_price, double coupon){
    if(callstore[k][no_of_divisions+i] != -1.0)
        return callstore[k][no_of_divisions+i];
    
    //check if current step is the dividend date;
    //if so, add coupon
    for(int i=0; i<dividend_date.size(); i++){
        if(k == dividend_date[i]){
            coupon += period_dividend_rate*face_value;
        }
    }
    
    //check if current step is the autocall date;
    //if so, check whether satisify the call condition
    for(int i=0; i<autocall_date.size(); i++){
        if(k == autocall_date[i]){
            if(current_stock_price >= initial_stock_price){
                return face_value+coupon;
            }else{
                break;
            }
        }
    }
    
    if(k == no_of_divisions){
        if(current_stock_price >= 32.78){
            callstore[k][no_of_divisions+i] = face_value + coupon;
            return callstore[k][no_of_divisions+i];
        }else{
            callstore[k][no_of_divisions+i] = 30.5064*current_stock_price + coupon;
            return callstore[k][no_of_divisions+i];
        }
    }
    
    callstore[k][no_of_divisions+i] =
    (uptick_prob * autocallable_yield_note(k+1,i+1,current_stock_price*up_factor, coupon)
     +downtick_prob * autocallable_yield_note(k+1,i-1,current_stock_price*down_factor, coupon))/R;
    
    return callstore[k][no_of_divisions+i];
}


int main(int argc, const char * argv[]) {
    initialize();
    double note_price = autocallable_yield_note(0, 0, initial_stock_price, 0);
    cout<<"Value of Autocallable Yield Notes is "<<note_price<<endl;
    return 0;
}

