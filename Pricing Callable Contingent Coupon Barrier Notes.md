Pricing Callable Contingent Coupon Barrier Notes 
====

# 1. Intro

* Full Version:
[Callable Contingent Coupon Barrier Notes Linked to the Common Stock of Bank of America Corporation](https://www.sec.gov/Archives/edgar/data/1000275/000114036120002500/form424b2.htm)

* Key Feature:

Valuation Date: February 4, 2022

Initial Stock Price:
$33.62, which was the closing price of the Reference Stock on the Trade Date.

Final Stock Price:
The closing price of the Reference Stock on the Valuation Date.

Trigger Price and
Coupon Barrier:
$25.22, which is 75% of the Initial Stock Price (rounded to two decimal places).

Contingent Coupon:
If the closing price of the Reference Stock is greater than or equal to the Coupon Barrier on the applicable Observation Date, we will pay the Contingent Coupon applicable to that Observation Date. You may not receive any Contingent Coupons during the term of the Notes.

Payment at Maturity(if held to maturity):
(1) If the Notes are not previously called, we will pay you at maturity an amount based on the Final Stock Price:
For each $1,000 in principal amount, $1,000 plus the Contingent Coupon at maturity, unless the Final Stock Price is less than the Trigger Price.
(2) If the Final Stock Price is less than the Trigger Price, then the investor will receive at maturity, for each $1,000 in principal amount, the number of shares of the Reference Stock equal to the Physical Delivery Amount, or at our election, the cash value of those shares.

Physical Delivery Amount:
For each $1,000 principal amount, a number of shares of the Reference Stock equal to the principal amount divided by the Initial Stock Price, subject to adjustment as described in the product prospectus supplement.

Call Feature:
The Notes may be called at our discretion on any Coupon Payment Date (other than the final Coupon Payment Date), if we send prior written notice, as described below.


Observation Dates:
Quarterly, on May 4, 2020, August 4, 2020, November 4, 2020, February 4, 2021, May 4, 2021, August 4, 2021, November 4, 2021 and the Valuation Date.

Coupon Payment Dates:
The Contingent Coupon, if payable, will be paid quarterly on May 7, 2020, August 7, 2020, November 9, 2020, February 9, 2021, May 7, 2021, August 9, 2021, November 9, 2021 and the Maturity Date.

Contingent Coupon:
We will pay you a Contingent Coupon during the term of the Notes, periodically in arrears on each Coupon Payment Date, under the conditions described below:
(1) If the closing price of the Reference Stock is greater than or equal to the Coupon Barrier on the applicable Observation Date, we will pay the Contingent Coupon applicable to that Observation Date.
(2) If the closing price of the Reference Stock is less than the Coupon Barrier on the applicable Observation Date, we will not pay you the Contingent Coupon applicable to that Observation Date.


Contingent Coupon Rate:
9.00% per annum (2.25% per quarter)




*The initial estimated value of the Notes as of the Trade Date is $975.37 per $1,000 in principal amount, which is less than the price to public. The actual value of the Notes at any time will reflect many factors, cannot be predicted with accuracy, and may be less than this amount. We describe our determination of the initial estimated value in more detail below.*




# 2. Pricing by Longstaff and Schwartz Algorithm(LSMC)

There is one specific point should we pay more attention to. 

When using Longstaff and Schwartz Algorithm to price American put option, we only choose in the money price paths to do regression, that's because the exercise boundary is clear, we never exercise it when S>K.

**But here, the boundary between call and non-call is unclear. We don't know the point at which the issuer switch from calling to not calling (the exercise boundary).**

**So actually we should choose all paths to do regression in this product.**


## 2.1 Using simple polynomial regression

```cpp
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
```

### 2.2 Using 4th order polynomial regression




```cpp
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


```

# 3. Outcome Analysis

I compared the simple regression and 4th order regression and different simulation.

Of course the senario with the 4th order regression and more simulation times is better.  Usually, using 4th order regression is more precise and have less variance.

```cpp
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
```

Outcome of one simulation:

```cpp
---Pricing by Longstaff&Schwartz---
Simulation times: 10,000
(1) Benchmark Sernario,simple polynomial regression: $982.9240456301
(2) Using 4th order polynomial regression: $974.0593603616

Simulation times: 1,000,000
(3) Using simple polynomial regression: $974.244949888
(4) Using 4th order polynomial regression: $973.5578147209
```
Remark: the considered real price is $975, pricing by Royal Bank of Canda, the issuer.

# 4. Appendix

[Full code](https://github.com/jieqian2/Numerical-Methods-in-Finance/blob/master/8.LSMC.cpp)

