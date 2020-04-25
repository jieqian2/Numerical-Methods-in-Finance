Pricing Auto Callable Contingent Interest Notes
====

> I am as smart as JPMorgan guys!

# 1. Intro 

Full details in https://www.sec.gov/Archives/edgar/data/19617/000161577420003465/s124044_424b2.htm

Brief introduction (not so brief, but..):

* Face Value: $1000
* Reference Stock: Amazon.com, Inc.
* Pricing Date: March 13, 2020. Amazon's closing price is $1785.00
* Initial Value: $1785.00
* Final Value: Closing price of Amazon on `2021.6.14`
* Contingent Interest Payments:(on review dates)
  * If closing stock price ≥ Interest barrier: get contingent interest payment of $19.125
  * If closing stock price < Interest barrier: no payment
  
* Automatic Call: 
on any review dates(**expect the first and final review dates**):
  * If closing stock price ≥ Initial value: notes will be automatically called for a cash payment, 
  
    cash payment = $1000 + Contingent Interest Payment applicable to that Review Date.
* Interest Barrier(Trigger Value): 70% of Inivial value, which is $1785 * 0.7 = $1249.50
* Review Dates:           `2020.6.15` `2020.9.14` `2020.12.14` `2021.3.15` `2021.6.14`
* Interest Payment Dates:  `2020.6.17` `2020.9.16` `2020.12.16` `2021.3.17` `2021.6.16`, two days after the Review Dates.
* Maturity Date:          `2021.6.16`
* Call Settlement Date: If the notes are automatically called on any Review Dates (other than the first and final Review Dates), the first Interest Payment Date immediately following that Review Date

 
* Payment at Maturity: 
  * If Final Value ≥ Initial Value, or trigger event not occurred:
  
    cash payment = $1000 + Contingent Interest Payment applicable to the final Review Date
 
  * If Final Value < Initial Value, and trigger event occurred: 
    
    cash payment = $1000 + $1000 * stock return, stock return = (Final Value - Initial Value)/Initial Value

* Trigger event: A Trigger Event occurs if the closing price ≤ the Trigger Value. **Remark: the trigger event is detected all the time, not just in callable time.**


  
![image](https://github.com/jieqian2/Numerical-Methods-in-Finance/blob/master/IMG/figure1.png)

![image](https://github.com/jieqian2/Numerical-Methods-in-Finance/blob/master/IMG/figure2.png)

![image](https://github.com/jieqian2/Numerical-Methods-in-Finance/blob/master/IMG/figure3.png)


# 2. Pricing using Crank-Nicolson Finite Difference Method 

## 2.1 准备工作

Suppose the stock price follow GBM

$dS/S = rdt+sigmadS $

formals,

value of the derivatives f(S,t)

BSPDE,

apply C-N FD,



## 2.2 Set Grid


We need to set up two grids. The first grid, named VT, is when the trigger event happened. The second grid, named V, is when it didn't. 

Denote:
```cpp
dS = small price steps
dt = small time steps
i = number of steps on time
j = number of steps on price
imax = T/dt
jmax = Smax/dS
IN = initial_price/dS
IB = interest_barrier/dS
review dates:
rd1 = (94.0/365.0)/dt
rd2 = (185.0/365.0)/dt
rd3 = (276.0/365.0)/dt
rd4 = (367.0/365.0)/dt
```
Terminal Boundary Condition(TBC)
```cpp
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
 ColumnVector V(M+1);
 V = 0.0;
 for(int i=0; i<M; i++){
    V.element(i) = (1019.125)*((double) i/(IN+MIN)) *exp(-r*2.0/365.0);
 }
```
Low Boundary Condition(LBC)
```cpp
VT.element(0) = 0.0;
V.element(0) = 0.0;
```



## 2.3 Algorithm

The key here is to maintain two grids at the same time. VT is for trigger event happened, V is not. 

Then go backward, whenever stock price hit the trigger(IB/interest_barrier), the value becomes the value we calculated in VT 
(i.e. set Vi,j become the VTi,j).

On coupon and autocall date, we have to change the UBC. The LBC doesn't change, for V, its LBC is actually always the IB.(Since we swithc to VT for price below that).

At last, we add on coupon when we finish our value for VT and V.


### 2.3.1. On non-coupon and non-autocall date:

At every step i, solve VT from j=0 to j=jmax.

Then, solve V on IB<j<jmax, 


### 2.3.2. On coupon and autocall date (rd2, rd3, rd4):

UBC: all j≥IN,

VTi,j = 1000 * exp(-r * days/365) 

Vi,j = 1000 * exp(-r * days/365) 

then solve for VT on 0<j<IN, 

solve V on IB<j<IN, i.e. set Vi,j = VTi,j when j≤ IB 

If j≥IB, add on discounted coupon on VTi,j and Vi,j


### 2.3.3. On coupon and non-autocall date (rd1):

solve VT,

set all j≤IB  Vi,j = VTi,j

solve V

If j≥IB, add on discounted coupon on VTi,j and Vi,j


# 3. Sensitivity Analysis

## 3.1 Sensitivity to dS/steps of stock price

## 3.2 Sensitivity to dt/steps of experiation time

## 3.3 Sensitivity to Volatility of underlying asset 






