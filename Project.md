Pricing Auto Callable Contingent Interest Notes
====

> I am as smart as JPMorgan guys!

# 1. Intro 

https://www.sec.gov/Archives/edgar/data/19617/000161577420003465/s124044_424b2.htm

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

## 2.1 Set Mesh


## 2.2 Set Boundary Condition

### 2.2.1 Up Boundary Condition

Vi,jmax = (1000 + coupon) * exp(-r*(Tnac - i * dt)
Tnac = review dates: 94;185;276;367;458(maturity)

### 2.2.2 Low Boundary Condition
Trigger event happened:

Vi,0 = 1000 + 1000 * stock return = 0

### 2.2.3 Terminal Boundary Condition

Terminal value depends on path.




 
## 2.3 Consider Callable Feature

we need two grids:

The first one is when the trigger event has occurred, the second is when it has not.

Then on the second grid, whenever you hit the trigger then value becomes the value that you calculated in the first grid. That is the lower boundary condition at S = trigger is the value already calculated in the first grid.


review dates: 94;185;276;367;458(maturity)


