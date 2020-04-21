Pricing Auto Callable Contingent Interest Notes
====

# Intro 

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

* Trigger event: A Trigger Event occurs if the closing price ≤ the Trigger Value
  
![image](https://github.com/jieqian2/Numerical-Methods-in-Finance/blob/master/IMG/figure1.png)

![image](https://github.com/jieqian2/Numerical-Methods-in-Finance/blob/master/IMG/figure2.png)

![image](https://github.com/jieqian2/Numerical-Methods-in-Finance/blob/master/IMG/figure3.png)

# Pricing Methods

## Finite Difference Method


