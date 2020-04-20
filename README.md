# Hull-White-model
Comparison between Monte-Carlo simulation and analytical solution of Zero-Coupen bond price under Hull-white model spot rate.

・This figure shows that some of sample paths of spot rate generated by Hull-white model.
![Hull_white](https://user-images.githubusercontent.com/54795218/79717075-8e1b6a00-8313-11ea-8d80-84211afd9948.png)

・In thsi repository, we show Zero-Coupen bond price by two different methods. One is analytical method, and the other one is Monte-Carlo simulation. We get the results:

Analytical Bond Price(*0.0,T*) = *0.825337*

Numerical Bond Price(*0.0,T*) = *0.825979*.

You can not compile this code unless you call nonapy2.hpp, random2.hpp and LibMK4.
