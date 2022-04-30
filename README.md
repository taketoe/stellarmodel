# Stellar Model Numerical Calculation Code
## Overview
 This package povides numerical calculation codes of stellar interia model. You can calculate and plot stellar interia structure which are temeperature, presurre, energy flux and so on.

## Installation
 The software is using following external packages.
### 1. matplotlib-cpp
  [Matplotlib-cpp](https://github.com/lava/matplotlib-cpp) pvovides to use matplotlib in C++.
  Installation is easy to place a heade file,matplotlibcpp.h, to your directory. The header file is already includeed this package. You can also use newer one through Github. 
### Eigen
  
## Sample codes and result
```c++
 #include "Stellar.h"
 int main(){
	double fact_star=1.;//Mass factor to solar mass;
	double X=0.70;//Hydrogen mass fraction
	double Y=0.28;//Helium mass fraction
	double Z=0.02;//Heavy element abundance
	double Ts=1.E3;//Boundary sondition;Temerature of stellar surface 
	double Ps=1.E8;//Boundary sondition;Temerature of stellar Pressure
	Stellar stellar(fact_star,X,Y,Z,Ts,Ps);
	stellar.setLog(true);
	stellar.calc();//Main calculation routine
	stellar.getResult();
	stellar.plot();
 }
```
<img src="https://github.com/taketoe/stellarmodel/blob/master/example/sun.png" width="600px">
<img src="https://github.com/taketoe/stellarmodel/blob/master/example/sun10_100.png" width="600px">
<img src="https://github.com/taketoe/stellarmodel/blob/master/example/energy_gen_KIP.png" width="600px">
(*) x-axis Temperature[K] vs y-axis energy generation [J/s/kg]

## References
 1. Python codes and PDF files.[SPA7023U,R.P.Nelson,QM Plus](https://2019.qmplus.qmul.ac.uk/course/view.php?id=9017),Queen Mary University of London.
 2. Stellar Structure and Evolution,R.Kippenhahan and A.Weigert.

