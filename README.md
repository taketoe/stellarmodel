# Stellar Model Numerical Calsulation Code
## Overview
 This package povides numerical calculation codes of stellar interia model. You can calculate 

## Installation
 The software is using following external packages.
### matplotlibcpp
  Matplotlibcpp pvovides to use matplotlib in C++.
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
	Stellar stellar(fact_star,X,Y,Z);
	stellar.setLog(true);
	stellar.calc();
	stellar.getResult();
	stellar.plot();
 }
```

## References
### 1. Python codes and PDF files.[SPA7023U,R.P.Nelson,QM Plus](https://2019.qmplus.qmul.ac.uk/course/view.php?id=9017),Queen Mary University of London.
### 2. Stellar Structure and Evolution,R.Kippenhahan and A.Weigert.

