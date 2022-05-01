# Stellar Model Numerical Calculation Code
## Overview
 This package povides numerical calculation codes of stellar interia model. You can calculate and plot stellar interia structure which are temeperature, presurre, energy flux and so on.

## Installation
 The software is using following external packages.
### 1. matplotlib-cpp
  [Matplotlib-cpp](https://github.com/lava/matplotlib-cpp) pvovides to use matplotlib in C++.
  Installation is easy to place a heade file,matplotlibcpp.h, to your directory. The header file is already includeed this package. You can also use newer one through Github. 
##### Modification
  The matplolibcpp.h are modifieded following three line of subplot routine,ll.1646~1648.
  From:
```c++
    PyTuple_SetItem(args, 0, PyFloat_FromDouble(nrows));
    PyTuple_SetItem(args, 1, PyFloat_FromDouble(ncols));
    PyTuple_SetItem(args, 2, PyFloat_FromDoubl(plot_number));
```
  To:
```c++
  PyTuple_SetItem(args, 0, PyLong_FromLong(nrows));
  PyTuple_SetItem(args, 1, PyLong_FromLong(ncols));
  PyTuple_SetItem(args, 2, PyLong_FromLong(plot_number));
```

### Eigen
  
## Sample codes and result
```c++
#include "Stellar.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace std;
// This program calculate stellar interia
// structure from mass factor 1 
// Boundary condition are surface temperature Ts and pressure Ps
// Ts:1E3[K], Ps1E8[K]
int main(){
	std::vector<double> M,R,T,P,L,rho;
	double fact_star=1.;//Mass factor to solar mass;M_{*}/M_{\circledcirc}
	double X=0.70;//Hydrogen mass fraction
	double Y=0.28;//Helium mass fraction
	double Z=0.02;//Heavy element abundance
	double Ts=1.E3;
	double Ps=1.E8;
	Stellar stellar(fact_star,X,Y,Z,Ts,Ps);
	stellar.setLog(true);
	stellar.calc();
	stellar.getResult();
	M = stellar.getM();
	R = stellar.getR();
	T = stellar.getT();
	P = stellar.getP();
	L = stellar.getL();
	rho = stellar.getRho();

	plt::figure_size(1000,600);
	plt::subplot(2,3,1);
	plt::plot(M,T,"r-");
	plt::ylim(*min_element(T.begin(),T.end()),*max_element(T.begin(),T.end()));
	plt::title("Temperature vs M");
	plt::xlabel("M [kg]");
	plt::ylabel("T [K]");
	
	plt::subplot(2,3,2);
	plt::plot(M,P,"r-");
	plt::ylim(*min_element(P.begin(),P.end()),*max_element(P.begin(),P.end()));
	plt::title("Pressure vs M");
	plt::xlabel("M [kg]");
	plt::ylabel("P [N/m$^2$]");

	plt::subplot(2,3,3);
	plt::plot(M,rho,"r-");
	plt::ylim(*min_element(rho.begin(),rho.end()),*max_element(rho.begin(),rho.end()));
	plt::title("Density vs M");
	plt::xlabel("M [kg]");
	plt::ylabel("$\\rho \\rm [kg/m^3]$");

	plt::subplot(2,3,4);
	plt::plot(M,L,"r-");
	plt::ylim(0.,1.1*(*max_element(L.begin(),L.end())));
	plt::title("Luminosity vs M");
	plt::xlabel("M [kg]");
	plt::ylabel("$L \\rm [W]$");

	plt::subplot(2,3,5);
	plt::plot(M,R,"r-");
	plt::ylim(*min_element(R.begin(),R.end()),*max_element(R.begin(),R.end()));
	plt::title("Radius vs M");
	plt::xlabel("M [kg]");
	plt::ylabel("$R \\rm[m]$");

	plt::subplot(2,3,6);
	plt::plot(R,T,"r-");
	plt::ylim(*min_element(T.begin(),T.end()),*max_element(T.begin(),T.end()));
	plt::title("Temperature vs Radius");
	plt::ylabel("T [K]");
	plt::xlabel("$R \\rm[m]$");
	plt::tight_layout();
	plt::show();
}

```
<img src="https://github.com/taketoe/stellarmodel/blob/master/example/sun.png" width="600px">
<img src="https://github.com/taketoe/stellarmodel/blob/master/example/sun10_100.png" width="600px">
<img src="https://github.com/taketoe/stellarmodel/blob/master/example/energy_gen_KIP.png" width="600px">

## References
 1. Python codes and PDF files.[SPA7023U,R.P.Nelson,QM Plus](https://2019.qmplus.qmul.ac.uk/course/view.php?id=9017),Queen Mary University of London.
 2. Stellar Structure and Evolution,R.Kippenhahan and A.Weigert.

