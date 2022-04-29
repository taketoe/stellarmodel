#include "Stellar.h"

// This program calculate stellar interia
// structure from mass factor 1 
// Boundary condition are surface temperature Ts and pressure Ps
// Ts:1E3[K], Ps1E8[K]
int main(){
	std::array<double,Ndim> M,R,T,P,L;
	double fact_star=1.;//Mass factor to solar mass;M_{*}/M_{\circledcirc}
	double X=0.70;//Hydrogen mass fraction
	double Y=0.28;//Helium mass fraction
	double Z=0.02;//Heavy element abundance
	Stellar stellar(fact_star,X,Y,Z);
	stellar.setLog(true);
	stellar.calc();
	stellar.getResult();
	M = stellar.getM();
	R = stellar.getR();
	P = stellar.getP();
	T = stellar.getT();


}
