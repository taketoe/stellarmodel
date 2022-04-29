#include "Stellar.h"

// This program calculate stellar interia
// structure from mass factor 1 
// Boundary condition are surface temperature Ts and pressure Ps
// Ts:1E3[K], Ps1E8[K]
int main(){
	double fact_star=1.;//Mass factor to solar mass;M_{*}/M_{\circledcirc}
	double X=0.70;//Hydrogen mass fraction
	double Y=0.28;//Helium mass fraction
	double Z=0.02;//Heavy element abundance
	Stellar stellar(fact_star,X,Y,Z);
	stellar.setLog(true);
	stellar.calc();
	stellar.getResult();
	stellar.plot();

}
