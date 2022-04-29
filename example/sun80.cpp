#include "Stellar.h"

int main(){
	double fact_star = 80.0;//Mass factor to solar mass;M_{*}/M_{\circledcirc}
	double X=0.70;//Hydrogen mass fraction
	double Y=0.28;//Helium mass fraction
	double Z=0.02;//Heavy element abundance
	double Pc = 8.54133e+15;//Initial guess value
	Stellar stellar(fact_star,X,Y,Z);
	stellar.setLog(true);
	stellar.setPc(Pc);
	stellar.calc();
	stellar.getResult();
	stellar.plot();
}
