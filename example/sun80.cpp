#include "Stellar.h"

int main(){
	double fact_star = 80.0;//Mass factor to solar mass;M_{*}/M_{\circledcirc}
	double X=0.70;//Hydrogen mass fraction
	double Y=0.28;//Helium mass fraction
	double Z=0.02;//Heavy element abundance
	double Ts=1.E3;
	double Ps=1.E8;
	double Pc = 8.54133e+15;//Initial guess value
//	double Pc = 1E15;
	Stellar stellar(fact_star,X,Y,Z,Ts,Ps);
	stellar.setLog(true);
	stellar.setPc(Pc);
	stellar.calc();
	stellar.getResult();
	stellar.plot();
}
