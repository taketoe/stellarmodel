#include "Stellar.h"

int main(){
	double fact_star = 60.0;//fact_star=<10 is now calculatable
	double x=0.70;
	double y=0.28;
	double z=0.02;
	Stellar stellar(fact_star,x,y,z);
	stellar.setLog(true);
	stellar.calc();
	stellar.getResult();
	stellar.plot();

}
