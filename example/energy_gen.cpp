#include <iostream>
#include "Stellar.h"

using namespace std;

int main(){
	double rho,T,P;
	double mu,X,Y,Z;
	X=0.70;
	Z=0.01;
	T=1E7;//[K],T_{7}=1
	rho=1E5;//[kg/m^3],100[g/cm^2]
//	rho=1E3;//[kg/m^3],1[g/cm^3]

	EnergyGeneration ene = EnergyGeneration();
	cout << "PP-chain" << endl;
	cout << " RPN:"<< scientific << ene.PP_RPN(rho,X,T) << "[J/s/kg]" << endl;
	cout << " KIP:"<< scientific << ene.PP_KIP(rho,X,T) << "[J/s/kg]" << endl;
	cout << " AML:"<< scientific << ene.PP_AML(rho,X,T) << "[J/s/kg]" << endl;
	cout << "    :"<< scientific << ene.PP_REDUCED(rho,X,T) << "[J/s/kg]" << endl;

	cout << "CNO-cycle" << endl;
	T=1.5E7;
	cout << " KIP:" << scientific << ene.CNO_KIP(rho,X,Z,T) << "[J/s/kg]" << endl;
}