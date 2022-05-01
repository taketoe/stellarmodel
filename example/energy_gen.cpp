#include <iostream>
#include "Stellar.h"
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

using namespace std;

int main(){
	double rho,T,P;
	double mu,X,Y,Z;
	X=0.70;
	Z=0.01;
//	T=1E7;//[K],T_{7}=1
	T=1E8;//[K],T_{7}=10
//	rho=1E5;//[kg/m^3],100[g/cm^2]
	rho=1E3;//[kg/m^3],1[g/cm^3]

	//EnergyGen ene = EnergyGen();
	cout << "PP-chain" << endl;
	cout << " RPN:"<< scientific << EnergyGen::PP_RPN(rho,X,T) << "[J/s/kg]" << endl;
	cout << " KIP:"<< scientific << EnergyGen::PP_KIP(rho,X,T) << "[J/s/kg]" << endl;
	cout << " AML:"<< scientific << EnergyGen::PP_AML(rho,X,T) << "[J/s/kg]" << endl;
	cout << "    :"<< scientific << EnergyGen::PP_REDUCED(rho,X,T) << "[J/s/kg]" << endl;

	cout << "CNO-cycle" << endl;
	T=1.5E7;
	cout << " KIP:" << scientific << EnergyGen::CNO_KIP(rho,X,Z,T) << "[J/s/kg]" << endl;
}