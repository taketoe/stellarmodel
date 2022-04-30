#include <iostream>
#include "Stellar.h"

using namespace std;
namespace plt = matplotlibcpp;
using namespace Eigen;

int main(){
	double rho,T,P;
	double mu,X,Y,Z;
	const int arr_size = 100;
	int n = 100;
	std::vector<double> T_,e_pp,e_CNO,e;
	X=0.70;
	Z=0.01;
	T=1E7;//[K],T_{7}=1
	rho=1E5;//[kg/m^3],100[g/cm^2]
//	rho=1E3;//[kg/m^3],1[g/cm^3]

	EnergyGen ene = EnergyGen();
	cout << "PP-chain" << endl;
	cout << " RPN:"<< scientific << ene.PP_RPN(rho,X,T) << "[J/s/kg]" << endl;
	cout << " KIP:"<< scientific << ene.PP_KIP(rho,X,T) << "[J/s/kg]" << endl;
	cout << " AML:"<< scientific << ene.PP_AML(rho,X,T) << "[J/s/kg]" << endl;
	cout << "    :"<< scientific << ene.PP_REDUCED(rho,X,T) << "[J/s/kg]" << endl;

	cout << "CNO-cycle" << endl;
	T=1.5E7;
	cout << " KIP:" << scientific << ene.CNO_KIP(rho,X,Z,T) << "[J/s/kg]" << endl;

	double Tl,Th,Tstep;
	Tl=1E6;Th=1E8;	
	Tstep=(Th-Tl)/arr_size;
	rho=1.E3;
	int i=0;
	for(T=1E6;T<1E8;T=T+Tstep){
		T_.push_back(T);
		e_pp.push_back(ene.PP_KIP(rho,X,T));
		e_CNO.push_back(ene.CNO_KIP(rho,X,Z,T));
		e.push_back(ene.PP_KIP(rho,X,T)+ene.CNO_KIP(rho,X,Z,T));
		i++;
	}

	plt::figure_size(900,700);
	plt::loglog(T_,e_pp);
	plt::loglog(T_,e_CNO);
	plt::loglog(T_,e);
	plt::xlim(1.E6,1.E8);
	plt::ylim(1.E-10,1.E8);
	plt::legend();
	plt::show();

}