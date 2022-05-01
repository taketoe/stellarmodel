#include <iostream>
#include "Stellar.h"
#include "matplotlibcppgithub.h"

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

//	EnergyGen ene = EnergyGen();
	cout << "PP-chain" << endl;
	cout << " RPN:"<< scientific << EnergyGen::PP_RPN(rho,X,T) << "[J/s/kg]" << endl;
	cout << " KIP:"<< scientific << EnergyGen::PP_KIP(rho,X,T) << "[J/s/kg]" << endl;
	cout << " AML:"<< scientific << EnergyGen::PP_AML(rho,X,T) << "[J/s/kg]" << endl;
	cout << "    :"<< scientific << EnergyGen::PP_REDUCED(rho,X,T) << "[J/s/kg]" << endl;

	cout << "CNO-cycle" << endl;
	T=1.5E7;
	cout << " KIP:" << scientific << EnergyGen::CNO_KIP(rho,X,Z,T) << "[J/s/kg]" << endl;

	double Tl,Th,Tstep;
	Tl=1E6;Th=1E8;	
	Tstep=(Th-Tl)/arr_size;
	rho=1.E3;
	int i=0;
	for(T=1E6;T<1E8;T=T+Tstep){
		T_.push_back(T);
		e_pp.push_back(EnergyGen::PP_KIP(rho,X,T));
		e_CNO.push_back(EnergyGen::CNO_KIP(rho,X,Z,T));
		e.push_back(EnergyGen::PP_KIP(rho,X,T)+EnergyGen::CNO_KIP(rho,X,Z,T));
		i++;
	}

    map<string, string> args_epp{
            {"label", "$\\epsilon_{pp}$"},
            {"c", "red"}
    };

    map<string, string> args_ecno{
            {"label", "$\\epsilon_{CNO}$"},
            {"c", "blue"}
    };

    map<string, string> args_e{
            {"label", "$\\epsilon_{pp}$+$\\epsilon_{CNO}$"},
            {"c", "green"}
    };

	plt::figure();
	//plt::figure_size(900,700);

	plt::loglog(T_,e_pp,args_epp);
	plt::loglog(T_,e_CNO,args_ecno);
	plt::loglog(T_,e,args_e);
	plt::xlim(1.E6,1.E8);
	plt::ylim(1.E-10,1.E8);
	plt::xlabel("T [K]");
	plt::ylabel("$\\epsilon$ [J/s/kg]");
	plt::title("Energy Generation of Hydrogen Burning in Stellar");
	plt::legend();
	plt::show();

}