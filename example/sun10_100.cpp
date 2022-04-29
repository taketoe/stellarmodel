#include <iostream>
#include "Stellar.h"

namespace plt = matplotlibcpp;
using namespace std;

// This program calculate stellar interia
// structure from mass factor 10 to 100 
// Boundary condition are surface temperature Ts and pressure Ps
// Ts:1E3[K], Ps1E8[K]
int main(){
	double fact_star=1.;//Mass factor to solar mass;M_{*}/M_{\circledcirc}
	double X=0.70;//Hydrogen mass fraction
	double Y=0.28;//Helium mass fraction
	double Z=0.02;//Heavy element abundance
	double Pc = 1.E15;//Initial value of central pressure
	Stellar stellar;//Stellar model class
	e_state state=notconverge;
	int i = 0;

	const int arr_size=9;
	Eigen::VectorXd _M(arr_size),_Pc(arr_size),_Tc(arr_size);
	Eigen::VectorXd _Ls(arr_size),_Rs(arr_size);
	std::array<double,arr_size> M_,Pc_,Tc_,Ls_,Rs_;

	int j=0;
	int step=10;
	for(fact_star=10;fact_star<int(step*arr_size)+1;fact_star+=step){
		cout << "Mass factor:" << fact_star;
		while(!(state==converge) && i<100000){
			stellar = Stellar(fact_star,X,Y,Z);
			Pc=Pc*1.001;
			stellar.setPc(Pc);
			stellar.setLog(false);
			state=stellar.calc();

			if(state == converge){
//				cout << " converge" << endl;
				_M(j)=fact_star;M_[j]=_M(j);
				_Pc(j)=stellar.getPhys(0).getP();Pc_[j]=_Pc(j);
				_Tc(j)=stellar.getPhys(0).getT();Tc_[j]=_Tc(j);
				_Ls(j)=stellar.getPhys(Ndim-2).getL();Ls_[j]=_Ls(j);
				_Rs(j)=stellar.getPhys(Ndim-2).getR();Rs_[j]=_Rs(j);
				cout << " Pc:"<< Pc_[j];
				cout << " Tc:"<< Tc_[j]<<endl;
				j++;
			}
			//if(i%100==0){cout << "  i:"<<i << endl;}
			i++;
		}if(state==notconverge){cout<<" not converge"<<endl;}
		state=notconverge;
	}

	cout << "M :" << _M.transpose() <<endl;
	cout << "Pc:" << _Pc.transpose() <<endl;
	cout << "Tc:" << _Tc.transpose() <<endl;
	cout << "Ls:" << _Ls.transpose() <<endl;
	cout << "Rs:" << _Rs.transpose() <<endl;

	plt::figure_size(900,700);
	plt::title("Stellar model;Numerical calculation");
	plt::subplot(2,2,1);
	plt::title("Mass factor vs Central Pressure");
	plt::xlabel("$\\rm M_{*}/M_{\\circledcirc} [-]$");
	plt::ylabel("$\\rm P_{c} [Pa]$");
	plt::plot(M_,Pc_,"or-");

	plt::subplot(2,2,2);
	plt::title("Mass factor vs Central Temperature");
	plt::xlabel("$\\rm M_{*}/M_{\\circledcirc} [-]$");
	plt::ylabel("$\\rm T_{c} [K]$");
	plt::plot(M_,Tc_,"or-");
	plt::tight_layout();

	plt::subplot(2,2,3);
	plt::title("Mass factor vs Surface Luminosity");
	plt::xlabel("$\\rm M_{*}/M_{\\circledcirc} [-]$");
	plt::ylabel("$\\rm L_{s} [W]$");
	plt::plot(M_,Ls_,"or-");
	plt::tight_layout();

	plt::subplot(2,2,4);
	plt::title("Mass factor vs Surface Radius");
	plt::xlabel("$\\rm M_{*}/M_{\\circledcirc} [-]$");
	plt::ylabel("$\\rm R_{s} [m]$");
	plt::plot(M_,Rs_,"or-");
	plt::tight_layout();

	plt::show();

}
