#include "Stellar.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace std;
// This program calculate stellar interia
// structure from mass factor 1 
// Boundary condition are surface temperature Ts and pressure Ps
// Ts:1E3[K], Ps1E8[K]
int main(){
	std::vector<double> M,R,T,P,L,rho;
	double fact_star=100.;//Mass factor to solar mass;M_{*}/M_{\circledcirc}
	double X=0.70;//Hydrogen mass fraction
	double Y=0.28;//Helium mass fraction
	double Z=0.02;//Heavy element abundance
	double Ts=1.E3;
	double Ps=1.E8;
	// Stellar stellar(fact_star,X,Y,Z,Ts,Ps);
	// stellar.setLog(true);
	// stellar.setPc(1.555e+17);
	// stellar.calc();
	// stellar.getResult();
	e_state state=notconverge;
	Stellar stellar;

	double Pc = 1.E16;//initial guess value
	int i=0;
	while(!(state==converge) && i<100000){
		stellar = Stellar(fact_star,X,Y,Z,Ts,Ps);
		Pc=Pc*1.001;
		stellar.setPc(Pc);
		stellar.setLog(false);
		state=stellar.calc();

		if(state == converge){
			cout << " converge" << endl;
			cout << " Pc:"<< Pc << endl;
		}
		i++;
	}if(state==notconverge){cout<<" not converge"<<endl;}

	stellar.getResult();

	M = stellar.getM();
	R = stellar.getR();
	T = stellar.getT();
	P = stellar.getP();
	L = stellar.getL();
	rho = stellar.getRho();

	plt::figure_size(1000,600);
	plt::subplot(2,3,1);
	plt::plot(M,T,"r-");
	plt::ylim(*min_element(T.begin(),T.end()),*max_element(T.begin(),T.end()));
	plt::title("Temperature vs M");
	plt::xlabel("M [kg]");
	plt::ylabel("T [K]");
	
	plt::subplot(2,3,2);
	plt::plot(M,P,"r-");
	plt::ylim(*min_element(P.begin(),P.end()),*max_element(P.begin(),P.end()));
	plt::title("Pressure vs M");
	plt::xlabel("M [kg]");
	plt::ylabel("P [N/m$^2$]");

	plt::subplot(2,3,3);
	plt::plot(M,rho,"r-");
	plt::ylim(*min_element(rho.begin(),rho.end()),*max_element(rho.begin(),rho.end()));
	plt::title("Density vs M");
	plt::xlabel("M [kg]");
	plt::ylabel("$\\rho \\rm [kg/m^3]$");

	plt::subplot(2,3,4);
	plt::plot(M,L,"r-");
	plt::ylim(0.,1.1*(*max_element(L.begin(),L.end())));
	plt::title("Luminosity vs M");
	plt::xlabel("M [kg]");
	plt::ylabel("$L \\rm [W]$");

	plt::subplot(2,3,5);
	plt::plot(M,R,"r-");
	plt::ylim(*min_element(R.begin(),R.end()),*max_element(R.begin(),R.end()));
	plt::title("Radius vs M");
	plt::xlabel("M [kg]");
	plt::ylabel("$R \\rm[m]$");

	plt::subplot(2,3,6);
	plt::plot(R,T,"r-");
	plt::ylim(*min_element(T.begin(),T.end()),*max_element(T.begin(),T.end()));
	plt::title("Temperature vs Radius");
	plt::ylabel("T [K]");
	plt::xlabel("$R \\rm[m]$");
	plt::tight_layout();

	plt::show();

}
