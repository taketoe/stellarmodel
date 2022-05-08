#include "Stellar.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace std;
// This program calculate stellar interia
// structure from mass factor 1 
// Boundary condition are surface temperature Ts and pressure Ps
// Ts:1E3[K], Ps1E8[K]
int main(){
	std::vector<double> M,R,T,P,L,rho,epsilon;
//	std::vector<double> M2,R2,T2,P2,L2,rho2,epsilon2;
	double fact_star;//Mass factor to solar mass;M_{*}/M_{\circledcirc}
//	double X=0.70;//Hydrogen mass fraction
//	double Y=0.28;//Helium mass fraction
//	double Z=0.02;//Heavy element abundance
	double X=0.685;//Hydrogen mass fraction
	double Y=0.294;//Helium mass fraction
	double Z=0.021;//Heavy element abundance
	double Ts,Ps,Pc;
//	fact_star=1.;Ts=Teff_sun;Ps=8.E8;Pc=1.0E16;//converge
	fact_star=1.;Ts=Teff_sun;Ps=1.E7;Pc=1.0E16;
	Stellar stellar1;

	double Ps_min = 1.E6;
	double Ps_max = 1.E15;
//	double dPs;
//	int i=0;
	Ps=Ps_min;
//	dPs = floor(log10f(Ps));

	e_state state;
	stringstream ss;
//	ss<<scientific<<setprecision(5);
	ss<<fixed<<setprecision(5);
	//for(Ps=Ps_min;Ps<Ps_max;Ps+=dPs){
	while(Ps<Ps_max){
		stellar1 = Stellar(fact_star,X,Y,Z,Ts,Ps);//Stellar model class
		stellar1.setPc(Pc);
		state=stellar1.calc();
		string msg;
		if(state == e_state::converge){
			P = stellar1.getPress();
			L = stellar1.getLuminoNormaSun();
			R = stellar1.getRadiusNormaSun();
			ss << "log10(Ps):" << log10f(P[Ndim-1]);
			ss << " Ls/Lsun:" << L[Ndim-1]; 
			ss << " Rs/Rsun:" <<  R[Ndim-1] << endl; 
		}else{
			cout << "Ps:"<< Ps << " : not converged"<< endl; 
		}
//		cout << "Ps:"<< Ps << endl; 
		Ps+=pow(10,floor(log10f(Ps)));
	}
	cout << "*** CONVERGED POINT ***"<< endl;
	cout << ss.str() << endl;

	/*

	M1 = stellar1.getMassNormaSrface();
	R1 = stellar1.getRadiusNormaSurface();
	T1 = stellar1.getTempLog();
	P1 = stellar1.getPress();
	L1 = stellar1.getLuminoNormaSurface();
	rho1 = stellar1.getRhoLog();
	epsilon1 = stellar1.getEpsilonLog();

	cout << " RESULT: Rstar/Rsun="<<R1[Ndim-1] << endl;
	cout << " RESULT:         Ts="<<T1[Ndim-1] << endl;
	cout << " RESULT: Lstar/Lsun="<<L1[Ndim-1] << endl;


    map<string, string> args1{
            {"label", "$1M_{\\circledcirc}$"},
            {"c", "red"}
    };

    map<string, string> args2{
            {"label", "$10M_{\\circledcirc}$"},
            {"c", "blue"}
    };

	plt::figure_size(600,650);
	plt::subplot(3,2,1);
	plt::plot(M1,rho1,args1);
	plt::plot(M2,rho2,args2);
	plt::xlim(0.,1.);
	plt::ylim(2.,5.);
	plt::title("Density");
	plt::xlabel("$\\rm m/M$");
	plt::ylabel("$\\rm log \\rho$");
	plt::grid(true);
	plt::legend();

	plt::subplot(3,2,2);
	plt::plot(R1,M1,"r");
	plt::plot(R2,M2,"b");
	plt::title("Mass");
	plt::xlabel("$\\rm r/R$");
	plt::ylabel("$\\rm m/M$");
	plt::grid(true);

	plt::subplot(3,2,3);
	plt::plot(M1,T1,"r");
	plt::plot(M2,T2,"b");
	plt::xlim(0.,1.);
	plt::ylim(6.4,8.0);
	plt::title("Temperature");
	plt::xlabel("$\\rm m/M$");
	plt::ylabel("log T");
	plt::grid(true);

	plt::subplot(3,2,4);
	plt::plot(M1,epsilon1,"r");
	plt::plot(M2,epsilon2,"b");
	plt::xlim(0.,0.7);
	plt::ylim(-6.,2.0);
	plt::title("Energy generation $\\epsilon$");
//	plt::xlabel("$\\rm m/M_{\\circledcirc}$");
	plt::xlabel("$\\rm m/M$");
	plt::ylabel("$\\rm log \\epsilon $");
	plt::grid(true);

	plt::subplot(3,2,5);
	plt::plot(M1,L1,"r");
	plt::plot(M2,L2,"b");
	plt::xlim(0.,0.7);
	plt::ylim(0.,1.);
	plt::title("Luminosity");
	plt::xlabel("$\\rm m/M$");
	plt::ylabel("$\\rm l/L$");
	plt::grid(true);

	plt::tight_layout();
	plt::show();
*/
}
