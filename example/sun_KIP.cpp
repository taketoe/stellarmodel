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
	double fact_star;//Mass factor to solar mass;M_{*}/M_{\circledcirc}
	double X=0.70;//Hydrogen mass fraction
	double Y=0.28;//Helium mass fraction
	double Z=0.02;//Heavy element abundance
	double Ts,Ps,Pc;
	fact_star=1.;Ts=Teff_sun;Ps=8.E8;Pc=1.0E16;
	Stellar stellar1 = Stellar(fact_star,X,Y,Z,Ts,Ps);//Stellar model class

	fact_star=10.;Ts=25119;Ps=Ps;Pc=1.0E16;//Converged
	Stellar stellar2 = Stellar(fact_star,X,Y,Z,Ts,Ps);


	stellar1.setPc(Pc);
	stellar1.setLog(true);
	stellar1.calc();
	stellar1.getResult();
	M = stellar1.getM();
	R = stellar1.getR();
	T = stellar1.getT();
	P = stellar1.getP();
	L = stellar1.getL();
	rho = stellar1.getRho();
	epsilon = stellar1.getEpsilon();

	cout << " RESULT: Rstar/Rsun="<<R[Ndim-1]/Rsun << endl;
	cout << " RESULT:         Ts="<<T[Ndim-1]<< endl;
	cout << " RESULT: Lstar/Lsun="<<L[Ndim-1]/Lsun<< endl;

	for(auto itr = M.begin(); itr != M.end(); ++itr) {
        *itr;
		*itr=*itr / Msun;
    }

	for(auto itr = rho.begin(); itr != rho.end(); ++itr) {
        *itr;
		*itr=log10(*itr);
    }

	for(auto itr = T.begin(); itr != T.end(); ++itr) {
        *itr;
		*itr=log10(*itr);
    }

	for(auto itr = epsilon.begin(); itr != epsilon.end(); ++itr) {
        *itr;
		*itr=log10(*itr);
    }

	double R_s=R[Ndim-1];
	for(auto itr = R.begin(); itr != R.end(); ++itr) {
        *itr;
		*itr=*itr / R_s;
    }

	double L_s=L[Ndim-1];
	for(auto itr = L.begin(); itr != L.end(); ++itr) {
        *itr;
		*itr=*itr / L_s;
    }

//	fact_star=10.;Ts=1.E3;Ps=1.E8;Pc=1.0E16;
//	fact_star=10.;Ts=1.E3;Ps=0.18E8;Pc=1.0E16;
//	fact_star=10.;Ts=1.5E3;Ps=1.E8;Pc=3.7E15;//target
//	fact_star=10.;Ts=1.5E3;Ps=0.248E8;Pc=1.0E16;//Converged
	stellar2.setPc(Pc);
//	stellar.setLog(true);
	stellar2.calc();
	stellar2.getResult();
	std::vector<double> M2,R2,T2,P2,L2,rho2,epsilon2;
	M2 = stellar2.getM();
	R2 = stellar2.getR();
	T2 = stellar2.getT();
	P2 = stellar2.getP();
	L2 = stellar2.getL();
	rho2 = stellar2.getRho();
	epsilon2 = stellar2.getEpsilon();

	cout << " RESULT: Rstar/Rsun="<<R2[Ndim-1]/Rsun << endl;
	cout << " RESULT:         Ts="<<T2[Ndim-1]<< endl;
	cout << " RESULT:         Ls="<<L2[Ndim-1]/Lsun<< endl;

	for(auto itr = M2.begin(); itr != M2.end(); ++itr) {
        *itr;
		*itr=*itr / (Msun*fact_star);
    }

	for(auto itr = rho2.begin(); itr != rho2.end(); ++itr) {
        *itr;
		*itr=log10(*itr);
    }

	for(auto itr = T2.begin(); itr != T2.end(); ++itr) {
        *itr;
		*itr=log10(*itr);
    }

	for(auto itr = epsilon2.begin(); itr != epsilon2.end(); ++itr) {
        *itr;
		*itr=log10(*itr);
    }

	double R2_s=R2[Ndim-1];
	for(auto itr = R2.begin(); itr != R2.end(); ++itr) {
        *itr;
		*itr=*itr / R2_s;
    }

	double L2_s=L2[Ndim-1];
	for(auto itr = L2.begin(); itr != L2.end(); ++itr) {
        *itr;
//		*itr=*itr / (Lsun*fact_star);
		*itr=*itr / L2_s;
    }

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
//	plt::loglog(R,rho,"r-");
	plt::plot(M,rho,"r-",args1);
	plt::plot(M2,rho2,"b-",args2);
//	plt::plot(M,rho,".");
	plt::xlim(0.,1.);
	plt::ylim(2.,5.);
	plt::title("Density");
//	plt::xlabel("$\\rm M/M_{\\circledcirc}[-]$");
	plt::xlabel("$\\rm m/M$");
	plt::ylabel("$\\rm log \\rho$");
	plt::grid(true);
	plt::legend();

	plt::subplot(3,2,2);
	plt::plot(R,M,"r-");
	plt::plot(R2,M2,"b-");
//	plt::plot(R,M,".");
//	plt::loglog(M,P,"r-");
//	plt::xlim(0.,0.7);
//	plt::ylim(0.,1.);
	plt::title("Mass");
	plt::xlabel("$\\rm r/R$");
	plt::ylabel("$\\rm m/M$");
	plt::grid(true);

	plt::subplot(3,2,3);
	plt::plot(M,T,"r-");
	plt::plot(M2,T2,"b-");
//	plt::plot(M,T,".");
	plt::xlim(0.,1.);
	plt::ylim(6.4,8.0);
	plt::title("Temperature");
	plt::xlabel("$\\rm m/M$");
	plt::ylabel("log T");
	plt::grid(true);

	plt::subplot(3,2,4);
	plt::plot(M,epsilon,"r-");
	plt::plot(M2,epsilon2,"b-");
//	plt::plot(M,epsilon,".");
	plt::xlim(0.,0.7);
	plt::ylim(-6.,2.0);
	plt::title("Energy generation $\\epsilon$");
//	plt::xlabel("$\\rm m/M_{\\circledcirc}$");
	plt::xlabel("$\\rm m/M$");
	plt::ylabel("$\\rm log \\epsilon $");
	plt::grid(true);

	plt::subplot(3,2,5);
	plt::plot(M,L,"r-");
	plt::plot(M2,L2,"b-");
//	plt::plot(M,L,".");
	plt::xlim(0.,0.7);
	plt::ylim(0.,1.);
	plt::title("Luminosity");
//	plt::xlabel("$\\rm m/M_{\\circledcirc}$");
	plt::xlabel("$\\rm m/M$");
//	plt::ylabel("$\\rm L/L_{\\circledcirc}$");
	plt::ylabel("$\\rm l/L$");
	plt::grid(true);

	plt::tight_layout();
	plt::show();

}
