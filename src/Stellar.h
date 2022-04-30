#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <array>
//#include "matplotlibcpp.h"

#include <Eigen/Dense>
#include <Eigen/LU>

const double Pi=3.141592;
const double G=6.67E-11;      // Gravitational constant
const double c=3.E8;          // Speed of light [m/s]
const double kB=1.38E-23;     // Boltzmann constant []
const double m_H=1.67E-27;    // Mass of hydrogen atom []
const double sigma_SB=5.67E-8;// Stefan-Boltzmann 
const double a=7.5657E-16;    // Radiation constant
const double gamma_c=5./3.;   // Ratio of specific heats
const double Msun=1.989E30;   // Solar mass [kg]
const double Rsun=6.957E8;    // Solar radius [m]
const double Lsun=3.828E26;   // Solar luminosity [J/s]

const double con_fact=0.3;	//????

const double tolerance=5.E-3;// convergence criteria
const long Ndim = 10000;// Number of grid points

const long maxIterate = 100;
const double pertubRatio = 0.01;

//Hear is sucesful parameters set
//mstar_fact:1, Ndim:10000,maxIterate:100,pertubRatio:0.01,tolerance5.E-3
//mstar_fact:10,Ndim:10000,maxIterate:100,pertubRatio:0.01,tolerance5.E-3

enum e_state{
	notconverge,
	converge,
	overshoot,
	overflow,
	steptoolarge,
};

class Phys
{
private:
	double M,R,P,T,L,Rho,Kappa,Epsilon;
public:
	double getM();
	double getR();
	double getP();
	double getT();
	double getL();
	double getRho();
	double getKappa();
	double getEpsilon();
	const Phys operator-(const Phys&) const;
	Phys();
	Phys(double M,double R,double P,double T,double L,double Rho,double Kappa,double Epsilon);
	~Phys();
	void Out();
};


class Stellar
{
private:
	double mstar_fact=1;
	double Rstar,Mstar;
	double X,Y,Z,mu;
	double epsilon_pp,epsilon_CNO;
	double kappa0;
	std::vector<double> M = std::vector<double>(Ndim);
	std::vector<double> R = std::vector<double>(Ndim);
	std::vector<double> T = std::vector<double>(Ndim);
	std::vector<double> P = std::vector<double>(Ndim);
	std::vector<double> L = std::vector<double>(Ndim);
	std::vector<double> rho = std::vector<double>(Ndim);
	std::vector<double> kappa = std::vector<double>(Ndim);
	std::vector<double> epsilon = std::vector<double>(Ndim);
	double dM,Pc,Tc,Ps,Ts,Rs,Ls;
	bool logOut;
	long numberOfIterate=0;
public:
	Stellar();
	//Stellar(double msun);
	Stellar(double msun,double X,double Y,double Z,double Ts,double Ps);
	~Stellar();
	e_state calc();
	void setLog(bool);
	void plot(){};
	Phys getPhys(long);
	void getResult();
	void setPc(double Pc){this->Pc=Pc;}
	std::vector<double> getM(){return M;}
	std::vector<double> getR(){return R;}
	std::vector<double> getP(){return P;}
	std::vector<double> getT(){return T;}
	std::vector<double> getL(){return L;}
	std::vector<double> getRho(){return rho;}
	std::vector<double> getKappa(){return kappa;}
	std::vector<double> getEpsilon(){return epsilon;}
private:
	double getPc(){return Pc;};
	double getTc(){return Tc;};
	double getPs(){return Ps;};
	double getTs(){return Ts;};
	double getRs(){return Rs;};
	double getLs(){return Ls;};

	void setTc(double Tc){this->Tc=Tc;};
	void setLs(double Ls){this->Ls=Ls;};
	void setRs(double Rs){this->Rs=Rs;};
	void setParameters(double msun,double X,double Y, double Z,double Ts,double Ps);
	void setInnerBoundary();
	void setOuterBoundary();
	Phys shootOut();
	Phys shootIn();
	void setPerturbPc();
	void setPerturbTc();
	void setPerturbLs();
	void setPerturbRs();
	bool checkConvergence(Phys,Phys);
	double eps_pp(double rho,double X,double T);
	double eps_CNO(double rho,double X1,double X_CNO,double T);

	//Utility methodes
	void outNumberOfIterate();
	void outBoundary();
	void outDifference(Eigen::MatrixXd,Eigen::VectorXd,Eigen::VectorXd);
	void checkOverflow(long);
};

class EnergyGen{
	public:
		static double PP_RPN(double rho,double X,double T);
		static double PP_KIP(double rho,double X,double T);
		static double PP_AML(double rho,double X,double T);
		static double PP_REDUCED(double rho,double X,double T);
		static double CNO_KIP(double rho,double X1,double X_CNO,double T);
	private:
		EnergyGen(){};
		~EnergyGen(){};
};