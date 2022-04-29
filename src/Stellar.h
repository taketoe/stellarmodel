#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <array>
#include "matplotlibcpp.h"

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
	std::array<double,Ndim> M,R,T,P,L,rho,kappa,epsilon;
	double dM,Pc,Tc,Ps,Ts,Rs,Ls;
	bool logOut;
	long numberOfIterate=0;
public:
	Stellar();
	Stellar(double msun);
	Stellar(double msun,double x,double y, double z);
	~Stellar();
	e_state calc();
	void setLog(bool);
	void plot();
	Phys getPhys(long);
	void getResult();
	void setPc(double);
private:
	double getPc();
	double getTc();
	double getPs();
	double getTs();
	double getRs();
	double getLs();

	void setTc(double);
	void setLs(double);
	void setRs(double);
	void setParameters(double msun,double x,double y, double z);
	void setInnerBoundary();
	void setOuterBoundary();
	Phys shootOut();
	Phys shootIn();
	void setPerturbPc();
	void setPerturbTc();
	void setPerturbLs();
	void setPerturbRs();
	bool checkConvergence(Phys,Phys);

	//Utility methodes
	void outNumberOfIterate();
	void outBoundary();
	void outDifference(Eigen::MatrixXd,Eigen::VectorXd,Eigen::VectorXd);
	void checkOverflow(long);
};

class EnergyGeneration{
	public:
		EnergyGeneration(){};
		~EnergyGeneration(){};
		double PP_RPN(double rho,double X,double T){
			//Ref. R.P.Nelson,https://2019.qmplus.qmul.ac.uk/course/;
			//SPA7023U - STELLAR STRUCTURE AND EVOLUTION
			//rho:[kg/m^3]
			//X:Hydrogen fraction[-]
			//T:Temperature[K]
			double e_;//[J/s/kg]
			double e_pp0 = 2.38E-37;
			double alpha_= 4;
			//3.5 ~ alpha ~ 6
			//T_{6}=10[K],P=1.E15[Pa] : alpha_=4.08
			//T_{6}= 1[K],P=1.E15[Pa] : alpha_=3.6
			e_ = e_pp0 * rho * pow(X,2.) * pow(T,alpha_); 
			return e_;
		}
		double PP_KIP(double rho,double X,double T){
			//Ref. R.Kippenhahn,p.163
			//rho:[kg/m^3]
			//X:Hydrogen fraction[-]
			//T:Temperature[K]
			double T6 = T/1E6;
			double e_;//[J/s/kg]
			double g_11,f_11,phi;
			rho = rho * 1E-3;//[kg/m^3]->[g/cm^3]
			g_11 = 1 + 0.0123*pow(T6,1./3.) 
				+ 0.0109*pow(T6,2./3.)
				+ 0.0009*T6;//around 1.
			f_11 = exp(5.92E-3*1*1* pow(1. * rho / pow(T6/10,3),1./2.));//around 1.
			phi  = 1.;
			e_ = 2.38E6 * phi * f_11 * g_11 * rho  * pow(X,2.)
				* pow(T6,-2./3.) * exp(-33.80 / pow(T6,1./3.)); 
			e_ = e_ * 1E-4;//[erg/s/g]->[J/s/kg]
			return e_;
		}
		double PP_AML(double rho,double X,double T){
			//Ref. Aller and McLaughlin, 1965
			//rho:[kg/m^3]
			//X:Hydrogen fraction[-]
			//T:Temperature[K]
			double T6 = T/1E6;
			double e_;//[J/s/kg]
			double g_11,f,phi;
			rho = rho * 1E-3;//[kg/m^3]->[g/cm^3]
			g_11 = 1 + 0.0123*pow(T6,1./3.) 
				+ 0.0109*pow(T6,2./3.)
				+ 0.0009*T6;//around 1.
			f = 1+0.25*pow(rho,1./2.)*pow(T6,-3./2.);
			phi  = 1.;
			e_ = 2.38E6 * phi * f * g_11 * rho  * pow(X,2.)
				* pow(T6,-2./3.) * exp(-33.80 / pow(T6,1./3.)); 
			e_ = e_ * 1E-4;//[erg/s/g]->[J/s/kg]
			return e_;
		}
		double PP_REDUCED(double rho,double X,double T){
			//rho:[kg/m^3]
			//X:Hydrogen fraction[-]
			//T:Temperature[K]
			double T6 = T/1E6;
			double e_;//[J/s/kg]
			double g_11,f_11,phi;
			rho = rho * 1E-3;//[kg/m^3]->[g/cm^3]
			g_11 = 1.;
			f_11 = 1.;
			phi  = 1.;
			e_ = 2.38E6 * phi * f_11 * g_11 * rho  * pow(X,2.)
				* pow(T6,-2./3.) * exp(-33.80 / pow(T6,1./3.)); 
			e_ = e_ * 1E-4;//[erg/s/g]->[J/s/kg]
			return e_;
		}
		double CNO_KIP(double rho,double X1,double X_CNO,double T){
			//Ref. R.Kippenhahn,p.163
			//rho:[kg/m^3]
			//X1:Hydrogen fraction[-]
			//X_CNO:,also Z
			//T:Temperature[K]
			double T6 = T/1E6;
			double e_;//[J/s/kg]
			double g_141;
			rho = rho * 1E-3;//[kg/m^3]->[g/cm^3]
			g_141 = 1 + 0.0027*pow(T6,1./3.) 
				- 0.00778*pow(T6,2./3.)
				- 0.000149*T6;//around 1.
			e_ = 8.67E27 * g_141 * X_CNO * pow(X1,2.) * rho  
				* pow(T6,-2./3.) * exp(-152.28 / pow(T6,1./3.)); 
			e_ = e_ * 1E-4;//[erg/s/g]->[J/s/kg]
			return e_;
		}
};