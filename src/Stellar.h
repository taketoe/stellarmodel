#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <array>

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
const double Teff_sun=5778;		  // Solar efective temperature [K]
const double con_fact=0.3;	// feedback ratio

const double tolerance=5.E-3;// convergence criteria
const long Ndim = 100000;// Number of grid points

const long maxIterate = 2000;
const double pertubRatio = 0.01;
//const double innerMassFract = 1.E-6;
//const double innerMassFract = 1.E-5;
const double innerMassFract = 1.E-6;

//Hear is sucesful parameters set

enum e_state{
	notconverge,
	converge,
	overshoot,
	overflow,
	steptoolarge,
	nonphysical,
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
	double dM,Mc,Pc,Tc,Ps,Ts,Rs,Ls;
	bool logOut;
	long numberOfIterate=0;
public:
	Stellar();
	Stellar(double msun,double X,double Y,double Z,double Ts,double Ps);
	~Stellar();
	e_state calc();
	e_state calcSO();
	void setLog(bool);
	void plot(){};
	Phys getPhys(long);
	void getResult();
	void setPc(double Pc){this->Pc=Pc;}
	void setRs(double Rs){this->Rs=Rs;};
	void setTc(double Tc){this->Tc=Tc;};
	void setLs(double Ls){this->Ls=Ls;};
	double getPc(){return Pc;};
	double getTc(){return Tc;};
	std::vector<double> getM(){return M;}
	std::vector<double> getR(){return R;}
	std::vector<double> getP(){return P;}
	std::vector<double> getT(){return T;}
	std::vector<double> getL(){return L;}
	std::vector<double> getRho(){return rho;}
	std::vector<double> getKappa(){return kappa;}
	std::vector<double> getEpsilon(){return epsilon;}
private:
	double getPs(){return Ps;};
	double getTs(){return Ts;};
	double getRs(){return Rs;};
	double getLs(){return Ls;};

	void init();
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
	double opacity(double rho,double T);
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

class MainSequence {
	public:
		MainSequence(){}
		MainSequence(double R_M1,double L_M1,double Teff_M1){
		 	this->R_M1 = R_M1;
		 	this->L_M1 = L_M1;
		 	this->Teff_M1 = Teff_M1;
		 };
		~MainSequence(){};
		double Radius(double mstar_fact){
			return(R_M1 * pow(mstar_fact,xsi));
		}
		double Luminosity(double mstar_fact){
			return(L_M1 * pow(mstar_fact,eta));
		}
		double Teff(double mstar_fact){
			return(Teff_M1 * pow(mstar_fact,eta/zeta));
		}
	private:
		double xsi  = 0.57;
		double eta  = 3.2;
		double zeta = 4./(1-2.*xsi/eta);
		double R_M1,L_M1,Teff_M1; 
		double X;
		double Y;
		double Z;
};

class ZAMS{
	public:
		static double Radius(double mstar_fact){
			MainSequence MS = MainSequence(Rsun,Lsun,Teff_sun);
			return(MS.Radius(mstar_fact));
		}
	private:
		ZAMS(){}
		~ZAMS(){}
};

class HeliumMS{
	public:
		HeliumMS(){
			/*
	 		this->R_M1 = R_M1;
		 	this->L_M1 = L_M1;
		 	this->Teff_M1 = Teff_M1;
			*/
		}
		~HeliumMS(){}
};

class CarbonMS{
	public:
		CarbonMS(){
			/*
	 		this->R_M1 = R_M1;
		 	this->L_M1 = L_M1;
		 	this->Teff_M1 = Teff_M1;
			*/
		}
		~CarbonMS(){}
};


class Opacity{
	public:
		static double KramApprox(double rho,double X,double Z,double T){
			//[m^2/kg]
			return 4.3E24*Z*(1+X)*rho*pow(T,-3.5);
		}
		static double Sum(double rho,double X,double Z,double T){
			double kappa = kappa_sc(X)
				+kappa_ff(rho,X,Z,T)
				+kappa_bf(rho,X,Z,T)
				+kappa_Hminus(rho,X,Z,T);
				return(kappa);
		}
		static double kappa_sc(double X){
			double kappa = 0.02*(1+X);//m^2/kg
			return(kappa);
		}
		static double kappa_ff(double rho,double X,double Z,double T){
			//double g_ff = 0.001;
			double kappa = 4.E21*(1.+X)*(1.-Z)*rho*pow(T,-7./2.); //m^2/kg
			return(kappa);
		}
		static double kappa_bf(double rho,double X,double Z,double T){
			double g_bf = 2.82;
			double t = 0.00001;
			double kappa = 4.E24*Z*(1+X)*rho*pow(T,-7./2.);
			return(kappa);
		}
		static double kappa_Hminus(double rho,double X,double Z,double T){
			double kappa = 0.;
			if(T>3E3 && T<6E3 && 1E-7>rho && 1E-2<rho
				&& 0.67>X && 0,73<X && 0.001>Z && 0.003<Z){
				kappa = 7.9E-34*Z/2E-2*pow(rho,1./2.)*pow(T,9);
			}
			return(kappa);
		}
	private:
		Opacity(){}
		~Opacity(){};
};