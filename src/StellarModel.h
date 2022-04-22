#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <array>
using namespace std;

#include <Eigen/Dense>
#include <Eigen/LU>
using namespace Eigen;

const double Pi=3.141592;
const double G=6.67E-11;      // Gravitational constant
const double c=3.E8;          // Speed of light
const double kB=1.38E-23;     // Boltzmann constant
const double m_H=1.67E-27;    // Mass of hydrogen atom
const double sigma_SB=5.67E-8;// Stefan-Boltzmann 
const double a=7.5657E-16;    // Radiation constant
const double gamma_c=5./3.;     // Ratio of specific heats
const double Msun=1.989E30;   // Solar mass
const double Rsun=6.957E8;    // Solar radius
const double Lsun=3.828E26;   // Solar luminosity

const double con_fact=0.3;	//????

const double tolerance=5.E-3;// convergence criteria
const int Ndim=10000;// Number of grid points

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
	const Phys operator/=(const Phys&) const;
	Phys();
	Phys(double M,double R,double P,double T,double L,double Rho,double Kappa,double Epsilon);
	~Phys();
	void Out();
};


class Stellar
{
private:
	double mstar_fact;
	double Rstar,Mstar;
	double X,Y,Z,mu;
	double epsilon_pp,epsilon_CNO;
	double kappa0;
	std::array<double,Ndim> M,R,T,P,L,rho,kappa,epsilon;
	double dM,Pc,Tc,Ps,Ts,Rs,Ls;
	bool logOut;
public:
	Stellar();
	Stellar(double msun);
	Stellar(double msun,double x,double y, double z);
	~Stellar();
	double getPc();
	double getTc();
	double getPs();
	double getTs();
	double getRs();
	double getLs();
	Phys getPhys(int);
	void setPc(double);
	void setTc(double);
	void setLs(double);
	void setRs(double);
	void Iterate();
	void Setup(double msun,double x,double y, double z);
	void SetInnerBoundary();
	Phys ShootIn();
	void SetOuterBoundary();
	Phys ShootOut();
	void SetPerturbPc();
	void SetPerturbTc();
	void SetPerturbLs();
	void SetPerturbRs();
	void SetLogOut(bool);
	bool CheckConverge(Phys,Phys);
	void Calc();
	void getResult();
	void Plot();
};

class StellarModel
{
private:
public:
	StellarModel(/* args */);
	~StellarModel();
	void Calc();
	bool CheckConverge(Phys,Phys);
	void Plot();
};
