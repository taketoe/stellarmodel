
#include "StellarModel.h"

Phys::Phys(){}

Phys::Phys(double m,double r,double p,double t,double l,double rho,double kappa,double epsilon)
{
	M=m;R=r;P=p;T=t;L=l;Rho=rho;Kappa=kappa;Epsilon=epsilon;
}

double Phys::getM(){return M;}
double Phys::getR(){return R;}
double Phys::getP(){return P;}
double Phys::getT(){return T;}
double Phys::getL(){return L;}
double Phys::getRho(){return Rho;}
double Phys::getKappa(){return Kappa;}
double Phys::getEpsilon(){return Epsilon;}

Phys::~Phys()
{
}

const Phys Phys::operator-(const Phys &phys) const
{
	double r = R - phys.R;
	double m = M - phys.M;
	double t = T - phys.T;
	double p = P - phys.P;
	double l = L - phys.L;
	double rho = Rho - phys.Rho;
	double kappa = Kappa - phys.Kappa;
	double epsilon = Epsilon - phys.Epsilon;
	Phys tmp(m,r,p,t,l,rho,kappa,epsilon);
	return tmp;
}

const Phys Phys::operator/=(const Phys &phys) const
{
	double r = R / phys.R;
	double m = M / phys.M;
	double t = T / phys.T;
	double p = P / phys.P;
	double l = L / phys.L;
	double rho = Rho / phys.Rho;
	double kappa = Kappa / phys.Kappa;
	double epsilon = Epsilon / phys.Epsilon;
	Phys tmp(m,r,t,p,l,rho,kappa,epsilon);
	return tmp;
}

void Phys::Out(){
	cout << "M=" << getM();
	cout << " R=" << getR();
	cout << " P=" << getP();
	cout << " T=" << getT();
	cout << " L=" << getL();
	cout << " rho=" << getRho();
	cout << " kappa=" << getKappa();
	cout << " epsilon=" << getEpsilon() << endl;
}

Stellar::Stellar():logOut(false){
	double mstar_fact = 1;
	X=0.70;
	Y=0.28;
	Z=0.02;
	Setup(mstar_fact,X,Y,Z);
}
Stellar::Stellar(double mstarFact):logOut(false){
	X=0.70;
	Y=0.28;
	Z=0.02;
	Setup(mstarFact,X,Y,Z);
}
Stellar::Stellar(double mstarFact,double x,double y, double z):logOut(false){
	Setup(mstarFact,x,y,z);
}

void Stellar::Setup(double mstarFact,double x,double y, double z){
	mstar_fact = mstarFact;
	Mstar=double(mstar_fact)*Msun;
//	cout << Msun << endl;
	if (Mstar < 1.66*Msun){
		Rstar=Rsun*1.06*pow((Mstar/Msun),0.945);
	}
	else{
		Rstar=Rsun*1.33*pow((Mstar/Msun),0.555);
	}
	X=x;Y=y;Z=z;
	//cout << X << endl;
	mu=1./(2.*X+0.75*Y+0.5*Z);
	epsilon_pp=2.6E-37*pow(X,2);
	//epsilon_CNO=7.9E-118*X*Z;//???
	//epsilon_CNO=8.24E-24*X*Z;//???
	epsilon_CNO=0.;
	kappa0=4.3E24*Z*(X+1);

	M[0]=1.E-06*Mstar;        		// Mass of inner cell
	dM=(Mstar-M[0])/(Ndim-1); 		// Mass of each mass shell
	Pc=1.E15;                 		// Initial guess for central pressure
	Tc=1.E7*pow((Mstar/Msun),0.5); 	// Initial guess for centraltemperature

	Ts=1000.;                 		// Surface temperature; input value
	Ps=1.E8;                  		// Surface pressure; input value
	Rs=1.*Rstar;              		// Stellar radius; input value
	Ls=Lsun*pow((Mstar/Msun),3.); 	// Luminosity; calculated value
}

double Stellar::getPc(){return Pc;}
double Stellar::getTc(){return Tc;}
double Stellar::getPs(){return Ps;}
double Stellar::getTs(){return Ts;}
double Stellar::getRs(){return Rs;}
double Stellar::getLs(){return Ls;}

Phys Stellar::getPhys(int index){
	return Phys(M[index],R[index],P[index],T[index],L[index],rho[index],kappa[index],epsilon[index]);
}

void Stellar::setPc(double pc){Pc=pc;}
void Stellar::setTc(double tc){Tc=tc;}
void Stellar::setLs(double ls){Ls=ls;}
void Stellar::setRs(double rs){Rs=rs;}



void Stellar::SetInnerBoundary(){
	P[0]=Pc;
	T[0]=Tc;
	rho[0]=P[0]*mu*m_H/(kB*T[0]);
	R[0]=pow(3.*M[0]/(4.*Pi*rho[0]),1./3.);
	epsilon[0]=epsilon_pp*rho[0]*pow(T[0],4.)+epsilon_CNO*rho[0]*pow(T[0]/1.E6,16);
	L[0]=epsilon[0]*M[0];
	kappa[0]=kappa0*rho[0]*pow(T[0],-3.5);
}

Phys Stellar::ShootIn(){
	double nabla_rad;
	for(int j=Ndim-2;j>=int(Ndim/2)-2;j--){
			R[j] = R[j+1] - dM/(4.*Pi*pow(R[j+1],2)*rho[j+1]);
			P[j] = P[j+1] + dM*G*M[j+1]/(4.*Pi*pow(R[j+1],4));
			L[j] = L[j+1] - dM*epsilon[j+1];
			//cout << "L[j]:" << L[j] << " " << L[j+1]<< endl;
			nabla_rad = (3.*kappa[j+1]*L[j+1]*P[j+1]/
				(16.*Pi*a*c*pow(T[j+1],4)*G*M[j+1]));
			if(nabla_rad < (gamma_c-1.)/gamma_c){
				T[j] = T[j+1] + (dM*3.*kappa[j+1]*L[j+1]/
				(16.*Pi*a*c*pow(R[j+1],2)*pow(T[j+1],3))/(4.*Pi*pow(R[j+1],2)));
			}
			else{
				T[j] = T[j+1] - (dM*(gamma_c-1.)/gamma_c*
				T[j+1]/P[j+1]*(P[j+1] - P[j])/dM);
			}
			//********************************************************
			// Check that we don't overshoot so that either the radius 
			// or luminosity become negative
			//********************************************************
			if(R[j] <= 0. or L[j] <= 0.){
				cout<<"R[j] <= 0. or L[j] <= 0."<<endl;
				cout << "j:"<<j<<" R[j]:"<<R[j]<<" L[j]:"<<L[j]<<endl;
				exit(0);
			}
			M[j] = M[j+1] - dM;
			rho[j] = P[j]*mu*m_H/(kB*T[j]);
			epsilon[j] = (epsilon_pp*rho[j]*pow(T[j],4.)+
					epsilon_CNO*rho[j]*pow((T[j]/1.E6),16));
			kappa[j] = kappa0*rho[j]*pow(T[j],-3.5);
	}
	int mid = int(Ndim/2)-1;
	return Phys(M[mid],R[mid],P[mid],T[mid],L[mid],rho[mid],kappa[mid],epsilon[mid]);
}

void Stellar::SetOuterBoundary(){
	M[Ndim-1]=Mstar;
	L[Ndim-1]=Ls;
	R[Ndim-1]=Rs;
	T[Ndim-1]=Ts;
	P[Ndim-1]=Ps;
	rho[Ndim-1]=P[Ndim-1]*mu*m_H/(kB*T[Ndim-1]);
	kappa[Ndim-1]=kappa0*rho[Ndim-1]*pow(T[Ndim-1],-3.5);
	epsilon[Ndim-1]=epsilon_pp*rho[Ndim-1]*pow(T[Ndim-1],4.)
		+epsilon_CNO*rho[Ndim-1]*pow(T[Ndim-1]/1.E6,16);
}

Phys Stellar::ShootOut(){
	double nabla_rad;
	for(int i=1;i<=int(Ndim/2);i++){
		R[i] = R[i-1] + dM/(4.*Pi*pow(R[i-1],2)*rho[i-1]);
		P[i] = P[i-1] - dM*G*M[i-1]/(4.*Pi*pow(R[i-1],4));
		L[i] = L[i-1] + dM*epsilon[i-1];
		nabla_rad = (3.*kappa[i-1]*L[i-1]*P[i-1]/
			(16.*Pi*a*c*pow(T[i-1],4)*G*M[i-1]));
		//********************************************
		// Radiative or adiabatic temperature gradient
		//********************************************
		if(nabla_rad < (gamma_c-1)/gamma_c){
			T[i] = T[i-1] - (dM*3.*kappa[i-1]*L[i-1]/
			(16.*Pi*a*c*pow(R[i-1],2)*pow(T[i-1],3))/(4.*Pi*pow(R[i-1],2)));
		}	
		else{
			T[i] = T[i-1] + (dM*(gamma_c-1.)/gamma_c*
					T[i-1]/P[i-1]*(P[i]-P[i-1])/dM);
		}
		//******************************************************
		// Check we don't overshoot the temperature or pressure 
		// such that one of them becomes negative
		//******************************************************
		if(T[i] <= 0. or P[i] <= 0.){
			cout<<"T[i] <= 0. or P[i] <= 0."<<endl;
			cout << "i:"<<i<<" T[i]:"<<T[i]<<" P[i]:"<<P[i]<<endl;
			exit(0);
		}
		M[i] = M[i-1] + dM;
		rho[i] = P[i]*mu*m_H/(kB*T[i]);
		epsilon[i] = (epsilon_pp*rho[i]*pow(T[i],4.)+
		    	      epsilon_CNO*rho[i]*pow((T[i]/1.E6),16));
		kappa[i] = kappa0*rho[i]*pow(T[i],(-3.5));
	}
	int mid = int(Ndim/2)-1;
	return Phys(M[mid],R[mid],P[mid],T[mid],L[mid],rho[mid],kappa[mid],epsilon[mid]);
}

void Stellar::SetPerturbPc(){
	P[0]=1.01*Pc;
	T[0]=Tc;
	rho[0]=P[0]*mu*m_H/(kB*T[0]);
	R[0]=pow(3.*M[0]/(4.*Pi*rho[0]),1./3.);
	epsilon[0]=epsilon_pp*rho[0]*pow(T[0],4.)+epsilon_CNO*rho[0]*pow(T[0]/1.E6,16);
	L[0]=epsilon[0]*M[0];
	kappa[0]=kappa0*rho[0]*pow(T[0],-3.5);	
}
void Stellar::SetPerturbTc(){
	P[0]=Pc;
	T[0]=1.01*Tc;
	rho[0]=P[0]*mu*m_H/(kB*T[0]);
	R[0]=pow(3.*M[0]/(4.*Pi*rho[0]),1./3.);
	epsilon[0]=epsilon_pp*rho[0]*pow(T[0],4.)+epsilon_CNO*rho[0]*pow(T[0]/1.E6,16);
	L[0]=epsilon[0]*M[0];
	kappa[0]=kappa0*rho[0]*pow(T[0],-3.5);
}
void Stellar::SetPerturbLs(){
	M[Ndim-1]=Mstar;
	L[Ndim-1]=Ls*1.01;
	R[Ndim-1]=Rs;
	T[Ndim-1]=Ts;
	P[Ndim-1]=Ps;
	rho[Ndim-1]=P[Ndim-1]*mu*m_H/(kB*T[Ndim-1]);
	kappa[Ndim-1]=kappa0*rho[Ndim-1]*pow(T[Ndim-1],-3.5);
	epsilon[Ndim-1]=epsilon_pp*rho[Ndim-1]*pow(T[Ndim-1],4.)
		+epsilon_CNO*rho[Ndim-1]*pow(T[Ndim-1]/1.E6,16);
}
void Stellar::SetPerturbRs(){
	M[Ndim-1]=Mstar;
	L[Ndim-1]=Ls;
	R[Ndim-1]=Rs*1.01;
	T[Ndim-1]=Ts;
	P[Ndim-1]=Ps;
	rho[Ndim-1]=P[Ndim-1]*mu*m_H/(kB*T[Ndim-1]);
	kappa[Ndim-1]=kappa0*rho[Ndim-1]*pow(T[Ndim-1],-3.5);
	epsilon[Ndim-1]=epsilon_pp*rho[Ndim-1]*pow(T[Ndim-1],4.)
		+epsilon_CNO*rho[Ndim-1]*pow(T[Ndim-1]/1.E6,16);


}

Stellar::~Stellar(){
}

void Stellar::Calc(){
	double Rdiff1,Pdiff1,Tdiff1,Ldiff1;
	double Rdiff2,Pdiff2,Tdiff2,Ldiff2;
	double Rdiff3,Pdiff3,Tdiff3,Ldiff3;
	double Rdiff4,Pdiff4,Tdiff4,Ldiff4;
	double d_deltaR_dPc,d_deltaP_dPc,d_deltaT_dPc,d_deltaL_dPc;
	double d_deltaR_dTc,d_deltaP_dTc,d_deltaT_dTc,d_deltaL_dTc;
	double d_deltaR_dLs,d_deltaP_dLs,d_deltaT_dLs,d_deltaL_dLs;
	double d_deltaR_dRs,d_deltaP_dRs,d_deltaT_dRs,d_deltaL_dRs;
	Phys physMidOp1,physMidOp2,physMidIp1,physMidIp2;
	Phys physDiff0,physDiff1,physDiff2,physDiff3,physDiff4;
	Phys phys_dPc,phys_dTc,phys_dRs,phys_dLs;
	bool converge = false;
	while(!converge){
		SetOuterBoundary();
		SetInnerBoundary();
		Phys physMidI=ShootIn();//OK
		Phys physMidO=ShootOut();//OK
		// physMidI.Out();
		// physMidO.Out();
		converge = CheckConverge(physMidI,physMidO);
		SetPerturbPc();
		physMidOp1=ShootOut();
		SetPerturbTc();
		physMidOp2=ShootOut();
		SetPerturbLs();
		physMidIp1=ShootIn();
		SetPerturbRs();
		physMidIp2=ShootIn();
		// physMidOp1.Out();	
		// physMidOp2.Out();	
		// physMidIp1.Out();	
		// physMidIp2.Out();	

		physDiff0 = physMidO - physMidI;
		physDiff1 = physMidOp1 - physMidI;//Rmid1-R2mid0
		physDiff2 = physMidOp2 - physMidI;//Rmid2-R2mid0
		physDiff3 = physMidO - physMidIp1;//Rmid0-R2mid1
		physDiff4 = physMidO - physMidIp2;//Rmid0-R2mid2

		phys_dPc = physDiff1 - physDiff0;
		phys_dTc = physDiff2 - physDiff0;
		phys_dLs = physDiff3 - physDiff0;
		phys_dRs = physDiff4 - physDiff0;

		d_deltaR_dPc = phys_dPc.getR() / (0.01*getPc());
		d_deltaR_dTc = phys_dTc.getR() / (0.01*getTc());
		d_deltaR_dLs = phys_dLs.getR() / (0.01*getLs());
		d_deltaR_dRs = phys_dRs.getR() / (0.01*getRs());

		d_deltaP_dPc = phys_dPc.getP() / (0.01*getPc());
		d_deltaP_dTc = phys_dTc.getP() / (0.01*getTc());
		d_deltaP_dLs = phys_dLs.getP() / (0.01*getLs());
		d_deltaP_dRs = phys_dRs.getP() / (0.01*getRs());

		d_deltaT_dPc = phys_dPc.getT() / (0.01*getPc());
		d_deltaT_dTc = phys_dTc.getT() / (0.01*getTc());
		d_deltaT_dLs = phys_dLs.getT() / (0.01*getLs());
		d_deltaT_dRs = phys_dRs.getT() / (0.01*getRs());

		d_deltaL_dPc = phys_dPc.getL() / (0.01*getPc());
		d_deltaL_dTc = phys_dTc.getL() / (0.01*getTc());
		d_deltaL_dLs = phys_dLs.getL() / (0.01*getLs());
		d_deltaL_dRs = phys_dRs.getL() / (0.01*getRs());

		MatrixXd A(4,4);
		A(0,0)=d_deltaR_dPc;
		A(0,1)=d_deltaR_dTc;
		A(0,2)=d_deltaR_dLs;
		A(0,3)=d_deltaR_dRs;
		A(1,0)=d_deltaP_dPc;
		A(1,1)=d_deltaP_dTc;
		A(1,2)=d_deltaP_dLs;
		A(1,3)=d_deltaP_dRs;
		A(2,0)=d_deltaT_dPc;
		A(2,1)=d_deltaT_dTc;
		A(2,2)=d_deltaT_dLs;
		A(2,3)=d_deltaT_dRs;
		A(3,0)=d_deltaL_dPc;
		A(3,1)=d_deltaL_dTc;
		A(3,2)=d_deltaL_dLs;
		A(3,3)=d_deltaL_dRs;

		VectorXd x(4),y(4);
		y(0)=con_fact*physDiff0.getR();
		y(1)=con_fact*physDiff0.getP();
		y(2)=con_fact*physDiff0.getT();
		y(3)=con_fact*physDiff0.getL();

		x = A.inverse()*y;

		// cout << "y" << endl;
		// cout << y.transpose() << endl;
		// cout << "A="<<endl;
		// cout << A << endl;
		// cout << "A^-1="<<endl;
		// cout << A.inverse() << endl;
		// cout << (A.inverse()*y).transpose() << endl;

		setPc(getPc()-x[0]);
		setTc(getTc()-x[1]);
		setLs(getLs()-x[2]);
		setRs(getRs()-x[3]);
	}

	Phys physSurface = getPhys(Ndim-2);
	Phys physCenter  = getPhys(0);
	//physSurface.Out();
	//physCenter.Out();
}

bool Stellar::CheckConverge(Phys phys1,Phys phys2){
	double Rdiff, Pdiff,Tdiff,Ldiff;
	double Rmidpoint, Pmidpoint,Tmidpoint,Lmidpoint;
	Rdiff = phys1.getR() - phys2.getR();
	Pdiff = phys1.getP() - phys2.getP();
	Tdiff = phys1.getT() - phys2.getT();
	Ldiff = phys1.getL() - phys2.getL();
	Rmidpoint = (phys1.getR() + phys2.getR())/2.;
	Pmidpoint = (phys1.getP() + phys2.getP())/2.;
	Tmidpoint = (phys1.getT() + phys2.getP())/2.;
	Lmidpoint = (phys1.getL() + phys2.getL())/2.;
	if(logOut){
		cout <<"**************************************"<<endl;
		cout <<"Convergence tolerance parameter="<<tolerance<<endl;
		cout <<"State of convergence shown below:"<<endl;
		cout <<" Rdiff/Rmidpoint="<<abs(Rdiff)/Rmidpoint;
		cout <<" Pdiff/Pmidpoint="<<abs(Pdiff)/Pmidpoint;
		cout <<" Tdiff/Tmidpoint="<<abs(Tdiff)/Tmidpoint;
		cout <<" Ldiff)/Lmidpoint="<<abs(Ldiff)/Lmidpoint<<endl;
		cout <<"**************************************"<<endl;
	}
	if(abs(Rdiff)/Rmidpoint < tolerance &&
		abs(Pdiff)/Pmidpoint < tolerance &&
		abs(Tdiff)/Tmidpoint < tolerance &&
		abs(Ldiff)/Lmidpoint < tolerance){
					if(logOut){cout << "Convergence achieved"<<endl;}
					return true;
	}
	return false;
}

void Stellar::SetLogOut(bool sw){logOut=sw;}


void Stellar::getResult(){
	Phys physSurface = getPhys(Ndim-2);
	Phys physCenter  = getPhys(0);

	cout << "Stellar Center" << endl;
	physCenter.Out(); 
	cout << "Stellar Surface" << endl;
	physSurface.Out(); 
}

void Stellar::Plot(){}