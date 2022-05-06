
#include "Stellar.h"
using namespace std;
using namespace Eigen;

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

void Phys::Out(){
	cout << " M=" << getM();
	cout << " R=" << getR();
	cout << " P=" << getP();
	cout << " T=" << getT();
	cout << " L=" << getL()<<endl;
	cout << " rho=" << getRho();
	cout << " kappa=" << getKappa();
	cout << " epsilon=" << getEpsilon() << endl;
}

Stellar::Stellar():logOut(false){
	double mstar_fact = 1;
	double X=0.70;
	double Y=0.28;
	double Z=0.02;
	double Ts=1.E3;
	double Ps=1.E8;
	setParameters(mstar_fact,X,Y,Z,Ts,Ps);
}

Stellar::Stellar(double mstarFact,double X,double Y, double Z, double Ts, double Ps):logOut(false){
	setParameters(mstarFact,X,Y,Z,Ts,Ps);
}

void Stellar::setParameters(double mstar_fact,double X,double Y, double Z, double Ts, double Ps){
	this->mstar_fact = mstar_fact;
	this->X=X;
	this->Y=Y;
	this->Z=Z;
	this->Ts=Ts;                 		// Surface temperature; input value
	this->Ps=Ps;                  		// Surface pressure; input value

	Mstar=double(mstar_fact)*Msun;
	if (Mstar < 1.66*Msun){
		Rstar=Rsun*1.06*pow((Mstar/Msun),0.945);
	}
	else{
		Rstar=Rsun*1.33*pow((Mstar/Msun),0.555);
	}
	mu=1./(2.*X+0.75*Y+0.5*Z);
	epsilon_pp=2.6E-37*pow(X,2);//Constant for PP-chain energy generation
	epsilon_CNO=0.;//Constant for CNO-cycle energy generation
	//epsilon_CNO=7.9E-118*X*Z;//
	//epsilon_CNO=8.24E-24*X*Z;//
	kappa0=4.3E24*Z*(X+1);

	Mc=innerMassFract*Mstar;        		// Mass of inner cell
	dM=(Mstar-Mc)/(double(Ndim-1)); 		// Mass of each mass shell
//	Pc=1.E15;                 		// Initial guess for central pressure;pp_RPN
	Pc=1.E16;                 		// Initial guess for central pressure;pp_KIP
	Tc=1.E7*pow((Mstar/Msun),0.5); 	// Initial guess for centraltemperature;RPN
//	Tc=1.6E7*pow((Mstar/Msun),0.5); 	// Initial guess for centraltemperature

	Rs=1.*Rstar;              		// Stellar radius; input value
	Ls=Lsun*pow((Mstar/Msun),3.2); 	// Luminosity; calculated value

}

Phys Stellar::getPhys(long index){
	return Phys(M[index],R[index],P[index],T[index],L[index],rho[index],kappa[index],epsilon[index]);
}

void Stellar::init(){
	for(int i=1;i<Ndim-1;i++){
		M[i]=0.;R[i]=0.;P[i]=0.;T[i]=0.;
		rho[i]=0.;kappa[i]=0;epsilon[i]=0.;
	}
}

void Stellar::setInnerBoundary(){
	M[0]=Mc;
	P[0]=Pc;
	T[0]=Tc;
	rho[0]=P[0]*mu*m_H/(kB*T[0]);
	R[0]=pow(3.*M[0]/(4.*Pi*rho[0]),1./3.);
	epsilon[0]=eps_pp(rho[0],X,T[0])+eps_CNO(rho[0],X,Z,T[0]);
	L[0]=epsilon[0]*M[0];
	kappa[0]=opacity(rho[0],T[0]);
}

void Stellar::setOuterBoundary(){
	M[Ndim-1]=Mstar;
	L[Ndim-1]=Ls;
	R[Ndim-1]=Rs;
	T[Ndim-1]=Ts;
	P[Ndim-1]=Ps;
	rho[Ndim-1]=P[Ndim-1]*mu*m_H/(kB*T[Ndim-1]);
	kappa[Ndim-1]=opacity(rho[Ndim-1],T[Ndim-1]);
	epsilon[Ndim-1]=eps_pp(rho[Ndim-1],X,T[Ndim-1])
		+eps_CNO(rho[Ndim-1],X,Z,T[Ndim-1]);
	if(!isfinite(epsilon[Ndim-1])){
		cout << " Exception occured: epsilon[Ndim-1] isn't finite" << endl;
		cout << "eps_pp:" << eps_pp(rho[Ndim-1],X,T[Ndim-1]) << endl;
		cout << "eps_CNO:" << eps_CNO(rho[Ndim-1],X,Z,T[Ndim-1]) << endl;
		cout << "rho[Ndim-1]:" << rho[Ndim-1] << " T[Ndim-1]:" << T[Ndim-1] << endl;
		throw e_state::overflow;
	}

	//*****************************************
	// Check the size of the first radial step
	// to make sure it isn't too large.
	//*****************************************
	double dr_Rstar = dM/(4.*Pi*pow(R[Ndim-1],3)*rho[Ndim-1]);
	if(dr_Rstar/Rs > 1.E-02){
		if(logOut){
			cout << "Ending run since radial step size is too large"<<endl;
			cout << "dr_Rstar/Rs="<< dr_Rstar/Rs << " " << dr_Rstar << "  " << Rs<< endl;
		}
		throw e_state::steptoolarge;
	}
}

Phys Stellar::shootIn(){
	double nabla_rad;
	for(long j=Ndim-2;j>=long(Ndim/2)-1;j--){
			R[j] = R[j+1] - dM/(4.*Pi*pow(R[j+1],2)*rho[j+1]);
			P[j] = P[j+1] + dM*G*M[j+1]/(4.*Pi*pow(R[j+1],4));
			L[j] = L[j+1] - dM*epsilon[j+1];
			nabla_rad = (3.*kappa[j+1]*L[j+1]*P[j+1]/
				(16.*Pi*a*c*pow(T[j+1],4)*G*M[j+1]));
			if(nabla_rad < (gamma_c-1.)/gamma_c){
				T[j] = T[j+1] + (dM*3.*kappa[j+1]*L[j+1]/
				(16.*Pi*a*c*pow(R[j+1],2)*pow(T[j+1],3))/(4.*Pi*pow(R[j+1],2)));
			}else{
				T[j] = T[j+1] - (dM*(gamma_c-1.)/gamma_c*
				T[j+1]/P[j+1]*(P[j+1] - P[j])/dM);
			}
			//
			// Check overshoot 
			//
			if(R[j] <= 0. || L[j] <= 0.){
				if(logOut){
					cout << "Error@shootIn" << endl;
					outNumberOfIterate();
					cout << " j:"<<j<<" R[j]:"<<R[j]<<" L[j]:"<<L[j]
					<< " epislon[j+1]:" << epsilon[j+1] << endl;
					cout << " estemated Ts:"<< Ts << " Rs:"<< Rs <<endl;
				}
				throw e_state::overshoot;
			}
			M[j] = M[j+1] - dM;
			rho[j] = P[j]*mu*m_H/(kB*T[j]);
			epsilon[j] = eps_pp(rho[j],X,T[j])
					+eps_CNO(rho[j],X,Z,T[j]);
			kappa[j]=opacity(rho[j],T[j]);
	}
	long mid = long(Ndim/2)-1;
	checkOverflow(mid);
	return Phys(M[mid],R[mid],P[mid],T[mid],L[mid],rho[mid],kappa[mid],epsilon[mid]);
}



Phys Stellar::shootOut(){
	double nabla_rad;
	for(long i=1;i<=long(Ndim/2)-1;i++){
		R[i] = R[i-1] + dM/(4.*Pi*pow(R[i-1],2)*rho[i-1]);
		P[i] = P[i-1] - dM*G*M[i-1]/(4.*Pi*pow(R[i-1],4));
		L[i] = L[i-1] + dM*epsilon[i-1];
		nabla_rad = (3.*kappa[i-1]*L[i-1]*P[i-1]/
			(16.*Pi*a*c*pow(T[i-1],4)*G*M[i-1]));
		//********************************************
		// Radiative or adiabatic temperature gradient
		//********************************************
		if(nabla_rad < (gamma_c-1.)/gamma_c){
			T[i] = T[i-1] - dM*3.*kappa[i-1]*L[i-1]/
			(4*a*c*pow(T[i-1],3.))/(16.*pow(Pi,2.)*pow(R[i-1],4.));
		}
		else{
			T[i] = T[i-1] + (gamma_c-1.)/gamma_c*
					T[i-1]/P[i-1]*(P[i]-P[i-1]);
		}
		//
		// Check overshoot
		// 
		if(T[i] <= 0. || P[i] <= 0.){
			if(logOut){
				cout << "Error@shootOut" << endl;
				if(nabla_rad < (gamma_c-1.)/gamma_c){cout << " radiative transport" << endl;}
				else{cout << " convective transport" << endl;}
				cout << " i:"<<i<<" T[i]:"<<T[i]<<" P[i]:"<<P[i]<<endl;
				cout << " estemated Tc:"<< Tc << " Pc:"<< Pc <<endl;
			}
			throw e_state::overshoot;
		}
		M[i] = M[i-1] + dM;
		rho[i] = P[i]*mu*m_H/(kB*T[i]);
		epsilon[i] = eps_pp(rho[i],X,T[i])
					+eps_CNO(rho[i],X,Z,T[i]);
		kappa[i]=opacity(rho[i],T[i]);
	}
	long mid = long(Ndim/2)-1;
	checkOverflow(mid);
	return Phys(M[mid],R[mid],P[mid],T[mid],L[mid],rho[mid],kappa[mid],epsilon[mid]);
}

void Stellar::setPerturbPc(){
	P[0]=Pc+pertubRatio*Pc;
	T[0]=Tc;
	rho[0]=P[0]*mu*m_H/(kB*T[0]);
	R[0]=pow(3.*M[0]/(4.*Pi*rho[0]),1./3.);
	epsilon[0]=eps_pp(rho[0],X,T[0])+eps_CNO(rho[0],X,Z,T[0]);
	L[0]=epsilon[0]*M[0];
	kappa[0]=opacity(rho[0],T[0]);
}
void Stellar::setPerturbTc(){
	P[0]=Pc;
	T[0]=Tc+pertubRatio*Tc;
	rho[0]=P[0]*mu*m_H/(kB*T[0]);
	R[0]=pow(3.*M[0]/(4.*Pi*rho[0]),1./3.);
	epsilon[0]=eps_pp(rho[0],X,T[0])+eps_CNO(rho[0],X,Z,T[0]);
	L[0]=epsilon[0]*M[0];
	kappa[0]=opacity(rho[0],T[0]);
}

void Stellar::setPerturbLs(){
	M[Ndim-1]=Mstar;
	L[Ndim-1]=Ls+pertubRatio*Ls;
	R[Ndim-1]=Rs;
	T[Ndim-1]=Ts;
	P[Ndim-1]=Ps;
	rho[Ndim-1]=P[Ndim-1]*mu*m_H/(kB*T[Ndim-1]);
	kappa[Ndim-1]=opacity(rho[Ndim-1],T[Ndim-1]);
	epsilon[Ndim-1]=eps_pp(rho[Ndim-1],X,T[Ndim-1])
		+eps_CNO(rho[Ndim-1],X,Z,T[Ndim-1]);
}
void Stellar::setPerturbRs(){
	M[Ndim-1]=Mstar;
	L[Ndim-1]=Ls;
	R[Ndim-1]=Rs+pertubRatio*Rs;
	T[Ndim-1]=Ts;
	P[Ndim-1]=Ps;
	rho[Ndim-1]=P[Ndim-1]*mu*m_H/(kB*T[Ndim-1]);
	kappa[Ndim-1]=opacity(rho[Ndim-1],T[Ndim-1]);
	epsilon[Ndim-1]=eps_pp(rho[Ndim-1],X,T[Ndim-1])
		+eps_CNO(rho[Ndim-1],X,Z,T[Ndim-1]);
}

Stellar::~Stellar(){
}

e_state Stellar::calc(){
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
	Eigen::MatrixXd A(4,4);
	Eigen::VectorXd x(4),y(4);

	double Ps = getPs();
	double Ps_min = Ps;
	double Ps_max = Ps*10.;
	double dPs = Ps;
	int i=0;
	Ps=Ps_min;
	while((!converge) && Ps<=Ps_max){
		cout <<"i:"<<i<<" Ps:" << Ps << endl;
		setPs(Ps);
		try{
			while(!converge && numberOfIterate<maxIterate){
				init();
				setOuterBoundary();
				setInnerBoundary();
				Phys physMidO=shootOut();
			//	if(physMidO.state == overshoot){return e_state::overshoot;}
				Phys physMidI=shootIn();
				converge = checkConvergence(physMidI,physMidO);
				if(!converge){
					setPerturbPc();
					physMidOp1=shootOut();
					setPerturbTc();
					physMidOp2=shootOut();
					setPerturbLs();
					physMidIp1=shootIn();
					setPerturbRs();
					physMidIp2=shootIn();

					x = guessNextPoint(physMidO,physMidI,physMidOp1,physMidOp2,physMidIp1,physMidIp2);

					setPc(getPc()-double(x[0]));
					setTc(getTc()-double(x[1]));
					setLs(getLs()-double(x[2]));
					setRs(getRs()-double(x[3]));

					//Migrade result
					init();
					setOuterBoundary();
					setInnerBoundary();
					shootOut();
					shootIn();

					//outDifference(A,x,y);
				}
				//outNumberOfIterate();
				//outBoundary();		
				numberOfIterate++;
			}
		}
		catch(e_state &e){
			if(logOut){
				if(e == overflow){cout << "Exception occured:overflow" << endl;}
				if(e == overshoot){cout << "Exception occured:overshoot" << endl;}	
			}
		}
		Ps+=dPs;
		i++;
	}
	if(converge == true){
//		cout <<"converged"<< endl;
		return(e_state::converge);
	}
	return(e_state::notconverge);
}

Eigen::VectorXd Stellar::guessNextPoint(Phys physMidO,Phys physMidI,Phys physMidOp1,Phys physMidOp2,Phys physMidIp1,Phys physMidIp2){
	double Rdiff1,Pdiff1,Tdiff1,Ldiff1;
	double Rdiff2,Pdiff2,Tdiff2,Ldiff2;
	double Rdiff3,Pdiff3,Tdiff3,Ldiff3;
	double Rdiff4,Pdiff4,Tdiff4,Ldiff4;
	double d_deltaR_dPc,d_deltaP_dPc,d_deltaT_dPc,d_deltaL_dPc;
	double d_deltaR_dTc,d_deltaP_dTc,d_deltaT_dTc,d_deltaL_dTc;
	double d_deltaR_dLs,d_deltaP_dLs,d_deltaT_dLs,d_deltaL_dLs;
	double d_deltaR_dRs,d_deltaP_dRs,d_deltaT_dRs,d_deltaL_dRs;
	Phys physDiff0,physDiff1,physDiff2,physDiff3,physDiff4;
	Phys phys_dPc,phys_dTc,phys_dRs,phys_dLs;
	Eigen::MatrixXd A(4,4);
	Eigen::VectorXd x(4),y(4);

	physDiff0 = physMidO - physMidI;
	physDiff1 = physMidOp1 - physMidI;
	physDiff2 = physMidOp2 - physMidI;
	physDiff3 = physMidO - physMidIp1;
	physDiff4 = physMidO - physMidIp2;

	phys_dPc = physDiff1 - physDiff0;
	phys_dTc = physDiff2 - physDiff0;
	phys_dLs = physDiff3 - physDiff0;
	phys_dRs = physDiff4 - physDiff0;

	d_deltaR_dPc = phys_dPc.getR() / (pertubRatio*getPc());
	d_deltaR_dTc = phys_dTc.getR() / (pertubRatio*getTc());
	d_deltaR_dLs = phys_dLs.getR() / (pertubRatio*getLs());
	d_deltaR_dRs = phys_dRs.getR() / (pertubRatio*getRs());

	d_deltaP_dPc = phys_dPc.getP() / (pertubRatio*getPc());
	d_deltaP_dTc = phys_dTc.getP() / (pertubRatio*getTc());
	d_deltaP_dLs = phys_dLs.getP() / (pertubRatio*getLs());
	d_deltaP_dRs = phys_dRs.getP() / (pertubRatio*getRs());

	d_deltaT_dPc = phys_dPc.getT() / (pertubRatio*getPc());
	d_deltaT_dTc = phys_dTc.getT() / (pertubRatio*getTc());
	d_deltaT_dLs = phys_dLs.getT() / (pertubRatio*getLs());
	d_deltaT_dRs = phys_dRs.getT() / (pertubRatio*getRs());

	d_deltaL_dPc = phys_dPc.getL() / (pertubRatio*getPc());
	d_deltaL_dTc = phys_dTc.getL() / (pertubRatio*getTc());
	d_deltaL_dLs = phys_dLs.getL() / (pertubRatio*getLs());
	d_deltaL_dRs = phys_dRs.getL() / (pertubRatio*getRs());

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

	y(0)=con_fact*physDiff0.getR();
	y(1)=con_fact*physDiff0.getP();
	y(2)=con_fact*physDiff0.getT();
	y(3)=con_fact*physDiff0.getL();

	x = A.inverse()*y;

	if(getPc()-double(x[0])<0 ||
		getTc()-double(x[1])<0 ||
		getLs()-double(x[2])<0 ||
		getRs()-double(x[3])<0)
	{
		y(0)=0.01*con_fact*physDiff0.getR();
		y(1)=0.01*con_fact*physDiff0.getP();
		y(2)=0.01*con_fact*physDiff0.getT();
		y(3)=0.01*con_fact*physDiff0.getL();
		x = A.inverse()*y;
	}

	if(getPc()-double(x[0])<0){
		cout << "Exception occured:";
		cout << "non-physical value guessed";
		cout << " Pc:" << getPc()-double(x[0]) << endl;
		throw e_state::nonphysical;
	}
	if(getTc()-double(x[1])<0){
		cout << "Exception occured:";
		cout << "non-physical value guessed";
		cout << " Tc:" << getTc()-double(x[1]) << endl;
		throw e_state::nonphysical;
	}
	return x;
};

bool Stellar::checkConvergence(Phys phys1,Phys phys2){
	double Rdiff, Pdiff,Tdiff,Ldiff;
	double Rmidpoint, Pmidpoint,Tmidpoint,Lmidpoint;
	if(!isfinite(phys1.getT()) == true || !isfinite(phys2.getT())){
		phys1.Out();phys2.Out();
		throw e_state::overflow;}
	Rdiff = phys1.getR() - phys2.getR();
	Pdiff = phys1.getP() - phys2.getP();
	Tdiff = phys1.getT() - phys2.getT();
	Ldiff = phys1.getL() - phys2.getL();
	Rmidpoint = (phys1.getR() + phys2.getR())/2.;
	Pmidpoint = (phys1.getP() + phys2.getP())/2.;
	Tmidpoint = (phys1.getT() + phys2.getT())/2.;
	Lmidpoint = (phys1.getL() + phys2.getL())/2.;
	if(logOut){
		cout <<"=== Convergence tolerance parameter="<<tolerance<<endl;
		cout <<"=== Diffrence of midpoint between SHoot-in and Shoot-out"<<endl;
		cout <<"=== Diff.R:"<<abs(Rdiff)/Rmidpoint<<endl;
		cout <<"=== Diff.P:"<<abs(Pdiff)/Pmidpoint<<endl;
		cout <<"=== Diff.T:"<<abs(Tdiff)/Tmidpoint<<endl;
		cout <<"=== Diff.L:"<<abs(Ldiff)/Lmidpoint<<endl;
		//phys1.Out();
		//phys2.Out();
	}
	if(abs(Rdiff)/Rmidpoint < tolerance &&
		abs(Pdiff)/Pmidpoint < tolerance &&
		abs(Tdiff)/Tmidpoint < tolerance &&
		abs(Ldiff)/Lmidpoint < tolerance){
					if(logOut){cout << "=== Convergence achieved"<<endl;}
					return true;
	}
	return false;
}

void Stellar::setLog(bool sw){logOut=sw;}


void Stellar::getResult(){
	Phys physSurface = getPhys(Ndim-1);
	Phys physCenter  = getPhys(0);

	cout << "Stellar Center" << endl;
	physCenter.Out(); 
	cout << "Stellar Surface" << endl;
	physSurface.Out(); 
}

void Stellar::outNumberOfIterate(){
		cout << " number of Iterate:" <<  numberOfIterate << endl;
}

void Stellar::outBoundary(){
		cout << " Pc:" << Pc;
		cout << " Tc:" << Tc;
		cout << " Ls:" << Ls;
		cout << " Rs:" << Rs  << endl;
}

void Stellar::outDifference(Eigen::MatrixXd A,Eigen::VectorXd x,Eigen::VectorXd y){
			cout << "A:\n" << A << endl;
			cout << "A^-1:\n" << A.inverse() << endl;
			cout << "x:" << x.transpose() << endl;
			cout << "y:" << y.transpose() << endl;
}

void Stellar::checkOverflow(long index){
	if(!isfinite(T[index])){throw e_state::overflow;}
}

double Stellar::eps_pp(double rho,double X,double T){
	//return(EnergyGen::PP_KIP(rho,X,T));
	//Ref.R.Kippenhahn,p.165
	//temperature region are sliced.
	// 
	double my_factor = 1.0;
	if(5.E6>T){
		return(0.);
	}else{
		return(EnergyGen::PP_KIP(rho,X,T)*my_factor);
	}
}

double Stellar::eps_CNO(double rho,double X1,double X_CNO,double T){
	double my_factor = 1.0;
	if(1.5E7>T){
		return(0.);
	}else{
		return(EnergyGen::CNO_KIP(rho,X1,X_CNO,T)*my_factor);
	}
}

double Stellar::opacity(double rho,double T){
	double my_factor = 1.;
	return(Opacity::KramApprox(rho,X,Z,T)*my_factor);
}

double EnergyGen::PP_RPN(double rho,double X,double T){
	//Ref. R.P.Nelson,https://2019.qmplus.qmul.ac.uk/course/;
	//SPA7023U - STELLAR STRUCTURE AND EVOLUTION
	// RPN and KIP are differed on low temperature.
	// RPN is survive at low temperature.
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
double EnergyGen::PP_KIP(double rho,double X,double T){
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
double EnergyGen::PP_AML(double rho,double X,double T){
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
double EnergyGen::PP_REDUCED(double rho,double X,double T){
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
double EnergyGen::CNO_KIP(double rho,double X1,double X_CNO,double T){
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
