#include "Stellar.h"

int main(){
	Stellar stellar(1.0);
	stellar.SetLogOut(true);
	stellar.Calc();
	stellar.getResult();
	stellar.Plot();

	Phys phys0 = stellar.getPhys((10000/2)-2);
	Phys phys1 = stellar.getPhys((10000/2)-1);
	Phys phys2 = stellar.getPhys((10000/2));

	phys0.Out();
	phys1.Out();
	phys2.Out();

	//Stellar stellar2(10.);
	//stellar2.SetLogOut(false);
	//stellar2.Calc();

}
