#include "Stellar.h"

int main(){
	Stellar stellar(1.0);
	stellar.SetLogOut(true);
	stellar.Calc();
	stellar.getResult();
	stellar.Plot();

	Stellar stellar2(10.);
	stellar2.SetLogOut(false);
	stellar2.Calc();

}
