#include "StellerModel.h"

int main(){
	Steller steller(1.0);
	steller.SetLogOut(true);
	steller.Calc();
	steller.getResult();
	steller.Plot();

	Steller steller2(10.);
	steller2.SetLogOut(true);
	steller2.Calc();

}
