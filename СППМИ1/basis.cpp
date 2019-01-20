#include "1.h"

basis::basis()
{
	function<double(double)> phi[2];
	phi[0] = [](double ksi) { return 1 - ksi; };
	phi[1] = [](double ksi) { return ksi; };

	psi[0] = [phi](double ksi, double etta) { return phi[0](ksi) * phi[0](etta); };
	psi[1] = [phi](double ksi, double etta) { return phi[1](ksi) * phi[0](etta); };
	psi[2] = [phi](double ksi, double etta) { return phi[0](ksi) * phi[1](etta); };
	psi[3] = [phi](double ksi, double etta) { return phi[1](ksi) * phi[1](etta); };

	dpsidx[0] = [phi](double etta, double hx) { return - phi[0](etta) / hx; };
	dpsidx[1] = [phi](double etta, double hx) { return phi[0](etta) / hx; };
	dpsidx[2] = [phi](double etta, double hx) { return - phi[1](etta) / hx; };
	dpsidx[3] = [phi](double etta, double hx) { return phi[1](etta) / hx; };
		
	dpsidy[0] = [phi](double ksi, double hy) { return - phi[0](ksi) / hy; };
	dpsidy[1] = [phi](double ksi, double hy) { return - phi[1](ksi) / hy; };
	dpsidy[2] = [phi](double ksi, double hy) { return phi[0](ksi) / hy; };
	dpsidy[3] = [phi](double ksi, double hy) { return phi[1](ksi) / hy; };
}