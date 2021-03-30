#include "lib.hpp"
#include "WignerFunction.hpp"
#include "poisson1D.hpp"

using namespace wigner;
using namespace poisson;


// TODO: solveWignerPoisson -> solved recursively ?

void WignerFunction::solveWignerPoisson(){

	Poisson1D p(nx_, dx_);
	p.dirichletL_ = 0, p.dirichletR_ = -uBias_;  // - bo obniżamy U w prawym kontakcie
	p.epsilonR_ = epsilonR_, p.temp_ = temp_;
	p.uNew_ = uStart_, p.uOld_ = uStart_;

	// Doping profile
	vec nD(nx_, fill::zeros);
	double s = 0.01;
	for (size_t i=0; i<nx_; ++i)
		nD(i) = cD_*(1+1/(1+exp((x_(i)-lC_)/s/l_))-1/(1+exp((x_(i)-l_+lC_)/s/l_)));

	vec j0(nx_, fill::zeros), j1(nx_, fill::zeros), nC_old(nx_, fill::zeros), nC_diff(nx_, fill::zeros);
	vec dj(nx_, fill::zeros);
	vec u_der(nx_, fill::zeros), nE(nx_, fill::zeros);
	vec nE_k(nx_, fill::zeros);

	// Convergence criteria
	size_t n_max = 10, n_it = 0, n_dj = 0, n_du = 0, n_conv = 1;
	double max_dj = 100/AU_Acm2, max_du = 1e-6/AU_eV, max_pFun = 1e-8/AU_eV;
	bool conv_J = false, conv_pot = false, pFun_zero = false;

	// Mixing parameters
	p.beta_ = 3e-2;  // Potential mixing parameter
	double alpha = 1;  // Density mixing parameter

	double curr = 0, nc = 0, nd = 0, q = 0;  // Current, carrier nr, dopant nr, total charge
	double dj_x = 0, du_x = 0, pFun_x = 0;  // Maximum and minimum values

	// Start electron concentration
	uC_ = p.uNew_;
	solveWignerEq();
	p.nC_ = nD - calcCD_X();
	// p.nC_ = (1.-alpha)*p.nC_ + alpha*(nD - p.nE_);

	std::ofstream poisson_step("test/poisson_step.out");
	// std::ofstream trChar;
	// trChar.open("wyniki/dane/poisson_trChar.out", std::ios::out);
	cout<<"# it.\tCurr. [Acm^-2]\tnE [cm^-2]\tnD [cm^-2]\tq [cm^-2]\tdj [Acm^-2]\tmax(du) [eV]\tmax(pFun) [ev]"<<endl;
	while ( !( (n_du > n_conv) && (n_dj > n_conv) ) && (n_it < n_max) ) {  // && pFun_zero

		// Solve Poisson eq.
		p.solve();
		u_der = calcDer(p.uNew_, dx_);

		// Solve Wigner/Boltzmann eq
		uC_ = p.uNew_;
		solveWignerEq();
		nE = calcCD_X();
		nE_k = calcCD_K();
		j0 = j1, j1 = calcCurrArr();
		nC_old = p.nC_;
		// Mixing old and new el. density
		p.nE_ = nE;
		p.nC_ = (1.-alpha)*p.nC_ + alpha*(nD - p.nE_);
		// p.nE_ = (1.-alpha)*p.nE_ + alpha*nE;
		// p.nC_ = nD - p.nE_;

		curr = calcCurr(), nc = calcNorm(), nd = calcInt(nD, dx_), q = calcInt(p.nC_, dx_);

		// Saving data
		p.uOld_.save("test/uOld.txt", arma_ascii);
		p.uNew_.save("test/uNew.txt", arma_ascii);
		u_der.save("test/u_der.txt", arma_ascii);
		p.du_.save("test/du.txt", arma_ascii);
		p.pFun_.save("test/pFun.txt", arma_ascii);
		// u_fit.save("test/u_fit.txt", arma_ascii);
		p.nC_.save("test/nC.txt", arma_ascii);
		nC_old.save("test/nC_old.txt", arma_ascii);
		j1.save("test/curr.txt", arma_ascii);
		nE_k.save("test/nE_k.txt", arma_ascii);

		saveWignerFun();

		//
		// Check current convergance
		//
		dj = j0 - j1;
		dj_x = max(abs(dj));
		conv_J = dj_x < max_dj ? true : false;
		n_dj = conv_J ? n_dj+1 : 0;

		//
		// Check potential convergance
		//
		du_x = max(abs(p.du_)), pFun_x = max(abs(p.pFun_));
		conv_pot = du_x < max_du ? true : false;
		n_du = conv_pot ? n_du + 1 : 0;
		pFun_zero = pFun_x < max_pFun ? true : false;
		p.uOld_ = p.uNew_;

		// TODO: rho change / charge change -> convergance

		cout<<n_it
			<<'\t'<<curr*AU_Acm2
			<<'\t'<<nc
			<<'\t'<<nd
			<<'\t'<<q
			<<'\t'<<dj_x*AU_Acm2
			<<'\t'<<du_x*AU_eV
			<<'\t'<<pFun_x*AU_eV
			<<'\t'<<n_du<<endl;

		//
		// Saving to file
		//
		poisson_step<<"## it  j [au]  n_el [au]  n_D [au]  q [au]  max(du)"<<endl;
		poisson_step<<"# "<<n_it
			// <<' '<<n_it*dt_*AU_s*1e12
			<<'\t'<<calcCurr()
			<<'\t'<<nc
			<<'\t'<<nd
			<<'\t'<<q
			<<'\t'<<du_x<<'\n';
		poisson_step<<"## it  x [au]  rho [au]  uNew [au]  du [au]  u_der [au]  j [au]\n";
		for (size_t i=0; i<nx_; ++i)
			poisson_step<<n_it<<'\t'<<x_(i)
				<<'\t'<<p.nC_(i)  // col. 3
				<<'\t'<<p.uNew_(i)  // col. 4
				<<'\t'<<u_der(i)  // col. 5
				<<'\t'<<p.du_(i)  // col. 6
				<<'\t'<<j1(i)<<'\n';  // col. 7
		poisson_step<<"\n";

		n_it += 1;
	}
	// trChar.close();
	poisson_step.close();

	uStart_ = p.uNew_;
	std::ofstream pot_out;
	pot_out.open("potentials/pot.out", std::ios::out);
	pot_out<<"# i u(i)\n";
	for (size_t i=0; i<nx_; ++i)
			pot_out<<x_(i)<<' '<<uStart_(i)<<'\n';
	pot_out.close();

}
