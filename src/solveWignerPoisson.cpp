#include "lib.hpp"
#include "WignerFunction.hpp"
#include "poisson1D.hpp"

using namespace wigner;
using namespace poisson;


// TODO: solveWignerPoisson -> solved recursively ?

void WignerFunction::solveWignerPoisson(){

	Poisson1D p(nx_, dx_);
	p.dirichletL_ = 0, p.dirichletR_ = -uB_;  // - bo obniżamy U w prawym kontakcie
	p.neumannL_ = 0, p.neumannR_ = 0;
	p.epsilonR_ = epsilonR_, p.temp_ = temp_;
	p.uNew_ = uStart_, p.uOld_ = uStart_;

	// Heterostructure energy profile
	vec band_off = u_;

	// Doping profile
	vec nD(nx_, fill::zeros);
	double s = 0.01;
	for (size_t i=0; i<nx_; ++i)
		nD(i) = cD_*(1+1/(1+exp((x_(i)-lC_)/s/l_))-1/(1+exp((x_(i)-l_+lC_)/s/l_)));

	vec j0(nx_, fill::zeros), j1(nx_, fill::zeros), rho_old(nx_, fill::zeros), rho_diff(nx_, fill::zeros);
	vec u_fit(nx_, fill::zeros), u_fit_params;

	// Convergence criteria
	size_t n_max = 10, n_min = 0, n_it = 0, n_dj = 0, n_du = 0, n_conv = 1;
	double max_dj = 1e-4, max_du = 1e-6/AU_eV, max_pFun = 1e-8/AU_eV;
	bool conv_J = false, conv_pot = false, pFun_zero = false;

	// Mixing parameters
	p.beta_ = 1e-2;  // Potential mixing parameter
	double alpha = 1e-2;  // Density mixing parameter

	double curr = 0, nc = 0, nd = 0, q = 0;  // Current, carrier nr, dopant nr, total charge
	double dj_x = 0, du_x = 0;  // Maximum and minimum values

	// Start electron concentration
	u_ = band_off + p.uNew_;
	solveWignerEq();
	p.nE_ = calcCD_X();
	for (size_t i=0; i<nx_; ++i)  // Mixing old and new el. density
		p.rho_(i) = (1.-alpha)*p.rho_(i) + alpha*(nD(i)-p.nE_(i));

	std::ofstream poisson_step("test/poisson_step.out");
	// std::ofstream trChar;
	// trChar.open("wyniki/dane/poisson_trChar.out", std::ios::out);
	cout<<"# it.\tCurr. [Acm^-2]\tnE [cm^-2]\tnD [cm^-2]\t(nD-nE)/nD [-]\tjMax-jMin\tmax(du) [eV]"<<endl;
	while ( !( (n_du > n_conv) ) && (n_it < n_max) ) {  //  && (n_dj > n_conv) && pFun_zero

		// Solve Poisson eq.
		p.solve();

		// Solve Wigner/Boltzmann eq.
		u_ = band_off + p.uNew_;
		// u_ = band_off + u_fit;
		solveWignerEq();
		p.nE_ = calcCD_X(), j0 = j1, j1 = calcCurrArr();
		rho_old = p.rho_;
		for (size_t i=0; i<nx_; ++i)  // Mixing old and new el. density
			p.rho_(i) = (1.-alpha)*p.rho_(i) + alpha*(nD(i)-p.nE_(i));

		// TODO: uNew -> polynomial fit
		// u_fit = polyval(polyfit(x_,p.uNew_,20), x_);
		// p.uNew_ = u_fit;

		curr = calcCurr(), nc = calcNorm(), nd = calcInt(nD, dx_), q = calcInt(p.rho_, dx_);

		// Saving data
		p.uOld_.save("test/uOld.txt", arma_ascii);
		p.uNew_.save("test/uNew.txt", arma_ascii);
		p.du_.save("test/du.txt", arma_ascii);
		// u_fit.save("test/u_fit.txt", arma_ascii);
		p.rho_.save("test/rho.txt", arma_ascii);
		rho_old.save("test/rho_old.txt", arma_ascii);
		p.nE_.save("test/nE.txt", arma_ascii);
		j1.save("test/curr.txt", arma_ascii);

		saveWignerFun();

		// TODO: approx_equal() <- use for checking convergance criterium

		//
		// Check current convergance
		//
		// j0 = j1, j1 = calcCurrArr();
		// dj = j0 - j1;
		// j1_x = max(j1), j1_n = min(j1), dj_x = max(dj);
		// conv_J = fabs(dj_x/j1_x) > max_dj ? false : true;
		// n_dj = conv_J ? n_dj+1 : 0;

		//
		// Check potential convergance
		//
		du_x = max(p.du_);
		conv_pot = fabs(du_x) < max_du ? true : false;
		// pFun_zero = fabs(max(p.pFun_)) > max_pFun ? false : true;
		n_du = conv_pot ? n_du + 1 : 0;
		p.uOld_ = p.uNew_;

		// TODO: rho change / charge change -> convergance

		cout<<n_it
			<<'\t'<<curr*AU_Acm2
			<<'\t'<<nc
			<<'\t'<<nd
			<<'\t'<<q
			<<'\t'<<range(j1)
			<<'\t'<<du_x*AU_eV
			<<'\t'<<approx_equal(p.uNew_, p.uOld_, "absdiff", 0)
			<<'\t'<<n_du<<endl;

		//
		// Saving to file
		//
		poisson_step<<"## it  j(it)  n_el(it)  n_D(it)  q(it)  av_el(it)  dj(it)  dv_av(it)  conv_J?  conv_pot?"<<endl;
		poisson_step<<"# "<<n_it
			// <<' '<<n_it*dt_*AU_s*1e12
			<<'\t'<<calcCurr()*AU_Acm2
			<<'\t'<<nc
			<<'\t'<<nd
			<<'\t'<<q
			<<'\t'<<du_x*AU_eV<<'\n';
		poisson_step<<"## it  x  rho  uNew  j\n";
		for (size_t i=0; i<nx_; ++i)
			poisson_step<<n_it<<'\t'<<x_(i)*AU_nm
				<<'\t'<<p.rho_(i)  // col. 3
				<<'\t'<<p.uNew_(i)  // col. 4
				<<'\t'<<j1(i)<<'\n';  // col. 5
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

	readPotential("potentials/pot.out");
	p.du_ = uStart_ - p.uNew_;
	uStart_.save("test/uStart.txt", arma_ascii), p.uNew_.save("test/uNew.txt", arma_ascii);
	p.du_.save("test/du.txt", arma_ascii);
	p.du_.print("du:");
	cout<<approx_equal(p.uNew_, uStart_, "absdiff", 0)<<endl;


}
