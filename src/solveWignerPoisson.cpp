#include "lib.hpp"
#include "WignerFunction.hpp"
#include "poisson1D.hpp"

using namespace wigner;
using namespace poisson;


// TODO: solveWignerPoisson -> solved recursively ?

void WignerFunction::solveWignerPoisson(){

	Poisson1D p(nx_, dx_);
	p.dirichletL_ = 0, p.dirichletR_ = -uB_;
	p.neumannL_ = 0, p.neumannR_ = 0;
	p.epsilonR_ = epsilonR_, p.temp_ = temp_;
	p.uNew_ = uStart_, p.uOld_ = uStart_;

	// Mixing parameters
	p.beta_ = 2e-2;  // Potential mixing parameter
	double alpha = 2e-2;  // Density mixing parameter

	// Structure potential
	// setLinPot(uB_);
	array<double> band_off = u_;

	// Dopng Array
	// TODO: Profil domieszkowania - osobna funkcja
	array<double> nD(nx_, 0.);
	// for (size_t i=0; i<nx_; ++i)
		// nD(i) = ( x_(i) <= lC_ || x_(i) >= l_-lC_ ) ? cD_ : 0;
	double s = 0.004;
	for (size_t i=0; i<nx_; ++i)
		nD(i) = cD_*(1+1/(1+exp((x_(i)-lC_)/s/l_))-1/(1+exp((x_(i)-l_+lC_)/s/l_)));
	// for (size_t i=0; i<nx_; ++i)
	//   cout<<x_(i)*AU_nm<<' '<<nD(i)/AU_cm3<<endl;

	array<double> nE(nx_, 0.), n_new(nx_, 0.), rho(nx_, 0.);
	double dj = 0, j0 = 0, j1 = 0;
	array<double> jt, nt, Q;  // Transient characteristics
	bool conv_J = false, conv_pot = false;
	size_t n_max = 1e4, n_min = 10, n_it = 0, n_dj = 0, n_du = 0, n_conv = 2;
	double max_dj = 1e-4, max_dv = 1e-6/AU_eV, dv_av = 0, min_j = 1e-2/AU_Acm2;

	double av_el = 0;

	dt_ = 1e-15/AU_s;

	std::ofstream poisson_step("wyniki/poisson_step_.out");  // "+std::to_string(uB_*AU_eV)+"
	while ( !(conv_pot && conv_J) && (n_it < n_max) ) {

		//
		// Solve Wigner/Boltzmann eq.
		//
		u_ = band_off + p.uNew_;
		solveWignerEq();
		// solveTimeEv();
		nE = calcCD_X();
		for (size_t i=0; i<nx_; ++i) // mixing old and new el. density
			rho(i) = (1.-alpha)*rho(i) + alpha*(nD(i)-nE(i));

		//
		// Solve Poisson eq.
		//
		p.rho_ = rho, p.nE_ = nE;
		p.solve();

		//
		// Check current convergance
		//
		j0 = j1, j1 = calcCurr(), dj = (j0-j1)/j1;
		if (n_it > n_min) n_dj = (fabs(dj) < max_dj || fabs(j1) < min_j) ? n_dj+1 : 0;
		conv_J = n_dj > n_conv ? true : false;

		//
		// Check potential convergance
		//
		dv_av = 0, conv_pot = true;
		for (size_t i=0; i<nx_; ++i)
			if ( fabs(p.uNew_(i)-p.uOld_(i)) > max_dv ) {
				conv_pot = false;
				break;
			}
		n_du = conv_pot ? n_du + 1 : 0;
		conv_pot = n_du > n_conv ? true : false;
		for (size_t i=0; i<nx_; ++i)
			dv_av += fabs(p.uNew_(i)-p.uOld_(i))/float(nx_);

		for (size_t i=0; i<nx_; ++i)
			p.uOld_(i) = p.uNew_(i);

		jt.add(j1), nt.add(calcInt(nE, dx_)), Q.add(calcInt(nD-nE, dx_));

		av_el = 0;
		for (size_t i=0; i<nx_; ++i)
			av_el = ( x_(i) <= lC_ || x_(i) >= l_-lC_ ) ? av_el + nE(i) : av_el;
		av_el = av_el * 2*lC_/l_;

		//
		// Check divergance
		//
		// if (calcInt(nE, dx_) < 1e-3*calcInt(nD, dx_)) {
		//   cout<<"ERROR DURING CONVERGENCE: Divergence criteria met, adjust parameters"<<endl;
		//   break;
		// }

		//
		// Saving to file
		//
		poisson_step<<"## it  j(it)  n_el(it)  n_D(it)  q(it)  av_el(it)  dj(it)  dv_av(it)  conv_J?  conv_pot?"<<endl;
		poisson_step<<"# "<<n_it
			// <<' '<<n_it*dt_*AU_s*1e12
			<<'\t'<<jt(n_it)*AU_Acm2
			<<'\t'<<nt(n_it)/AU_cm2
			<<'\t'<<calcInt(nD, dx_)/AU_cm2
			<<'\t'<<Q(n_it)/AU_cm2
			<<'\t'<<av_el/AU_cm3
			<<'\t'<<dj<<'\t'<<dv_av
			<<'\t'<<conv_J<<'\t'<<conv_pot<<'\n';
		poisson_step<<"## x  rho  nD  nE  band_off  phi_hartree\n";
		for (size_t i=0; i<nx_; ++i)
				poisson_step<<x_(i)*AU_nm<<'\t'<<p.rho_(i)/AU_cm3<<'\t'<<
					nD(i)/AU_cm3<<'\t'<<nE(i)/AU_cm3<<'\t'<<
					band_off(i)*AU_eV<<'\t'<<p.uNew_(i)*AU_eV<<'\n';
		poisson_step<<"\n\n";
		cout<<n_it
			<<'\t'<<j1*AU_Acm2
			<<'\t'<<nt(n_it)/AU_cm2
			<<'\t'<<calcInt(nD, dx_)/AU_cm2
			<<'\t'<<Q(n_it)/AU_cm2
			<<'\t'<<dj
			<<'\t'<<dv_av
			<<'\t'<<conv_J<<'\t'<<conv_pot<<endl;

		n_it += 1;
	}
	poisson_step.close();

	double rho_sum = 0;
	for (size_t i=0; i<nx_; ++i)
		rho_sum += p.rho_(i);

	uStart_ = p.uNew_;
	std::ofstream pot_out;
	pot_out.open("potentials/pot.out", std::ios::out);
	pot_out<<"# i u(i)\n";
	for (size_t i=0; i<nx_; ++i)
			pot_out<<x_(i)<<' '<<uStart_(i)<<'\n';
	pot_out.close();

	// cout<<"# rho_sum = "<<rho_sum<<endl;

	std::ofstream file;
	file.open("wyniki/poisson_tranChar.out", std::ios::out);
	file<<"# i j(i) n(i)\n";
	for (size_t i=0; i<jt.size(); ++i)
			file<<i<<' '<<jt(i)<<' '<<nt(i)<<'\n';
	file.close();

	saveWignerFun();

}
