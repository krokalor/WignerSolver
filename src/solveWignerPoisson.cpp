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
	p.beta_ = 1e-2;  // Potential mixing parameter
	double alpha = 1e-2;  // Density mixing parameter

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

	array<double> nE(nx_, 0.), n_new(nx_, 0.), rho(nx_, 0.), nE_prev(nx_, 0.), dnE(nx_, 0.);
	array<double> dj(nx_, 0.), j0(nx_, 0.), j1(nx_, 0.), du(nx_, 0);
	// array<double> jt, nt, Q;  // Transient characteristics
	bool conv_J = false, conv_pot = false;
	size_t n_max = 1e4, n_min = 2, n_it = 0, n_dj = 0, n_du = 0, n_conv = 2;
	double max_dj = 1e-4, max_du = 1e-6/AU_eV;

	double curr = 0, nc = 0, nd = 0, q = 0, rmse_j = 0;
	double j1_x = 0, j1_n = 0, dj_x = 0, du_x = 0, nE_x = 0, dnE_x = 0;

	std::ofstream poisson_step("wyniki/dane/poisson_step.out");  // "+std::to_string(uB_*AU_eV)+"
	std::ofstream out;
	out.open("wyniki/dane/poisson_trChar.out", std::ios::out);
	// out<<"# i j(i) n(i)\n";
	out<<"# it. \tt [ps] \tCurr. [Acm^-2] \tn_E [cm^-2] \tn_D [cm^-2] \tQ. [C*cm^-2]"
		<<"\tdJ [Acm^-2] \tmax(J) \tmin(J) \t|max(J) - min(J)| \tdJ/max(J) \tRMSE \tdn [cm^-3] \tdn/max(n) \tmax(du) [eV] \tCONV?"<<endl;
	cout<<"# it. \tt [ps] \tCurr. [Acm^-2] \tn_E [cm^-2] \tn_D [cm^-2] \tQ. [C*cm^-2]"
		<<"\tdJ [Acm^-2] \tmax(J) \tmin(J) \t|max(J) - min(J)| \tdJ/max(J) \tRMSE \tdn [cm^-3] \tdn/max(n) \tmax(du) [eV] \tCONV?"<<endl;
	while ( !( (n_du > n_conv) && (n_dj > n_conv) ) && (n_it < n_max) ) {

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

		if (n_it > n_min) {
			//
			// Check current convergance
			//
			// j0 = j1, j1 = calcCurr(), dj = (j0-j1)/j1;
			// if (n_it > n_min) n_dj = (fabs(dj) < max_dj || fabs(j1) < min_j) ? n_dj+1 : 0;
			// conv_J = n_dj > n_conv ? true : false;
			j0 = j1, j1 = calcCurrArr();
			dj = j0 - j1;
			j1_x = j1.max(), j1_n = j1.min(), dj_x = dj.max();
			conv_J = fabs(dj_x/j1_x) > max_dj ? false : true;
			n_dj = conv_J ? n_dj+1 : 0;

			//
			// Check potential convergance
			//
			// dv_av = 0, conv_pot = true;
			// for (size_t i=0; i<nx_; ++i)
			// 	if ( fabs(p.uNew_(i)-p.uOld_(i)) > max_du ) {
			// 		conv_pot = false;
			// 		break;
			// 	}
			du = p.uNew_-p.uOld_;
			du_x = du.max();
			conv_pot = fabs(du_x) > max_du ? false : true;
			n_du = conv_pot ? n_du + 1 : 0;
			// for (size_t i=0; i<nx_; ++i)
			// 	dv_av += fabs(p.uNew_(i)-p.uOld_(i))/float(nx_);
			// out<<i<<' '<<jt(i)<<' '<<nt(i)<<'\n';
			curr = calcCurr()*AU_Acm2, nc = calcNorm()/AU_cm2, nd = calcInt(nD, dx_)/AU_cm2, q = calcInt(nD-nE, dx_)*E0/AU_cm2;
			rmse_j = sqrt(((j1-j0)*(j1-j0)).sum())/nx_*AU_Acm2;

			//
			// Carrier density convergance
			//
			dnE = nE-nE_prev;
			nE_x = nE.max(), dnE_x = dnE.max();

			out<<n_it
				<<'\t'<<n_it*dt_*AU_s*1e12
				<<'\t'<<curr
				<<'\t'<<nc
				<<'\t'<<nd
				<<'\t'<<q
				<<'\t'<<dj_x*AU_Acm2
				<<'\t'<<j1_x*AU_Acm2
				<<'\t'<<j1_n*AU_Acm2
				<<'\t'<<fabs(j1_x-j1_n)*AU_Acm2
				<<'\t'<<dj_x/j1_x
				<<'\t'<<rmse_j
				<<'\t'<<dnE_x/AU_cm3
				<<'\t'<<dnE_x/nE_x
				<<'\t'<<du_x*AU_eV
				<<'\t'<<conv_J
				<<'\t'<<conv_pot<<endl;
			cout<<n_it
				<<'\t'<<n_it*dt_*AU_s*1e12
				<<'\t'<<curr
				<<'\t'<<nc
				<<'\t'<<nd
				<<'\t'<<q
				<<'\t'<<dj_x*AU_Acm2
				<<'\t'<<j1_x*AU_Acm2
				<<'\t'<<j1_n*AU_Acm2
				<<'\t'<<fabs(j1_x-j1_n)*AU_Acm2
				<<'\t'<<dj_x/j1_x
				<<'\t'<<rmse_j
				<<'\t'<<dnE_x/AU_cm3
				<<'\t'<<dnE_x/nE_x
				<<'\t'<<du_x*AU_eV
				<<'\t'<<conv_J
				<<'\t'<<conv_pot<<endl;

			saveWignerFun();
		}

		p.uOld_ = p.uNew_;
		nE_prev = nE;

		// av_el = 0;
		// for (size_t i=0; i<nx_; ++i)
		// 	av_el = ( x_(i) <= lC_ || x_(i) >= l_-lC_ ) ? av_el + nE(i) : av_el;
		// av_el = av_el * 2*lC_/l_;

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
			<<'\t'<<curr*AU_Acm2
			<<'\t'<<nc/AU_cm2
			<<'\t'<<calcInt(nD, dx_)/AU_cm2
			<<'\t'<<q*E0/AU_cm2
			// <<'\t'<<av_el/AU_cm3
			<<'\t'<<dj.max()/j1.max()<<'\t'<<du.max()
			<<'\t'<<conv_J<<'\t'<<conv_pot<<'\n';
		poisson_step<<"## x  rho  nD  nE  band_off  phi_hartree\n";
		for (size_t i=0; i<nx_; ++i)
				poisson_step<<x_(i)*AU_nm
					<<'\t'<<p.rho_(i)/AU_cm3
					<<'\t'<<nD(i)/AU_cm3
					<<'\t'<<nE(i)/AU_cm3
					<<'\t'<<band_off(i)*AU_eV
					<<'\t'<<p.uNew_(i)*AU_eV
					<<'\t'<<j1(i)*AU_Acm2
					<<'\t'<<dj(i)*AU_Acm2
					<<'\t'<<du(i)*AU_eV<<'\n';
		poisson_step<<"\n\n";

		n_it += 1;
	}
	out.close();
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

}
