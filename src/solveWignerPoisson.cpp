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

	// Mixing parameters
	p.beta_ = 1e-3;  // Potential mixing parameter
	double alpha = 1e-3;  // Density mixing parameter

	// Structure potential
	// setLinPot(uB_);
	vec band_off = u_;

	// Doping profile
	vec nD(nx_, fill::zeros);
	// for (size_t i=0; i<nx_; ++i)
		// nD(i) = ( x_(i) <= lC_ || x_(i) >= l_-lC_ ) ? cD_ : 0;
	double s = 0.01;
	for (size_t i=0; i<nx_; ++i)
		nD(i) = cD_*(1+1/(1+exp((x_(i)-lC_)/s/l_))-1/(1+exp((x_(i)-l_+lC_)/s/l_)));
	// for (size_t i=0; i<nx_; ++i)
	//   cout<<x_(i)*AU_nm<<' '<<nD(i)/AU_cm3<<endl;

	vec nE(nx_, fill::zeros), n_new(nx_, fill::zeros), rho(nx_, fill::zeros), nE_prev(nx_, fill::zeros), dnE(nx_, fill::zeros);
	vec dj(nx_, fill::zeros), j0(nx_, fill::zeros), j1(nx_, fill::zeros), du(nx_, fill::zeros);
	// array<double> jt, nt, Q;  // Transient characteristics
	bool conv_J = false, conv_pot = false, pFun_zero = false;
	size_t n_max = 100, n_min = 0, n_it = 0, n_dj = 0, n_du = 0, n_conv = 2;
	double max_dj = 1e-4, max_du = 1e-6/AU_eV, max_pFun = 1e-8/AU_eV;

	double curr = 0, nc = 0, nd = 0, q = 0, rmse_j = 0;
	double j1_x = 0, j1_n = 0, dj_x = 0, du_x = 0, nE_x = 0, dnE_x = 0;

	std::ofstream poisson_step("wyniki/dane/poisson_step.out");  // "+std::to_string(uB_*AU_eV)+"
	std::ofstream out;
	out.open("wyniki/dane/poisson_trChar.out", std::ios::out);
	// out<<"# i j(i) n(i)\n";
	// out<<"# it.\tt [ps]\tCurr. [Acm^-2]\tn_E [cm^-2]\tn_D [cm^-2]\tQ. [C*cm^-2]"
	// 	<<"\tdJ/max(J)\t|max(J) - min(J)|\tmax(du) [eV]\tpFun [eV]\tCONV?"<<endl;
	// cout<<"# it.\tt [ps]\tCurr. [Acm^-2]\tn_E [cm^-2]\tn_D [cm^-2]\tQ. [C*cm^-2]"
	// 	<<"\tdJ/max(J)\t|max(J) - min(J)|\tmax(du) [eV]\tpFun [eV]\tCONV?"<<endl;
	cout<<"# it.\tCurr. [Acm^-2]\tn_E [cm^-2]\tn_D [cm^-2]\tQ. [C*cm^-2]\tmax(du) [eV]"<<endl;
	while ( (n_it < n_max) ) {  //  && (n_dj > n_conv) && pFun_zero !( (n_du > n_conv) && pFun_zero ) &&

		//
		// Solve Wigner/Boltzmann eq.
		//
		u_ = band_off + p.uNew_;
		solveWignerEq();
		nE = calcCD_X();

		for (size_t i=0; i<nx_; ++i) // mixing old and new el. density
			rho(i) = (1.-alpha)*rho(i) + alpha*(nD(i)-nE(i));
			// nE(i) = (1.-alpha)*nE(i) + alpha*nE(i);
		// rho = nD-nE;

		//
		// Solve Poisson eq.
		//
		p.rho_ = rho, p.nE_ = nE;
		p.solve();

		p.uNew_.save("uNew.txt", arma_ascii);
		rho.save("rho.txt", arma_ascii);
		nE.save("nE.txt", arma_ascii);

		cout<<n_it
			<<'\t'<<calcCurr()*AU_Acm2
			<<'\t'<<calcNorm()
			<<'\t'<<calcInt(nD, dx_)
			<<'\t'<<calcInt(nD-nE, dx_)
			<<'\t'<<max(p.uNew_-p.uOld_)*AU_eV<<endl;

		saveWignerFun();

		// TODO: approx_equal() <- use for checking convergance criterium

		// if (n_it > n_min) {
		// 	//
		// 	// Check current convergance
		// 	//
		// 	j0 = j1, j1 = calcCurrArr();
		// 	dj = j0 - j1;
		// 	j1_x = max(j1), j1_n = min(j1), dj_x = max(dj);
		// 	conv_J = fabs(dj_x/j1_x) > max_dj ? false : true;
		// 	n_dj = conv_J ? n_dj+1 : 0;
		//
		// 	//
		// 	// Check potential convergance
		// 	//
		// 	du = p.uNew_-p.uOld_;
		// 	du_x = max(du);
		// 	conv_pot = fabs(du_x) > max_du ? false : true;
		// 	pFun_zero = fabs(max(p.pFun_)) > max_pFun ? false : true;
		// 	n_du = conv_pot ? n_du + 1 : 0;
		// 	curr = calcCurr()*AU_A, nc = calcNorm(), nd = calcInt(nD, dx_), q = calcInt(nD-nE, dx_);
		//
		// 	//
		// 	// Carrier density convergance
		// 	//
		// 	dnE = nE-nE_prev;
		// 	nE_x = max(nE), dnE_x = max(dnE);
		//
		// 	out<<n_it
		// 		<<'\t'<<n_it*dt_*AU_s*1e12
		// 		<<'\t'<<curr
		// 		<<'\t'<<nc
		// 		<<'\t'<<nd
		// 		<<'\t'<<q
		// 		<<'\t'<<dj_x/j1_x
		// 		<<'\t'<<fabs(j1_x-j1_n)/j1_x
		// 		<<'\t'<<du_x*AU_eV
		// 		<<'\t'<<max(p.pFun_)
		// 		<<'\t'<<conv_J
		// 		<<'\t'<<conv_pot<<endl;
		// 	cout<<n_it
		// 		<<'\t'<<n_it*dt_*AU_s*1e12
		// 		<<'\t'<<curr
		// 		<<'\t'<<nc
		// 		<<'\t'<<nd
		// 		<<'\t'<<q
		// 		<<'\t'<<dj_x/j1_x
		// 		<<'\t'<<fabs(j1_x-j1_n)/j1_x
		// 		<<'\t'<<du_x*AU_eV
		// 		<<'\t'<<max(p.pFun_)
		// 		<<'\t'<<conv_J
		// 		<<'\t'<<conv_pot<<endl;
		//
		// 	saveWignerFun();
		// }

		p.uOld_ = p.uNew_;

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
		// poisson_step<<"## it  j(it)  n_el(it)  n_D(it)  q(it)  av_el(it)  dj(it)  dv_av(it)  conv_J?  conv_pot?"<<endl;
		// poisson_step<<"# "<<n_it
		// 	// <<' '<<n_it*dt_*AU_s*1e12
		// 	<<'\t'<<curr*AU_Acm2
		// 	<<'\t'<<nc/AU_cm2
		// 	<<'\t'<<calcInt(nD, dx_)/AU_cm2
		// 	<<'\t'<<q*E0/AU_cm2
		// 	// <<'\t'<<av_el/AU_cm3
		// 	<<'\t'<<dj.max()/j1.max()<<'\t'<<du.max()
		// 	<<'\t'<<conv_J<<'\t'<<conv_pot<<'\n';
		// poisson_step<<"## x  rho  nD  nE  band_off  phi_hartree\n";
		// for (size_t i=0; i<nx_; ++i)
		// 		poisson_step<<x_(i)*AU_nm
		// 			<<'\t'<<p.rho_(i)/AU_cm3
		// 			<<'\t'<<nD(i)/AU_cm3
		// 			<<'\t'<<nE(i)/AU_cm3
		// 			<<'\t'<<band_off(i)*AU_eV
		// 			<<'\t'<<p.uNew_(i)*AU_eV
		// 			<<'\t'<<j1(i)*AU_Acm2
		// 			<<'\t'<<dj(i)*AU_Acm2
		// 			<<'\t'<<du(i)*AU_eV<<'\n';
		// poisson_step<<"\n\n";

		n_it += 1;
	}
	out.close();
	poisson_step.close();

	// nE = calcCD_X();
	// nE.print("RIGHT AFTER while loop:");

	uStart_ = p.uNew_;
	std::ofstream pot_out;
	pot_out.open("potentials/pot.out", std::ios::out);
	pot_out<<"# i u(i)\n";
	for (size_t i=0; i<nx_; ++i)
			pot_out<<x_(i)<<' '<<uStart_(i)<<'\n';
	pot_out.close();

}
