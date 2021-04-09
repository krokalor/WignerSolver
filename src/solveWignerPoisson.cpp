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
	vec u_der(nx_, fill::zeros);
	mat out_data;

	// Convergence criteria
	size_t n_max = 5000, n_it = 0, n_dj = 0, n_du = 0, n_conv = 1;
	double max_dj = 100/AU_Acm2, max_du = 1e-6/AU_eV, max_pFun = 1e-8/AU_eV;
	bool conv_J = false, conv_pot = false, pFun_zero = false;

	// Mixing parameters
	p.beta_ = 4e-4;  // Potential mixing parameter
	double alpha = 1;  // Density mixing parameter

	double curr = 0, nc = 0, nd = 0, q = 0;  // Current, carrier nr, dopant nr, total charge
	double dj_x = 0, du_x = 0, pFun_x = 0;  // Maximum and minimum values

	// Start electron concentration
	uC_ = p.uNew_;
	solveWignerEq();
	p.nC_ = nD - calcCD_X();
	// p.nC_ = (1.-alpha)*p.nC_ + alpha*(nD - p.nE_);

	std::ofstream poisson_step("out_data/poisson_step.out");
	// std::ofstream trChar;
	// trChar.open("out_data/poisson_trChar.out", std::ios::out);
	cout<<"# it.\tCurr. [Acm^-2]\tnE [cm^-2]\tnD [cm^-2]\tq [cm^-2]\tdj [Acm^-2]\tmax(du) [eV]\tmax(pFun) [ev]"<<endl;
	while ( !( (n_du > n_conv) && (n_dj > n_conv) ) && (n_it < n_max) ) {  // && pFun_zero

		// p.beta_ = 3e-4*exp(n_it*2.3e-3);

		// Solve Poisson eq.
		p.solve();
		u_der = calcDer(p.uNew_, dx_);

		// Solve Wigner/Boltzmann eq
		uC_ = p.uNew_;
		solveWignerEq();
		calcCD_X();
		calcCD_K();
		curr = calcCurr();
		j0 = j1, j1 = currD_;
		nC_old = p.nC_;
		// Mixing old and new el. density
		p.nE_ = cdX_;
		p.nC_ = (1.-alpha)*p.nC_ + alpha*(nD - p.nE_);

		nc = calcNorm(), nd = calcInt(nD, dx_), q = calcInt(p.nC_, dx_);

		// Saving data
		out_data.resize(0, 0);
		out_data.insert_cols(0, x_*AU_nm);
		out_data.insert_cols(1, p.uNew_*AU_eV);
		out_data.insert_cols(2, p.du_*AU_eV);
		out_data.insert_cols(3, p.nC_*AU_eV);
		out_data.insert_cols(4, j1*AU_Acm2);

		out_data.save("out_data/poisson_test.txt", arma_ascii);

		saveWignerFun();

		//
		// Check current convergance
		//
		dj = j0 - j1;  // Vectors
		dj_x = max(abs(dj));
		conv_J = dj_x < max_dj ? true : false;
		n_dj = conv_J ? n_dj+1 : 0;  // Ile razy został spełniony warunek

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
			<<'\t'<<n_du
			<<'\t'<<p.beta_<<endl;

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
		// for (size_t j=0; j<nk_; ++j)
		// 	poisson_step<<n_it<<'\t'<<k_(j)
		// 		<<'\t'<<nE_k(j)<<'\n';
		poisson_step<<"\n";

		n_it += 1;
	}
	// trChar.close();
	poisson_step.close();

	uStart_ = p.uNew_;
	std::ofstream pot_out;
	pot_out.open("out_data/poisson_potential.out", std::ios::out);
	pot_out<<"# i u(i)\n";
	for (size_t i=0; i<nx_; ++i)
			pot_out<<x_(i)<<' '<<uStart_(i)<<'\n';
	pot_out.close();

}
