#include "lib.hpp"
#include "WignerFunction.hpp"
#include "poisson1D.hpp"

using namespace wigner;
using namespace poisson;


void WignerFunction::solveWignerPoisson
	(double u_bias, double i_alpha, double i_beta, size_t i_n_max, bool timeDependent){

	Poisson1D p(nx_, dx_);
	p.dirichletL_ = u_bias/2., p.dirichletR_ = -u_bias/2.;  // - bo obniżamy U w prawym kontakcie
	p.epsilonR_ = epsilonR_, p.temp_ = temp_;
	p.uNew_ = uStart_, p.uOld_ = uStart_;

	// Doping profile
	arma::vec nD(nx_, arma::fill::zeros);
	double s = 0.005;
	for (size_t i=0; i<nx_; ++i)
		nD(i) = cD_*(1+1/(1+exp((x_(i)-lC_)/s/l_))-1/(1+exp((x_(i)-l_+lC_)/s/l_)));

	arma::vec j0(nx_, arma::fill::zeros), j1(nx_, arma::fill::zeros);
	arma::vec dj(nx_, arma::fill::zeros);
	arma::vec u_der(nx_, arma::fill::zeros);
	arma::vec rho_old(nx_, arma::fill::zeros), rho_new(nx_, arma::fill::zeros);
	arma::vec currD_a, nc_a, q_a, dJ_x_a, dU_x_a, dRho_x_a;
	arma::vec rho_test(nx_, arma::fill::zeros);

	// Convergence criteria
	size_t n_max = i_n_max, n_min = 10, n_it = 0;
	size_t n_dJ = 0, n_dU = 0, n_dRho = 0, n_conv = 10;
	double max_dJ = 1e-4, max_dU = 1e-3, max_pFun = 1e-8/AU_eV;
	double max_dRho = 1e-3;
	bool conv = false, pFun_zero = false;

	// Mixing parameters
	double alpha = i_alpha;  // Density mixing parameter
	p.beta_ = i_beta;  // Potential mixing parameter

	double curr = 0, nc = 0, nd = 0, q = 0;  // Current, carrier nr, dopant nr, total charge
	double dJ_x = 0, dU_x = 0, pFun_x = 0, dRho_x;  // Maximum and minimum values

	// TODO: Pointers to functions
	// void (WignerFunction::*solve)();
	// solve = timeDependent ? &WignerFunction::solveWignerEq : &WignerFunction::solveTimeEv;

	// Start electron concentration
	// uC_ = uStart_;
	// solveWignerEq();
	// calcCD_X();

	std::ofstream poisson_step("out_data/poisson_step.out");
	poisson_step<<"it,x [nm],rho [cm^{-3}],uNew [eV],{/Symbol d}u [eV],du/dx [au],J [Acm^{-2}]\n";

	cout.width( 5 );
	cout<<"# it.\tCurr. [Acm^-2]\tnE [cm^-2]\tnD [cm^-2]\tq [cm^-2]\tmax(dJ/J)\tmax(dU/U)\tmax(pFun) [ev]\tmax(dRho/rho)"<<endl;
	while ( !(n_dJ >= n_conv) && n_it < n_max ) {
		// !(n_dJ > n_conv && n_dU > n_conv) &&

		n_it += 1;
		// TODO: Change order (BTE -> PE)

		//
		// Solve Wigner/Boltzmann eq
		uC_ = p.uNew_;
		timeDependent ? solveTimeEv() : solveWignerEq();
		calcCD_X(), calcCD_K();
		curr = calcCurr();
		j0 = j1, j1 = currD_;

		//
		// Mixing old and new el. density
		if (n_it > 1) rho_new = (1.-alpha)*rho_old + alpha*(nD - cdX_);
		else rho_new = nD - cdX_;
		p.rho_ = rho_new, p.nE_ = cdX_;

		//
		// Solve Poisson eq.
		p.solve();

		//
		// Check current convergance
		dj = j0 - j1;  // Vectors
		dJ_x = max(abs(dj/j1));
		conv = dJ_x < max_dJ ? true : false;
		n_dJ = conv && n_it > n_min ? n_dJ+1 : 0;  // Ile razy został spełniony warunek
		// Check potential convergance
		dU_x = max(abs(p.du_/p.uNew_)), pFun_x = max(abs(p.pFun_));
		conv = dU_x < max_dU ? true : false;
		// n_dU = conv && n_it > n_min  ? n_dU + 1 : 0;
		// pFun_zero = pFun_x < max_pFun ? true : false;
		// Check charge den. convergance
		dRho_x = max(abs((rho_new-rho_old)/rho_new));
		conv = dRho_x < max_dRho ? true : false;
		// n_dRho = conv && n_it > n_min  ? n_dRho + 1 : 0;

		nc = calcNorm(), nd = calcInt(nD, dx_), q = calcInt(p.rho_, dx_);
		u_der = calcFirstDer(p.uNew_, dx_);
		rho_test = calcSecondDer(p.uNew_, dx_)*epsilonR_/4./M_PI;

		cout.width( 5 );
		cout<<n_it;
		cout.width( 5 );
		cout.setf( ios::scientific );
		cout.precision( 3 );
		cout<<'\t'<<curr*AU_Acm2
			<<'\t'<<nc
			<<'\t'<<nd
			<<'\t'<<q
			<<'\t'<<dJ_x
			<<'\t'<<dU_x
			<<'\t'<<pFun_x*AU_eV
			<<'\t'<<dRho_x
			<<'\t'<<n_dJ<<endl;

		//
		// Saving data
		field<std::string> header(10);
		arma::mat out_data;
		out_data.insert_cols(0, x_*AU_nm), header(0) = "x [nm]";
		out_data.insert_cols(1, p.uNew_*AU_eV), header(1) = "U_{new} [eV]";
		out_data.insert_cols(2, p.uOld_*AU_eV), header(2) = "U_{old} [eV]";
		out_data.insert_cols(3, p.du_*AU_eV), header(3) = "{/Symbol d}U [eV]";
		out_data.insert_cols(4, rho_new/AU_cm3), header(4) = "{/Symbol r}_{new} [cm^{-3}]";
		out_data.insert_cols(5, rho_old/AU_cm3), header(5) = "{/Symbol r}_{old} [cm^{-3}]";
		out_data.insert_cols(6, nD/AU_cm3), header(6) = "n_{D} [cm^{-3}]";
		out_data.insert_cols(7, cdX_/AU_cm3), header(7) = "n_{E} [cm^{-3}]";
		out_data.insert_cols(8, j1*AU_Acm2), header(8) = "J [Acm^{-2}]";
		out_data.insert_cols(9, rho_test/AU_cm3), header(9) = "{/Symbol r}_{test} [cm^{-3}]";
		out_data.save(csv_name("out_data/poisson_test.csv", header));

		saveWignerFun();

		poisson_step<<"## it  j [au]  n_el [au]  n_D [au]  q [au]  max(du)"<<endl;
		poisson_step<<"# "<<n_it
			// <<' '<<n_it*dt_*AU_s*1e12
			<<'\t'<<curr
			<<'\t'<<nc
			<<'\t'<<nd
			<<'\t'<<q
			<<'\t'<<dU_x<<'\n';
		poisson_step<<"## it  x [au]  rho [au]  uNew [au]  du [au]  u_der [au]  j [au]\n";
		for (size_t i=0; i<nx_; ++i)
			poisson_step<<n_it<<','<<x_(i)
				<<','<<p.rho_(i)/AU_cm3  // col. 3
				<<','<<p.uNew_(i)*AU_eV  // col. 4
				<<','<<p.du_(i)*AU_eV  // col. 5
				<<','<<u_der(i)  // col. 6
				<<','<<j1(i)*AU_Acm2<<'\n';  // col. 7
		poisson_step<<"\n";
		currD_a.resize(n_it), nc_a.resize(n_it), q_a.resize(n_it);
		dJ_x_a.resize(n_it), dU_x_a.resize(n_it), dRho_x_a.resize(n_it);
		currD_a(n_it-1) = curr, nc_a(n_it-1) = nc, q_a(n_it-1) = q;
		dJ_x_a(n_it-1) = dJ_x, dU_x_a(n_it-1) = dU_x, dRho_x_a(n_it-1) = dRho_x;

		rho_old = rho_new;
		p.uOld_ = p.uNew_;

	}
	poisson_step.close();

	std::ofstream tr_char("out_data/tr_char.csv");
	tr_char<<"it.,Curr. den. [au],n_E [au],q [au],dJ/J,max(dU/U),max(dRho/rho)"<<endl;
	for (size_t i=0; i<n_it; ++i)
		tr_char<<i
		<<','<<currD_a(i)
		<<','<<nc_a(i)
		<<','<<q_a(i)
		<<','<<dJ_x_a(i)
		<<','<<dU_x_a(i)
		<<','<<dRho_x_a(i)<<endl;
	tr_char.close();

	uStart_ = p.uNew_;
	// std::ofstream pot_out;
	// pot_out.open("out_data/poisson_pot/poisson_pot.out", std::ios::out);
	// pot_out<<"# i u(i)\n";
	// for (size_t i=0; i<nx_; ++i)
	// 		pot_out<<x_(i)<<' '<<uStart_(i)<<'\n';
	// pot_out.close();
	uStart_.save("out_data/poisson_pot.bin");

}
