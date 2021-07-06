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
	double s = 0.01;
	for (size_t i=0; i<nx_; ++i)
		nD(i) = cD_*(1+1/(1+exp((x_(i)-lC_)/s/l_))-1/(1+exp((x_(i)-l_+lC_)/s/l_)));

	arma::vec j0(nx_, arma::fill::zeros), j1(nx_, arma::fill::zeros);
	arma::vec dj(nx_, arma::fill::zeros);
	arma::vec u_der(nx_, arma::fill::zeros);
	arma::vec rho_old(nx_, arma::fill::zeros), rho_new(nx_, arma::fill::zeros);
	arma::vec currD_a, nc_a, q_a, dJ_x_a, dU_x_a, dRho_x_a;
	arma::vec rho_test(nx_, arma::fill::zeros);

	// Convergence criteria
	size_t n_max = i_n_max, n_min = 10;
	size_t n_dJ = 0, n_dU = 0, n_conv = 10;
	double max_dJ = 1e-6, max_dU = 1e-3;
	bool conv = false;

	// Mixing parameters
	double alpha = i_alpha;  // Density mixing parameter
	p.beta_ = i_beta;  // Potential mixing parameter

	double curr = 0, nc = 0, nd = 0, q = 0;  // Current, carrier nr, dopant nr, total charge
	double dJ_x = 0, dU_x = 0, pFun_x = 0;  // Maximum and minimum values

	// TODO: Pointers to functions
	// void (WignerFunction::*solve)();
	// solve = timeDependent ? &WignerFunction::solveWignerEq : &WignerFunction::solveTimeEv;

	// Start electron concentration
	// uC_ = uStart_;
	// solveWignerEq();
	// calcCD_X();

	// std::ofstream poisson_step("out_data/poisson_step.out");
	// poisson_step<<"it,x [nm],rho [cm^{-3}],uNew [eV],{/Symbol d}u [eV],du/dx [au],J [Acm^{-2}]\n";

	std::ofstream tr_char("out_data/tr_char.csv");
	tr_char<<"it.,J [au],n_E [au],q [au],max({/Symbol d}J/J),max({/Symbol d}U/U),max({/Symbol d}Rho/rho),(J_1-J_0)/J_1,range(J)"<<endl;
	cout.width( 5 );
	cout<<"# it.\tCurr. [Acm^-2]\tnE [cm^-2]\tnD [cm^-2]\tq [cm^-2]\tmax(dJ/J)\tmax(dU/U)\tmax(pFun) [ev]"<<endl;
	for ( size_t n_it = 0; n_it < n_max; ++n_it ) {

		// TODO: Change order (BTE -> PE)

		//
		// Solve Poisson eq.
		p.solve();

		//
		// Solve Wigner/Boltzmann eq
		uC_ = p.uNew_;
		timeDependent ? solveTimeEv() : solveWignerEq();
		calcCD_X();
		// calcCD_K();
		curr = calcCurr();
		j0 = j1, j1 = currD_;

		//
		// Mixing old and new el. density
		rho_new = (1.-alpha)*rho_old + alpha*(nD - cdX_);
		p.rho_ = rho_new, p.nE_ = cdX_;

		//
		// Maximum values
		dj = j0 - j1;  // Vectors
		dJ_x = arma::max(arma::abs(dj/j1));
		dU_x = arma::max(arma::abs(p.du_/p.uNew_));
		pFun_x = arma::max(arma::abs(p.pFun_));

		nc = calcNorm(), nd = calcInt(nD, dx_), q = nd - nc;
		// u_der = calcFirstDer(p.uNew_, dx_);
		// rho_test = calcSecondDer(p.uNew_, dx_)*epsilonR_/4./M_PI;

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
			<<'\t'<<pFun_x*AU_eV<<endl;

		tr_char<<n_it
		<<','<<curr
		<<','<<nc
		<<','<<q
		<<','<<dJ_x
		<<','<<dU_x
		<<','<<calcInt(dj,dx_)/l_/curr<<endl;

		//
		// Saving data

		std::ofstream test("out_data/poisson_test.csv");
		test << "x [nm],U_{new} [eV],U_{old} [eV],{/Symbol d}U [eV],"
		"{/Symbol r}_{new} [C/cm^{3}],{/Symbol r}_{old} [C/cm^{3}],"
		"n_{E} [cm^{-3}],n_{D} [cm^{-3}],J [Acm^{-2}]"<<endl;
		for (size_t i=0; i<nx_; ++i)
			test<<x_(i)*AU_nm<<','<<p.uNew_(i)*AU_eV<<','<<p.uOld_(i)*AU_eV<<','
			<<p.du_(i)*AU_eV<<','<<rho_new(i)*E0/AU_cm3<<','<<rho_old(i)*E0/AU_cm3<<','
			<<cdX_(i)/AU_cm3<<','<<nD(i)/AU_cm3<<','<<j1(i)*AU_Acm2<<endl;
		test.close();

		// field<std::string> header(9);
		// arma::mat out_data;
		// out_data.insert_cols(0, x_*AU_nm), header(0) = "x [nm]";
		// out_data.insert_cols(1, p.uNew_*AU_eV), header(1) = "U_{new} [eV]";
		// out_data.insert_cols(2, p.uOld_*AU_eV), header(2) = "U_{old} [eV]";
		// out_data.insert_cols(3, p.du_*AU_eV), header(3) = "{/Symbol d}U [eV]";
		// out_data.insert_cols(4, rho_new/AU_cm3), header(4) = "{/Symbol r}_{new} [cm^{-3}]";
		// out_data.insert_cols(5, rho_old/AU_cm3), header(5) = "{/Symbol r}_{old} [cm^{-3}]";
		// out_data.insert_cols(6, cdX_/AU_cm3), header(6) = "n_{E} [cm^{-3}]";
		// out_data.insert_cols(7, nD/AU_cm3), header(7) = "n_{D} [cm^{-3}]";
		// out_data.insert_cols(8, j1*AU_Acm2), header(8) = "J [Acm^{-2}]";
		// // out_data.insert_cols(9, rho_test/AU_cm3), header(9) = "{/Symbol r}_{test} [cm^{-3}]";
		// out_data.save(csv_name("out_data/poisson_test.csv", header));

		saveWignerFun();

		// poisson_step<<"## it  j [au]  n_el [au]  n_D [au]  q [au]  max(du)"<<endl;
		// poisson_step<<"# "<<n_it
		// 	// <<' '<<n_it*dt_*AU_s*1e12
		// 	<<'\t'<<curr
		// 	<<'\t'<<nc
		// 	<<'\t'<<nd
		// 	<<'\t'<<q
		// 	<<'\t'<<dU_x<<'\n';
		// poisson_step<<"## it  x [au]  rho [au]  uNew [au]  du [au]  u_der [au]  j [au]\n";
		// for (size_t i=0; i<nx_; ++i)
		// 	poisson_step<<n_it<<','<<x_(i)
		// 		<<','<<p.rho_(i)/AU_cm3  // col. 3
		// 		<<','<<p.uNew_(i)*AU_eV  // col. 4
		// 		<<','<<p.du_(i)*AU_eV  // col. 5
		// 		<<','<<u_der(i)  // col. 6
		// 		<<','<<j1(i)*AU_Acm2<<'\n';  // col. 7
		// poisson_step<<"\n";
		// currD_a.resize(n_it), nc_a.resize(n_it), q_a.resize(n_it);
		// dJ_x_a.resize(n_it), dU_x_a.resize(n_it), dRho_x_a.resize(n_it);
		// currD_a(n_it-1) = curr, nc_a(n_it-1) = nc, q_a(n_it-1) = q;
		// dJ_x_a(n_it-1) = dJ_x, dU_x_a(n_it-1) = dU_x, dRho_x_a(n_it-1) = dRho_x;

		//
		// Check convergance
		conv = dJ_x < max_dJ ? true : false;
		n_dJ = conv && n_it > n_min ? n_dJ+1 : 0;  // Ile razy został spełniony warunek
		//
		conv = dU_x < max_dU ? true : false;
		n_dU = conv && n_it > n_min  ? n_dU + 1 : 0;
		//
		if (n_dJ > n_conv) break;  //  && n_dU > n_conv

		rho_old = rho_new;
		p.uOld_ = p.uNew_;

	}
	// poisson_step.close();

	p.uNew_.save("out_data/poisson_pot.bin");
	p.rho_.save("out_data/poisson_rho.bin");

	uStart_ = p.uNew_;
	// std::ofstream pot_out;
	// pot_out.open("out_data/poisson_pot/poisson_pot.out", std::ios::out);
	// pot_out<<"# i u(i)\n";
	// for (size_t i=0; i<nx_; ++i)
	// 		pot_out<<x_(i)<<' '<<uStart_(i)<<'\n';
	// pot_out.close();

}
