#include "lib.hpp"
#include "WignerFunction.hpp"
#include "poisson1D.hpp"

using namespace wigner;
using namespace poisson;


void WignerFunction::solveWignerPoisson
	(double u_bias, double i_alpha, double i_beta, size_t i_n_max){

	Poisson1D p(nx_, dx_);
	p.dirichletL_ = u_bias/2., p.dirichletR_ = -u_bias/2.;  // - bo obniżamy U w prawym kontakcie
	p.epsilonR_ = epsilonR_, p.temp_ = temp_;
	p.uNew_ = uStart_, p.uOld_ = uStart_;

	// Doping profile
	vec nD(nx_, fill::zeros);
	double s = 0.005;
	for (size_t i=0; i<nx_; ++i)
		nD(i) = cD_*(1+1/(1+exp((x_(i)-lC_)/s/l_))-1/(1+exp((x_(i)-l_+lC_)/s/l_)));

	vec j0(nx_, fill::zeros), j1(nx_, fill::zeros);
	vec dj(nx_, fill::zeros);
	vec u_der(nx_, fill::zeros);
	vec rho_old(nx_, fill::zeros), rho_new(nx_, fill::zeros);

	// Convergence criteria
	size_t n_max = i_n_max, n_it = 0;
	size_t n_dJ = 0, n_dU = 0, n_dRho = 0, n_conv = 1;
	double max_dJ = 100/AU_Acm2, max_dU = 1e-6/AU_eV, max_pFun = 1e-8/AU_eV;
	double max_dRho = 2e-6*AU_cm3/E0;
	bool conv = false, pFun_zero = false;

	// Mixing parameters
	double alpha = i_alpha;  // Density mixing parameter
	p.beta_ = i_beta;  // Potential mixing parameter

	double curr = 0, nc = 0, nd = 0, q = 0;  // Current, carrier nr, dopant nr, total charge
	double dJ_x = 0, dU_x = 0, pFun_x = 0, dRho_x;  // Maximum and minimum values

	// Start electron concentration
	// uC_ = uStart_;
	// solveWignerEq();
	// calcCD_X();
	// p.rho_ = (nD - cdX_);
	// p.nE_ = cdX_;
	// rho_new = p.rho_, rho_old = p.rho_;

	std::ofstream poisson_step("out_data/poisson_step.out");
	std::ofstream tr_char("out_data/tr_char.csv");
	poisson_step<<"it,x [nm],rho [cm^{-3}],uNew [eV],{/Symbol d}u [eV],du/dx [au],J [Acm^{-2}]\n";
	tr_char<<"it.,Curr. [Acm^{-2}],nE [cm^{-2}],q [cm^{-2}],dj [Acm^{-2}],max(du) [eV]";

	cout.width( 5 );
	cout<<"# it.\tCurr. [Acm^-2]\tnE [cm^-2]\tnD [cm^-2]\tq [cm^-2]\tdj [Acm^-2]\tmax(du) [eV]\tmax(pFun) [ev]\tmax(rho) [C/cm^3]\n"<<endl;
	while ( !( n_dU > n_conv ) && n_it < n_max ) {
		// && pFun_zero  !( n_dRho > n_conv &&  (n_dU > n_conv) && (n_dJ > n_conv) )

		// TODO: Change order (BTE -> PE)

		//
		// Solve Poisson eq.
		p.solve();
		u_der = calcDer(p.uNew_, dx_);

		//
		// Solve Wigner/Boltzmann eq
		uC_ = p.uNew_;
		// solveWignerEq();
		solveTimeEv();
		calcCD_X(), calcCD_K();
		curr = calcCurr();
		j0 = j1, j1 = currD_;

		//
		// Mixing old and new el. density
		rho_new = (1.-alpha)*rho_old + alpha*(nD - cdX_);
		p.rho_ = rho_new, p.nE_ = cdX_;

		//
		// Check current convergance
		dj = j0 - j1;  // Vectors
		dJ_x = max(abs(dj));
		conv = dJ_x < max_dJ ? true : false;
		n_dJ = conv ? n_dJ+1 : 0;  // Ile razy został spełniony warunek
		// Check potential convergance
		dU_x = max(abs(p.du_)), pFun_x = max(abs(p.pFun_));
		conv = dU_x < max_dU ? true : false;
		n_dU = conv ? n_dU + 1 : 0;
		pFun_zero = pFun_x < max_pFun ? true : false;
		// Check charge den. convergance
		dRho_x = max(abs(rho_new-rho_old));
		conv = dRho_x < max_dRho ? true : false;
		n_dRho = conv ? n_dRho + 1 : 0;

		nc = calcNorm(), nd = calcInt(nD, dx_), q = calcInt(p.rho_, dx_);

		cout.width( 5 );
		cout<<n_it;
		cout.width( 5 );
		cout.setf( ios::scientific );
		cout.precision( 2 );
		cout<<'\t'<<curr*AU_Acm2
			<<'\t'<<nc
			<<'\t'<<nd
			<<'\t'<<q
			<<'\t'<<dJ_x*AU_Acm2
			<<'\t'<<dU_x*AU_eV
			<<'\t'<<pFun_x*AU_eV
			<<'\t'<<dRho_x*E0/AU_cm3
			<<'\t'<<n_dU<<endl;

		//
		// Saving data
		field<std::string> header(9);
		mat out_data;
		out_data.insert_cols(0, x_*AU_nm), header(0) = "x [nm]";
		out_data.insert_cols(1, p.uNew_*AU_eV), header(1) = "U_{new} [eV]";
		out_data.insert_cols(2, p.uOld_*AU_eV), header(2) = "U_{old} [eV]";
		out_data.insert_cols(3, p.du_*AU_eV), header(3) = "{/Symbol d}U [eV]";
		out_data.insert_cols(4, rho_new/AU_cm3), header(4) = "{/Symbol r}_{new} [cm^{-3}]";
		out_data.insert_cols(5, rho_old/AU_cm3), header(5) = "{/Symbol r}_{old} [cm^{-3}]";
		out_data.insert_cols(6, nD/AU_cm3), header(6) = "n_{D} [cm^{-3}]";
		out_data.insert_cols(7, cdX_/AU_cm3), header(7) = "n_{E} [cm^{-3}]";
		out_data.insert_cols(8, j1*AU_Acm2), header(8) = "J [Acm^{-2}]";
		out_data.save(csv_name("out_data/poisson_test.csv", header));

		saveWignerFun();

		poisson_step<<"## it  j [au]  n_el [au]  n_D [au]  q [au]  max(du)"<<endl;
		poisson_step<<"# "<<n_it
			// <<' '<<n_it*dt_*AU_s*1e12
			<<'\t'<<calcCurr()
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
		tr_char<<n_it
			<<','<<curr*AU_Acm2
			<<','<<nc
			<<','<<q
			<<','<<dJ_x*AU_Acm2
			<<','<<dU_x*AU_eV<<endl;
		// for (size_t j=0; j<nk_; ++j)
		// 	poisson_step<<n_it<<'\t'<<k_(j)
		// 		<<'\t'<<nE_k(j)<<'\n';
		poisson_step<<"\n";

		rho_old = rho_new;
		p.uOld_ = p.uNew_;

		n_it += 1;
	}
	poisson_step.close();
	tr_char.close();

	uStart_ = p.uNew_;
	// std::ofstream pot_out;
	// pot_out.open("out_data/poisson_pot/poisson_pot.out", std::ios::out);
	// pot_out<<"# i u(i)\n";
	// for (size_t i=0; i<nx_; ++i)
	// 		pot_out<<x_(i)<<' '<<uStart_(i)<<'\n';
	// pot_out.close();
	uStart_.save("out_data/poisson_pot.bin");
	// out_data.reset();
	// out_out.save("out_data/poisson_pot/poissson_pot.out", arma_ascii);
	// out_data.insert_cols(x_(i)), out_data.insert_cols(uStart_(i));
	// out_data.save();

}
