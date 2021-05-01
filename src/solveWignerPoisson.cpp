#include "lib.hpp"
#include "WignerFunction.hpp"
#include "poisson1D.hpp"

using namespace wigner;
using namespace poisson;


// TODO: solveWignerPoisson -> solved recursively ?

void WignerFunction::solveWignerPoisson(){

	Poisson1D p(nx_, dx_);
	p.dirichletL_ = uBias_/2., p.dirichletR_ = -uBias_/2.;  // - bo obniżamy U w prawym kontakcie
	p.epsilonR_ = epsilonR_, p.temp_ = temp_;
	p.uNew_ = uStart_, p.uOld_ = uStart_;

	// Doping profile
	vec nD(nx_, fill::zeros);
	double s = 0.001;
	for (size_t i=0; i<nx_; ++i)
		nD(i) = cD_*(1+1/(1+exp((x_(i)-lC_)/s/l_))-1/(1+exp((x_(i)-l_+lC_)/s/l_)));

	vec j0(nx_, fill::zeros), j1(nx_, fill::zeros);
	vec dj(nx_, fill::zeros);
	vec u_der(nx_, fill::zeros);
	vec rho_old(nx_, fill::zeros), rho_new(nx_, fill::zeros);

	// Convergence criteria
	size_t n_max = 200, n_it = 0, n_dj = 0, n_du = 0, n_conv = 1;
	double max_dj = 100/AU_Acm2, max_du = 1e-6/AU_eV, max_pFun = 1e-8/AU_eV;
	bool conv_J = false, conv_pot = false, pFun_zero = false;

	// Mixing parameters
	p.beta_ = .1;  // Potential mixing parameter
	double alpha = 1;  // Density mixing parameter

	double curr = 0, nc = 0, nd = 0, q = 0;  // Current, carrier nr, dopant nr, total charge
	double dj_x = 0, du_x = 0, pFun_x = 0;  // Maximum and minimum values

	std::ofstream poisson_step("out_data/poisson_step.out");
	poisson_step<<"it,x [nm],rho [cm^{-3}],uNew [eV],{/Symbol d}u [eV],du/dx [au],J [Acm^{-2}]\n";
	cout<<"# it.\tCurr. [Acm^-2]\tnE [cm^-2]\tnD [cm^-2]\tq [cm^-2]\tdj [Acm^-2]\tmax(du) [eV]\tmax(pFun) [ev]"<<endl;
	while ( !( (n_du > n_conv) && (n_dj > n_conv) ) && (n_it < n_max) ) {  // && pFun_zero

		//
		// Solve Poisson eq.
		p.solve();
		u_der = calcDer(p.uNew_, dx_);

		//
		// Solve Wigner/Boltzmann eq
		uC_ = p.uNew_;
		solveWignerEq();
		calcCD_X(), calcCD_K();
		curr = calcCurr();
		j0 = j1, j1 = currD_;

		//
		// Mixing old and new el. density
		p.nE_ = cdX_;
		rho_new = (1.-alpha)*rho_old + alpha*(nD - cdX_);
		p.rho_ = rho_new;

		//
		// Check current convergance
		dj = j0 - j1;  // Vectors
		dj_x = max(abs(dj));
		conv_J = dj_x < max_dj ? true : false;
		n_dj = conv_J ? n_dj+1 : 0;  // Ile razy został spełniony warunek
		// Check potential convergance
		du_x = max(abs(p.du_)), pFun_x = max(abs(p.pFun_));
		conv_pot = du_x < max_du ? true : false;
		n_du = conv_pot ? n_du + 1 : 0;
		pFun_zero = pFun_x < max_pFun ? true : false;
		// TODO: rho change / charge change -> convergance

		nc = calcNorm(), nd = calcInt(nD, dx_), q = calcInt(p.rho_, dx_);

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
			<<'\t'<<du_x<<'\n';
		poisson_step<<"## it  x [au]  rho [au]  uNew [au]  du [au]  u_der [au]  j [au]\n";
		for (size_t i=0; i<nx_; ++i)
			poisson_step<<n_it<<','<<x_(i)
				<<','<<p.rho_(i)/AU_cm3  // col. 3
				<<','<<p.uNew_(i)*AU_eV  // col. 4
				<<','<<p.du_(i)*AU_eV  // col. 5
				<<','<<u_der(i)  // col. 6
				<<','<<j1(i)*AU_Acm2<<'\n';  // col. 7
		// for (size_t j=0; j<nk_; ++j)
		// 	poisson_step<<n_it<<'\t'<<k_(j)
		// 		<<'\t'<<nE_k(j)<<'\n';
		poisson_step<<"\n";

		rho_old = rho_new;
		p.uOld_ = p.uNew_;

		n_it += 1;
	}
	poisson_step.close();

	uStart_ = p.uNew_;
	std::ofstream pot_out;
	pot_out.open("out_data/poisson_pot/poisson_pot.out", std::ios::out);
	pot_out<<"# i u(i)\n";
	for (size_t i=0; i<nx_; ++i)
			pot_out<<x_(i)<<' '<<uStart_(i)<<'\n';
	pot_out.close();
	// uStart_ = p.uNew_;
	// out_data.reset();
	// out_out.save("out_data/poisson_pot/poissson_pot.out", arma_ascii);
	// out_data.insert_cols(x_(i)), out_data.insert_cols(uStart_(i));
	// out_data.save();

}
