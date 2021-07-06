#include "src/lib.hpp"
#include "src/WignerFunction.hpp"
#include "src/poisson1D.hpp"

#include <chrono>


using namespace wigner;
using namespace poisson;
using namespace std::chrono;


void dissDecoh(WignerFunction&, double, double, size_t, double, double, size_t);
void calcTimeChar(WignerFunction&);
void gaussBarrier(WignerFunction&);


int main(){

	// TODO: ALL OUTPUT AND INPUT TO .CSV FORMAT, USE SAVE/LOAD FUNCTIONS

	// omp_set_dynamic(0);
	// omp_set_num_threads(1);
	// cout<<"omp_get_max_threads(): "<<omp_get_max_threads() <<endl;
	// cout<<"omp_get_num_threads(): "<<omp_get_num_threads()<<endl;

	high_resolution_clock::time_point t_start, t_end;
	duration<double> t_elapsed;
	t_start = high_resolution_clock::now();

	size_t nx = 150, nk = 150;
	double lD = 1000/AU_nm, lC = 1500/AU_nm;
	double k_max = 0.15;  // -1, 0.15

	WignerFunction f(nx, lD, lC, nk, k_max);

	arma::vec x_val = f.get_x_arr(), k_val = f.get_k_arr();
	// double dx = f.get_dx(), dk = f.get_dk();

	f.set_m(0.067);
	f.set_temp(4.2);
	f.set_epsilonR(13.1);
	f.set_cD(2e18*AU_cm3);
	f.set_uL( calcFermiEn(f.get_cD(), f.get_m(), f.get_temp()) );
	f.set_uR( calcFermiEn(f.get_cD(), f.get_m(), f.get_temp()) ); // calcFermiEn(f.get_cD(), f.get_m(), f.get_temp())
	f.set_dt(5*1e-15/AU_s);  // .005*1e-15/AU_s

	f.set_rR(0), f.set_rM(0); // 1./(1e-12/AU_s)
	f.set_rG(0), f.set_rF(0), f.set_lambda(0);

	f.set_useQC(false);
	f.set_useNLP(false);
	f.set_uBias_BC(true);

	// "UDS1", "UDS2", "UDS3", "HDS22"
	f.set_diffSch_K("HDS22");
	f.set_diffSch_P("HDS22");
	f.set_diffSch_J("HDS22");

	// Warunek brzegowy
	// 0 -> 0, 1 -> SF, 2:4 -> splot with SF, -1 -> Gauss, -2:-4 splot with Gauss
	cout<<"# Setting up BC"<<endl;
	f.set_bcType(1);

	//
	// Potencjał
	cout<<"# Setting up potential"<<endl;
	f.set_uBias(0./AU_eV);
	// f.setPotBias(0.1/AU_eV);
	//
	// f.addRectBarr(0.3/AU_eV, 17.5/AU_nm, 2/AU_nm, 10);
	// f.addRectBarr(0.3/AU_eV, 22.5/AU_nm, 2/AU_nm, 10);
	//
	// f.addRectBarr(0.3/AU_eV, 1750/AU_nm, 200/AU_nm, 10);
	// f.addRectBarr(0.3/AU_eV, 2250/AU_nm, 200/AU_nm, 10);
	// f.addGaussBarr(0.3/AU_eV, 2000/AU_nm, 500/AU_nm);
	//
	// f.load_poisson_pot("poisson_pot_100meV_4e4it.bin");
	// f.set_uC(f.get_uStart());
	//
	// Read potential from poisson_test
	// mat a; a.load("IV_char/out_data_10meV_2e5it/poisson_test.csv", csv_ascii);
	// f.set_uC(a.col(1)/AU_eV);

	//
	// FUNKCJA RÓWNOWAGOWA
	// f.setEquilibriumFunction("out_data/wf_feq_BP.bin", true);

	f.printParam();

	//
	// Iterator test
	// mat X(5, 5, fill::randn);
	// X.print();
	// for (size_t j=5; j--;) {
	// 	for(arma::mat::col_iterator i = X.begin_col(j)+1; i != X.end_col(j); ++i)
	// 		cout << *i + *(i-1) << ' ' ;
	// 	cout << endl;
	// }

	//
	// Boltzmann
	// cout<<"# Solving BTE"<<endl;
	// f.solveWignerEq();
	// f.saveWignerFun();

	// f.addWavePacket(500/AU_nm, 100/AU_nm, 0.05, 0.005);  // sqrt(2*f.get_m()*f.get_uL())
	// f.addWavePacket(3500/AU_nm, 100/AU_nm, -0.05, 0.005);  // sqrt(2*f.get_m()*f.get_uL())
	// double t_total = 50e-15/AU_s, t = 0, dt = f.get_dt();
	// while (t <= t_total) {
	// 	t += dt;
	// 	f.solveTimeEv();
	// 	f.saveWignerFun();
	// 	cout<<t*AU_s*1e15<<' '<<f.calcEK()<<' '<<sqrt(f.calcEK2())<<' '<<f.calcEX()*AU_nm<<endl;
	// }

	//
	// Boltzmann-Poisson
	cout<<"# Solving B-P set of equations"<<endl;
	// (uBias, alpha, beta, n_max, timeDependent)
	f.solveWignerPoisson(0.01/AU_eV, 2e-5, 1, 2e5, false);
	f.saveWignerFun();

	// arma::mat wf1, wf2, wf;
	// wf = f.get_wf();
	// wf1.load("out_data/wf_BP_100meV.bin");
	// wf2.load("out_data/wf_feq_BP.bin");
	// wf = f.get_wf();
	// for (size_t i=0; i<nx; ++i)
	// 	for (size_t j=0; j<nk; ++j)
	// 		wf(i,j) = ( wf1(i,j) - wf2(i,j) );  //  * k_val(j)
	// f.set_wf(wf);
	// f.saveWignerFun();

	//
	// Poisson test
	// Read potential from poisson_test
	// mat a; a.load("IV_char/out_data_10meV_2e5it/poisson_test.csv", csv_ascii);
	// Poisson1D p(nx, dx);
	// p.dirichletL_ = 0.01/AU_eV/2., p.dirichletR_ = -0.01/AU_eV/2.;  // - bo obniżamy U w prawym kontakcie
	// p.epsilonR_ = f.get_epsilonR(), p.temp_ = f.get_temp();
	// p.rho_ = a.col(4)*AU_cm3;
	// p.uOld_ = a.col(2)/AU_eV;
	// // p.rho_ = -1*normpdf(linspace(0, 200, 200), 100, 10);
	// p.solve();
	// cout.setf( ios::scientific ), cout.precision( 5 );
	// for (size_t i=0; i<nx; ++i)
	// 	cout<<x_val(i)*AU_nm<<'\t'<<p.rho_(i)/AU_cm3<<'\t'<<p.uOld_(i)*AU_eV<<'\t'<<p.uNew_(i)*AU_eV<<endl;
	// p.testPoisson();

	//
	// Barriera gaussowska
	// gaussBarrier(f);

	//
	// Zalezność funkcji rozkładu od dyssypacji
	// dissDecoh(f, 1e-13, 1e-10, 5, 0.0/AU_eV, 0.001/AU_eV, 10);

	//
	// Charakterystyka J-V
	// f.calc_IVchar(0.0/AU_eV, 0.1/AU_eV, 4);
	// f.saveWignerFun();

	//
	// FINAL RESULTS
	//

	cout.setf( ios::scientific ), cout.precision( 5 );

	double curr = f.calcCurr();
	f.calcCD_X(), f.calcCD_K();
	arma::vec cdX = f.get_cdX(), cdK = f.get_cdK();

	// Doping profile
	arma::vec nD(nx, arma::fill::zeros), rho(nx, arma::fill::zeros);
	double s = 0.005, l = 2*lC+lD;
	for (size_t i=0; i<nx; ++i)
		nD(i) = f.get_cD()*(1+1/(1+exp((x_val(i)-lC)/s/l))-1/(1+exp((x_val(i)-l+lC)/s/l)));
	rho = nD - cdX;

	double n = f.calcNorm();
	cout<<"# <p> = "<<f.calcEK()<<", sqrt(<p^2>) = "<<sqrt(f.calcEK2())
		<<", J(<p>) = "<<f.calcEK()/f.get_m() * n/f.get_l() * AU_Acm2
		<<", J(sqrt(<p^2>)) = "<<sqrt(f.calcEK2())/f.get_m() * n/f.get_l() * AU_Acm2<<endl;

	//  exK (nx, arma::fill::zeros), exK2 (nx, arma::fill::zeros);
	// for (size_t i=0; i<nx; ++i) {
	// 	for (size_t j=1; j<nk; ++j) {
	// 		exK(i) += ( wf(i,j) + wf(i,j-1) )/n * k_val(j) * dk/2.;
	// 		exK2(i) += ( wf(i,j) + wf(i,j-1) )/n * k_val(j)*k_val(j) * dk/2.;
	// 	}
	// }
	// for (size_t i=0; i<nx; ++i)
	// 	cout<<x_val(i)<<' '<<exK(i)<<' '<<sqrt(exK2(i))
	// 		<<' '<<exK(i)*cdX(i)/f.get_m()*AU_Acm2
	// 		<<' '<<sqrt(exK2(i))*cdX(i)/f.get_m()*AU_Acm2<<endl;

	//  nx_1 (f.nx_, arma::fill::zeros), nx_2 (f.nx_, arma::fill::zeros);
	//  jx_1 (f.nx_, arma::fill::zeros), jx_2 (f.nx_, arma::fill::zeros);
	// for (size_t i=0; i<f.nx_; ++i) {
	// 	for (size_t j=0; j<f.nk2_; ++j) {
	// 		nx_1(i) += f.f_(i,j)*f.dk_;
	// 		jx_1(i) += f.k_(j)/f.m_*f.f_(i,j)*f.dk_;
	// 	}
	// 	for (size_t j=f.nk2_; j<f.nk_; ++j) {
	// 		nx_2(i) += f.f_(i,j)*f.dk_;
	// 		jx_2(i) += f.k_(j)/f.m_*f.f_(i,j)*f.dk_;
	// 	}
	// }

	arma::vec an_pot (nx, arma::fill::zeros);
	double sig = 20/AU_nm, x0 = 35/AU_nm;
	for (size_t i=0; i<nx; ++i)
		an_pot(i) = 0.03/AU_eV*exp(-(x_val(i)-x0)*(x_val(i)-x0)/sig/sig)
		*(-8/pow(sig,6)*pow(x_val(i)-x0,3)+12/pow(sig,4)*(x_val(i)-x0));

	field<std::string> header(8);
	arma::mat out_data;
	out_data.insert_cols(0, x_val*AU_nm), header(0) = "x [nm]";
	out_data.insert_cols(1, f.get_u()*AU_eV), header(1) = "U [eV]";  // col. 2
	out_data.insert_cols(2, f.get_currD()*AU_Acm2), header(2) = "J(x) [Acm^{-2}]";  // col. 3
	out_data.insert_cols(3, cdX/AU_cm3), header(3) = "n [cm^{-3}]";  // col. 4
	out_data.insert_cols(4, f.get_du()*AU_eV/AU_nm), header(4) = "U' [eV/nm]";  // col. 5
	out_data.insert_cols(5, f.get_d3u()), header(5) = "U''' [au]";  // col. 6
	out_data.insert_cols(6, f.get_uB()*AU_eV), header(6) = "U^B [eV]";  // col. 7
	out_data.insert_cols(7, f.get_uC()*AU_eV), header(7) = "U^C [eV]";  // col. 7
	// out_data.insert_cols(5, nx_1/AU_cm3), header(5) = "nx+[cm^-3]";  // col. 6
	// out_data.insert_cols(6, nx_2/AU_cm3), header(6) = "nx-[cm^-3]";  // col. 7
	// out_data.insert_cols(7, (nx_2-nx_1)/AU_cm3), header(7) = "(nx+)-(nx-)[cm^-3]";  // col. 8
	// out_data.insert_cols(8, jx_1*AU_Acm2), header(8) = "J+[Acm^-2]";  // col. 9
	// out_data.insert_cols(9, jx_2*AU_Acm2), header(9) = "J-[Acm^-2]";  // col. 10
	// out_data.insert_cols(10, (jx_2+jx_1)*AU_Acm2), header(10) = "(J+)+(J-)[Acm^-2]";  // col. 11
	out_data.save( csv_name("out_data/test.csv", header) );

	std::ofstream file;
	file.open("out_data/cdX.out", std::ios::out);
	file<<"# Carrier density in 'x' space\n";
	file<<"# x [au]  n(x) [au]\n";
	for (size_t i=0; i<nx; ++i)
		file<<x_val(i)<<'\t'<<cdX(i)<<'\n';
	file.close();

	file.open("out_data/cdK.out", std::ios::out);
	file<<"# Carrier density in 'k' space\n";
	file<<"# p [au]  n(k) [au]\n";
	for (size_t j=0; j<nk; ++j)
		file<<k_val(j)<<'\t'<<cdK(j)<<'\n';
	file.close();

	cout<<"# Final current = "<<curr*AU_Acm2<<" [Acm^-2]"<<endl;
	cout<<"# N = "<<f.calcNorm()/AU_cm2<<" [cm^-2]"<<endl;
	cout<<"# int{dp} f_BC = "<<calcInt(f.get_bc(), f.get_dk())/2./M_PI/AU_cm3<<" [a.u.]"<<endl;

	t_end = high_resolution_clock::now();
	t_elapsed =  duration_cast<duration<double>>(t_end - t_start);
	if (t_elapsed.count() < 60)
		cout<<"# RUN TIME: "<<t_elapsed.count()<<" s, "<<endl;
	else
		cout<<"# RUN TIME: "<<t_elapsed.count()/60.<<" min, "<<endl;

	return 0;

}


/*
void gaussBarrier(WignerFunction& f) {

	size_t nx = f.get_nx(), nk = f.get_nk();
	double m = f.get_m();
	arma::mat f;

	arma::vec nx_1 (nx, arma::fill::zeros), nx_2 (nx, arma::fill::zeros);
	arma::vec jx_1 (nx, arma::fill::zeros), jx_2 (nx, arma::fill::zeros);


	// f.setPotBias(0.02/AU_eV);
	f.uB_.zeros();
	f.addGaussBarr(0.1/AU_eV, 500/AU_nm, 100/AU_nm);

	size_t n = 30;
	double u0_start = 0, u0_end = 0.15/AU_eV;
	double du0 = (u0_end-u0_start)/(float(n)-1.);

	std::ofstream out_data("out_data/gaussBarrier.out");
	for (size_t l=0; l<n; ++l) {
		// f.uB_.zeros();
		// f.addGaussBarr(l*du0, 500/AU_nm, 100/AU_nm);
		f.setPotBias(l*du0);
		f.solveWignerEq();
		f.calcCurr(), f.calcCD_X();
		nx_1.zeros(), nx_2.zeros();
		jx_1.zeros(), jx_2.zeros();
		for (size_t i=0; i<f.nx_; ++i) {
			for (size_t j=0; j<f.nk2_; ++j) {
				nx_1(i) += f.f_(i,j)*f.dk_;
				jx_1(i) += f.k_(j)/f.m_*f.f_(i,j)*f.dk_;
			}
			for (size_t j=f.nk2_; j<f.nk_; ++j) {
				nx_2(i) += f.f_(i,j)*f.dk_;
				jx_2(i) += f.k_(j)/f.m_*f.f_(i,j)*f.dk_;
			}
		}
		for (size_t i=0; i<f.nx_; ++i) {
			out_data<<l*du0*AU_eV  // col. 1
				<<'\t'<<f.x_(i)*AU_nm  // col. 2
				<<'\t'<<f.u_(i)*AU_eV  // col. 3
				<<'\t'<<f.currD_(i)*AU_Acm2  // col. 4
				<<'\t'<<f.cdX_(i)/AU_cm3  // col. 5
				<<'\t'<<(jx_2(i)+jx_1(i))*AU_Acm2  // col. 6
				<<'\t'<<(nx_2(i)-nx_1(i))/AU_cm3<<endl;  // col. 7
		}
		out_data<<endl;
	}
	out_data.close();

	// std::ofstream uBarr_uBias_map("out_data/uBarr_uBias_map.out");
	// for (size_t l=0; l<n; ++l) {
	// 	f.uB_.zeros();
	// 	f.addGaussBarr(l*du0, 500/AU_nm, 100/AU_nm);
	// 	f.calc_IVchar(u0_start, u0_end, n);
	// 	for (size_t i=0; i<n; ++i)
	// 		uBarr_uBias_map<<l*du0*AU_eV
	// 			<<'\t'<<f.iv_v_(i)*AU_eV
	// 			<<'\t'<<f.iv_i_(i)*AU_Acm2
	// 			<<'\t'<<f.iv_iRange_(i)*AU_Acm2
	// 			<<'\t'<<f.iv_n_(i)*AU_Acm2<<endl;
	// 	uBarr_uBias_map<<endl;
	// }
	// uBarr_uBias_map.close();
}


void dissDecoh(WignerFunction& f, double t_start, double t_end, size_t n_step, double v_min, double v_max, size_t nv) {

	double t_step = (log10(t_end) - log10(t_start))/float(n_step-1);

	// f.solveWignerEq();
	f.solveWignerPoisson();
	double j0 = f.calcCurr(), n0 = f.calcNorm();
	double j=0, n=0, t=0;
	arma::vec fit, nE_k;
	arma::mat tpMap(n_step, f.nk_, arma::fill::zeros);
	arma::mat ivMap(n_step, nv, arma::fill::zeros);

	bool calcIVchar = false;
	double r0 = 0, r = 0;

	if (calcIVchar) {
		f.calc_IVchar(v_min, v_max, nv);
		fit = polyfit(f.iv_v_, f.iv_i_, 1);
		r0 = 1/fit(0);
	}

	std::ofstream out;
	out.open("out_data/diss_tau.out", std::ios::out);
	out<<"## j0 = "<<j0*AU_Acm2<<"  n0 = "<<n0/AU_cm2<<endl;
	cout<<"## j0 = "<<j0*AU_Acm2<<"  n0 = "<<n0/AU_cm2<<endl;
	out<<"# i  tS [s]  j/j0  n/n0  R\n";
	cout<<"# i  tS [s]  j/j0  n/n0  R\n";
	for (size_t i = 0; i < n_step; ++i) {
		t = t_start*pow(10,(n_step-1-i)*t_step);
		// f.rR_ = 1/(t/AU_s);
		f.rM_ = 1/(t/AU_s);
		// f.solveWignerEq();
		f.solveWignerPoisson();
		nE_k = f.calcCD_K();
		for (size_t k=0; k<f.nk_; ++k)
			tpMap(i, k) = nE_k(k);
		j = f.calcCurr(), n = f.calcNorm();
		//
		// I-V char
		if (calcIVchar) {
			f.calc_IVchar(v_min, v_max, nv);
			fit = polyfit(f.iv_v_, f.iv_i_, 1);
			r = 1/fit(0);
			for (size_t l=0; l<nv; ++l)
				ivMap(i, l) = f.iv_i_(l);
		}
		// out<<i<<'\t'<<t<<'\t'<<j/j0<<'\t'<<1./fit(0)<<endl;
		// cout<<i<<'\t'<<t<<'\t'<<j/j0<<'\t'<<1./fit(0)<<endl;
		out<<i<<'\t'<<t<<'\t'<<j/j0<<'\t'<<n/n0;
		if (calcIVchar) out<<'\t'<<r/r0;
		out<<endl;
		cout<<i<<'\t'<<t<<'\t'<<j/j0<<'\t'<<n/n0;
		if (calcIVchar) cout<<'\t'<<r/r0;
		cout<<endl;
		f.saveWignerFun();
	}
	out.close();

	std:: ofstream tpMap_out;
	tpMap_out.open("out_data/tpMap.out", std::ios::out);
	tpMap_out<<"# tau [s]  p [au]  f(p) [au]\n";
	for (size_t i = 0; i < n_step; ++i) {
		t = t_start*pow(10,(n_step-1-i)*t_step);
		for (size_t k=0; k<f.nk_; ++k) {
			tpMap_out<<t<<' '<<f.k_(k)
				// <<'\t'<<f.f_(size_t(f.nx_/2.),k)
				<<'\t'<<tpMap(i,k)
				<<'\n';
		}
		tpMap_out<<'\n';
	}
	tpMap_out.close();

	if (calcIVchar) {
		std:: ofstream ivMap_out;
		ivMap_out.open("out_data/ivMap.out", std::ios::out);
		ivMap_out<<"# tau [s]  p [au]  f(p) [au]\n";
		for (size_t i = 0; i < n_step; ++i) {
			t = t_start*pow(10,(n_step-1-i)*t_step);
			for (size_t l=0; l<nv; ++l) {
				ivMap_out<<t<<' '<<f.iv_v_(l)
					// <<'\t'<<f.f_(size_t(f.nx_/2.),k)
					<<'\t'<<ivMap(i,l)
					<<'\n';
			}
			ivMap_out<<'\n';
		}
		ivMap_out.close();
	}

}


void calcTimeChar(WignerFunction& f) {

	//
	// Pakiet gaussowski
	f.gwp_x0_ = 100/AU_nm;
	f.gwp_dx_ = 25/AU_nm;
	f.gwp_dp_ = 0.005; // 1./(2*f.gwp_dx_);  // a.u.
	f.gwp_p0_ = sqrt(2*f.m_*f.uL_);
	f.addWavePacket();
	f.gwp_x0_ = f.l_-100/AU_nm;
	f.gwp_p0_ = -sqrt(2*f.m_*f.uR_);
	f.addWavePacket();

	//
	// Time evolution
	double t_total = 2.05e-12/AU_s, t = 0, max_dj = 1e-4;
	arma::vec dj(f.nx_, arma::fill::zeros), j0(f.nx_, arma::fill::zeros), j1(f.nx_, arma::fill::zeros);
	int n_min = 0, n_conv = 2, n_dj = 0, n_it = 0;
	double nrmse_j = 0;
	double j = 0, nc = 0, nc_prev = 0, dj_x = 0, j1_x = 0, j1_n = 0;
	bool conv_J = false;
	std::ofstream out;
	out.open("out_data/tev.out", std::ios::out);
	cout<<"# t [ps] \tCurr. [Acm^-2] \tN [cm^-2] \tdN \tdJ [Acm^-2] \tdJ/J \tCONV? \tn_dj"<<endl;
	out<<"# t [ps] \t Curr. [Acm^-2] \tN [cm^-2] \tdN \tdJ [Acm^-2] \tdJ/J \tCONV? \tn_dj"<<endl;
	while (t<t_total && n_dj < n_conv) {  //
		f.solveTimeEv();
		if (n_it > n_min) {
			j = f.calcCurr();
			nc_prev = nc, nc = f.calcNorm();
			j0 = j1, j1 = f.currD_;
			dj = j1 - j0;
			dj_x = max(dj);
			conv_J = fabs(dj_x/j) > max_dj ? false : true;
			n_dj = conv_J ? n_dj+1 : 0;
			out<<t*AU_s*1e12
				<<'\t'<<j
				<<'\t'<<nc
				<<'\t'<<nc-nc_prev
				<<'\t'<<dj_x*AU_Acm2
				<<'\t'<<dj_x/j
				<<'\t'<<conv_J
				<<'\t'<<n_dj<<endl;
			cout<<t*AU_s*1e12
				<<'\t'<<j
				<<'\t'<<nc
				<<'\t'<<nc-nc_prev
				<<'\t'<<dj_x*AU_Acm2
				<<'\t'<<dj_x/j
				<<'\t'<<conv_J
				<<'\t'<<n_dj<<endl;
		}
		f.saveWignerFun();
		t += f.dt_, n_it++;
	}
	out.close();

}
*/
