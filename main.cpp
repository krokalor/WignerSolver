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

	double nx = 150, nk = 150;
	double lD = 1000/AU_nm, lC = 1500/AU_nm;
	double k_max = 0.15;  // -1

	double m = 0.067;
	double temp = 300;
	double cD = 2e18*AU_cm3;
	double uF = calcFermiEn(cD, m, temp);  // calcFermiEn(cD, m, temp) 0.087/AU_eV
	double dt = 5*1e-15/AU_s;

	double rR = 1/(1e-12/AU_s);  // 1/(1e-12/AU_s)
	double rM = 1/(1e-12/AU_s);
	double rF = 0;
	double rG = 0;
	double lambda = 0;

	WignerFunction f(nx, lD, lC, nk, k_max);

	vec x_arr = f.get_x_arr(), k_arr = f.get_k_arr();

	f.set_m(m);
	f.set_temp(temp);
	f.set_cD(cD);
	f.set_uF(uF);
	f.set_dt(dt);
	f.set_rR(rR), f.set_rM(rM), f.set_rF(rF), f.set_rG(rG), f.set_lambda(lambda);

	f.set_useQC(false);
	f.set_useNLP(false);
	f.set_uBias_BC(true);

	// Warunek brzegowy
	// 0 -> 0, 1 -> SF, 2:4 -> splot with SF, -1 -> Gauss, -2:-4 splot with Gauss
	cout<<"Set up BC"<<endl;
	f.set_bcType(1);

	//
	// Potencjał
	cout<<"# Setting up potential"<<endl;
	f.set_uBias(0./AU_eV);
	// f.setPotBias(0.1/AU_eV);
	// f.addRectBarr(0.3/AU_eV, 17.5/AU_nm, 2/AU_nm, 10);
	// f.addRectBarr(0.3/AU_eV, 22.5/AU_nm, 2/AU_nm, 10);
	// f.addGaussBarr(0.3/AU_eV, 2000/AU_nm, 500/AU_nm);
	// f.readPotential("out_data/poisson_pot/poisson_pot_2um_0eV.in");
	// f.set_uC(f.get_uStart());

	//
	// FUNKCJA RÓWNOWAGOWA
	f.setEquilibriumFunction("out_data/wf_feq_p.bin", true);

	f.printParam();

	//
	// Boltzmann
	// cout<<"# Solving BTE"<<endl;
	// f.solveWignerEq();
	// f.saveWignerFun();

	// // #pragma omp parallel for schedule(dynamic)
	// for (size_t i=0; i<50; ++i) {
	// 	cout<<"solveWigner: "<<i<<' '<<omp_get_thread_num()<<endl;
	// 	// f.set_uBias(i*0.01);
	// 	f.solveWignerEq();
	// }

	//
	// Boltzmann-Poisson
	cout<<"# Solving BTE+PE"<<endl;
	f.solveWignerPoisson(0.2/AU_eV, 4e-5, 1, 2000);  // uBias, alpha, beta, n_max
	f.saveWignerFun();

	//
	// Poisson test
	// Poisson1D p(200, 1/AU_nm);
	// p.rho_ = -1*normpdf(linspace(0, 200, 200), 100, 10);
	// p.solve();
	// p.rho_.print("# rho:"), p.uNew_.print("\n\n # U^H:");
	// p.testPoisson();

	//
	// Barriera gaussowska
	// gaussBarrier(f);

	//
	// Zalezność funkcji rozkładu od dyssypacji
	// dissDecoh(f, 1e-13, 1e-10, 5, 0.0/AU_eV, 0.001/AU_eV, 10);

	//
	// Charakterystyka J-V
	// f.calc_IVchar(0.0/AU_eV, 0.4/AU_eV, 5);
	// f.saveWignerFun();

	double curr = f.calcCurr();
	f.calcCD_X(), f.calcCD_K();
	vec cdX = f.get_cdX(), cdK = f.get_cdK();

	// vec nx_1 (f.nx_, fill::zeros), nx_2 (f.nx_, fill::zeros);
	// vec jx_1 (f.nx_, fill::zeros), jx_2 (f.nx_, fill::zeros);
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

	vec an_pot (nx, fill::zeros);
	double sig = 20/AU_nm, x0 = 35/AU_nm;
	for (size_t i=0; i<nx; ++i)
		an_pot(i) = 0.03/AU_eV*exp(-(x_arr(i)-x0)*(x_arr(i)-x0)/sig/sig)
		*(-8/pow(sig,6)*pow(x_arr(i)-x0,3)+12/pow(sig,4)*(x_arr(i)-x0));

	field<std::string> header(8);
	mat out_data;
	out_data.insert_cols(0, x_arr*AU_nm), header(0) = "x [nm]";
	out_data.insert_cols(1, f.get_u()*AU_eV), header(1) = "U [eV]";  // col. 2
	out_data.insert_cols(2, f.get_currD()*AU_Acm2), header(2) = "J(x) [Acm^{-2}]";  // col. 3
	out_data.insert_cols(3, cdX/AU_cm3), header(3) = "n [cm^{-3}]";  // col. 4
	out_data.insert_cols(4, f.get_du()*AU_eV/AU_nm), header(4) = "U' [eV/nm]";  // col. 5
	out_data.insert_cols(5, f.get_dddu()), header(5) = "U''' [au]";  // col. 6
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
		file<<x_arr<<'\t'<<cdX(i)<<'\n';
	file.close();

	file.open("out_data/cdK.out", std::ios::out);
	file<<"# Carrier density in 'k' space\n";
	file<<"# p [au]  n(k) [au]\n";
	for (size_t j=0; j<nk; ++j)
		file<<k_arr(j)<<'\t'<<cdK(j)<<'\n';
	file.close();

	// file.open("out_data/pot.out", std::ios::out);
	// file<<"# x [au]  U [au]\n";
	// for (size_t i=0; i<f.nx_; ++i)
	// 	file<<f.x_(i)<<' '<<f.u_(i)<<'\n';
	// file.close();

	cout<<"# Final current = "<<curr*AU_Acm2<<" [Acm^-2]"<<endl;
	cout<<"# N = "<<f.calcNorm()/AU_cm2<<" [cm^-2]"<<endl;
	// cout<<"# int{dp} f_BC = "<<calcInt(f.bc_, f.dk_)/AU_cm3<<" [a.u.]"<<endl;

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
	mat f;

	vec nx_1 (nx, fill::zeros), nx_2 (nx, fill::zeros);
	vec jx_1 (nx, fill::zeros), jx_2 (nx, fill::zeros);


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
	vec fit, nE_k;
	mat tpMap(n_step, f.nk_, fill::zeros);
	mat ivMap(n_step, nv, fill::zeros);

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
	f.gwp_p0_ = sqrt(2*f.m_*f.uF_);
	f.addWavePacket();
	f.gwp_x0_ = f.l_-100/AU_nm;
	f.gwp_p0_ = -sqrt(2*f.m_*f.uF_);
	f.addWavePacket();

	//
	// Time evolution
	double t_total = 2.05e-12/AU_s, t = 0, max_dj = 1e-4;
	vec dj(f.nx_, fill::zeros), j0(f.nx_, fill::zeros), j1(f.nx_, fill::zeros);
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
