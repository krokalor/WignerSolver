#include "src/lib.hpp"
#include "src/WignerFunction.hpp"
#include "src/poisson1D.hpp"

#include <chrono>


using namespace wigner;
using namespace poisson;
using namespace std::chrono;


void dissDecoh(WignerFunction&, double, double, size_t, double, double, size_t);

void simGWP(WignerFunction&);
void simGWP(WignerFunction&, WignerFunction&, WignerFunction&);


int main(){

	// omp_set_dynamic(0);
	// omp_set_num_threads(1);
	// cout<<"omp_get_max_threads(): "<<omp_get_max_threads() <<endl;
	// cout<<"omp_get_num_threads(): "<<omp_get_num_threads()<<endl;

	high_resolution_clock::time_point t_start, t_end;
	duration<double> t_elapsed;
	t_start = high_resolution_clock::now();

	WignerFunction f;

	f.set_m(0.067);
	f.set_temp(77);
	f.set_uF(0.087/AU_eV);  // 0.087/AU_eV
	f.set_nx(150), f.set_nk(150);
	f.set_lD(600/AU_nm), f.set_lC(200/AU_nm);
	f.set_kmax(0.04);
	f.set_dt(5*1e-15/AU_s);
	f.set_cD(2e18*AU_cm3); // 2e18*AU_cm3
	// f.set_lYZ(1e-8/AU_cm2);
	// f.set_part_num(1);

	f.update();

	//
	// Warunek brzegowy
	// -1 -> 0, 0 -> SF, 1-3 -> splot, -2 -> MB, -3 -> Gauss
	//
	f.bcType_ = -3, f.rG_ = 0;  // 1./(1e-15/AU_s);

	//
	// FUNKCJA RÓWNOWAGOWA
	//
	f.setEquilibriumFunction("potentials/pot_0meV_100nm_cl.in", false);

	//
	// Pakiet gaussowski
	//
	// f.gwp_x0_ = 100/AU_nm;
	// f.gwp_dx_ = 25/AU_nm;
	// f.gwp_dp_ = 0.005; // 1./(2*f.gwp_dx_);  // a.u.
	// f.gwp_p0_ = sqrt(2*f.m_*f.uF_);
	// f.addWavePacket();
	// f.gwp_x0_ = f.l_-100/AU_nm;
	// f.gwp_p0_ = -sqrt(2*f.m_*f.uF_);
	// f.addWavePacket();

	//
	// Dyssypacja
	//
	// f.rR_ = 1./(1e-13/AU_s);  // 1./(1e-12/AU_s);
	// f.rM_ = 1./(1e-12/AU_s);  // 1./(1e-12/AU_s);
	// f.lambda_ = 0*AU_nm*AU_nm*AU_s;  // [nm^-2*s^-1]

	//
	// Potencjał
	//
	// cout<<"# Setting up potential"<<endl;
	f.setPotBias(0.02/AU_eV);
	// f.addGaussPot(0.01/AU_eV, 50/AU_nm, 10/AU_nm);
	// f.readPotential("potentials/pot_20meV_1um_cl.in");
	// f.u_ = f.u_ + f.uStart_;
	f.uStart_ = f.uC_;
	// f.set_useNLP(false);

	f.printParam();

	//
	// Boltzmann
	//
	// cout<<"# Solving BTE"<<endl;
	// f.solveWignerEq();
	// // // calc_fd(f);
	// // // f.set_time_dependent(true);
	// // // for (size_t i=0; i < 100; ++i) f.solveTimeEv();
	// f.saveWignerFun();

	//
	// Poisson test
	//
	// Poisson1D p(200, 1/AU_nm);
	// p.nC_ = -1*normpdf(linspace(0, 200, 200), 100, 10), p.q_ = -1;
	// p.solve();
	// p.nC_.print("# rho:"), p.uNew_.print("\n\n # U^H:");
	// p.testPoisson();

	//
	// Boltzmann-Poisson
	//
	cout<<"# Solving BTE+PE"<<endl;
	f.solveWignerPoisson();

	//
	// Zalezność funkcji rozkładu od dyssypacji
	// diss_type = 'R', 'M', 'F'
	//
	// dissDecoh(f, 1e-13, 1e-10, 20, 0.0/AU_eV, 0.001/AU_eV, 20);

	//
	// Time evolution
	//
	// double t_total = 2.05e-12/AU_s, t = 0, max_dj = 1e-4;
	// vec dj(f.nx_, fill::zeros), j0(f.nx_, fill::zeros), j1(f.nx_, fill::zeros);
	// int n_min = 0, n_conv = 2, n_dj = 0, n_it = 0;
	// double nrmse_j = 0;
	// double j = 0, nc = 0, nc_prev = 0, dj_x = 0, j1_x = 0, j1_n = 0;
	// bool conv_J = false;
	// std::ofstream out;
	// out.open("wyniki/dane/tev.out", std::ios::out);
	// cout<<"# t [ps] \tCurr. [Acm^-2] \tN [cm^-2] \tdN \tdJ [Acm^-2] \tdJ/J \tCONV? \tn_dj"<<endl;
	// out<<"# t [ps] \t Curr. [Acm^-2] \tN [cm^-2] \tdN \tdJ [Acm^-2] \tdJ/J \tCONV? \tn_dj"<<endl;
	// while (t<t_total && n_dj < n_conv) {  //
	// 	f.solveTimeEv();
	// 	if (n_it > n_min) {
	// 		j = f.calcCurr();
	// 		nc_prev = nc, nc = f.calcNorm();
	// 		j0 = j1, j1 = f.calcCurrArr();
	// 		dj = j1 - j0;
	// 		dj_x = max(dj);
	// 		conv_J = fabs(dj_x/j) > max_dj ? false : true;
	// 		n_dj = conv_J ? n_dj+1 : 0;
	// 		out<<t*AU_s*1e12
	// 			<<'\t'<<j
	// 			<<'\t'<<nc
	// 			<<'\t'<<nc-nc_prev
	// 			<<'\t'<<dj_x*AU_Acm2
	// 			<<'\t'<<dj_x/j
	// 			<<'\t'<<conv_J
	// 			<<'\t'<<n_dj<<endl;
	// 		cout<<t*AU_s*1e12
	// 			<<'\t'<<j
	// 			<<'\t'<<nc
	// 			<<'\t'<<nc-nc_prev
	// 			<<'\t'<<dj_x*AU_Acm2
	// 			<<'\t'<<dj_x/j
	// 			<<'\t'<<conv_J
	// 			<<'\t'<<n_dj<<endl;
	// 	}
	// 	f.saveWignerFun();
	// 	t += f.dt_, n_it++;
	// }
	// out.close();

	//
	// Charakterystyka J-V
	//
	// f.calc_IVchar(0.0/AU_eV, 0.02/AU_eV, 10);
	// vec fit = polyfit(f.iv_v_, f.iv_i_, 1);
	// fit.print();
	// f.saveWignerFun();

	cout<<"# Calulating electron density"<<endl;
	vec nE_x = f.calcCD_X();
	vec nE_k = f.calcCD_K();

	std::ofstream file;
	file.open("wyniki/dane/nE_x.out", std::ios::out);
	file<<"# x [au]  n(x) [au]\n";
	for (size_t i=0; i<f.nx_; ++i)
		file<<f.x_(i)<<' '<<nE_x(i)<<'\n';
	file.close();

	file.open("wyniki/dane/nE_k.out", std::ios::out);
	file<<"# p [au]  n(k) [au]\n";
	for (size_t i=0; i<f.nk_; ++i)
		file<<f.k_(i)<<' '<<nE_k(i)<<'\n';
	file.close();

	file.open("wyniki/dane/pot.out", std::ios::out);
	file<<"# x [au]  U [au]\n";
	for (size_t i=0; i<f.nx_; ++i)
		file<<f.x_(i)<<' '<<f.u_(i)<<'\n';
	file.close();

	cout<<"# Final current = "<<f.calcCurr()*AU_Acm2<<" [Acm^-2]"<<endl;
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


void dissDecoh(WignerFunction& f, double t_start, double t_end, size_t n_step, double v_min, double v_max, size_t nv) {

	double t_step = (log10(t_end) - log10(t_start))/float(n_step-1);

	f.solveWignerEq();
	// f.solveWignerPoisson();
	double j0 = f.calcCurr(), n0 = f.calcNorm();
	double j=0, n=0, t=0;
	vec fit, nE_k;
	mat tpMap(n_step, f.nk_, fill::zeros);

	std::ofstream out, tpMap_out;
	out.open("wyniki/dane/diss_tau.out", std::ios::out);
	out<<"## j0 = "<<j0*AU_Acm2<<"  n0 = "<<n0/AU_cm2<<endl;
	cout<<"## j0 = "<<j0*AU_Acm2<<"  n0 = "<<n0/AU_cm2<<endl;
	out<<"# i  tS [s]  j/j0  R\n";
	cout<<"# i  tS [s]  j/j0  R\n";
	for (size_t i = 0; i < n_step; ++i) {
		t = t_start*pow(10,(n_step-1-i)*t_step);
		// f.rR_ = 1/(t/AU_s);
		f.rM_ = 1/(t/AU_s);
		f.solveWignerEq();
		// f.solveWignerPoisson();
		nE_k = f.calcCD_K();
		for (size_t k=0; k<f.nk_; ++k)
			tpMap(i, k) = nE_k(k);
		j = f.calcCurr(), n = f.calcNorm();
		// I-V char
		// f.calc_IVchar(v_min, v_max, nv);
		// fit = polyfit(f.iv_v_, f.iv_i_, 1);
		// out<<i<<'\t'<<t<<'\t'<<j/j0<<'\t'<<1./fit(0)<<endl;
		// cout<<i<<'\t'<<t<<'\t'<<j/j0<<'\t'<<1./fit(0)<<endl;
		out<<i<<'\t'<<t<<'\t'<<j/j0<<endl;
		cout<<i<<'\t'<<t<<'\t'<<j/j0<<endl;
		f.saveWignerFun();
	}
	out.close();

	tpMap_out.open("wyniki/dane/tpMap.out", std::ios::out);
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

}


/*
void simGWP(WignerFunction& f, WignerFunction& f2, WignerFunction& f3) {
	// Wave packet simulation

	size_t nt=int(40);
	double t;
	double dp0, dx0;

	f.f_.zeros(), f.addWavePacket();
	f2.f_.zeros(), f2.addWavePacket();
	f3.f_.zeros(), f3.addWavePacket();

	cout<<"# t x0 p0 dx dp"<<endl;
	for (size_t l = 0; l < nt+1; ++l) {
		// f.solveWignerEq();
		f.solveTimeEv(), f2.solveTimeEv(), f3.solveTimeEv();
		t = l*f.dt_;
		dp0 = f.gwp_dp_, dx0 = f.gwp_dx_;
		cout<<t*AU_s*1e15<<' '<<f2.calcSDX()/dx0<<' '<<f2.calcSDK()/dp0;
		cout<<' '<<f3.calcSDX()/dx0<<' '<<f3.calcSDK()/dp0;
		cout<<' '<<dx0*dp0<<' '<<f2.calcSDX()*f2.calcSDK()<<' '<<f3.calcSDX()*f3.calcSDK()<<endl;
	}
	cout<<"\n\n";

}
*/


/*
void simGWP(WignerFunction& f) {
	// Wave packet simulation

	size_t nt=int(40);
	double t, vt;

	f.f_.zeros();
	f.addWavePacket();

	cout<<"# t x0 p0 dx dp"<<endl;
	for (size_t l = 0; l < nt+1; ++l) {
		// f.solveWignerEq();
		f.solveTimeEv();
		t = l*f.dt_;
		vt = f.gwp_p0_*t/f.m_;
		cout<<t*AU_s*1e15<<' '<<f.calcSDX()/f.gwp_dx_<<' '<<f.calcSDK()/f.gwp_dp_<<' '<<f.calcSDX()*f.calcSDK()<<endl;
		// cout<<t*AU_s*1e15<<' '<<f.calcEX()<<' '<<f.calcEK()<<' '<<f.calcSDX()<<' '<<f.calcSDK()<<endl;
		// if (l%1 == 0) {
		//   cout<<"# t = "<<t<<endl;
		//   cout<<"# N = "<<f.calcNorm()<<endl;
		//   cout<<"# xt "<<f.gwp_x0_+vt<<"  E[x] "<<f.calcEX()<<endl;
		//   cout<<"# p0 "<<f.gwp_p0_<<"  E[p] "<<f.calcEK()<<endl;
		//   cout<<"# dp "<<1/f.gwp_dx_<<"  SD[p] "<<f.calcSDK()<<endl;
		//   cout<<"# dx "<<2/f.gwp_dx_*t/f.m_<<"  SD[x] "<<f.calcSDX()<<endl;
		//   // cout<<"# x f(x)"<<endl;
		//   // for (size_t i=0; i<f.nx_; ++i)
		//   //   cout<<f.x_(i)<<' '<<f.f_(i, f.nk2_+int(0.02/f.dk_))<<endl;
		//   // cout<<"\n\n";
		// }
	}
	cout<<"\n\n";

	//
	// TEST
	//
	cout<<"# t x0 p0 dx dp"<<endl;
	for (size_t l = 0; l < nt+1; ++l) {
		// f.solveWignerEq();
		t = l*f.dt_;
		vt = f.gwp_p0_*t/f.m_;

		double s2 = 2*f.gwp_dp_*f.gwp_dp_;
		for (size_t i=0; i<f.nx_; ++i) {
			for (size_t j=0; j<f.nk_; ++j){
				double vt = f.k_(j)*t/f.m_;
				f.f_(i,j) = exp(
					- (f.k_(j)-f.gwp_p0_)*(f.k_(j)-f.gwp_p0_)/s2
					- (f.x_(i)-f.gwp_x0_-vt)*(f.x_(i)-f.gwp_x0_-vt)*s2  )/M_PI; // * gwp_A_;
			}
		}

		if (l%1 == 0) {
			cout<<t*AU_s*1e15<<' '<<f.calcEX()<<' '<<f.calcEK()<<' '<<f.calcSDX()<<' '<<f.calcSDK()<<endl;
			// cout<<"# t = "<<t<<endl;
			// cout<<"# N = "<<f.calcNorm()<<endl;
			// cout<<"# xt "<<f.gwp_x0_+vt<<"  E[x] "<<f.calcEX()<<endl;
			// cout<<"# p0 "<<f.gwp_p0_<<"  E[p] "<<f.calcEK()<<endl;
			// cout<<"# dp "<<f.gwp_dp_<<"  SD[p] "<<f.calcSDK()<<endl;
			// cout<<"# dx "<<2/f.gwp_dx_*t/f.m_<<"  SD[x] "<<f.calcSDX()<<endl;
			// cout<<"# x f(x)"<<endl;
			// for (size_t i=0; i<f.nx_; ++i)
			//   cout<<f.x_(i)<<' '<<f.f_(i, f.nk2_+int(0.02/f.dk_))<<endl;
			// cout<<"\n\n";
		}
	}
	//
	//
	//

}
*/
