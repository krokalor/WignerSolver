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
	f.bcType_ = -3, f.rG_ = 0;  // 1./(1e-15/AU_s);

	//
	// Potencjał
	// cout<<"# Setting up potential"<<endl;
	// f.setPotBias(0.02/AU_eV);
	// f.addGaussPot(0.1/AU_eV, 500/AU_nm, 100/AU_nm);
	// f.readPotential("pot.in");
	// f.u_ = f.u_ + f.uStart_;
	// f.uStart_ = f.uC_;
	// f.set_useNLP(false);

	//
	// FUNKCJA RÓWNOWAGOWA
	// f.setEquilibriumFunction("pot.in", false);

	//
	// Dyssypacja
	// f.rR_ = 1./(1e-12/AU_s);  // 1./(1e-12/AU_s);
	// f.rM_ = 1./(1e-11/AU_s);  // 1./(1e-12/AU_s);
	// f.lambda_ = 0*AU_nm*AU_nm*AU_s;  // [nm^-2*s^-1]

	f.printParam();

	//
	// Boltzmann
	// cout<<"# Solving BTE"<<endl;
	// f.solveWignerEq();
	// f.saveWignerFun();

	//
	// Poisson test
	// Poisson1D p(200, 1/AU_nm);
	// p.nC_ = -1*normpdf(linspace(0, 200, 200), 100, 10), p.q_ = -1;
	// p.solve();
	// p.nC_.print("# rho:"), p.uNew_.print("\n\n # U^H:");
	// p.testPoisson();

	//
	// Boltzmann-Poisson
	// cout<<"# Solving BTE+PE"<<endl;
	// f.solveWignerPoisson();

	//
	// Barriera gaussowska
	gaussBarrier(f);

	//
	// Zalezność funkcji rozkładu od dyssypacji
	// dissDecoh(f, 1e-13, 1e-10, 20, 0.0/AU_eV, 0.001/AU_eV, 20);

	//
	// Charakterystyka J-V
	// f.calc_IVchar(0.0/AU_eV, 0.04/AU_eV, 40);
	// vec fit = polyfit(f.iv_v_, f.iv_i_, 1);
	// fit.print();

	f.calcCurr(), f.calcCD_X(), f.calcCD_K();

	vec nx_1 (f.nx_, fill::zeros), nx_2 (f.nx_, fill::zeros);
	vec jx_1 (f.nx_, fill::zeros), jx_2 (f.nx_, fill::zeros);
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

	mat out_data;
	out_data.insert_cols(0, f.x_*AU_nm);
	out_data.insert_cols(1, f.u_*AU_eV);  // col. 2
	out_data.insert_cols(2, f.currD_*AU_Acm2);  // col. 3
	out_data.insert_cols(3, f.cdX_/AU_cm3);  // col. 4
	out_data.insert_cols(4, f.du_*AU_eV/AU_nm);  // col. 5
	out_data.insert_cols(5, nx_1/AU_cm3);  // col. 6
	out_data.insert_cols(6, nx_2/AU_cm3);  // col. 7
	out_data.insert_cols(7, (nx_2-nx_1)/AU_cm3);  // col. 8
	out_data.insert_cols(8, jx_1*AU_Acm2);  // col. 9
	out_data.insert_cols(9, jx_2*AU_Acm2);  // col. 10
	out_data.insert_cols(10, (jx_2+jx_1)*AU_Acm2);  // col. 11
	out_data.print("##  x [nm]  U [eV]  j [Acm^-2]  n(x) [cm^-3]  du/dx [eV/nm]"
		"  nx(k<0) [cm^-3]  nx(k>0) [cm^-3]  nx_diff [cm^3]  jx(k<0) [Acm^-2]  jx(k>0) [Acm^-2]  jx [Acm^-2]");

	std::ofstream file;
	file.open("out_data/nE_x.out", std::ios::out);
	file<<"# x [au]  n(x) [au]\n";
	for (size_t i=0; i<f.nx_; ++i)
		file<<f.x_(i)<<' '<<f.cdX_(i)<<'\n';
	file.close();

	file.open("out_data/nE_k.out", std::ios::out);
	file<<"# p [au]  n(k) [au]\n";
	for (size_t j=0; j<f.nk_; ++j)
		file<<f.k_(j)<<' '<<f.cdK_(j)<<'\n';
	file.close();

	file.open("out_data/pot.out", std::ios::out);
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

void gaussBarrier(WignerFunction& f) {
	vec nx_1 (f.nx_, fill::zeros), nx_2 (f.nx_, fill::zeros);
	vec jx_1 (f.nx_, fill::zeros), jx_2 (f.nx_, fill::zeros);

	// f.setPotBias(0.02/AU_eV);
	f.uB_.zeros();
	f.addGaussPot(0.1/AU_eV, 500/AU_nm, 100/AU_nm);

	size_t n = 40;
	double u0_start = 0, u0_end = 0.15/AU_eV;
	double du0 = (u0_end-u0_start)/(float(n)-1.);

	std::ofstream out_data("out_data/gaussBarrier.out");
	for (size_t l=0; l<n; ++l) {
		// f.uB_.zeros();
		// f.addGaussPot(l*du0, 500/AU_nm, 100/AU_nm);
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
	out.open("out_data/diss_tau.out", std::ios::out);
	out<<"## j0 = "<<j0*AU_Acm2<<"  n0 = "<<n0/AU_cm2<<endl;
	cout<<"## j0 = "<<j0*AU_Acm2<<"  n0 = "<<n0/AU_cm2<<endl;
	out<<"# i  tS [s]  j/j0  R\n";
	cout<<"# i  tS [s]  j/j0  R\n";
	for (size_t i = 0; i < n_step; ++i) {
		t = t_start*pow(10,(n_step-1-i)*t_step);
		f.rR_ = 1/(t/AU_s);
		// f.rM_ = 1/(t/AU_s);
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
		out<<i<<'\t'<<t<<'\t'<<j/j0<<'\t'<<n/n0<<endl;
		cout<<i<<'\t'<<t<<'\t'<<j/j0<<'\t'<<n/n0<<endl;
		f.saveWignerFun();
	}
	out.close();

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
