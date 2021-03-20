#include "src/lib.hpp"
#include "src/WignerFunction.hpp"
#include "src/poisson1D.hpp"

#include <chrono>


using namespace wigner;
using namespace poisson;
using namespace std::chrono;


void simGWP(WignerFunction&);
void dissDecoh(WignerFunction&, double, double, size_t);
void dissDecoh_BC(WignerFunction&, double, double, size_t);
void simGWP(WignerFunction&, WignerFunction&, WignerFunction&);
void convolution(WignerFunction&);
double sfx(double, double, WignerFunction&);
void calc_fd(WignerFunction&);
template <class T>
void linear_fit(vec, vec, size_t);


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
	f.set_nx(150), f.set_nk(150);
	f.set_lD(600/AU_nm), f.set_lC(200/AU_nm);
	f.set_lYZ(1e-8/AU_cm2);
	f.set_part_num(1);
	f.set_cD(2e18*AU_cm3); // 2e18*AU_cm3
	f.set_kmax(0.04);
	f.set_temp(77);
	f.set_dt(.5*1e-15/AU_s);
	f.set_uF(0.087/AU_eV);  // 0.087/AU_eV
	f.set_time_dependent(false);

	f.update();

	//
	// Warunek brzegowy
	// -1 -> 0, 0 -> SF, 1-3 -> splot, -2 -> MB, -3 -> Gauss
	//
	f.bcType_ = -3, f.rG_ = 0;  // 1./(1e-15/AU_s);
	// convolution(f);

	//
	// Dyssypacja
	//
	// f.rR_ = 1./(1e-12/AU_s);  // 1./(1e-12/AU_s);
	// f.rM_ = 1./(1e-12/AU_s);  // 1./(1e-12/AU_s);
	// f.lambda_ = 0*AU_nm*AU_nm*AU_s;  // [nm^-2*s^-1]

	//
	// FUNKCJA RÓWNOWAGOWA
	//
	// f.setLinPot(0.);
	// f.readPotential("potentials/pot_0V_1um_cl.in");
	// f.u_ = f.u_ + f.uStart_;
	// f.solveWignerEq();
	// f.fe_ = f.f_;

	//
	// Pakiet gaussowski
	//
	// f.gwp_x0_ = 10/AU_nm;
	// f.gwp_dx_ = 2.5/AU_nm;
	// f.gwp_dp_ = 0.005; // 1./(2*f.gwp_dx_);  // a.u.
	// f.gwp_p0_ = sqrt(2*f.m_*f.uF_);
	// f.setWavePacket();
	// f.gwp_x0_ = f.l_-10/AU_nm;
	// f.gwp_p0_ = -sqrt(2*f.m_*f.uF_);
	// f.setWavePacket();

	//
	// Potencjał
	//
	// cout<<"# Setting up potential"<<endl;
	// f.setGaussPot(0.1/AU_eV, 800, 50);
	f.setLinPot(0.0/AU_eV);
	// f.uStart_.zeros();
	f.readPotential("potentials/pot_0V_1um_cl.in");
	// f.u_ = f.u_ + f.uStart_;
	// f.set_useNLP(false);

	f.printParam();

	//
	// Boltzmann
	//
	// cout<<"# Solving BTE"<<endl;
	// f.solveWignerEq();
	// calc_fd(f);
	// f.set_time_dependent(true);
	// for (size_t i=0; i < 100; ++i) f.solveWignerEq();
	// f.saveWignerFun();

	//
	// Poisson test
	//
	// Poisson1D p(200, 1./AU_nm);
	// p.testPoisson();

	//
	// Boltzmann-Poisson
	//
	cout<<"# Solving BTE+PE"<<endl;
	f.solveWignerPoisson();

	cout<<"# Calulating electron density"<<endl;
	vec nE_x = f.calcCD_X();
	vec nE_k = f.calcCD_K();

	// for (size_t i=0; i<f.nx_; ++i)
	// 	cout<<f.x_(i)*AU_nm<<'\t'<<nE_x(i)/AU_cm3<<endl;

	//
	// Zalezność gęstości prądu od dyssypacji
	// diss_type = 'R', 'M', 'F'
	//
	// dissDecoh(f, 1e-12, 1e-10, 10);
	// dissDecoh_BC(f, 1e-14, 1e-11, 10);

	//
	// Time evolution
	//
	// f.set_time_dependent(true);
	// double t_total = 2e-12/AU_s, t = 0, max_dj = 1e-4;
	// vec dj(f.nx_, fill::zeros), j0(f.nx_, fill::zeros), j1(f.nx_, fill::zeros);
	// int n_min = 2, n_conv = 2, n_dj = 0, n_it = 0;
	// double nrmse_j = 0;
	// double j = 0, nc = 0, nc_prev = 0, dj_x = 0, j1_x = 0, j1_n = 0;
	// bool conv_J = false;
	// std::ofstream out;
	// out.open("wyniki/dane/tev.out", std::ios::out);
	// cout<<"# t [ps] \tCurr. [Acm^-2] \tNorm. [cm^-2] \tdN \tdJ [Acm^-2] \tmax(J) [Acm^-2] \t|max(J) - min(J)| [Acm^-2] \tdJ/max(J) \tCONV? \tn_dj"<<endl;
	// out<<"# t [ps] \t Curr. [Acm^-2] \tNorm. [cm^-2] \tdN \tdJ [Acm^-2] \tmax(J) [Acm^-2] \t|max(J) - min(J)| [Acm^-2] \tdJ/max(J) \tCONV? \tn_dj"<<endl;
	// while (t<t_total && n_dj < n_conv) {  //
	// 	f.solveWignerEq();
	// 	if (n_it > n_min) {
	// 		j = f.calcCurr();
	// 		nc_prev = nc, nc = f.calcNorm();
	// 		j0 = j1, j1 = f.calcCurrArr();
	// 		dj = j1 - j0;
	// 		dj_x = max(dj), j1_x = max(j1), j1_n = min(j1);
	// 		conv_J = fabs(dj_x/j1_x) > max_dj ? false : true;
	// 		n_dj = conv_J ? n_dj+1 : 0;
	// 		out<<t*AU_s*1e12
	// 			<<'\t'<<j
	// 			<<'\t'<<nc
	// 			<<'\t'<<(nc-nc_prev)/nc
	// 			<<'\t'<<dj_x*AU_Acm2
	// 			<<'\t'<<j1_x*AU_Acm2
	// 			<<'\t'<<range(j1)*AU_Acm2
	// 			<<'\t'<<dj_x/j1_x
	// 			<<'\t'<<conv_J
	// 			<<'\t'<<n_dj<<endl;
	// 		cout<<t*AU_s*1e12
	// 			<<'\t'<<j
	// 			<<'\t'<<nc
	// 			<<'\t'<<(nc-nc_prev)/nc
	// 			<<'\t'<<dj_x*AU_Acm2
	// 			<<'\t'<<j1_x*AU_Acm2
	// 			<<'\t'<<range(j1)*AU_Acm2
	// 			<<'\t'<<dj_x/j1_x
	// 			<<'\t'<<conv_J
	// 			<<'\t'<<n_dj<<endl;
	// 	}
	// 	t += f.dt_, n_it++;
	// }
	// out.close();

	//
	// Charakterystyka J-V
	//

	// f.v_max_ = 0.001/AU_eV;
	// f.v_min_ = 0.0/AU_eV;
	// f.nv_ = 20;
	// f.calc_IVchar();
	// vec fit = polyfit(f.iv_v_, f.iv_i_, 1);
	// fit.print();

	std::ofstream file;
	file.open("wyniki/dane/nE_x.out", std::ios::out);
	file<<"# x [um]  n(x) [cm^-3]\n";
	for (size_t i=0; i<f.nx_; ++i)
		file<<f.x_(i)*AU_nm*1e-3<<' '<<nE_x(i)/AU_cm3<<'\n';
	file.close();

	file.open("wyniki/dane/nE_p.out", std::ios::out);
	file<<"# p [a.u.]  n(k) [cm^-1]\n";
	for (size_t i=0; i<f.nk_; ++i)
		file<<f.k_(i)<<' '<<nE_k(i)/AU_cm<<'\n';
	file.close();

	file.open("wyniki/dane/pot.out", std::ios::out);
	file<<"# x [um]  u(x) [eV]\n";
	for (size_t i=0; i<f.nx_; ++i)
		file<<f.x_(i)*AU_nm*1e-3<<' '<<f.u_(i)*AU_eV<<'\n';
	file.close();

	// std::ofstream pot_out;
	// pot_out.open("pot_02V.in", std::ios::out);
	// pot_out<<"# Input potential\n";
	// pot_out<<"# v_b = 0.2 V; L_D = 3 um, L_C = 0.5 um; T = 300 K; m*/m0 = 0.067; n_D = 2e18 cm^-3; N = 200\n";
	// pot_out<<"# x [a.u.]  u [a.u.]\n";
	// for (size_t i=0; i<f.nx_; ++i)
	//     pot_out<<f.x_(i)<<' '<<f.uStart_(i)<<'\n';
	// pot_out.close();

	cout<<"# l_YZ: "<<f.lYZ_<<" [a.u.]"<<endl;
	cout<<"# Final current = "<<f.calcCurr()*AU_Acm2<<" [Acm^-2]"<<endl;
	cout<<"# Norma: = "<<f.calcNorm()<<" [1]"<<endl;
	cout<<"# Int. BC: "<<calcInt(f.bc_, f.dk_)/AU_cm3<<" [a.u.]"<<endl;
	cout<<"# Fermi energy:\t"<<f.calcFermiEn(f.cD_, f.m_, f.temp_)*AU_eV<<'\t'
		<<f.calcFermiEn_MB(f.cD_, f.m_, f.temp_)*AU_eV<<endl;

	t_end = high_resolution_clock::now();
	t_elapsed =  duration_cast<duration<double>>(t_end - t_start);
	if (t_elapsed.count() < 60)
		cout<<"# RUN TIME: "<<t_elapsed.count()<<" s"<<endl;
	else
		cout<<"# RUN TIME: "<<t_elapsed.count()/60.<<" min"<<endl;

	return 0;

}


void calc_fd(WignerFunction& f) {

	double k = 0;
	for (size_t i=0; i<f.nk_; ++i) {
		k = f.k_(i);
		cout<<k<<'\t'
			<<1/(exp((k*k/2/f.m_-f.uF_)/(KB/AU_eV*f.temp_))+1)<<'\t'
			<<f.supplyFunction(k)<<'\t'
			<<f.maxwell_boltzmann(k)<<'\t'
			<<f.gaussian_bc(k)<<'\t'<<endl;
	}
	cout<<"# pF = "<<sqrt(f.uF_*f.m_*2.)<<endl;

}


void convolution(WignerFunction& f) {
	// size_t N1 = 100, N2 = 100;
	// cout<<"# Funkcja licząca splot"<<endl;
	// vec a(N1, fill::zeros), b(N2, fill::zeros), c(N1, fill::zeros);
	// for (size_t i=N1; i>0; --i) {
	// 	a(N1-i) = i/float(N1)*4;
	// }
	// for (size_t i=0; i<N2; ++i) b(i) = 1;
	// c = conv(a, b, "same");
	// for (size_t i=0; i<N1; ++i)
	// 	cout<<i/float(N1)*4<<'\t'<<a(i)<<'\t'<<c(i)<<endl;

	// double g = .5/(1e-12/AU_s);
	// double k = f.dk_;
	// double u = k*k/f.m_/2.;
	// auto fun = [g](double x) { return g/(x*x+g*g)/M_PI; };
	// size_t N = 1e2;
	// double h = (40/f.beta_+f.uF_)/float(N);
	// vec sf = vec(N, fill::zeros), p = vec(N, fill::zeros), r = vec(N, fill::zeros);
	// for (size_t i=0; i<N; ++i) {
	// 	p(i) = fun(i*h-u), sf(i) = sfx(f.uF_, i*h, f);
	// }
	// r = conv(sf, p);
	//
	// for (size_t i=0; i<N; ++i)
	// 	cout<<i*h<<'\t'<<sf(i)<<'\t'<<p(i)<<'\t'<<r(i)<<endl;

}


// Supply function (x)
double sfx(double mu, double x, WignerFunction& f){
	double m = f.m_;
	double c = m/M_PI*KB/AU_eV*f.temp_, ex = -(x-mu)/(KB/AU_eV*f.temp_);	 // [au]
	if (ex < 700)
		return c * log(exp(ex)+1);
	else{
		if (ex > 0)
			return c * ex;
		else return 0;
	}
}


void dissDecoh(WignerFunction& f, double t_start, double t_end, size_t n_step) {

	double t_step = (log10(t_end) - log10(t_start))/float(n_step-1);

	f.solveWignerEq();
	// f.solveWignerPoisson();
	double j0 = f.calcCurr(), n0 = f.calcNorm();
	double j=0, n=0, t=0;
	vec fit;

	cout.fill(' ');
	std::ofstream out, tpMap;
	out.open("wyniki/dane/diss_tau.out", std::ios::out);
	tpMap.open("wyniki/dane/tpMap.out", std::ios::out);
	out<<"## j0 = "<<j0*AU_Acm2<<"  n0 = "<<n0/AU_cm2<<endl;
	cout<<"## j0 = "<<j0*AU_Acm2<<"  n0 = "<<n0/AU_cm2<<endl;
	out<<"# i  tS [s]  J/J0  N/N0  R\n";
	cout<<"# i  tS [s]  J/J0  N/N0  R\n";
	tpMap<<"# tau [s]  p [a.u]  1/4  1/2  3/4\n";
	for (size_t i = 0; i < n_step; ++i) {
		t = t_start*pow(10,(n_step-1-i)*t_step);
		f.rR_ = 1/(t/AU_s);
		// f.rM_ = 1/(t/AU_s);
		f.solveWignerEq();
		// f.solveWignerPoisson();
		j = f.calcCurr(), n = f.calcNorm();
		// I-V char
		f.v_max_ = 0.001/AU_eV;
		f.v_min_ = 0.0/AU_eV;
		f.nv_ = 5;
		f.calc_IVchar();
		fit = polyfit(f.iv_v_, f.iv_i_, 1);
		//
		out<<i<<'\t'<<t<<'\t'<<j/j0<<'\t'<<1./fit(0)<<endl;
		cout<<i<<'\t'<<t<<'\t'<<j/j0<<'\t'<<1./fit(0)<<endl;
		for (size_t l=0; l<f.nk_; ++l) {
			tpMap<<t<<' '<<f.k_(l)
				<<'\t'<<f.f_(size_t(f.nx_/4.),l)
				<<'\t'<<f.f_(size_t(f.nx_/2.),l)
				<<'\t'<<f.f_(size_t(f.nx_*3/4.),l)<<'\n';
		}
		tpMap<<'\n';
	}
	out.close();
	tpMap.close();

	// cout<<"# int SF "<<calcInt(f.bc_, f.dk_)/2./M_PI/AU_cm3<<endl;
	// cout<<"# dx*dk "<<f.dx_*f.dk_<<" pi/2/dx "<<M_PI/2./f.dx_<<endl;

	// f.calc_IVchar();
	// cout<<"\n\n";

}


void dissDecoh_BC(WignerFunction& f, double t_start, double t_end, size_t n_step) {

	double t_step = (log10(t_end) - log10(t_start))/float(n_step-1);

	// f.solveWignerEq();
	f.solveWignerPoisson();
	double j0 = f.calcCurr(), n0 = f.calcNorm();
	double j=0, n=0, t=0;

	cout.fill(' ');
	std::ofstream out, tpMap;
	out.open("wyniki/dane/diss_tau.out", std::ios::out);
	tpMap.open("wyniki/dane/tpMap.out", std::ios::out);
	out<<"## j0 = "<<j0*AU_Acm2<<"  n0 = "<<n0/AU_cm2<<endl;
	cout<<"## j0 = "<<j0*AU_Acm2<<"  n0 = "<<n0/AU_cm2<<endl;
	out<<"# i  tS [s]  J/J0  N/N0  J [Acm^-2]  N [cm^-2]\n";
	cout<<"# i  tS [s]  J/J0  N/N0  J [Acm^-2]  N [cm^-2]\n";
	tpMap<<"# tau [s]  p [a.u]  1/4  1/2  3/4\n";
	for (size_t i = 0; i < n_step; ++i) {
		t = t_start*pow(10,(n_step-1-i)*t_step);
		// -1 -> 0, 0 -> SF, 1-3 -> splot, -2 -> MB, -3 -> Gauss
		f.bcType_ = 1;
		f.rG_ = 1./(t/AU_s);
		// f.solveWignerEq();
		f.solveWignerPoisson();
		j = f.calcCurr(), n = f.calcNorm();
		out<<i<<'\t'<<t<<'\t'<<j/j0<<'\t'<<n/n0<<'\t'<<j*AU_Acm2<<'\t'<<n/AU_cm2<<endl;
		cout<<i<<'\t'<<t<<'\t'<<j/j0<<'\t'<<n/n0<<'\t'<<j*AU_Acm2<<'\t'<<n/AU_cm2<<endl;
		for (size_t l=0; l<f.nk_; ++l) {
			tpMap<<t<<'\t'<<f.k_(l)
				<<'\t'<<f.f_(size_t(f.nx_/4.),l)
				<<'\t'<<f.f_(size_t(f.nx_/2.),l)
				<<'\t'<<f.f_(size_t(f.nx_*3/4.),l)<<'\n';
		}
		tpMap<<'\n';
	}
	out.close();
	tpMap.close();

	// cout<<"# int SF "<<calcInt(f.bc_, f.dk_)/2./M_PI/AU_cm3<<endl;
	// cout<<"# dx*dk "<<f.dx_*f.dk_<<" pi/2/dx "<<M_PI/2./f.dx_<<endl;

	// f.calc_IVchar();
	// cout<<"\n\n";

}


/*
void simGWP(WignerFunction& f, WignerFunction& f2, WignerFunction& f3) {
	// Wave packet simulation

	size_t nt=int(40);
	double t;
	double dp0, dx0;

	f.f_.zeros(), f.setWavePacket();
	f2.f_.zeros(), f2.setWavePacket();
	f3.f_.zeros(), f3.setWavePacket();

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
	f.setWavePacket();

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
