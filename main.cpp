#include "src/lib.hpp"
#include "src/WignerFunction.hpp"
#include "src/poisson1D.hpp"

#include <chrono>


using namespace wigner;
using namespace poisson;
using namespace std::chrono;


void simGWP(WignerFunction&);
void dissDecoh(WignerFunction&, double, double, size_t, char);
void simGWP(WignerFunction&, WignerFunction&, WignerFunction&);
void convolution(WignerFunction&);
double sfx(double, double, WignerFunction&);


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
	f.set_nx(200), f.set_nk(200);
	f.set_lD(3000/AU_nm), f.set_lC(500/AU_nm);
	f.set_kmax(0.1);
	f.set_temp(300);
	f.set_cD(2e18*AU_cm3);
	f.set_dt(5*1e-15/AU_s);
	f.set_time_dependent(true);
	f.update();

	f.printParam();

	//
	// Pakiet gaussowski
	//
	f.gwp_x0_ = 250/AU_nm;
	f.gwp_dx_ = 50/AU_nm;
	f.gwp_dp_ = 0.005; // 1./(2*f.gwp_dx_);  // a.u.
	f.gwp_p0_ = sqrt(2*f.m_*f.uF_);
	f.setWavePacket();
	f.gwp_x0_ = f.l_-250/AU_nm;
	f.gwp_p0_ = -sqrt(2*f.m_*f.uF_);
	f.setWavePacket();

	//
	// Warunek brzegowy
	// -1 -> 0, 0 -> SF, 1-3 -> splot
	//
	f.bcType_ = 0, f.rG_ = 0;  // 1./(1e-15/AU_s);
	// convolution(f);

	//
	// FUNKCJA RÓWNOWAGOWA
	//
	// f.setLinPot(0.);
	// f.readPotential("potentials/pot_0V.in");
	// f.u_ = f.u_ + f.uStart_;
	// f.solveWignerEq();
	// f.solveWignerPoisson();
	// f.fe_ = f.f_;

	//
	// Potencjał
	//
	cout<<"# Setting up potential"<<endl;
	f.setLinPot(0.0/AU_eV);
	// f.setGaussPot(0.1/AU_eV, 800, 50);
	// f.readPotential("potentials/pot_0V.in");
	// f.u_ = f.u_ + f.uStart_;
	// f.uStart_.zero();
	f.set_useNLP(false);

	//
	// Dyssypacja
	//
	// f.rR_ = 1./(1e-13/AU_s);  // 1./(1e-12/AU_s);
	// f.lambda_ = 0*AU_nm*AU_nm*AU_s;  // [nm^-2*s^-1]
	// f.rR_ = 1./(1e-11/AU_s);

	//
	// Rozwiazanie równania
	//
	// cout<<"# Solving WTE"<<endl;
	//
	// Boltzmann
	//
	// f.solveWignerEq();
	//
	// Boltzmann-Poisson
	//
	// f.solveWignerPoisson();
	//
	// Time evolution
	//
	double t_total = 100e-12/AU_s, t = 0, max_dj = 1e-5;
	array<double> dj(f.nx_, 0.), j0(f.nx_, 0.), j1(f.nx_, 0.);
	int n_min = 2, n_conv = 2, n_dj = 0, n_it = 0;
	double nrmse_j = 0;
	double j = 0, nc = 0, nc_prev = 0, dj_x = 0, j1_x = 0, j1_n = 0;
	bool conv_J = false;
	std::ofstream out;
	out.open("wyniki/dane/tev.out", std::ios::out);
	cout<<"# t [ps] \tCurr. [Acm^-2] \tNorm. [cm^-2] \tdN \tdJ [Acm^-2] \tmax(J) [Acm^-2] \t|max(J) - min(J)| [Acm^-2] \tdJ/max(J) \tRMSE \tCONV? \tn_dj"<<endl;
	out<<"# t [ps] \t Curr. [Acm^-2] \tNorm. [cm^-2] \tdN \tdJ [Acm^-2] \tmax(J) [Acm^-2] \t|max(J) - min(J)| [Acm^-2] \tdJ/max(J) \tRMSE \tCONV? \tn_dj"<<endl;
	while (t<t_total && n_dj < n_conv) {
		f.solveWignerEq();
		if (n_it > n_min) {
			j = f.calcCurr();
			nc_prev = nc, nc = f.calcNorm();
			j0 = j1, j1 = f.calcCurrArr();
			dj = j1 - j0;
			dj_x = dj.max(), j1_x = j1.max(), j1_n = j1.min();
			nrmse_j = sqrt(((j1-j0)*(j1-j0)).sum())/fabs(j1_x-j1_n);
			conv_J = fabs(dj_x/j1_x) > max_dj ? false : true;
			n_dj = conv_J ? n_dj+1 : 0;
			out<<t*AU_s*1e12
				<<'\t'<<j*AU_Acm2
				<<'\t'<<nc/AU_cm2
				<<'\t'<<(nc-nc_prev)/nc
				<<'\t'<<dj_x*AU_Acm2
				<<'\t'<<j1_x*AU_Acm2
				<<'\t'<<fabs(j1_x-j1_n)*AU_Acm2
				<<'\t'<<dj_x/j1_x
				<<'\t'<<nrmse_j
				<<'\t'<<conv_J
				<<'\t'<<n_dj<<endl;
			cout<<t*AU_s*1e12
				<<'\t'<<j*AU_Acm2
				<<'\t'<<nc/AU_cm2
				<<'\t'<<(nc-nc_prev)/nc
				<<'\t'<<dj_x*AU_Acm2
				<<'\t'<<j1_x*AU_Acm2
				<<'\t'<<fabs(j1_x-j1_n)*AU_Acm2
				<<'\t'<<dj_x/j1_x
				<<'\t'<<nrmse_j
				<<'\t'<<conv_J
				<<'\t'<<n_dj<<endl;
		}
		t += f.dt_, n_it++;
	}
	out.close();

	//
	// Charakterystyka J-V
	//
	// f.v_max_ = 0.2/AU_eV;
	// f.v_min_ = 0.0/AU_eV;
	// f.nv_ = 20;
	// f.calc_IVchar();

	//
	// Zalezność gęstości prądu od dyssypacji
	// diss_type = 'R', 'M', 'F'
	//
	// dissDecoh(f, 1e-14, 1e-10, 20, 'R');

	cout<<"# Calulating electron density"<<endl;
	array<double> nE_x =  f.calcCD_X();
	array<double> nE_k =  f.calcCD_K();

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

	// std::ofstream pot_out;
	// pot_out.open("pot_02V.in", std::ios::out);
	// pot_out<<"# Input potential\n";
	// pot_out<<"# v_b = 0.2 V; L_D = 3 um, L_C = 0.5 um; T = 300 K; m*/m0 = 0.067; n_D = 2e18 cm^-3; N = 200\n";
	// pot_out<<"# x [a.u.]  u [a.u.]\n";
	// for (size_t i=0; i<f.nx_; ++i)
	//     pot_out<<f.x_(i)<<' '<<f.uStart_(i)<<'\n';
	// pot_out.close();

	cout<<"# Saving wigner function"<<endl;
	f.saveWignerFun();

	cout<<"# Final current = "<<f.calcCurr()*AU_Acm2<<" [Acm^-2]"<<endl;
	cout<<"# Electron number = "<<f.calcNorm()/AU_cm2<<" [cm^-2]"<<endl;
	cout<<"# Int. SF "<<calcInt(f.bc_, f.dk_)/2./M_PI/AU_cm3<<" [cm^-3]"<<endl;

	t_end = high_resolution_clock::now();
	t_elapsed =  duration_cast<duration<double>>(t_end - t_start);
	if (t_elapsed.count() < 60)
		cout<<"# RUN TIME: "<<t_elapsed.count()<<" s"<<endl;
	else
		cout<<"# RUN TIME: "<<t_elapsed.count()/60.<<" min"<<endl;

	return 0;

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


void dissDecoh(WignerFunction& f, double t_start, double t_end, size_t n_step, char diss_type) {

	double t_step = (log10(t_end) - log10(t_start))/float(n_step-1);

	f.solveWignerEq();
	// f.solveWignerPoisson();
	double j0 = f.calcCurr(), n0 = f.calcNorm();
	double j=0, n=0, t=0;

	cout.fill(' ');
	std::ofstream out, tpMap;
	out.open("wyniki/dane/diss_tR.out", std::ios::out);
	tpMap.open("wyniki/dane/tpMap.out", std::ios::out);
	out<<"## j0 = "<<j0*AU_Acm2<<"  n0 = "<<n0/AU_cm2<<endl;
	cout<<"## j0 = "<<j0*AU_Acm2<<"  n0 = "<<n0/AU_cm2<<endl;
	out<<"# i  tS [s]  J/J0  N/N0  J [Acm^-2]  N [cm^-2]\n";
	cout<<"# i  tS [s]  J/J0  N/N0  J [Acm^-2]  N [cm^-2]\n";
	tpMap<<"# tau [s]  p [a.u]  1/4  1/2  3/4\n";
	for (size_t i = 0; i < n_step; ++i) {
		t = t_start*pow(10,i*t_step);
		switch(diss_type) {
			case 'R' : f.rR_ = 1/(t/AU_s); break;
			case 'M' : f.rM_ = 1/(t/AU_s); break;
			case 'F' : f.rF_ = 1/(t/AU_s); break;
			default : f.rR_ = 1/(t/AU_s);
		}

		f.solveWignerEq();
		// f.solveWignerPoisson();
		j = f.calcCurr(), n = f.calcNorm();
		out<<i<<' '<<t<<' '<<j/j0<<' '<<n/n0<<' '<<j*AU_Acm2<<' '<<n/AU_cm2<<endl;
		cout<<i<<' '<<t<<' '<<j/j0<<' '<<n/n0<<' '<<j*AU_Acm2<<' '<<n/AU_cm2<<endl;
		for (size_t l=0; l<f.nk_; ++l) {
			tpMap<<t<<' '<<f.k_(j)
				<<' '<<f.f_(size_t(f.nx_/4.),l)
				<<' '<<f.f_(size_t(f.nx_/2.),l)
				<<' '<<f.f_(size_t(f.nx_*3/4.),l)<<'\n';
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

	f.f_.zero(), f.setWavePacket();
	f2.f_.zero(), f2.setWavePacket();
	f3.f_.zero(), f3.setWavePacket();

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

	f.f_.zero();
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
