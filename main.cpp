#include "src/lib.hpp"
#include "src/WignerFunction.hpp"
#include "src/poisson1D.hpp"

#include <chrono>

using namespace wigner;
using namespace poisson;
using namespace std::chrono;

void simGWP(WignerFunction&);
void dissDecoh(WignerFunction&);
void simGWP(WignerFunction&, WignerFunction&, WignerFunction&);


int main(){

	high_resolution_clock::time_point t_start, t_end;
	duration<double> t_elapsed;

	t_start = high_resolution_clock::now();

	WignerFunction f;

	f.set_m(0.067);
	f.set_nx(200), f.set_nk(200);
	f.set_lD(60/AU_nm), f.set_lC(20/AU_nm);
	f.set_kmax(0.1);
	f.set_temp(300);
	f.set_cD(2e18*AU_cm3);
	f.set_useNLP(false);
	f.update();

	f.printParam();
	// // f.setGaussPot(0.1/AU_eV, 800, 50);
	// // f.setGaussPot(0.1/AU_eV, 1200, 50);
	// f.driftTermType_ = 0; // 1 - NLP, else - Classic
	// f.bcType_ = 0, f.rG_ = 0;  // 1./(1e15/AU_s);
	// f.gwp_x0_ = 150;
	// f.gwp_dx_ = 50;
	// f.gwp_p0_ = 0.04;
	// // f.gwp_dp_ = 0.01;  // a.u.
	// f.dt_ = 1e-15/AU_s;   // [au]
	// // f.rR_ = 1./(1e-13/AU_s);  // 1./(1e-12/AU_s);
	// // f.lambda_ = 0*AU_nm*AU_nm*AU_s;  // [nm^-2*s^-1]
	//
	// f.v_max_ = 0.4/AU_eV;
	// f.v_min_ = 0.0/AU_eV;
	// f.nv_ = 10;

	// FUNKCJA RÓWNOWAGOWA
	/*
	f.setLinPot(0);
	f.readPotential("potentials/pot_0V.in");
	f.u_ = f.u_ + f.uStart_;
	f.solveWignerEq();
	// f.solveWignerPoisson();
	f.fe_ = f.f_;
	*/

	cout<<"Setting up potential"<<endl;
	// f.rR_ = 1./(1e-12/AU_s);
	f.setLinPot(.2/AU_eV);
	// f.readPotential("potentials/pot_02V_tR_1e-12.in");
	// f.u_ = f.u_ + f.uStart_;
	// f.uStart_.zero();

	cout<<"Solving WTE"<<endl;
	// f.solveWignerEq();
	// f.solveWignerPoisson();

	f.v_max_ = 0.1/AU_eV;
	f.v_min_ = 0.0/AU_eV;
	f.nv_ = 10;
	f.calc_IVchar();

	// dissDecoh(f);

	cout<<"Calulating electron density"<<endl;
	array<double> nE_p =  f.calcCD_K();

	std::ofstream file;
	file.open("wyniki/nE_p.out", std::ios::out);
	file<<"# p  n(p)\n";
	for (size_t i=0; i<f.nk_; ++i)
		file<<f.k_(i)<<' '<<nE_p(i)/AU_cm2<<'\n';
	file.close();

	// std::ofstream pot_out;
	// pot_out.open("pot_02V.in", std::ios::out);
	// pot_out<<"# Input potential\n";
	// pot_out<<"# v_b = 0.2 V; L_D = 3 um, L_C = 0.5 um; T = 300 K; m*/m0 = 0.067; n_D = 2e18 cm^-3; N = 200\n";
	// pot_out<<"# x [a.u.]  u [a.u.]\n";
	// for (size_t i=0; i<f.nx_; ++i)
	//     pot_out<<f.x_(i)<<' '<<f.uStart_(i)<<'\n';
	// pot_out.close();

	cout<<"Saving wigner function"<<endl;
	f.saveWignerFun();

	cout<<"# Final current = "<<f.calcCurr()*AU_Acm2<<endl;
	cout<<"# Electron number = "<<f.calcNorm()/AU_cm2<<endl;
	cout<<"# Int. SF "<<calcInt(f.bc_, f.dk_)/2./M_PI/AU_cm3<<endl;

	t_end = high_resolution_clock::now();
	t_elapsed =  duration_cast<duration<double>>(t_end - t_start);
	cout<<"# RUN TIME: "<<t_elapsed.count()<<" s"<<endl;

	return 0;

}


void dissDecoh(WignerFunction& f) {

	size_t nt = 40;  // 1.69
	double tS = 1e-10;  // *pow(dt,nt-1);
	double dt = pow(1e-14/tS, 1./(nt-1));

	f.solveWignerEq();
	// f.solveWignerPoisson();
	double j0 = f.calcCurr(), n0 = f.calcNorm();
	double j=0, n=0;

	std::ofstream out, tpMap;
	out.open("wyniki/diss_tF.out", std::ios::out);
	tpMap.open("wyniki/tpMap.out", std::ios::out);
	out<<"## j0 = "<<j0*AU_Acm2<<"  n0 = "<<n0/AU_cm2<<endl;
	cout<<"## j0 = "<<j0*AU_Acm2<<"  n0 = "<<n0/AU_cm2<<endl;
	out<<"# i  tS [s]  J [Acm^-2]  N [cm^-2]\n";
	cout<<"# i  tS [s]  J [Acm^-2]  N [cm^-2]\n";
	tpMap<<"# tau [s]  p [a.u]  1/4  1/2  3/4\n";
	for (size_t i = 0; i < nt; ++i) {
		f.rF_ = 1./(tS/AU_s);
		f.solveWignerEq();
		// f.solveWignerPoisson();
		j = f.calcCurr(), n = f.calcNorm();
		out<<i<<' '<<tS<<' '<<j/j0<<' '<<n/n0<<' '<<j<<' '<<n<<endl;
		cout<<i<<' '<<tS<<' '<<j/j0<<' '<<n/n0<<' '<<j<<' '<<n<<endl;
		for (size_t l=0; l<f.nk_; ++l) {
			tpMap<<tS<<' '<<f.k_(j)
				<<' '<<f.f_(size_t(f.nx_/4.),l)
				<<' '<<f.f_(size_t(f.nx_/2.),l)
				<<' '<<f.f_(size_t(f.nx_*3/4.),l)<<'\n';
		}
		tpMap<<'\n';
		tS *= dt;  // tS *= 1.3
	}
	out.close();
	tpMap.close();

	// cout<<"# int SF "<<calcInt(f.bc_, f.dk_)/2./M_PI/AU_cm3<<endl;
	// cout<<"# dx*dk "<<f.dx_*f.dk_<<" pi/2/dx "<<M_PI/2./f.dx_<<endl;

	// f.calc_IVchar();
	// cout<<"\n\n";

}


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
