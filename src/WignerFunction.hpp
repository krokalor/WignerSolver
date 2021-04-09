#ifndef WIGNERFUNCTION_HPP
#define WIGNERFUNCTION_HPP

#include "lib.hpp"

std::map<std::string, double> readParam(std::string);

namespace wigner
{

// #################### klasa WignerFunction ####################
class WignerFunction{

public:

	/*
	int calc_mode_;
	double v_max_;
	int nv_;
	void setSysContacts(double irG, char idist)
			{ rG_ = irG*t0_, Gamma_ = rG_*.5, bcType_ = idist; }

	// #################### Set system dissipation ####################
	void setSysDissipation(double irR, double irM, double ilambda)
			{ rR_ = irR*t0_, rM_ = irM*t0_, lambda_ = ilambda*a0_*a0_*t0_; }
	*/

	// Default constructor
	WignerFunction() :
		m_ (0.067),
		lD_ (60./AU_nm),
		lC_ (20./AU_nm),
		l_ (lD_ + 2*lC_),
		lYZ_ (1),
		nx_ (100),
		nk_ (100),
		nk2_ (size_t(nk_/2.)),
		nxk_ (nx_*nk_),
		dx_ (l_/float(nx_-1)),
		kmax_ (M_PI/2./dx_),
		dk_ (2.*kmax_/float(nk_-1)),
		dt_ (1e-15/AU_s),
		courant_num_ (dt_*kmax_/m_/dx_),  // Courant-Friedricks-Lewy
		uF_ (0.087/AU_eV), uR_ (0.087/AU_eV), uL_ (0.087/AU_eV),
		set_uF_ (false),
		cD_(2e18*AU_cm3),
		temp_ (300), beta_ (1./KB/temp_*AU_eV),
		epsilonR_(13.1),
		uBias_(0),
		part_num_(1e3),
		gwp_tl_ (sqrt(2*M_PI*beta_/m_)),
		gwp_dx_ (0.28 * gwp_tl_),
		gwp_dp_ (sqrt(2.*m_/beta_)),
		gwp_x0_ (3*gwp_dx_),
		gwp_p0_ (sqrt(2*m_*uL_)),
		gwp_A_(part_num_ / (2.*M_PI*gwp_dx_*gwp_dp_)),
		rR_(0),
		rM_(0),
		rG_(0), Gamma_ (rG_*.5), rF_(0),
		lambda_(0),
		bcType_(0),
		useNLP_(false),
		f_(mat(nx_, nk_)),
		fe_(mat(nx_, nk_)),
		f0_(mat(nx_, nk_)),
		fL_(mat(nx_, nk_)),
		fR_(mat(nx_, nk_)),
		u_(vec(nx_, fill::zeros)),
		uC_(vec(nx_, fill::zeros)),
		uB_(vec(nx_, fill::zeros)),
		du_(vec(nx_, fill::zeros)),
		uStart_(vec(nx_, fill::zeros)),
		bc_(vec(nk_, fill::zeros)),
		x_(vec(nx_, fill::zeros)),
		k_(vec(nk_, fill::zeros)),
		sin_(mat(nk_,nk_*nk2_)),
		cdX_(vec(nx_, fill::zeros)),
		cdK_(vec(nk_, fill::zeros)),
		currD_(vec(nx_, fill::zeros))
		// cR_(iv["dconc_right"]),
		{

		// ########## Configuration space array values ##########
		for (size_t i=0; i<nx_; ++i) x_(i) = i*dx_;

		// ########## Wave vector space array values ##########
		for (size_t j=0; j<nk_; ++j) k_(j) = dk_*(j-(nk_-1)*.5);

		// ########## NLP sinus values ##########
		// #pragma omp parallel for collapse(3)
		for (size_t j=0; j<nk_; ++j)
				for (size_t g=0; g<nk_; g++)
						for (size_t h=0; h<nk2_; h++)
								sin_(j,g*nk2_+h) = sin(2*M_PI/nk_*h*(j-g));

	}  // End of constructor

	~WignerFunction(){}

	double get_dk() { return dk_; }
	double get_dx() { return dx_; }
	double get_nk() { return nk_; }
	double get_nx() { return nx_; }
	double get_l() { return l_; }
	double get_m_eff() { return m_; }
	double get_uF() { return uF_; }
	double get_m() { return m_; }
	double get_dt(){ return dt_; }
	double get_rR() { return rR_; }
	double get_rM() { return rM_; }
	double get_rG() { return rG_; }
	double get_lambda() { return lambda_; }

	void set_m(double m) { m_ = m; }
	void set_nx(size_t nx) { nx_ = nx; }
	void set_nk(size_t nk) { nk_ = nk; }
	void set_lD(double lD) { lD_ = lD; }
	void set_lC(double lC) { lC_ = lC; }
	void set_lYZ(double lYZ) { lYZ_ = lYZ; }
	void set_part_num(double part_num) { part_num_ = part_num; }
	void set_kmax(double kmax) { kmax_ = kmax; }
	void set_dt(double dt) { dt_ = dt; }
	void set_cD(double cD) { cD_ = cD; }
	void set_uF(double uF) { uF_ = uF, set_uF_ = true; }
	void set_temp(double temp) { temp_ = temp; }
	void set_epsilonR(double epsilonR) { epsilonR_ = epsilonR; }
	void set_uBias(double uBias) { uBias_ = uBias; }
	void set_useNLP(bool useNLP) { useNLP_ = useNLP; }

	void update() {
		nxk_ = nx_*nk_, nk2_ = size_t(nk_/2.);
		l_ = lD_ + 2*lC_;
		dx_ = l_/float(nx_-1);
		kmax_ < 0 ? M_PI/2./dx_ : kmax_;
		dk_ = 2.*kmax_/float(nk_-1);
		courant_num_ = dt_*kmax_/m_/dx_;
		beta_ = 1./KB/temp_*AU_eV;
		uF_ = set_uF_ ? uF_ : calcFermiEn(cD_, m_, temp_);  // Fermi energy
		f_ = mat(nx_, nk_, fill::zeros);
		fe_ = mat(nx_, nk_, fill::zeros);
		f0_ = mat(nx_, nk_, fill::zeros);
		fL_ = mat(nx_, nk_, fill::zeros);
		fR_ = mat(nx_, nk_, fill::zeros);
		u_ = vec(nx_, fill::zeros);
		uC_ = vec(nx_, fill::zeros);
		uB_ = vec(nx_, fill::zeros);
		du_ = vec(nx_, fill::zeros);
		uStart_ = vec(nx_, fill::zeros);
		bc_ = vec(nk_, fill::zeros);
		x_ = vec(nx_, fill::zeros);
		k_ = vec(nk_, fill::zeros);
		sin_ = mat(nk_,nk_*nk2_);
		cdX_ = vec(nx_, fill::zeros), cdK_ = vec(nk_, fill::zeros);
		currD_ = vec(nx_, fill::zeros);
		for (size_t i=0; i<nx_; ++i) x_(i) = i*dx_;
		for (size_t j=0; j<nk_; ++j) k_(j) = dk_*(j-(nk_-1)*.5);
		// #pragma omp parallel for collapse(3)
		for (size_t j=0; j<nk_; ++j)
				for (size_t g=0; g<nk_; g++)
						for (size_t h=0; h<nk2_; h++)
								sin_(j,g*nk2_+h) = sin(2*M_PI/nk_*h*(j-g));
	}

	double calcNorm();
	double calcEX();
	double calcEK();
	double calcSDK();
	double calcSDX();
	vec calcCD_X();
	vec calcCD_K();

	void readPotential(std::string);                       // Read potential from file pot.in
	void printParam();

	void solveWignerEq();
	void solveWignerEq_A();
	void solveTimeEv();
	double calcCurr();					// Current density
	vec calcCurrArr();				// Current density array
	void saveWignerFun();
	void clearWignerFun();

	void initEq();

	void setBoundCond();						// Boundary conditions
	void setEquilibriumFunction(std::string, bool);				// Calculating equilibrium function

	void diffusionTerm(size_t, size_t, double);       // Filling matrice with drift term
	void driftTerm(size_t, size_t, double);       // Filling matrice with drift term
	void scatteringTerm(size_t, size_t, double);       // Filling matrice with drift term

	void solveRec();							// Wigner equation solved recursively
	void solveMatrixEq();						// Wigner equation solved by solving matrix equation

	void solveWignerPoisson();
	void calc_IVchar(double, double, size_t);
	void calcMobility();

	double supplyFunction(double);		// Supply function
	double sf_x(double, double);					// Supply function as function of x

	double maxwell_boltzmann(double);
	double gaussian_bc(double);
	double gaussian_x(double, double);

	double lorentz(double);	 // Lorentz profile
	double gauss(double);	 // Gauss
	double voigt(double);	 // Voigt profile

	// wignerTools.cpp
	void setPotBias(double);
	void addGaussPot(double, double, double);
	void setRTD(double, double, double, double, double);
	void setSysContacts(double, char);
	void setSysDissipation(double, double, double);
	void addWavePacket();
	double WavePacket_TE(double, double);
	double nC(double, double);
	double fermiInt(double, double);
	double calcFermiEn(double, double, double);
	double calcFermiEn_MB(double, double, double);

	// PUBLIC VARIABLES

	double m_;									// effective mass [1]
	double lD_, lC_, l_;					// total lenght, devie len., contacts len.
	double lYZ_;
	size_t nx_, nk_, nk2_, nxk_;   				// x/k-space nr of steps
	double dx_;									// x-space step size (lattice constant)
	double kmax_;  // k-space range
	double dk_;									// k-space step size (Brillouin zone / Nk)
	double dt_;                                 // time step size
	double courant_num_;                        // Courant number
	double uF_, uR_, uL_;						// Fermi energy in left/right contact
	double set_uF_;
	double cD_;  // cL_, cR_;
	double temp_;								// Contacts temperature
	double beta_;								// beta = 1/kb/T
	double epsilonR_;

	double uBias_;

	double part_num_, gwp_tl_;
	double gwp_dx_, gwp_dp_;  			// Wave packet parameters
	double gwp_x0_, gwp_p0_;
	double gwp_A_;
	// double ex_, ek_;                            // Expected values

	double rR_, rM_, rG_, Gamma_, rF_;				// Scattering rate (1/tau)
	double lambda_;								// Localization rate
	int bcType_;									// Distribution used as bc
	bool useNLP_;

	// private:

	mat f_;  // Wigner function
	mat fe_;  // Equilibrium Wigner function
	mat f0_;  // Wigner function before time evolution
	mat fL_;  // Wigner function for el. from LEFT contact
	mat fR_;  // Wigner function for el. from RIGHT contact
	vec u_;  // Potential energy
	vec uC_; // Hartree potential / bias potential
	vec uB_;  // Band offset
	vec du_;  // Potential derivative
	vec uStart_;  // Potential
	vec bc_;  // Boundary condition
	vec x_;  // Position values
	vec k_;  // Wave vector values
	mat sin_;  // Sine function values
	vec cdX_;  // Carrier density in x / k
	vec cdK_;
	vec currD_;  // Current density

	sp_mat a_;
	vec b_;

	vec iv_i_, iv_v_;
};


}


#endif


/*

*/
