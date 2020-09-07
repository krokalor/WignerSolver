#ifndef WIGNERFUNCTION_HPP
#define WIGNERFUNCTION_HPP

#include "lib.hpp"
using namespace arma;

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
    cD_(2e18*AU_cm3),
    temp_ (300), beta_ (1./KB/temp_*AU_eV),
		epsilonR_(13.1),
    uB_(0),
    part_num_(1e3),
		gwp_tl_ (sqrt(2*M_PI*beta_/m_)),
    gwp_dx_ (0.28 * gwp_tl_),
    gwp_dp_ (sqrt(2.*m_/beta_)),
    gwp_x0_ (3*gwp_dx_),
    gwp_p0_ (sqrt(2*m_*uL_)),
    gwp_A_(part_num_ / (2.*M_PI*gwp_dx_*gwp_dp_)),
    v_max_ (0),
    v_min_ (0.1/AU_eV),
    nv_ (10),
    rR_(0),
    rM_(0),
    rG_(0), Gamma_ (rG_*.5), rF_(0),
    lambda_(0),
    bcType_(0),
		useNLP_(false),
		f_(matrix<double>(nx_, nk_)),
		fe_(f_),
    f0_(f_),  fL_(f_),  fR_(f_),
    u_(array<double>(nx_)), du_(u_), uStart_(u_),
		bc_(array<double>(nk_)),
		x_(array<double>(nx_)),
		k_(array<double>(nk_)),
		sin_(matrix<double>(nk_,nk_*nk2_))
    // cR_(iv["dconc_right"]),
    {

    // ########## Configuration space array values ##########
    for (size_t i=0; i<nx_; ++i) x_(i) = i*dx_;

    // ########## Wave vector space array values ##########
    for (size_t j=0; j<nk_; ++j) k_(j) = dk_*(j-(nk_-1)*.5);

    // ########## NLP sinus values ##########
		// num_threads(3)
    #pragma omp parallel for collapse(3) num_threads(N_THREADS)
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
  double get_v_max() { return v_max_; }
  double get_nv() { return nv_; }
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
	void set_kmax(double kmax) { kmax_ = kmax; }
	void set_dt(double dt) { dt_ = dt; }
	void set_cD(double cD) { cD_ = cD; }
	void set_temp(double temp) { temp_ = temp; }
	void set_epsilonR(double epsilonR) { epsilonR_ = epsilonR; }
	void set_uBias(double uB) { uB_ = uB; }
	void set_useNLP(bool useNLP) { useNLP_ = useNLP; }

	void update() {
		nxk_ = nx_*nk_, nk2_ = size_t(nk_/2.);
		l_ = lD_ + 2*lC_;
		dx_ = l_/float(nx_-1);
		kmax_ < 0 ? M_PI/2./dx_ : kmax_;
		dk_ = 2.*kmax_/float(nk_-1);
		courant_num_ = dt_*kmax_/m_/dx_;
		beta_ = 1./KB/temp_*AU_eV;
		f_ = matrix<double>(nx_, nk_);
		fe_ = f_, f0_ = f_,  fL_ = f_,  fR_ = f_;
		u_ = array<double>(nx_), du_ = u_, uStart_ = u_;
		bc_ = array<double>(nk_);
		x_ = array<double>(nx_);
		k_ = array<double>(nk_);
		sin_ = matrix<double>(nk_,nk_*nk2_);
    for (size_t i=0; i<nx_; ++i) x_(i) = i*dx_;
		for (size_t j=0; j<nk_; ++j) k_(j) = dk_*(j-(nk_-1)*.5);
    #pragma omp parallel for collapse(3) num_threads(N_THREADS)
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
  array<double> calcCD_X();
  array<double> calcCD_K();

  void readPotential(std::string);                       // Read potential from file pot.in
  void printParam();

	void solveWignerEq();
	void solveWignerEq_A();
	void solveTimeEv();
	double calcCurr();					// Current density
	array<double> calcCurrArr();				// Current density array
	void saveWignerFun();
	void clearWignerFun();

  void initEq();

  void setBoundCond();						// Boundary conditions
  void setEquilibriumFunction();				// Calculating equilibrium function

  void diffusionTerm(size_t, size_t, double);       // Filling matrice with drift term
  void driftTerm(size_t, size_t, double);       // Filling matrice with drift term
  void scatteringTerm(size_t, size_t, double);       // Filling matrice with drift term

	void solveRec();							// Wigner equation solved recursively
	void solveMatrixEq();						// Wigner equation solved by solving matrix equation

	void solveWignerPoisson();
	void calc_IVchar();
	void calcMobility();

	double supplyFunction(double);		// Supply function
	double sf_x(double, double);					// Supply function as function of x

	double lorentz(double);	 // Lorentz profile
	double gauss(double);	 // Gauss
	double voigt(double);	 // Voigt profile

	// wignerTools.cpp
	void setLinPot(double);
	void setGaussPot(double, double, double);
	void setRTD(double, double, double, double, double);
	void setSysContacts(double, char);
	void setSysDissipation(double, double, double);
	void setWavePacket();
	double WavePacket_TE(double, double);
	double nC(double, double);
	double fermiInt(double, double);
	double calcFermiEn(double, double, double);

  // PUBLIC VARIABLES

	double m_;									// effective mass [1]
	double lD_, lC_, l_;					// total lenght, devie len., contacts len.
	size_t nx_, nk_, nk2_, nxk_;   				// x/k-space nr of steps
	double dx_;									// x-space step size (lattice constant)
	double kmax_;  // k-space range
	double dk_;									// k-space step size (Brillouin zone / Nk)
	double dt_;                                 // time step size
  double courant_num_;                        // Courant number
	double uF_, uR_, uL_;						// potential in left/right contact
  double cD_;  // cL_, cR_;
	double temp_;								// Contacts temperature
	double beta_;								// beta = 1/kb/T
	double epsilonR_;

  double uB_;

	double part_num_, gwp_tl_;
	double gwp_dx_, gwp_dp_;  			// Wave packet parameters
	double gwp_x0_, gwp_p0_;
	double gwp_A_;
  // double ex_, ek_;                            // Expected values

  double v_max_, v_min_;
	double nv_;

	double rR_, rM_, rG_, Gamma_, rF_;				// Scattering rate (1/tau)
	double lambda_;								// Localization rate
	int bcType_;									// Distribution used as bc
  bool useNLP_;

	// private:

  /*
	std::vector<double> bc_;					// Boundary condition
	std::vector<double> fe_;					// Equilibrium Wigner function
	std::vector< std::vector<double> > f_;		// Wigner function
	std::vector< std::vector<double> > f0_;		// Wigner function before time evolution
	std::vector< std::vector<double> > fL_;		// Wigner function for el. from LEFT contact
	std::vector< std::vector<double> > fR_;		// Wigner function for el. from RIGHT contact
	std::vector< std::vector<double> > sin_;	// Sine function values
	std::vector<double> k_;			            // Wave vector values
	std::vector<double> x_;			            // Position values
	std::vector<double> u_;						// Potential
	std::vector<double> du_;					// Potential derivative
  */

  matrix<double> f_;		                    // Wigner function
	matrix<double> fe_;					        // Equilibrium Wigner function
	matrix<double> f0_;		                    // Wigner function before time evolution
	matrix<double> fL_;		                    // Wigner function for el. from LEFT contact
	matrix<double> fR_;		                    // Wigner function for el. from RIGHT contact
	array<double> u_;						    // Potential
	array<double> du_;					        // Potential derivative
	array<double> uStart_;						    // Potential
  array<double> bc_;					        // Boundary condition
	array<double> x_;		                    // Position values
	array<double> k_;	                        // Wave vector values
	matrix<double> sin_;	                    // Sine function values
	sp_mat a_;
  vec b_;
};


}


#endif


/*

*/
