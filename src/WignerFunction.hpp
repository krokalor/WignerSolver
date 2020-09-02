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

	WignerFunction(std::map<std::string, double> iv) :
    m_ (iv["effective_mass"]),
    l_ (0),
    lD_ (iv["device_lenght"]),
    lC_ (iv["contact_lenght"]),
    kmax_ (iv["max_k"]),
    nx_ ((unsigned int)iv["xspace_step_nr"]),
    nk_ ((unsigned int)iv["kspace_step_nr"]),
    nk2_ ((unsigned int)(nk_/2.)),
    nxk_ (nx_*nk_),
    dx_ (0), dk_ (0), dt_ (0),  // To po to żeby uniknąć warningów
    courant_num_ (iv["courant_num"]),
    uF_ (iv["fermi_level"]), uR_ (0), uL_ (0),
    cD_(iv["dconc"]),
    temp_ (iv["contact_temp"]), beta_ (0),
		epsilonR_(iv["rel_diel_per"]),
    uB_(iv["v_bias"]),
    gwp_x0_ (iv["gwp_x0"]),
    gwp_dx_ (iv["gwp_dx"]),
    gwp_p0_ (iv["gwp_p0"]),
    gwp_dp_ (iv["gwp_dp"]),
    gwp_tl_ (0), gwp_A_(0),
    v_max_ (iv["max_voltage"]),
    v_min_ (iv["min_voltage"]),
    nv_ ((unsigned int)iv["voltage_step_nr"]),
    part_num_((unsigned long)iv["part_num"]),
    rR_(iv["inelastic_sc"]),
    rM_(iv["elastic_sc"]),
    rG_(iv["contact_diss"]), Gamma_ (rG_*.5), rF_(0),
    lambda_(iv["spatial_decoh"]),
    bcType_(iv["contact_dist"]),
		driftTermType_(0)
    // cR_(iv["dconc_right"]),
    {

    // ########## nm, eV, s -> a.u. ##########
    // uF_ = uF_/ha_;
    // kmax_ = kmax_*a0_;
    // lC_ = lC_/a0_;
    // lD_ = lD_/a0_;
    l_ = lD_ + 2*lC_;
    // gwp_x0_ = gwp_x0_/a0_;
    // gwp_dx_ = gwp_dx_/a0_;
    // rR_ = rR_*t0_;
    // rM_ = rM_*t0_;
    // rG_ = rG_*t0_, Gamma_ = rG_*.5;
    // lambda_ = lambda_*a0_*a0_*t0_;

    // ########## setting up parameters ##########
    dx_ = l_/(float)(nx_-1);
 		kmax_ = kmax_ < 0 ? M_PI/2./dx_ : kmax_;
    dk_ = 2.*kmax_/(float)(nk_-1);
    nk2_ = (unsigned int)(nk_/2.);
    nxk_ = (unsigned int)(nx_*nk_);
    beta_ = 1./KB/temp_*AU_eV;

		if (gwp_dx_ < 0)
		{
			gwp_tl_ = sqrt(2*M_PI*beta_/m_);
			gwp_dx_ = 0.28 * gwp_tl_;  // 0.28
		}
		if (gwp_dp_ < 0)
			gwp_dp_ = sqrt(2.*m_/beta_);
		if (gwp_x0_ < 0)
			gwp_x0_ = 3*gwp_dx_;
		if (gwp_p0_ < 0)
			gwp_p0_ = sqrt(2*m_*uL_);

		gwp_A_ = part_num_ / (2.*M_PI*gwp_dx_*gwp_dp_);

		driftTermType_ = 0;

    // ########## Courant-Friedricks-Lewy ##########
    dt_ = dx_/kmax_*m_*courant_num_;

    // ########## Arrays initialization ##########
    bc_ = array<double>(nk_);
    fe_ = matrix<double>(nx_, nk_);		// Equilibrium WF depends only on k
    f_ = matrix<double>(nx_, nk_);
    f0_ = f_,  fL_ = f_,  fR_ = f_;
    u_ = array<double>(nx_), du_ = u_, uStart_ = u_;

    // ########## Wave vector array values ##########
    k_ = array<double>(nk_);
    for (size_t j=0; j<nk_; ++j) k_(j) = dk_*(j-(nk_-1)*.5);

    // ########## Wave vector array values ##########
    x_ = array<double>(nx_);
    for (size_t i=0; i<nx_; ++i) x_(i) = i*dx_;

    // ########## NLP sinus values ##########
    sin_ = matrix<double>(nk_,nk_*nk2_);
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

	void set_kmax(double kmax_i) {
		kmax_ = kmax_i;
		dk_ = 2.*kmax_/(float)(nk_-1);
		for (size_t j=0; j<nk_; ++j) k_(j) = dk_*(j-(nk_-1)*.5);
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
	double l_, lD_, lC_, kmax_;					// total lenght, devie len., contacts len., k-space range
	size_t nx_, nk_, nk2_, nxk_;   				// x/k-space nr of steps
	double dx_;									// x-space step size (lattice constant)
	double dk_;									// k-space step size (Brillouin zone / Nk)
	double dt_;                                 // time step size
  double courant_num_;                        // Courant number
	double uF_, uR_, uL_;						// potential in left/right contact
  double cD_;  // cL_, cR_;
	double temp_;								// Contacts temperature
	double beta_;								// beta = 1/kb/T
	double epsilonR_;

  double uB_;

	double gwp_x0_, gwp_dx_;
	double gwp_p0_, gwp_dp_, gwp_tl_;  			// Wave packet parameters
	double gwp_A_;
  // double ex_, ek_;                            // Expected values

  double v_max_, v_min_;
	double nv_;
	double part_num_;

	double rR_, rM_, rG_, Gamma_, rF_;				// Scattering rate (1/tau)
	double lambda_;								// Localization rate
	int bcType_;									// Distribution used as bc
  int driftTermType_;

  matrix<double> f_;		                    // Wigner function

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

  array<double> bc_;					        // Boundary condition
	matrix<double> fe_;					        // Equilibrium Wigner function
	matrix<double> f0_;		                    // Wigner function before time evolution
	matrix<double> fL_;		                    // Wigner function for el. from LEFT contact
	matrix<double> fR_;		                    // Wigner function for el. from RIGHT contact
	matrix<double> sin_;	                    // Sine function values
	array<double> k_;	                        // Wave vector values
	array<double> x_;		                    // Position values
	array<double> u_;						    // Potential
	array<double> du_;					        // Potential derivative
	array<double> uStart_;						    // Potential


	sp_mat a_;
  vec b_;
};


}


#endif


/*

*/
