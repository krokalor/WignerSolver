#ifndef WIGNERFUNCTION_HPP
#define WIGNERFUNCTION_HPP

#include "lib.hpp"

namespace wigner
{

// #################### klasa WignerFunction ####################
class WignerFunction{

	size_t nx_;
	double lD_, lC_;  // total lenght, devie len., contacts len.
	double l_;
	double dx_;  // x-space step size (lattice constant)
	double nk_;
	double kmax_;  // k-space range
	double dk_;  // k-space step size (Brillouin zone / Nk)
	size_t nk2_, nxk_;  // x/k-space nr of steps

	double dt_ = 1e-15/AU_s;  // time step size (1 fs)
	double courant_num_;  // Courant number

	double m_ = 0.067;
	double temp_ = 77;  // Contacts temperature [K]
	double beta_ = 1./KB/temp_*AU_eV;  // beta = 1/kb/T
	double uF_ = 0.087/AU_eV;
	double uR_ = uF_, uL_ = uF_;  // Fermi energy in left/right contact
	double cD_ = 2e18/AU_cm3;  // cL_, cR_;
	double epsilonR_ = 13.1;  // relative permitivitty (for GaAs)

	double uBias_ = 0;  // bias voltage [eV]
	double rR_ = 0, rM_ = 0, rG_ = 0, Gamma_ = 0, rF_ = 0;  // Scattering rate (1/tau)
	double lambda_ = 0;  // Localization rate
	int bcType_ = 1;  // boundary condition type

	// device cross section (1um x 1um)
	double lYZ_ = 1e-6*1e-6 * 1e18/AU_nm/AU_nm;
	double part_num_;  // particle number
	double gwp_dx_, gwp_dp_; // Wave packet parameters
	double gwp_x0_, gwp_p0_;
	double gwp_A_;

	bool useNLP_, useQC_, uBias_BC_;  // bool variables

	mat f_;  // Wigner function
	mat fe_;  // Equilibrium Wigner function
	mat f0_;  // Wigner function before time evolution
	mat fL_;  // Wigner function for el. from LEFT contact
	mat fR_;  // Wigner function for el. from RIGHT contact
	vec u_;  // Potential energy
	vec uC_; // Hartree potential / bias potential
	vec uB_;  // Band offset
	vec du_;  // Potential derivative
	vec dddu_;  // Potential third derivative
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

	vec iv_i_, iv_v_, iv_iRange_, iv_n_;

public:

	// Default constructor
	WignerFunction() :
		nx_ (100),
		lD_ (60./AU_nm),
		lC_ (20./AU_nm),
		l_ (lD_ + 2*lC_),
		dx_ (l_/float(nx_)),
		nk_ (100),
		kmax_ (M_PI/2./dx_),
		dk_ (2.*kmax_/float(nk_)),
		nk2_ (size_t(nk_/2.)),
		nxk_ (nx_*nk_),
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
		currD_(vec(nx_, fill::zeros)),
		a_(sp_mat(nxk_, nxk_)),
		b_(vec(nxk_, fill::zeros))
		{
		cout<<"## Start: WignerFunction default constructor"<<endl;
		// ########## Configuration space array values ##########
		cout<<"# Setting up configuration space array values"<<endl;
		for (size_t i=0; i<nx_; ++i) x_(i) = i*dx_;
		// ########## Wave vector space array values ##########
		cout<<"# Setting up wave vector space array values"<<endl;
		for (size_t j=0; j<nk_; ++j) k_(j) = dk_*(j-(nk_-1)*.5);
		// ########## NLP sinus values ##########
		// #pragma omp parallel for collapse(3)
		cout<<"# Setting up NLP sine values"<<endl;
		for (size_t j=0; j<nk_; ++j)
				for (size_t g=0; g<nk_; g++)
						for (size_t h=0; h<nk2_; h++)
								sin_(j,g*nk2_+h) = sin(2*M_PI/nk_*h*(j-g));
		cout<<"## End: WignerFunction default constructor"<<endl;
	}  // End of constructor

	WignerFunction(size_t i_nx, double i_lD, double i_lC,
		size_t i_nk, double i_kmax) :
		nx_ (i_nx),
		lD_ (i_lD),
		lC_ (i_lC),
		l_ (lD_ + 2*lC_),
		dx_ (l_/float(nx_)),
		nk_ (i_nk),
		kmax_ (i_kmax > 0 ? i_kmax : M_PI/2./dx_),
		dk_ (2.*kmax_/float(nk_)),
		nk2_ (size_t(nk_/2.)),
		nxk_ (nx_*nk_),
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
		currD_(vec(nx_, fill::zeros)),
		a_(sp_mat(nxk_, nxk_)),
		b_(vec(nxk_, fill::zeros))
		{
		cout<<"## Start: WignerFunction constructor"<<endl;
		// ########## Configuration space array values ##########
		cout<<"# Setting up configuration space array values"<<endl;
		for (size_t i=0; i<nx_; ++i) x_(i) = i*dx_;
		// ########## Wave vector space array values ##########
		cout<<"# Setting up wave vector space array values"<<endl;
		for (size_t j=0; j<nk_; ++j) k_(j) = dk_*(j-(nk_-1)*.5);
		// ########## NLP sinus values ##########
		// #pragma omp parallel for collapse(3)
		cout<<"# Setting up NLP sine values"<<endl;
		for (size_t j=0; j<nk_; ++j)
				for (size_t g=0; g<nk_; g++)
						for (size_t h=0; h<nk2_; h++)
								sin_(j,g*nk2_+h) = sin(2*M_PI/nk_*h*(j-g));
		cout<<"## End: WignerFunction constructor"<<endl;
	}  // End of constructor

	~WignerFunction(){}

	double get_nx() { return nx_; }
	double get_nk() { return nk_; }
	double get_dk() { return dk_; }
	double get_dx() { return dx_; }
	double get_l() { return l_; }
	double get_lD() { return lD_; }
	double get_lC() { return lC_; }
	double get_m() { return m_; }
	double get_temp() { return temp_; }
	double get_cD() { return cD_; }
	double get_uF() { return uF_; }
	double get_dt(){ return dt_; }
	double get_rR() { return rR_; }
	double get_rM() { return rM_; }
	double get_rF() { return rF_; }
	double get_rG() { return rG_; }
	double get_lambda() { return lambda_; }
	vec get_x_arr() { return x_; }
	vec get_k_arr() { return k_; }
	vec get_u() { return u_; }
	vec get_uB() { return uB_; }
	vec get_uC() { return uC_; }
	vec get_uStart() { return uStart_; }
	vec get_du() { return du_; }
	vec get_dddu() { return dddu_; }
	vec get_currD() { return currD_; }
	vec get_cdX() { return cdX_; }
	vec get_cdK() { return cdK_; }

	void set_m(double m) { m_ = m; }
	void set_temp(double temp) {
		temp_ = temp;
		beta_ = 1./KB/temp_*AU_eV;
	}
	void set_cD(double cD) { cD_ = cD; }
	void set_uF(double uF) {
		uF_ = uF;
		uL_ = uF_, uR_ = uF_;
	}
	void set_dt(double dt) {
		dt_ = dt;
		courant_num_ = dt_*kmax_/m_/dx_;
	}

	void set_rR(double rR) { rR_ = rR; }
	void set_rM(double rM) { rM_ = rM; }
	void set_rF(double rF) { rF_ = rF; }
	void set_rG(double rG) { rG_ = rG; }
	void set_lambda(double lambda) { lambda_ = lambda; }

	void set_epsilonR(double epsilonR) { epsilonR_ = epsilonR; }
	void set_uBias(double uBias) { uBias_ = uBias; }

	void set_useNLP(bool useNLP) { useNLP_ = useNLP; }
	void set_useQC(bool useQC) { useQC_ = useQC; }
	void set_uBias_BC(bool uBias_BC) { uBias_BC_ = uBias_BC; }

	void set_bcType(int bcType) { bcType_ = bcType; }

	void set_lYZ(double lYZ) { lYZ_ = lYZ; }
	void set_part_num(double part_num) { part_num_ = part_num; }

	void set_uC(vec uC) { uC_ = uC; }

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
	void quantumCorrTerm(size_t, size_t, double);

	void solveRec();							// Wigner equation solved recursively
	void solveMatrixEq();						// Wigner equation solved by solving matrix equation

	void solveWignerPoisson(double, double, double, size_t);
	void calc_IVchar(double, double, size_t);
	void calcMobility();

	double fermiDirac(double);
	double supplyFunction(double);		// Supply function
	double sf_x(double, double);					// Supply function as function of x

	double maxwell_boltzmann(double);
	double gaussian_bc(double);
	double gaussian_x(double, double);

	double eqFun_x(double, double);
	double lorentz(double);	 // Lorentz profile
	double gauss(double);	 // Gauss
	double voigt(double);	 // Voigt profile

	// wignerTools.cpp
	void setPotBias(double);
	void addGaussBarr(double, double, double);
	void addRectBarr(double, double, double, double);
	void addWavePacket();
	double WavePacket_TEV(double, double);
	double nC(double, double);
	double fermiInt(double, double);
	double calcFermiEn(double, double, double);
	double calcFermiEn_MB(double, double, double);
};

}

#endif
