// Copyright 2008-2009,2016-2017 Maciej Woloszyn

#ifndef WIGNER_EQUATION_SOLVER_H
#define WIGNER_EQUATION_SOLVER_H 1

//# include <vector>

namespace mw {

/**
 * Wigner-Poisson self-consistent solver. All values in atomic units [a.u.]!
 * See notes in MW/QM1/37.
 *
 * @author Maciej Woloszyn
 *
 * @todo electric perm. eps = eps(x)
 * @todo better way of choosing dk and k_max ... --> affects discr. of potential term ???
 * @todo time-dependent simulation
 * @todo add method to apply extra voltage only in some given range (e.g. gate)
 * @todo relaxation time as a function of k
 * @todo allow Dirac delta terms in potential ?
 * @todo effective mass m = m(x)
 * @todo other discretisation schemes (UDS3?)
 *
 * @todo inline functions for xi, ki, fW etc. (and use them in code!)
 * @todo better documentation (arguments of functions, returned values, ...)
 * @todo cleanup #include's
 */
class WignerEquationSolver {
//  protected:
public:
    enum BoundaryConditionType
    {
       SUPPLYFUNCTION = 0,
       FERMIDIRAC
    };

    // parameters used while solving the Wigner Eq.
    const double m; ///< effective mass of electron [a.u. = m_e]

    const size_t Nx; ///< \f$ N_x; i = 0,..,(N_x - 1) \f$
    const size_t Nk; ///< \f$ N_k; j = 0,..,(N_k - 1) \f$
    const size_t Nxk; ///< \f$ N_{xk} = N_x  N_k \f$
    const size_t K; ///< K = 2*N_k; number of super/subdiagonals in A, for band storage
    const size_t LDAB; ///< leading dim. of AB; LDAB = 3*K+1
    const size_t sizeAB; ///< size of the array AB; sizeAB = (LDAB*N_xk)
    const double dx; ///< \f$ \Delta_x \f$ [a.u.]
    const double dk; ///< \f$ \Delta_k \f$ [a.u.] = \f$ ( \pi / (N_k*dx) ) \f$
    const double EF; ///< Fermi energy of contacts [a.u.]

    // other parameters:
    const std::vector<double>& U_B; ///< confining potential [a.u.], as loaded from the input file (el.field is included inly in U[] !)
    const size_t Nx_contact; ///< number of points before "spacers" (left and right) [a.u.]
    //const double x_space_left; ///< size of the left contact "spacer"
    //const double x_space_right; ///< size of the right contact "spacer"
    const std::string info; ///< description, e.g. about potential
    //
    double* U; ///< potential U(x) (including electric field, if present!) [a.u.]
    double* AB; ///< band matrix for storing A
    double* B; ///< matrix for Ax=B; size of A: Nxk
    double T; ///< temperature [K]
    double beta; ///< 1/(k_B*T) [a.u.]
    double V; ///< voltage [a.u.]
    double* Sin; ///< array of sin((2.0*M_PI/(double)Nk)*(double)(q*(j-l))) values; size: 2+(Nk+2)(Nk-1)/2
    double tau; ///< relaxation time [a.u.] (enables scattering only if tau>0)
    double tauGamma; ///< relaxation time [a.u.] (enables scattering only if tauGamma>0)
    double Gamma; ///< 1/(2 tauGamma) [a.u.]

    double* B0; ///< solution for zero voltage, used in scattering 'equillibrium' term
    double* C0; ///< solution(charge) for zero voltage, used in scattering 'equillibrium' term
    double* BC0; ///< boundary condition bc(k,V=0), tabelarized

  public:
    size_t self_consistent_max_steps = 0; ///< max. no. of self-consistent iterations
    double self_consistent_epsilon = 1e-3; ///< convergence parameter
    double self_consistent_omega = 0.05; ///< self-consistent underrelaxation parameter

    double p = 0.004; ///< N_D(x) "smoothness"
    double epsilon = 1.0; ///< relative permittivity

    bool save_relax_steps = false; ///< save rho, phi rtc in all iterations

    double sc_diff = 0.; ///< "diff" after the self-cons calc
    size_t sc_iter = 0; ///< no. of performed iterations

    BoundaryConditionType bctype = SUPPLYFUNCTION; ///< type of boundary condition

  public:
    /**
     * Default constructor.
     */
    WignerEquationSolver(size_t _Nk, double _T, double _m, double _dx,
                         double _EF, const std::vector<double>& _potential,
                         size_t _Nx_space, double _tau,  double _tauGamma,
                         std::string _info,
                         BoundaryConditionType _bctype = SUPPLYFUNCTION);

    /**
     * Default destructor.
     */
    ~WignerEquationSolver();

    /**
     * Set temp. [K]
     */
    void setT(double _T);

    /**
     * Set relaxation time [a.u.]
     */
    void setTau(double _tau) { tau = _tau; }

    /**
     * Set relaxation time for contacts [a.u.]
     */
    void setTauGamma(double _tauGamma) { tauGamma = _tauGamma; Gamma = 1./(2*tauGamma); }

    /**
     * Set applied voltage [a.u.]
     * @todo TODO call fillMatrices() from here ?
     */
    void setV(double _V);

    /**
     * Apllied voltage [a.u.]
     */
    double getV() const { return V; }

    /**
     * List of parameters + date and hostname
     */
    std::string parameters() const;

    /**
     * Prepare A and B by calling discretising procedures..
     */
    void fillMatrices(bool scattering = true);

    /**
     * Prepare A and B by calling discretising procedures..
     */
    void fillBC0();

    /**
     * one iteration of self-consist procedure
     */
    void solve_one_iter();

    /**
     * Solve eq. using LAPACK
     */
    void solve(double V_bias, bool scattering = true);

    /**
     * Average density of current (for UDS2 discr.).
     */
    double currentDensity() const;

    /**
     * Density of charge (for UDS2 discr.) at x=i*dx
     * @param i number of node along x-axis, i=1..Nx
     */
    double chargeDensity(size_t i) const;

    /**
     * Density of charge (for UDS2 discr.)
     */
    std::vector<double> chargeDensity() const;

    /**
     * Wigner function f_w(x,k,t) values.
     */
    void printWignerFunction(std::ostream& _out) const;

    /**
     * Charge density n(x,t) values (for UDS2 discr.)
     */
    void printChargeDensity(std::ostream& _out) const;

    /**
     * Current j(x,t) values (for UDS2 discr.)
     * @todo TODO merge with @c currentDensity()
     */
    void printCurrentDensity(std::ostream& _out) const;

    /**
     * Potential U(x)
     */
    void printU(std::ostream& _out) const;

    /**
     * b.c. BC0(k)
     */
    void printBC0(std::ostream& _out) const;
    
    /**
     * Non-classicality parameter nu
     */
    double nonclParam() const;

    /**
     * Expectation vlaues of <x>, <x^2>, <k>, <k^2>, <xk>
     */
    std::vector<double> expectVals() const;

//  protected:
    /**
     * Boundary condition using the Fermi-Dirac distribution.
     * @param k wave number [a.u.]
     * @param V_shift voltage value in the contact [a.u.]
     */
    double bc(double k, double V_shift) const;

    /**
     * Convenience function calling bc(k, 0.0)
     * (boundary condition at the left edge)
     */
    double bc0(double k) const  { return bc(k, 0.); }

    /**
     * Convenience function calling bc(k, V);
     * (boundary condition at the right edge)
     */
    //double bcV(double k) const { return bc(k, V); }
    
    /**
     * Distr. function for multi-layer 1D system f0(E, mu);
     */
    double f0(double E, double mu) const;

    /**
      * N_D(x) - density of dopants (n-doped contacts)
      */
    double nD(size_t i) const;

    /**
     * Lorentzian
     */
    double gGamma(double E) const { return Gamma/(M_PI*(Gamma*Gamma+E*E)); }
    
    /**
     * Discretise diffusion term (with d/dx) using UDS2.
     * \$ ... \$
     * @param i 1..Nx
     * @param j 1..Nk
     */
    void discretiseDiffusion(size_t i, size_t j);

    /**
     * Discretise drift term (with Wigner potential).
     * @param i 1..Nx
     * @param j 1..Nk
     */
    void discretiseDrift(size_t i, size_t j);

    /**
     * Discretise scattering term (with Wigner potential).
     * @param i 1..Nx
     * @param j 1..Nk
     */
    void discretiseScattering(size_t i, size_t j);


    inline double& _AB(size_t i, size_t j) {
#ifdef DEBUG
      // sizeAB = (LDAB*N_xk); LDAB rows, N_xk columns
      if(i>=LDAB) {
        std::cerr << "\nERR: i="<<i<<" j="<<j<<std::endl;
        throw std::out_of_range(std::string("_AB() i>=LDAB"));
      }
      if(j>=Nxk) {
        std::cerr << "\nERR: i="<<i<<" j="<<j<<std::endl;
        throw std::out_of_range(std::string("_AB() j>=Nxk"));
      }
#endif
      return AB[j*LDAB+i]; // for explanation see below (in _A(i,j) )
    }


    inline double& _A(size_t i, size_t j) {
#ifdef DEBUG
      try {
        if(i>=Nxk)  throw std::out_of_range(std::string("_A() i>=Nxk"));
        if(j>=Nxk)  throw std::out_of_range(std::string("_A() j>=Nxk"));
      // if(i < max(0,j-K) ) ...
      // if(i > min(Nxk-1, j+k) ) ...
      // max(0,j-K)<=i<=min(Nxk-1,j+K)
        size_t min = 0;
        if(j>K) min=j-K;
        if(i < min ) throw std::out_of_range(std::string("_A() i<max(0,j-K)"));
        size_t max = j+K;
        if(max>Nxk-1) max=Nxk-1;
        if(i > max ) throw std::out_of_range(std::string("_A() i>min(Nxk-1, j+k)"));
      } catch(std::out_of_range e) {
        //std::cerr << "\nERR: i="<<i<<" j="<<j<<" K="<<K<<std::endl;
        throw;
      }
#endif
      return _AB(2*K+i-j, j);

      // return A[i*Nxk+j];
      // return A[j*Nxk+i]; // avoid transposing
      
      // http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html
      // Order of multi dimensional arrays in C/C++ is the opposite of FORTRAN.
      // Native FORTRAN layout (collumn-major order): INTEGER A(2,3)
      //      a(1,1) a(2,1) a(1,2) a(2,2) a(1,3) a(2,3)
      // or starting from 0:
      //      a(0,0) a(1,0) a(0,1) a(1,1) a(0,2) a(1,2)
      // Native C layout (row-major order) is NOT equivalent to the FORTRAN layout: int a[2][3];
      //      a[0][0] a[0][1] a[0][2] a[1][0] a[1][1] a[1][2]
      //
      // For A[N][M] or A(N,M) : (N rows, M columns)
      // - Fortran-style: A[i][j] <--> A[j*N+i]
      // - C-style:       A[i][j] <--> A[i*M+j]
    }

    inline double& _B(size_t i) {
#ifdef DEBUG
      if(i>=Nxk)  throw std::out_of_range(std::string("_B()"));
#endif
      return B[i];
    }

    inline double& _B0(size_t i) {
#ifdef DEBUG
      if(i>=Nxk)  throw std::out_of_range(std::string("_B0()"));
#endif
      return B0[i];
    }

    inline double& _C0(size_t i) {
#ifdef DEBUG
      if(i>=Nx)  throw std::out_of_range(std::string("_C0()"));
#endif
      return C0[i];
    }

    inline double& _B(size_t i) const {
#ifdef DEBUG
      if(i>=Nxk)  throw std::out_of_range(std::string("_B()"));
#endif
      return B[i];
    }

    inline double& _U(size_t i) {
#ifdef DEBUG
      if(i>=Nx) throw std::out_of_range(std::string("_U()"));
#endif
      return U[i];
    }

    inline double& _Sin(size_t q, size_t j, size_t l) {
      size_t i = q*(j-l);
#ifdef DEBUG
      if(j<l) {
        std::cerr << "\nERR: i="<<i<<" q="<<q <<" j="<<j <<" l="<<l <<std::endl;
        throw std::out_of_range(std::string("_Sin() j<l "));
      }
      if(i>(Nk*(Nk-1)/2) ) throw std::out_of_range(std::string("_Sin() i>(Nk*(Nk-1)/2)  "));
#endif
      return Sin[i];
    }


}; // class WignerEquationSolver

} // namespace wignerequation

#endif // WignerEquationSolver.h
