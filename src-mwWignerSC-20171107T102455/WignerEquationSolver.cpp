// File WignerEquationSolver.cpp
// Copyright 2008,2009,2016 Maciej Woloszyn, Pawel Wojcik

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <omp.h>
#include <acml.h>

#include "mwAtomicUnits.h"
#include "mwTools.h"
#include "WignerEquationSolver.h"
#include "Poisson1D.h"

namespace mw {

WignerEquationSolver::WignerEquationSolver(
        size_t _Nk, double _T, double _m, double _dx, double _EF,
        const std::vector<double>& _potential,
        size_t _Nx_space, double _tau,  double _tauGamma,
        std::string _info, BoundaryConditionType _bctype
        )
    : m(_m),
      Nx(_potential.size()),
      Nk(_Nk),
      Nxk(Nx*Nk),
      K(2*Nk),
      LDAB(3*K+1),
      sizeAB(LDAB*Nxk),
      dx(_dx),
      dk(M_PI/((int)Nk*dx)),
      EF(_EF),
      U_B(_potential),
      Nx_contact(_Nx_space),
      info(_info),
      V(0.0),
      tau(_tau),
      tauGamma(_tauGamma),
      bctype(_bctype)
{
#ifdef DEBUG
    std::cerr<< "#E> WignerEquationSolver(), sizeAB = "<<LDAB<<"*"<<Nxk<<" = "<<sizeAB<<"\n";
#endif
    B = new double[Nxk];
    AB = new double[sizeAB];
    U = new double[Nx];

    // prepare potential en. values etc.
    for(size_t i=0; i<Nx; ++i) {
        U[i] = U_B[i];
    }

    setT(_T); // sets beta as well
    setTauGamma(_tauGamma); // sets Gamma as well

    BC0 = new double[Nk];
    fillBC0();

    // TMP
    std::ofstream fu("BC0.out");
    //fu << infoProgram << " " << w.parameters() << endl;
    printBC0(fu);
    fu.close();

    // prepare array of sin(...) values --> optimization
    Sin = new double[1+Nk*(Nk-1)/2];
    for(size_t i=0; i<=Nk*(Nk-1)/2; ++i) {
        Sin[i] = sin((2.0*M_PI/(double)Nk)*(double)(i));
    }

#ifdef DEBUG
    std::cerr<< "#E> WignerEquationSolver() V=0\n";
#endif
    /// @todo TODO do it only if scattering turned on
    /// @todo TODO do it only once if the result for V=0 needed anyway
    B0 = new double[Nxk];
    C0 = new double[Nx];
    // solve for V=0.0 --> needed for scattering term
    if(tau>0.0) {
        //throw std::logic_error(std::string("equil solution not ready!! --> WignerEquationSolver()"));
        //TODO--
        //fillMatrices(false); //scattering=false; here we want 'equilibrium' solution
        //solve();
        solve(0.0,false); // check!!!
        //--

        // arma::fill in B0 and C0 :
#pragma omp parallel for
        for(size_t i=0; i<Nxk; ++i) B0[i] = B[i];
#pragma omp parallel for
        for(size_t i=0; i<Nx; ++i) C0[i] = chargeDensity(i+1);
#ifdef DEBUG
        std::cerr << "#E> WignerEquationSolver() (V, j)= " << getV() << " " <<
                     currentDensity() << std::endl;
#endif
    }
}


WignerEquationSolver::~WignerEquationSolver()
{
#ifdef DEBUG
    std::cerr << "#E> ~WignerEquationSolver()\n";
#endif
    delete[] U;
    delete[] B;
    delete[] AB;
    delete[] B0;
    delete[] C0;
    delete[] BC0;
    delete[] Sin;
}


void WignerEquationSolver::setT(double _T) {
    T = _T;
    beta = 1.0/(AU_kB*T);
}


void WignerEquationSolver::setV(double _V) {
#ifdef DEBUG
    std::cerr << "#E> WignerEquationSolver::setV("<<_V<<")\n";
#endif
    V = _V;
    // modify the potential accordingly to V_bias
    const double step = V/(int)(Nx-2*Nx_contact);
    for (size_t i = Nx_contact; i < (Nx-Nx_contact); i++) {
        // linear change inside the device
        U[i] = U_B[i] - step*(int)(i-Nx_contact);
    }
    for (size_t i = (Nx-Nx_contact); i < Nx; i++) {
        // for contacts: only right one is shifted (as a whole)
        U[i] = U_B[i] - V;
    }
}


std::string WignerEquationSolver::parameters() const
{
    std::stringstream s;
    s << "\n#[atomic units]";
    s << "\n#  T= " << T;
    s << "\n#  M_EFF= " << (m);
    s << "\n#  E_F= " << (EF);
    s << "\n#  TAU= " << (tau);
    s << "\n#  TAUGAMMA= " << (tauGamma);
    s << "\n#  NX= " << (Nx);
    s << "\n#  NK= " << (Nk);
    s << "\n#  DX= " << (dx);
    s << "\n#  DK= " << (dk);
    s << "\n#  NX_CONTACT= " << (Nx_contact);
    s << "\n#  ND_SMOOTHNESS_P= " << (p);
    s << "\n#  SELF_CONS_MAX_STEPS= " << (self_consistent_max_steps);
    s << "\n#  SELF_CONS_OMEGA= " << (self_consistent_omega);
    s << "\n#  SELF_CONS_EPS= " << (self_consistent_epsilon);
    s << "\n#  BC_TYPE= " << (bctype);
    s << "\n#  V= " << (V);
    //s << "  X_R= " + std::to_string(x_space_right);
    s << "\n#  potential_info= " << info << "\n#\n";
    s << "## Start (" << mw::host() << ") " << mw::time() ;
    return s.str();
}

// TODO: rename this --> f0(k,V)
double WignerEquationSolver::f0(double E, double mu) const
{
    double result = 0.0;

    if(bctype == SUPPLYFUNCTION) {
        //const double C = -beta*((k*k)/(2.0*m)-(EF-V_shift));
        const double C = -beta*(E-mu);
        result = ( m / (M_PI*beta) );
        // for C>>1 ln(1+exp(C)) --> C
        // if C=30 relative difference is approx. 10^-15
        // (see e.g. gnuplot: plot [29:31] (x-log(1+exp(x)))/x )
        /// @todo TODO replace 30. with some general "epsilon"...
        if (C > 30.) {
            result *= C;
        } else {
            result *= log(1.0 + exp(C));
        }
    } else if(bctype == FERMIDIRAC) {
        result = f_FD(E, mu, T);
    }

    return result;
}

double WignerEquationSolver::nD(size_t i) const
{
    // not ready for 3D -- other formula required
    if(bctype != SUPPLYFUNCTION) {
        throw std::logic_error(std::string("nD() : implemented only for 1D (i.e. supply function as b.c.)"));
    }

    // charge density [au] for uniform distribution function in layered systems
    const double ND = pow(2*m*EF,1.5)/(3*M_PI*M_PI);
    // N_D(x) for uniform distribution function in layered systems, exp-smoothed
    return ND*(1.0
               + 1.0/(1.0+std::exp((int)(i-Nx_contact)/(p*(double)Nx)))
               - 1.0/(1.0+std::exp((int)(i-Nx+Nx_contact)/(p*(double)Nx)))
               );
}


double WignerEquationSolver::bc(double k, double V_shift) const
{
    // we should never have k=0
    /// @todo TODO replace 1.e-8 with some general "epsilon"... e.g. MIN_DOUBLE ?
    if(fabs(k) < 1.e-8) {
        std::cerr << ">> bc: k=" << k << "\n";
        throw std::logic_error(std::string("bc() : |k| < 1.e-8"));
    }

    double mu = EF - V_shift;
    double Ek = k*k/(2.*m);
    double result = 0.;
    // TODO Gamma --> new function bc(k,V) which calls f0(k,V)
    if (tauGamma > 0) { // integrate gGamma*f0
        //result = f0(Ek,mu); // TMP !!!
        // integrate_here
        double dE = Gamma/100; // TMP
        double E = 0;
        int iE = 0;
        double F = 1;
        while( E<Ek || F>1e-15 ) { // TMP
            E = iE*dE;
            F = gGamma(Ek-E)*f0(E,mu);
            result += F;
            iE++;
        }
        result *= dE;
        std::cerr << "#E> iE = " << iE <<std::endl;
    } else {
        result = f0(Ek,mu);
    }
    return result;
}

void WignerEquationSolver::fillBC0() {
    std::cerr << "#E> fillBC0()  : " << mw::time() <<std::endl;
    /// @todo TODO OpenMP
#pragma omp parallel for
    for(size_t j=1;j<=Nk;j++) {
        const double k = (M_PI/((int)Nk*dx))*((int)(j-1)-0.5*((int)Nk-1));
        BC0[j-1] = bc0(k);
    }
}


void WignerEquationSolver::discretiseDiffusion(size_t i, size_t j) {
    const double k = (M_PI/((int)Nk*dx))*((int)(j-1)-0.5*((int)Nk-1));
    const double C = k/(2.*m*dx);
    //Since we start from i=1 and j=1 we have to substract 1
    const size_t row=(i-1)*Nk+(j-1);

    // we should never have k=0
    /// @todo TODO replace 1.e-8 with some general "epsilon"... e.g. MIN_DOUBLE ?
    if(fabs(k) < 1.e-8) {
        std::cerr << ">> discretiseDiffusion: row=" << row << " k=" << k << " Nxk="
                  << Nxk << "\n";
        throw std::logic_error(std::string("discretiseDiffusion() : |k| < 1.e-8"));
    }

    /// @todo TODO cleanup: e.g. don't use 'i*Nk+j-1' but 'row-Nk' or similar
    if(k<=0.0)
    {
        /// @todo check why we use bc(k,0) and not bc(k,Vdrain)?
        /// seems OK: E(k) is also shifted when U(R)=E-e*Vdrain (but check it!)
        if(i==Nx) {
            _A(row, row) += -3.0*C;
            //_B(row) += -3.*C*bc0(k);
            _B(row) += -3.*C*BC0[j-1];
        }
        else if(i==(Nx-1)) {
            _A(row, row) += -3.0*C;
            _A(row, i*Nk+j-1) += 4.0*C;
            //_B(row) += C*bc0(k);
            _B(row) += C*BC0[j-1];
        }
        else {
            _A(row, row) += -3.0*C;
            _A(row, i*Nk+j-1) += 4.0*C;
            _A(row, (i+1)*Nk+j-1) += -C;
        }
    }
    else { // <=> if(k>0)
        if(i==1) {
            _A(row, row) += 3.0*C;
            //_B(row) += 3.*C*bc0(k);
            _B(row) += 3.*C*BC0[j-1];
        }
        else if(i==2) {
            _A(row, row) += 3.0*C;
            _A(row, (i-2)*Nk+j-1) += -4.0*C;
            //_B(row) += -C*bc0(k);
            _B(row) += -C*BC0[j-1];
        }
        else {
            _A(row, row) += 3.0*C;
            _A(row, (i-2)*Nk+j-1) += -4.0*C;
            _A(row, (i-3)*Nk+j-1) += C;
        }
    }
}


void WignerEquationSolver::discretiseDrift(size_t i, size_t j) {
    size_t l,q;
    double sum,v1,v2;

    /// @todo TODO cleanup: row-->(row-1), to use _A(row,..), not _A(row-1,..)
    const size_t row=(i-1)*Nk+j;
    //double k=(M_PI/(Nk*dx))*((j-1)-0.5*(Nk-1));

    for(l=1;l<=Nk;l++) {
        const size_t vrs=(i-1)*Nk+l;
        // If we want to use parity of matrix
        //if( ((i-1)*Nk+l)<=row) {  // same as if(l<=j) ?
        if( vrs <= row ) {
            /// @todo TODO is it possible to have this integral tabelarised?
            sum=0.0;
            for(q=1;q<=Nk/2;q++) {
                if ((i+q)>Nx)
                    v1=_U(Nx-1);
                else
                    v1=_U(i+q-1);
                //
                if (i < q+1) //((i-q)<1)
                    v2=_U(0);
                else
                    v2=_U(i-q-1);
                //
                sum += (v1-v2)*_Sin(q,j,l);
                /// @todo TODO in case of j==l Sin is always zero ! optimize, see below
#ifdef DEBUG
                // j and l are unsigned !
                //   ( "j>=l" should be true beacause of "if(((i-1)*Nk+l)<=row)" above )
                if(j<l) throw std::logic_error(std::string("j<l in discretiseDrift()"));
#endif
            } // end for(q=...)

            const double C = 2.*sum/(int)Nk;
            // to avoid taking diagonal terms twice...
            if( vrs < row ) {
                _A( (row-1), (vrs-1) ) += C;
                _A( (vrs-1), (row-1) ) += -C;
            }
            else { // ((i-1)*Nk+l) == row)
                _A( (row-1), (row-1) ) += C;
                /// @todo TODO in this case C is always zero
                //  std::cerr<<"### C="<<C;
            }

        }
    } // end for(l=...)
}


void WignerEquationSolver::discretiseScattering(size_t i, size_t j) {
#ifdef DEBUG
    std::cerr << "#E> discretiseScattering i="<<i<<" j="<<j<< "\n";
#endif
    const size_t row=(i-1)*Nk+(j-1);
    size_t q;
    double alfa=dk/(2.*M_PI);
    for(q=1; q<=Nk; ++q) {
        const size_t vrs=(i-1)*Nk+(q-1);
        if(q==j) {
            _A(row,row) += -( alfa*_B0(row)/_C0(i-1)-1.0 )/tau; // TODO CHECK sign!
        }
        else {  // q!=j
            _A(row,vrs) += -( alfa*_B0(row)/_C0(i-1) )/tau; // TODO CHECK sign!
        }
    }
}


void WignerEquationSolver::fillMatrices(bool scattering) {

    std::cerr << "#E> fillMatrices() - init with zeros  : " << mw::time()
              <<std::endl;
    // initialise A and B with zeros
#pragma omp parallel for
    for(size_t i=0; i<Nxk; ++i) B[i] = 0.0;
#pragma omp parallel for
    for(size_t i=0; i<sizeAB; ++i) AB[i] = 0.0;

    std::cerr << "#E> fillMatrices() - discretiseDiffusion  : "<< mw::time()
              << std::endl;
    for(size_t i=1; i<=Nx; i++)
        for(size_t j=1; j<=Nk; j++)
            discretiseDiffusion(i,j);

    std::cerr << "#E> fillMatrices() - discretiseDrift  : " << mw::time()
              << std::endl;
    {
        size_t i, j;
#pragma omp parallel for default(shared) private(i, j) schedule(static)
        for(i=1;i<=Nx;i++) {
            for(j=1;j<=Nk;j++) {
                discretiseDrift(i,j);
            }
        }
    }

    /// @todo TODO clean up this condition (?)
    if(scattering && tau>0.0) {
        std::cerr << "#E> fillMatrices() - discretiseScattering  : " <<
                     mw::time() << std::endl;
        {
            size_t i, j;
#pragma omp parallel for default(shared) private(i, j) schedule(static)
            for(i=1;i<=Nx;i++) {
                for(j=1;j<=Nk;j++) {
                    discretiseScattering(i,j);
                }
            }
        }
    }

}


void WignerEquationSolver::solve_one_iter() {
    // http://www.netlib.org/lapack/double/dgbsv.f
    // F77:  DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
    // ACML: dgbsv(int n, int nsub, int nsuper, int nrhs, double *a, int lda,
    //             int *ipiv, double *b, int ldb, int *info);

    int info = 0;
    int* ipiv = new int[Nxk];
    for(size_t i=0; i<Nxk; ++i) ipiv[i] = 0;

#ifdef DEBUG
    {
        std::cerr << "#E> Saving matrix A to file 'A.out' ...\n";
        std::ofstream of("A.out");
        of << "# matrix A  [ from WignerEquationSolver::solve_one_iter() ]";
        of << " Nx="<<Nx;
        of << " Nk="<<Nk;
        of << " Nxk="<<Nxk;
        of << "\n";
        of << "# '1' denotes non-zero array element!\n";
        for(size_t j=0; j<Nxk; ++j) {
            for(size_t i=0; i<Nxk; ++i) {
                double aij = 0.0;
                try {
                    aij = _A(j,i);
                } catch (std::out_of_range e) {
                }
                of << " " << ( aij != 0.0 ? 1 : 0) ;
                //of << " " <<  aij ;
            }
            of << "\n";
        }
        of.close();
    }
    {
        std::cerr << "#E> Saving matrix AB (A-band-storage) to file 'AB.out' ...\n";
        std::ofstream of("AB.out");
        of << "# matrix AB  [ from WignerEquationSolver::solve_one_iter() ]";
        of << " LDAB="<<LDAB;
        of << "\n";
        of << "# '1' denotes non-zero array element!\n";
        for(size_t j=0; j<LDAB; ++j) {
            for(size_t i=0; i<Nxk; ++i) {
                of << " " << (_AB(j,i) != 0.0 ? 1 : 0) ;
            }
            of << "\n";
        }
        of.close();
    }
#endif

    std::cerr << "#E> solve_one_iter() : DGBSV start " << mw::time() << std::endl;
    //dgbsv(n,nsub,nsuper,nrhs,a,lda, ipiv, b, ldb, info);
    dgbsv((int)Nxk, (int)K, (int)K, 1, AB, (int)LDAB, ipiv, B, (int)Nxk, &info);
    std::cerr << "#E> solve_one_iter() : DGBSV stop  " << mw::time() << std::endl;

    if(info < 0)
        throw std::logic_error(
                std::string("dgesv/dbgsv: the i-th argument had an illegal value, i=")
                + std::to_string(info)
                );
    if(info > 0)
        throw std::logic_error(
                std::string("dgesv/dbgsv: U(i,i) is exactly zero.  The factorization has been completed, but the factor U is exactly singular, so the solution could not be computed. i=")
                + std::to_string(info)
                );

#ifdef DEBUG
    {
        std::cerr << "#E> Saving solution matrix B to 'B.out' ...\n";
        std::ofstream of("B.out");
        of << "# matrix B  [ from WignerEquationSolver::solve_one_iter() and DGBSV ]\n";
        for(size_t i=0; i<Nxk; ++i) {
            of << _B(i) << "\n";
        }
        of.close();
    }
#endif

    delete[] ipiv;
}

void WignerEquationSolver::solve(double V_bias, bool scattering)
{
    // TODO use n(x) from previous V_bias as starting n


    // a) without self-cons
    if(self_consistent_max_steps<=0) {
        setV(V_bias); // set profile for spacers
        fillMatrices(scattering);
        solve_one_iter();
        return;
    }

    // b) self-cons calc
    V = V_bias;
    double w = self_consistent_omega;

    // solve Poisson to find U_H(x)
    std::vector<double> rho(Nx,0.0);
    std::vector<double> rho_prev(Nx,0.0);
    std::vector<double> n(Nx,0.0);
    std::vector<double> n_prev(Nx,0.0);
    double n0 = 0.;

    // TODO check Phi.out ---> what is written to that file?
    std::ofstream fphi("Phi.out");
    fphi << "#\n#x phi(x) rho(x) n(x) U(x)\n";

    for(size_t iter=0; iter<self_consistent_max_steps; iter++) {

        Poisson1D poisson(rho);
        poisson.epsilonc = epsilon*AU_epsilon0;
        // TODO check points at i=0 and i=Nx --> sth may be wrongly shifted by one dx...?
        poisson.phiL = 0.;
        poisson.phiR =  V_bias;
        poisson.dx = dx;
        poisson.solve();

        // U(x) = U_B(x) + U_H(x)
        // where U_H(x) = -e phi(x)
        for (size_t i = 0; i < Nx; i++) {
            U[i] = U_B[i] - poisson.phiAt(i);
        }

        // fillMatrices(scattering); // ???
        fillMatrices();
        solve_one_iter();

        n = chargeDensity();
        if(iter==0) //n0 = n;
            n0 = n.at(0);

        double diff = 0.;
        for(size_t i=0; i<Nx; ++i) {
            // under-relax:
            rho.at(i) = (1-w)*rho_prev.at(i) + w*(nD(i) - n.at(i));
            //
            diff += (n[i]-n_prev[i])*(n[i]-n_prev[i])/(n0*n0);
            if(save_relax_steps) {
                fphi << (int)i*dx << "\t" << poisson.phiAt(i) << "\t"
                     << rho.at(i) << "\t" << n.at(i) << "\t" << U[i] << "\n";
            }
        }
        fphi << "#diff,j= " << diff << " " << currentDensity() << "\n" << "\n\n";
        std::cerr << "#E> solve() : diff= " << diff << std::endl;

        sc_diff = diff;
        sc_iter = iter;

        if(iter>0 && diff<self_consistent_epsilon) {
            std::cerr << "#E> solve() : iter= " << iter << std::endl;
            break;
        }

        rho_prev = rho;
        n_prev = n;
    }
    fphi.close();
}

double WignerEquationSolver::nonclParam() const {
    size_t i;
    double w = 0.;
    double Ip = 0.;
    double Im = 0.;
    double I = 0.;
    // TODO OpenMP
    for (i=0;i<Nxk;i++) {
        w = _B(i);
        I += w;
        if(w>0) Ip += w;
        else Im += (-w);
    }
    return 1. - (Ip-Im)/(Ip+Im);
}

std::vector<double> WignerEquationSolver::expectVals() const {
    std::vector<double> result;
    size_t i,j;
    size_t site=0;
    double x,k,fW;
    double ex=0.0, ex2=0.0, ek=0.0, ek2=0.0, exk=0.0, norm=0.0;
    for (i=1;i<=Nx;i++)
    {
        x=((int)i-1)*dx;
        for (j=1;j<=Nk;j++)
        {
            k=dk*(((int)j-1)-0.5*((int)Nk-1));
            fW = _B(site); // Wigner function fW(x_i,k_j)
            ex += fW*x;
            ex2 += fW*x*x;
            ek += fW*k;
            ek2 += fW*k*k;
            exk += fW*x*k;
            norm += fW;
            site++;
        }
    }
    result.push_back(ex*dx*dk);
    result.push_back(ex2*dx*dk);
    result.push_back(ek*dx*dk);
    result.push_back(ek2*dx*dk);
    result.push_back(exk*dx*dk);
    result.push_back(norm*dx*dk);
    return result;
}

void WignerEquationSolver::printWignerFunction(std::ostream& _out) const {
    size_t i,j;
    size_t site=0;
    double x,k;
    _out << "## Wigner function [begin]\n";
    _out << "#  x[a.u.]   k[a.u.]   f_W(x,k)[a.u.]\n";
    for (i=1;i<=Nx;i++)
    {
        x=((int)i-1)*dx;
        for (j=1;j<=Nk;j++)
        {
            k=dk*(((int)j-1)-0.5*((int)Nk-1));
            _out << x << " " << k << " " << _B(site) << "\n";
            site++;
        }
        _out << "\n";
    }
    _out << "## Wigner function [end]\n";
}


void WignerEquationSolver::printChargeDensity(std::ostream& _out) const {

    _out << "## Charge density [begin]\n";
    _out << "#  x[a.u.]   n(x)[a.u.]\n";
    for (size_t i=1;i<=Nx;i++) {
        double x=((int)i-1)*dx;
        _out << x << " " << chargeDensity(i) << "\n";
    }
    _out << "## Charge density [end]\n";
}


double WignerEquationSolver::chargeDensity(size_t i) const {
#ifdef DEBUG
    if(i<1 || i>Nx) {
        std::cerr << ">> chargeDensity: i="<<i<<"\n";
        throw std::logic_error(std::string("chargeDensity: i out of range"));
    }
#endif
    size_t site;
    double charge=0.0;
    /*
  for(size_t j=1;j<=Nk-1;j++) {
    site=(i-1)*Nk+j;
    charge += (_B(site-1)+_B(site))/2.0; //trapezoid method
  }
*/
    for(size_t j=1;j<=Nk-1;j++) {
        site=(i-1)*Nk+j;
        charge += _B(site-1);
    }
    return  charge*dk/(2.*M_PI) ;
}

std::vector<double> WignerEquationSolver::chargeDensity() const
{
    std::vector<double> n(Nx);
    for (size_t i=0;i<Nx;i++) {
        n.at(i)=chargeDensity(i+1); //
    }
    return n;
}


void WignerEquationSolver::printCurrentDensity(std::ostream& _out) const {
    double k, x;
    double cur = 0.;
    _out << "## Current density [begin]\n";
    _out << "#  x[a.u.]   j(x)[a.u.]\n";

    // see SQAUDS for algorithm

    /// @todo TODO use currentDensity() instead of havin the same code in 2 places
    //if(uds==2)
    for (size_t i=2; i<=Nx-2; i++) {
        x=((int)i+0.5-1)*dx;
        cur = 0.;
        for (size_t j=1;j<=Nk;j++) {
            k=dk*(((int)j-1)-0.5*((int)Nk-1));
            if (k<=0.0)
                cur += k*(-_B((i+1)*Nk+j-1)+3.0*_B(i*Nk+j-1));
            else //if (k>0.0)
                cur += k*(3.0*_B((i-1)*Nk+j-1)-_B((i-2)*Nk+j-1));
        }
        _out << x << " " << cur*dk/(4.*M_PI*m) << "\n";
    }
    _out << "## Current density [end]\n";
}


double WignerEquationSolver::currentDensity() const {
    size_t i,j;
    double k;
    //double* current = new double[Nx];
    //for(i=0; i<Nx; i++) current[i]=0.0;
    double cur = 0.;

    // see SQAUDS for algorithm

    // TODO OpenMP
    //if(uds==2)
    for (i=2; i<=Nx-2; i++) {
        for (j=1;j<=Nk;j++) {
            k=dk*(((int)j-1)-0.5*((int)Nk-1)); //*(RTD_au2k/1e-9);
            if (k<=0.0)
                //current[i-1] += k*(-_B((i+1)*Nk+j-1)+3.0*_B(i*Nk+j-1));
                cur += k*(-_B((i+1)*Nk+j-1)+3.0*_B(i*Nk+j-1));
            else //if (k>0.0)
                //current[i-1] += k*(3.0*_B((i-1)*Nk+j-1)-_B((i-2)*Nk+j-1));
                cur += k*(3.0*_B((i-1)*Nk+j-1)-_B((i-2)*Nk+j-1));
        }
    }
    return cur*dk/(4.*M_PI*m)/((int)Nx-3.);
}


void WignerEquationSolver::printU(std::ostream& _out) const {
    _out << "## Potential [begin]\n";
    _out << "#  x[a.u.]   U(x)[a.u.]\n";
    for (size_t i=0;i<Nx;i++) {
        double x=(int)i*dx;
        _out << x << " " << U[i] << "\n";
    }
    _out << "## Potential [end]\n";
}

void WignerEquationSolver::printBC0(std::ostream& _out) const {
    _out << "## Boundary condition [begin]\n";
    _out << "#  k[a.u.]   BC0(k)[a.u.]\n";
    for (size_t j=1;j<=Nk;j++) {
        double k=dk*(((int)j-1)-0.5*((int)Nk-1));
        _out << k << " " << BC0[j] << "\n";
    }
    _out << "## Boundary condition [end]\n";
}

} // namespace wignerequation
