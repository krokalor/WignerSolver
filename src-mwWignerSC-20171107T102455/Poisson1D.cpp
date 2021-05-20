// File Poisson1D.cpp
//
// DO NOT edit outside SYNC/Sci/src/Poisson/
//
// Copyright 2014 Maciej Woloszyn

#include <iostream>
#include <cmath>
#include <stdexcept>
//#include <typeinfo>
//#include <iomanip>
//#include <sstream>
//#include <omp.h>

#include <acml.h>
// TODO : MKL

#include "mwTools.h"
#include "Poisson1D.h"


namespace mw {

// makes my "<<" operator available for vectors
//#ifdef DEBUG
//using namespace tools;
//#endif

void Poisson1D::solve()
{
    if(solved)
        throw std::logic_error( "Poisson1D::solve() : already solved" );

    DU = v1D(N-1, 1.);
    D = v1D(N, -2.);
    DL = v1D(N-1, 1.);

    if(const_epsilon) {
        const double minus_dx2e = -dx*dx / epsilonc;
        for(int i=0; i<N; ++i) {
            phi[i] *=  minus_dx2e; // -dx^2/epsilon * rho[i]
        }
        phi[0]   += -phiL;
        phi[N-1] += -phiR;
    }
    else { // epsilon=epsilon(x)

        /*   // OLD = QM1/61
      const double minus_dx2 = -dx*dx;
        double ksi = 0.;
        // i=0
        //ksi = ksiAt(0); //
        ksi = 0.25*(epsilon[2]-epsilon[0])/epsilon[0]; // ksi_1
        DU[0] += ksi;
        phi[0] = minus_dx2*phi[0]/epsilon[0] - phiL*(1-ksi);
        // i=N-1
        //ksi = ksiAt(N-2); //
        ksi = 0.25*(epsilon[N-1]-epsilon[N-3])/epsilon[N-1]; // ksi_N
        DL[N-2] += -ksi;
        phi[N-1] = minus_dx2*phi[N-1]/epsilon[N-1] - phiR*(1+ksi);
        // i= 1,...,N-2
        for(int i=1; i<N-1; ++i) {
            //ksi = ksiAt(i); //
            ksi = 0.25*(epsilon[i+1]-epsilon[i-1])/epsilon[i]; // ksi_i
            //DEBUG_PRINT(ksi) ;
            DU[i]   +=  ksi;
            DL[i-1] += -ksi;
            phi[i]  =  minus_dx2*phi[i]/epsilon[i]; // -dx^2/epsilon * rho[i]
        }
        */

        // NEW = QM1/142
        double dx2 = 2*dx*dx;
        for(int i=0; i<N-1; ++i) {
            double t = -(epsilon[i]+epsilon[i+1]);
            DU[i]=t;
            DL[i]=t;
        }
        D[0]   = 3*epsilon[0] + epsilon[1];
        D[N-1] = epsilon[N-2] + 3*epsilon[N-1];
        phi[0] = dx2*phi[0]+2*epsilon[0]*phiL;
        phi[N-1] = dx2*phi[N-1]+2*epsilon[N-1]*phiR;
        for(int i=1; i<N-1; ++i) {
            D[i] = epsilon[i-1]+2*epsilon[i]+epsilon[i+1];
            phi[i] = dx2*phi[i];
        }
        //

        
    }

    int info = 0;
    //DEBUG_PRINT(DU) ;
    //DEBUG_PRINT(D) ;
    //DEBUG_PRINT(DL) ;
    //DEBUG_PRINT(phi) ;
    //DEBUG_PRINT2("Poisson1D::solve() : DGTSV start") ;
    //DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
    // output: B=phi=solution, DL,D,DU also overwritten
    //dgtsv(N, 1, &DL[0], &D[0], &DU[0], &phi[0], N, &info);
    dgtsv(N, 1, DL.data(), D.data(), DU.data(), phi.data(), N, &info); // use c++11 ! vector::data() !
    //DEBUG_PRINT2("Poisson1D::solve() : DGTSV stop") ;
    if(info==0) {
        solved=true;
    }
    else {
        if(info < 0)
            throw std::logic_error(
                    "dgtsv: the i-th argument had an illegal value, i=" + std::to_string(info) );
        if(info > 0)
            throw std::logic_error(
                    "dgtsv: U(i,i) is exactly zero, and the solution has not been computed.  The factorization has not been completed unless i = N. i="
                    + std::to_string(info) );
    }
}


////////////////////////////////////////////////////////////////////////////////


} // end of namespace mw
