// File Poisson1D.h
//
// DO NOT edit outside SYNC/Sci/src/Poisson/
//
// Copyright 2014 Maciej Woloszyn

#ifndef POISSON1D_H
#define POISSON1D_H

#include "mwTools.h"
#include "Poisson.h"

namespace mw {


/**
 * @brief 1D Poisson solver
 *
 * Simple 1D finite difference solver using
 * the Dirichlet (or first-type) boundary condition.
 *
 * @todo TODO replace Poisson1D::N with Poisson::Nx
 *
 * See notes in MW/QM1/59
 */
class Poisson1D : public Poisson
{
public:
    Poisson1D(const v1D& rho) :
        N((int)rho.size()),
        phi(rho)
    {
        if(rho.size()<3)
            throw std::logic_error( "Poisson1D() : rho-array size must be > 2" );
    }

    Poisson1D(const v1D& rho, const v1D& vepsilon) :
        N((int)rho.size()),
        phi(rho),
        epsilon(vepsilon)
    {
        if(rho.size()<3)
            throw std::logic_error( "Poisson1D() : rho-array size must be > 2" );
        if(epsilon.size() != rho.size())
            throw std::logic_error( "Poisson1D() : epsilon-array and rho-array sizes not equal" );
        const_epsilon = false;
    }

    /// @brief solve eq.
    void solve();

    /// @brief phi(x_i)
    /// @todo TODO rename
    inline double phiAt(size_t i) const {
        return phi[i];
    }

//    inline double ksiAt(size_t i) const {
//        //return 0.25*(epsilon[i+1]-epsilon[i-1])/epsilon[i]; // uses 2nd order eps'(x)
//        return 0.5*(epsilon[i+1]-epsilon[i])/epsilon[i]; // uses 1st order eps'(x)
//    }

    double phiL = 0.; ///< boundary condition fi(xMin)
    double phiR = 0.; ///< boundary condition fi(xMax)

//protected:
    const int N; ///< rho.size()
    v1D phi; ///< electrostatic potential
    v1D epsilon; ///< dielectric constant epsilon(x)
    v1D DL, DU, D; ///< sub-, super-, and diagonal elements
    int info = 0;  ///< return code from LAPACK's routines

//    double dx = 1.0; // already in Poisson.h

//    bool epsilon_const_in_ranges = false;

}; // end of class Poisson1D

} // end of namespace mw

#endif // POISSON1D_H
