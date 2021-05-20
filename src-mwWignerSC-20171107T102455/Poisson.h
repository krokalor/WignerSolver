// File Poisson.h
//
// DO NOT edit outside SYNC/Sci/src/Poisson/
//
// Copyright 2014 Maciej Woloszyn


#ifndef MW_POISSON_H
#define MW_POISSON_H

#include <stdexcept>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>
//#include "mwTools.h"

#include <limits> // temp

//#define DEBUG

namespace mw {

/**
 * @brief Base class for Poisson equation solvers
 *
 * All values in atomic units [a.u.]!
 *
 * @todo TODO solve(), init() etc as virtual functions
 * @todo TODO multigrid -->
 http://www.earth.lsa.umich.edu/~keken/534/poisson/poisson.pdf
 http://www.physics.buffalo.edu/phy410-505/2011/topic3/index.html    (src files!!!)
 https://nanohub.org/resources/9130/download/Multigrid_method.pdf
 http://cs.ucsb.edu/~koc/docs/c13.pdf
 http://vuduc.org/teaching/cse8803-pna-sp08/slides/cse8803-pna-sp08-12.pdf
 http://hod4.net/~hod/papers/SCF/Roberts.pdf
 http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.39.3428
 * @todo TODO FEM
 *
 * @author Maciej Woloszyn
 */
class Poisson
{
public:
    Poisson(size_t Nx = 0, size_t Ny = 0, size_t Nz = 0)
        : Nx(Nx), Ny(Ny), Nz(Nz)  { }

    const size_t Nx = 0; ///< mesh size along the x-axis
    const size_t Ny = 0; ///< mesh size along the y-axis, used only for 2D and 3D
    const size_t Nz = 0; ///< mesh size along the z-axis, used only for 3D

    double dx = 1.0; ///< dx=dy=dz in the current version
    double dx2e = dx*dx / epsilonc; ///< dx^2/epsilon
    double dx2 = dx*dx;
    bool use_openmp = false; ///< OpenMP usage allowed (default false)
    bool use_sor_red_black = false; ///< use SOR red-black scheme (default false)
    bool use_sor_omega_auto = true; ///< use SOR's omega automatic choice (default true)
    bool use_multigrid = true; ///< use simple multigrid (default true)
    bool const_epsilon = true; ///< use constant el. permittivity (default true)
    bool solved = false;

    double epsilonc = 0.25/M_PI; ///< absolute permittivity [au] (defualt: vacuum permittivity = 0.25/pi [au])

    double sor_omega = 2./(1.+2*M_PI/double(Nx+Ny)); ///< relaxation param. omega for SOR
    size_t sor_max_iter = 10000; ///< max number of SOR iterations
    double sor_tol = 1e-3; ///< max residual value
    double sor_action_tol = 1e-9;  ///< max decrease of a functional

protected:
    double sor_residual = 0.;

}; // class Poisson

////////////////////////////////////////////////////////////////////////////////

} // end of namespace mw

#endif // MW_POISSON_H
