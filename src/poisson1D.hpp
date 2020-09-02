#ifndef POISSON1D_HPP
#define POISSON1D_HPP

#include "lib.hpp"

using namespace arma;

namespace poisson
{

class Poisson1D {

  public:

    Poisson1D(size_t nx, double h) : nx_ (nx), h_ (h) {
      rho_ = array<double>(nx_);
      nE_ = array<double>(nx_);
      uOld_ = array<double>(nx_);
      uNew_ = array<double>(nx_);
      dPu_ = sp_mat(nx_, nx_);
      pFun_ = vec(nx_);
      epsilonR_ = 1, temp_ = 300, beta_ = 1;
      dirichletL_ = 0, dirichletR_ = 0;
      neumannL_ = 0, neumannR_ = 0;
    }
    ~Poisson1D(){};

    void solve();
    void solve_gummel();
    void solve_tridiag();

    void testPoisson();

    size_t nx_;
    double h_, l_;
    double dirichletL_, dirichletR_;
    double neumannL_, neumannR_;
    double epsilonR_, temp_;
    double beta_;  // potential mixing

    array<double> rho_;
    array<double> nE_;
    array<double> uOld_;
    array<double> uNew_;

    vec pFun_;
    sp_mat dPu_;

};


}


# endif
