#ifndef POISSON1D_HPP
#define POISSON1D_HPP

#include "lib.hpp"

using namespace arma;

namespace poisson
{

class Poisson1D {

	public:

		Poisson1D(size_t nx, double h) : nx_ (nx), h_ (h) {
			rho_ = vec(nx_, fill::zeros);
			nE_ = vec(nx_, fill::zeros);
			uOld_ = vec(nx_, fill::zeros);
			uNew_ = vec(nx_, fill::zeros);
			du_ = vec(nx_, fill::zeros);
			dPu_ = sp_mat(nx_, nx_);
			pFun_ = vec(nx_, fill::zeros);
			epsilonR_ = 1, temp_ = 300, beta_ = 1;
			dirichletL_ = 0, dirichletR_ = 0;
		}
		~Poisson1D(){};

		void solve();
		void solve_gummel();
		void solve_tridiag();

		void testPoisson();

		size_t nx_;
		double h_, l_;
		double dirichletL_, dirichletR_;
		double epsilonR_, temp_;
		double beta_;  // potential mixing

		vec rho_;
		vec nE_;
		vec uOld_;
		vec uNew_;
		vec du_;

		vec pFun_;
		sp_mat dPu_;

};


}


# endif
