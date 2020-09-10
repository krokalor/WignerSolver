#include "lib.hpp"
#include "poisson1D.hpp"

using namespace poisson;

void Poisson1D::solve() { solve_gummel (); }  // solve_gummel solve_tridiag

void Poisson1D::solve_gummel() {
	// Solving Poisson equation using Gummel algorithm

	double epsilon = epsilonR_/4./M_PI;
	double c = h_*h_/epsilon;

	// P_i

	// TODO: add von Neumann b.c.
	double phi_L = -uOld_(0)+dirichletL_+h_*neumannL_;
	double phi_R = -uOld_(nx_-1)+dirichletR_-h_*neumannR_;

	array<double> uD(nx_);
	for (size_t i=1; i<nx_-1; ++i)
		uD(i) = uOld_(i+1)-2.*uOld_(i)+uOld_(i-1);
	uD(0) = phi_L, uD(nx_-1) = phi_R;

	pFun_ = vec(nx_);
	for (size_t i=0; i<nx_; ++i)
		pFun_(i) = -(uD(i) - c*rho_(i));  // '-' because in equation -pFun

	// dP_i/du_j

	dPu_ = sp_mat(nx_, nx_);
	for (size_t i=0; i<nx_; ++i) {
		for (size_t j=0; j<nx_; ++j) {
			if (i-1 == j || i+1 == j)
				dPu_(i, j) = 1;
			else if (j == i)
				dPu_(i, j) = -2;  // - c*nE_(i)/temp_;
		}
	}

	/*
	for (size_t i=0; i<nx_; ++i) {
		for (size_t j=0; j<nx_; ++j) {
			cout<<dPu_(i, j)<<' ';  // + c*nE_(i)/temp_;
		}
		cout<<endl;
	}
	*/

	vec x(nx_);
	superlu_opts settings;
	settings.symmetric = true;
	//    settings.refine = superlu_opts::REF_DOUBLE;
	spsolve(x, dPu_, pFun_, "superlu", settings);

	for (size_t i=0; i<nx_; ++i)
		uNew_(i) = uOld_(i) + beta_*x(i);

}


void Poisson1D::solve_tridiag() {
	// Solving Poisson equation using tridiagonal matrix

	double epsilon = epsilonR_/4./M_PI;

	array<double> b(nx_, -2.);  // Diagonal
	array<double> c(nx_, 1.);  // Upper diagonal
	array<double> a(nx_, 1.);  // Lower diagonal

	array<double> d(nx_), x(nx_);  // A*x = d

	for (size_t i=0; i<nx_; ++i)
		d(i) = h_*h_/epsilon*rho_(i);

	// Dirichlet BC
	d(0) -= dirichletL_, d(nx_-1) -= dirichletR_;
	// von Neumann BC
	d(0) -= h_*neumannL_, d(nx_-1) += h_*neumannR_;
	b(0) = -1, b(nx_-1) = -1;
	a(nx_-1) = 0, c(0) = 0;

	c(0) /= b(0);
	d(0) /= b(0);

	double m;
	for (size_t i=1; i<nx_; ++i) {
		m = 1./(b(i) - a(i)*c(i-1));
		c(i) *= m;
		d(i) = (d(i) - a(i)*d(i-1)) * m;
	}

	x(nx_-1) = d(nx_-1);
	for (size_t i=nx_-2; i>0; --i)
		x(i) = d(i) - c(i)*x(i+1);
	x(0) = d(0) - c(0)*x(1);

	for (size_t i=0; i<nx_; ++i)  // mixing old and new potential
		uNew_(i) = (1-beta_)*uOld_(i) + beta_*x(i);

}


void Poisson1D::testPoisson() {

	double sigma = 1e-3; // 10^-3 C/m^2

	rho_(nx_/2) = -sigma/(h_*AU_nm*1e-9);
	rho_(nx_/2) *= 1e-6*AU_cm3/E0;

	dirichletL_ = -5.64717/AU_eV;
	dirichletR_ = -5.64717/AU_eV;
	epsilonR_ = 1.;  // 8.854e-12;

	solve();

	cout<<"# sigma = "<<sigma<<" dirichletL "<<dirichletL_<<" dirichletR "<<dirichletR_<<" n "<<nx_<<" h "<<h_<<endl;
	for (size_t i = 0; i < nx_; ++i)
		cout<<i*h_*AU_nm<<' '<<uNew_(i)*AU_eV<<' '<<rho_(i)*E0/AU_cm3<<endl;

	// Run example:
	// Poisson1D p(200, 1./AU_nm);
	// p.testPoisson();

}
