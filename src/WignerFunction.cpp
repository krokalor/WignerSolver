#include "lib.hpp"
#include "WignerFunction.hpp"

// #include "Eigen/Dense"
// using namespace Eigen;

using namespace wigner;


// ############################## Equilibrium WF ##############################

void WignerFunction::initEq() {
	setBoundCond();
	// setEquilibriumFunction();
	u_ = uB_ + uC_;
	du_ = calcFirstDer(u_, dx_);
	d3u_ = calcThirdDer(u_, dx_);
	a_.zeros(), b_.zeros();
}

// ############################## Numerical WF solution ##############################

void WignerFunction::solveWignerEq(){

	initEq();
	// readPotential();

	// setBoundCond();
	// setEquilibriumFunction();

	// TODO: OpenMP parallel calculations
	// #pragma omp parallel for collapse(2) shared(a_, b_)
	for (size_t i=0; i<nx_; ++i) {
		for (size_t j=0; j<nk_; ++j) {
			diffusionTerm(i, j, -1);
			driftTerm(i, j, -1);
			scatteringTerm(i, j, -1);
			if (useQC_) quantumCorrTerm(i, j, -1);
		}  // end j loop
	}  // end i loop

	// cout<<"# min(A) = "<<min(min(a_))<<endl;
	// cout<<"# max(A) = "<<max(max(a_))<<endl;

	// Setting up solver options
	arma::superlu_opts opts;
	opts.symmetric = false;
	opts.equilibrate = true;
	opts.permutation = arma::superlu_opts::NATURAL;
	opts.refine = arma::superlu_opts::REF_EXTRA;  // 	iterative refinement in extra precision
	opts.allow_ugly  = true;
	opts.pivot_thresh = 0.01;

	arma::vec x(nxk_, arma::fill::zeros);
	arma::spsolve(x, a_, b_, "superlu", opts);  // use SuperLU solver

	//
	// EIGEN
	//
	// MatrixXd a((int)nxk_,(int)nxk_);
	// for (size_t i=0; i<nxk_; ++i)
	// 	for (size_t j=0; j<nxk_; ++j)
	// 		a(i,j) = a_(i,j);
	// VectorXd b((int)nxk_);
	// for (size_t i=0; i<nxk_; ++i)
	// 	b(i) = b_(i);
	// VectorXd x =  a.lu().solve(b);
	//
	// Matrix3d A;
	// Vector3d B;
	// A << 1,2,3,  4,5,6,  7,8,10;
	// B << 3, 3, 4;
	// cout << "Here is the matrix A:\n" << A << endl;
	// cout << "Here is the vector b:\n" << B << endl;
	// Vector3d X = A.lu().solve(B);
	// cout << "The solution is:\n" << X << endl;

	f_.zeros();
	for (size_t i=0; i<nx_; ++i)
		for (size_t j=0; j<nk_; ++j)
			f_(i,j) = x(i*nk_+j);

	/*
	// Checking A matrix
	for (i=0; i<nxk_; ++i) {
		for (j=0; j<nxk_; ++j) {
			char ch = a_(i, j) != 0 ? '0' : '-';
			cout<<ch<<' ';
		}
		cout<<endl;
	}
	*/

}


// ############################## WF time evolution ##############################
void WignerFunction::solveTimeEv(){

	initEq();

	// #pragma omp parallel for collapse(2)
	for (size_t i=0; i<nx_; ++i) {
		for (size_t j=0; j<nk_; ++j) {
			a_(i*nk_+j, i*nk_+j) += 1.;
			b_(i*nk_+j) = 2.*f_(i, j);
			diffusionTerm(i, j, dt_);
			driftTerm(i, j, dt_);
			scatteringTerm(i, j, dt_);
			if (useQC_) quantumCorrTerm(i, j, dt_);
		}  // end j loop
	}  // end i loop

	arma::superlu_opts opts;
	opts.symmetric = false;
	opts.equilibrate = true;
	// opts.refine = superlu_opts::REF_NONE;
	// opts.refine = superlu_opts::REF_EXTRA;  // 	iterative refinement in extra precision
    // opts.allow_ugly  = true;
	arma::vec x(nxk_, arma::fill::zeros);
	arma::spsolve(x, a_, b_, "superlu", opts);  // use SuperLU solver
	// spsolve(x, a, b, "superlu");

	for (size_t i=0; i<nx_; ++i)
		for (size_t j=0; j<nk_; ++j)
			f_(i,j) = x(i*nk_+j) - f_(i,j);
}


// #################### diffusionTerm ####################
void WignerFunction::diffusionTerm(size_t i, size_t j, double dt){
	// Fills Boltzmann equation matrix with diffusion term
	size_t r = i*nk_ + j;
	if (diffSch_K_ == "UDS1") {
		// UDS1
		double C = k_(j)/m_/dx_;
		double B = bc_(j)*C;
		if (dt > 0) C *= dt/2., B *= dt;
		if (k_(j)<0.) {
			if (i==nx_-1) {
				a_(r, r) += -C;
				b_(r) += -bc_(j)*C;
			}
			else {
				a_(r, r) += -C;
				a_(r, (i+1)*nk_+j) += C;
			}
		}
		if (k_(j)>0.) {
			if (i==0) {
				a_(r, r) += C;
				b_(r) += bc_(j)*C;
			}
			else {
				a_(r, r) += C;
				a_(r, (i-1)*nk_+j) += -C;
			}
		}
	}
	else if (diffSch_K_ == "UDS2") {
		// UDS1
		double C = k_(j)/m_/dx_/2.;
		double B = bc_(j)*C;
		if (dt > 0) C *= dt/2., B *= dt;
		if (k_(j)<0.) {
			if (i==nx_-1) {
				a_(r, r) += -3*C;
				b_(r) += -3*B;
			}
			else if (i==nx_-2) {
				a_(r, r) += -3*C;
				a_(r, (i+1)*nk_+j) += 4*C;
				b_(r) += B;
			}
			else {
				a_(r, r) += -3*C;
				a_(r, (i+1)*nk_+j) += 4*C;
				a_(r, (i+2)*nk_+j) += -C;
			}
		}
		else if (k_(j)>0.) {
			if (i==0) {
				a_(r, r) += 3*C;
				b_(r) += 3*B;
			}
			else if (i==1) {
				a_(r, r) += 3*C;
				a_(r, (i-1)*nk_+j) += -4*C;
				b_(r) += -B;
			}
			else {
				a_(r, r) += 3*C;
				a_(r, (i-1)*nk_+j) += -4*C;
				a_(r, (i-2)*nk_+j) += C;
			}
		}
	}
	else if (diffSch_K_ == "UDS3") {
		// UDS1
		double C = k_(j)/m_/dx_/6.;
		double B = bc_(j)*C;
		if (dt > 0) C *= dt/2., B *= dt;
		if (k_(j)<0.) {
			if (i==nx_-1) {
				a_(r, r) += -11*C;
				b_(r) += -11*B;
			}
			else if (i==nx_-2) {
				a_(r, r) += -11*C;
				a_(r, (i+1)*nk_+j) += 18*C;
				b_(r) += 7*B;
			}
			else if (i==nx_-3) {
				a_(r, r) += -11*C;
				a_(r, (i+1)*nk_+j) += 18*C;
				a_(r, (i+2)*nk_+j) += -9*C;
				b_(r) += -2*B;
			}
			else {
				a_(r, r) += -11*C;
				a_(r, (i+1)*nk_+j) += 18*C;
				a_(r, (i+2)*nk_+j) += -9*C;
				a_(r, (i+3)*nk_+j) += 2*C;
			}
		}
		else if (k_(j)>0.) {
			if (i==0) {
				a_(r, r) += 11*C;
				b_(r) += 11*B;
			}
			else if (i==1) {
				a_(r, r) += 11*C;
				a_(r, (i-1)*nk_+j) += -18*C;
				b_(r) += -7*B;
			}
			else if (i==2) {
				a_(r, r) += 11*C;
				a_(r, (i-1)*nk_+j) += -18*C;
				a_(r, (i-2)*nk_+j) += 9*C;
				b_(r) += 2*B;
			}
			else {
				a_(r, r) += 11*C;
				a_(r, (i-1)*nk_+j) += -18*C;
				a_(r, (i-2)*nk_+j) += 9*C;
				a_(r, (i-3)*nk_+j) += -2*C;
			}
		}
	}
	else if (diffSch_K_ == "UDS4") {
		// UDS4
		double C = k_(j)/m_/dx_/12.;
		double B = bc_(j)*C;
		if (dt > 0) C *= dt/2., B *= dt;
		if (k_(j)<0.) {
			if (i==nx_-1) {
				a_(r, r) += -25*C;
				b_(r) += -25*B;
			}
			else if (i==nx_-2) {
				a_(r, r) += -25*C;
				a_(r, (i+1)*nk_+j) += 48*C;
				b_(r) += 23*B;
			}
			else if (i==nx_-3) {
				a_(r, r) += -25*C;
				a_(r, (i+1)*nk_+j) += 48*C;
				a_(r, (i+2)*nk_+j) += -36*C;
				b_(r) += -13*B;
			}
			else if (i==nx_-4) {
				a_(r, r) += -25*C;
				a_(r, (i+1)*nk_+j) += 48*C;
				a_(r, (i+2)*nk_+j) += -36*C;
				a_(r, (i+3)*nk_+j) += 16*C;
				b_(r) += 3*B;
			}
			else {
				a_(r, r) += -25*C;
				a_(r, (i+1)*nk_+j) += 48*C;
				a_(r, (i+2)*nk_+j) += -36*C;
				a_(r, (i+3)*nk_+j) += 16*C;
				a_(r, (i+4)*nk_+j) += -3*C;
			}
		}
		else if (k_(j)>0.) {
			if (i==0) {
				a_(r, r) += 25*C;
				b_(r) += 25*C;
			}
			else if (i==1) {
				a_(r, r) += 25*C;
				a_(r, (i-1)*nk_+j) += -48*C;
				b_(r) += -23*B;
			}
			else if (i==2) {
				a_(r, r) += 25*C;
				a_(r, (i-1)*nk_+j) += -48*C;
				a_(r, (i-2)*nk_+j) += 36*C;
				b_(r) += 13*B;
			}
			else if (i==3) {
				a_(r, r) += 25*C;
				a_(r, (i-1)*nk_+j) += -48*C;
				a_(r, (i-2)*nk_+j) += 36*C;
				a_(r, (i-3)*nk_+j) += -16*C;
				b_(r) += -3*B;
			}
			else {
				a_(r, r) += 25*C;
				a_(r, (i-1)*nk_+j) += -48*C;
				a_(r, (i-2)*nk_+j) += 36*C;
				a_(r, (i-3)*nk_+j) += -16*C;
				a_(r, (i-4)*nk_+j) += 3*C;
			}
		}
	}
	else if (diffSch_K_ == "HDS22") {
		// HDS22
		double alpha = 2., beta = 1.;
		double C = k_(j)/m_/dx_/2.;
		double D = C/(alpha+beta);
		double B = bc_(j)*D;
		if (dt > 0) C *= dt/2., D *= dt/2., B *= dt;
		if (k_(j)<0.) {
			if (i==0) {  // UDS2 is used at outgoing boundary
				// b_(r) += B*alpha;
				// a_(r, r) += -3.*beta*D;
				// a_(r, (i+1)*nk_ + j) += (alpha+4.*beta)*D;
				// a_(r, (i+2)*nk_ + j) += -beta*D;
				a_(r, r) += -3.*C;
				a_(r, (i+1)*nk_ + j) += 4.*C;
				a_(r, (i+2)*nk_ + j) += -C;
			}
			else if (i==nx_-1) {
				a_(r, (i-1)*nk_ + j) += -alpha*D;
				a_(r, r) += -3.*beta*D;
				b_(r) += -B*(alpha+3.*beta);
				// a_(r, r) += -3.*C;
				// b_(r) += -3.*B;
			}
			else if (i==nx_-2) {
				a_(r, (i-1)*nk_ + j) += -alpha*D;
				a_(r, r) += -3.*beta*D;
				a_(r, (i+1)*nk_ + j) += (alpha+4.*beta)*D;
				b_(r) += B*beta;
				// a_(r, r) += -3.*C;
				// a_(r, (i+1)*nk_ + j) += 4.*C;
				// b_(r) += B;
			}
			else {
				a_(r, (i-1)*nk_ + j) += -alpha*D;
				a_(r, r) += -3.*beta*D;
				a_(r, (i+1)*nk_ + j) += (alpha+4.*beta)*D;
				a_(r, (i+2)*nk_ + j) += -beta*D;
			}
		}
		if (k_(j)>0.) {
			if (i==nx_-1) {  // UDS2 is used at outgoing boundary
				// b_(r) += -alpha*B;
				// a_(r, r) += 3.*beta*D;
				// a_(r, (i-1)*nk_ + j) += -(alpha+4.*beta)*D;
				// a_(r, (i-2)*nk_ + j) += beta*D;
				a_(r, r) += 3.*C;
				a_(r, (i-1)*nk_ + j) += -4.*C;
				a_(r, (i-2)*nk_ + j) += C;
			}
			else if (i==0) {
				a_(r, (i+1)*nk_ + j) += alpha*D;
				a_(r, r) += 3.*beta*D;
				b_(r) += (alpha+3.*beta)*B;
				// a_(r, r) += 3.*C;
				// b_(r) += 3.*B;
			}
			else if (i==1) {
				a_(r, (i+1)*nk_ + j) += alpha*D;
				a_(r, r) += 3.*beta*D;
				a_(r, (i-1)*nk_ + j) += -(alpha+4.*beta)*D;
				b_(r) += -beta*B;
				// a_(r, r) += 3.*C;
				// a_(r, (i-1)*nk_ + j) += -4.*C;
				// b_(r) += -B;
			}
			else {
				a_(r, (i+1)*nk_ + j) += alpha*D;
				a_(r, r) += 3.*beta*D;
				a_(r, (i-1)*nk_ + j) += -(alpha+4.*beta)*D;
				a_(r, (i-2)*nk_ + j) += beta*D;
			}
		}
	}
	else {
		cout<<"ERROR: WRONG DIFFERENTIATION SCHEME IN DIFFUSION TERM"
		"PLEASE CHOOSE: UDS1, UDS2 OR HDS22"<<endl;
		exit(0);
	} // End DS ifs
}


// #################### Drift term ####################
void WignerFunction::driftTerm(size_t i, size_t j, double dt){
	// Fills Boltzmann equation matrix with drift terms
	size_t r = i*nk_ + j;
	if (useNLP_) {
		// Non-local potential
		// TODO: Armadillo FFT
		size_t v;
		double sum, u1, u2;
		double C;
		// #pragma omp parallel for
		for (size_t l=0; l<nk_; l++){  // TODO: 0 -> 1 ?
			v = i*nk_+l;
			if (v <= r) {
				sum = 0;
				for (size_t q=0; q<nk2_; q++){
					if (i+q > nx_-1)
						u1 = u_(nx_-1);
					else
						u1 = u_(i+q);
					if ( int(i-q) < 0)
						u2 = u_(0);
					else
						u2 = u_(i-q);
					sum += sin_(j,l*nk2_+q) * (u1 - u2);  // sin(2*M_PI/nk_*q*(j-l))
				}  // end q loop
				C = 2./float(nk_) * sum;
				if (v < r){  // To avoid taking diagonal term twice
					a_(r, v) += C;
					a_(v, r) += -C;
				}
				else
					a_(r, r) += C;
			}  // end if
		}  // end l loop
	}  // end if
	else {
		// Classical force
		if (diffSch_P_ == "UDS1") {
			// UDS1
			double F = -du_(i);  // klasyczna siła równa -du/dx
			double C = F/dk_;
			if (dt > 0) C *= dt/2.;
			if (F <= 0) {
				if (j == nk_-1) {
					a_(r, r) += -C;
					b_(r) += -fermiDirac(kmax_)*C;
				}
				else{
					a_(r, r) += -C;
					a_(r, r+1) += C;
				}
			}
			else if (F > 0) {
				if (j == 0) {
					a_(r, r) += C;
					b_(r) += fermiDirac(-kmax_)*C;
				}
				else {
					a_(r, r) += C;
					a_(r, r-1) += -C;
				}
			}
		}
		else if (diffSch_P_ == "UDS2") {
			// UDS2
			double F = -du_(i);  // klasyczna siła równa -du/dx
			double C = F/dk_/2.;
			if (dt > 0) C *= dt/2.;
			if (F <= 0) {
				if (j==nk_-1) {
					a_(r, r) += -3*C;
					b_(r) += -3*fermiDirac(kmax_)*C;
				}
				else if (j==nk_-2) {
					a_(r, r) += -3*C;
					a_(r, r+1) += 4*C;
					b_(r) += fermiDirac(kmax_)*C;
				}
				else {
					a_(r, r) += -3*C;
					a_(r, r+1) += 4*C;
					a_(r, r+2) += -C;
				}
			}
			else if (F > 0) {
				if (j==0) {
					a_(r, r) += 3*C;
					b_(r) += 3*fermiDirac(-kmax_)*C;
				}
				else if (j==1) {
					a_(r, r) += 3*C;
					a_(r, r-1) += -4*C;
					b_(r) += -fermiDirac(-kmax_)*C;
				}
				else {
					a_(r, r) += 3*C;
					a_(r, r-1) += -4*C;
					a_(r, r-2) += C;
				}
			}
		}
		else if (diffSch_P_ == "UDS3") {
			// UDS3
			double F = -du_(i);  // klasyczna siła równa -du/dx
			double C = F/dk_/6.;
			if (dt > 0) C *= dt/2.;
			if (F <= 0) {
				if (j==nk_-1) {
					a_(r, r) += -11*C;
					b_(r) += -11*fermiDirac(kmax_)*C;
				}
				else if (j==nk_-2) {
					a_(r, r) += -11*C;
					a_(r, r+1) += 18*C;
					b_(r) += 7*fermiDirac(kmax_)*C;
				}
				else if (j==nk_-3) {
					a_(r, r) += -11*C;
					a_(r, r+1) += 18*C;
					a_(r, r+2) += -9*C;
					b_(r) += -2*fermiDirac(kmax_)*C;
				}
				else {
					a_(r, r) += -11*C;
					a_(r, r+1) += 18*C;
					a_(r, r+2) += -9*C;
					a_(r, r+3) += 2*C;
				}
			}
			else if (F > 0) {
				if (j==0) {
					a_(r, r) += 11*C;
					b_(r) += 11*fermiDirac(-kmax_)*C;
				}
				else if (j==1) {
					a_(r, r) += 11*C;
					a_(r, r-1) += -18*C;
					b_(r) += -7*fermiDirac(-kmax_)*C;
				}
				else if (j==2) {
					a_(r, r) += 11*C;
					a_(r, r-1) += -18*C;
					a_(r, r-1) += 9*C;
					b_(r) += 2*fermiDirac(-kmax_)*C;
				}
				else {
					a_(r, r) += 11*C;
					a_(r, r-1) += -18*C;
					a_(r, r-2) += 9*C;
					a_(r, r-3) += -2*C;
				}
			}
		}
		else if (diffSch_P_ == "HDS22") {
			// HDS22
			double F = -du_(i);  // klasyczna siła równa -du/dx
			double C = F/dk_/2.;
			if (dt > 0) C *= dt/2.;
			double alpha = 2., beta = 1.;
			double D = C/(alpha+beta);
			if (F <= 0) {
		  		if (j==0) {
				    b_(r) += fermiDirac(kmax_)*alpha*D;
				    a_(r, r) += -3.*beta*D;
				    a_(r, r+1) += (alpha+4.*beta)*D;
				    a_(r, r+2) += -beta*D;
				    // a_(r, r) += -3.*C;
				    // a_(r, r+1) += 4.*C;
				    // a_(r, r+2) += -C;
		  		}
		  		else if (j==nk_-1) {
				    a_(r, r-1) += -alpha*D;
				    a_(r, r) += -3.*beta*D;
				    b_(r) += -fermiDirac(kmax_)*(alpha+3.*beta)*D;
				    // a_(r, r) += -3.*C;
				    // b_(r) += -3.*fermiDirac(kmax_)*C;
		  		}
		  		else if (j==nk_-2) {
				    a_(r, r-1) += -alpha*D;
				    a_(r, r) += -3.*beta*D;
				    a_(r, r+1) += (alpha+4.*beta)*D;
				    b_(r) += fermiDirac(kmax_)*beta*D;
				    // a_(r, r) += -3.*C;
				    // a_(r, r+1) += 4.*C;
				    // b_(r) += fermiDirac(kmax_)*C;
		  		}
		  		else {
				    a_(r, r-1) += -alpha*D;
				    a_(r, r) += -3.*beta*D;
				    a_(r, r+1) += (alpha+4.*beta)*D;
				    a_(r, r+2) += -beta*D;
		  		}
			}
			else if (F > 0) {
		  		if (j==nk_-1) {
				    b_(r) += -fermiDirac(-kmax_)*alpha*D;
				    a_(r, r) += 3.*beta*D;
				    a_(r, r-1) += -(alpha+4.*beta)*D;
				    a_(r, r-2) += beta*D;
				    // a_(r, r) += 3.*C;
				    // a_(r, r-1) += -4.*C;
				    // a_(r, r-2) += C;
		  		}
	  			else if (j==0) {
				    a_(r, r+1) += alpha*D;
				    a_(r, r) += 3.*beta*D;
				    b_(r) += fermiDirac(-kmax_)*(alpha+3.*beta)*D;
				    // a_(r, r) += 3.*C;
				    // b_(r) += 3.*fermiDirac(-kmax_)*C;
		  		}
		  		else if (j==1) {
				    a_(r, r+1) += alpha*D;
				    a_(r, r) += 3.*beta*D;
				    a_(r, r-1) += -(alpha+4.*beta)*D;
				    b_(r) += -fermiDirac(-kmax_)*beta*D;
				    // a_(r, r) += 3.*C;
				    // a_(r, r-1) += -4.*C;
				    // b_(r) += -fermiDirac(-kmax_)*C;
		  		}
		  		else {
				    a_(r, r+1) += alpha*D;
				    a_(r, r) += 3.*beta*D;
				    a_(r, r-1) += -(alpha+4.*beta)*D;
				    a_(r, r-2) += beta*D;
			  	}
			}
		}
		else {
			cout<<"ERROR: WRONG DIFFERENTIATION SCHEME IN DRIFT TERM"
			"PLEASE CHOOSE: UDS1, UDS2 OR HDS22"<<endl;
			exit(0);
		} // End DS ifs
	}
}

// #################### scatteringTerm ####################
inline void WignerFunction::scatteringTerm(size_t i, size_t j, double dt){
	size_t r = i*nk_ + j;
	double cR = rR_, cM = rM_, cL = lambda_/dk_/dk_;
	if (dt > 0) cR *= dt/2., cM *= dt/2., cL *= dt/2.;
	// #################### rR term ####################
	a_(r, r) += cR;
	b_(r) += fEq_(i, j)*cR;
	// #################### rM term ####################
	a_(r, r) += cM;
	a_(r, i*nk_+(nk_-j-1)) += -cM;
	// #################### Lambda term ####################
	a_(r, r) += 2*cL;
	if (j==0)
		a_(r, r+1) += -cL;
	else if (j==nk_-1)
		a_(r, r-1) += -cL;
	else{
		a_(r, r-1) += -cL;
		a_(r, r+1) += -cL;
	}
	// #################### gamma term ####################
	double C = -rF_;
	if (dt > 0) C *= dt/2.;
	double F = -du_(i);  // klasyczna siła równa -du/dx
	// UDS1
	// C *= k_(j)/dk_;
	// a_(r,r) += C;
	// if (F > 0){
	// 	if (j == 0)
	// 		a_(r, r) += C;
	// 	else{
	// 		a_(r, r) += C;
	// 		a_(r, r-1) += -C;
	// 	}
	// }
	// else if (F <= 0){
	// 	if (j == nk_-1)
	// 		a_(r, r) += -C;
	// 	else{
	// 		a_(r, r) += -C;
	// 		a_(r, r+1) += C;
	// 	}
	// }
	// HDS22
	C *= k_(j)/dk_;
	a_(r,r) += C;
	double alpha = 2., beta = 1.;
	double D = C/(alpha+beta);
	if (F <= 0) {
		if (j==0) {
			b_(r) += fermiDirac(kmax_)*alpha*D;
			a_(r, r) += -3.*beta*D;
			a_(r, r+1) += (alpha+4.*beta)*D;
			a_(r, r+2) += -beta*D;
			// a_(r, r) += -3.*C;
			// a_(r, r+1) += 4.*C;
			// a_(r, r+2) += -C;
		}
		else if (j==nk_-1) {
			a_(r, r-1) += -alpha*D;
			a_(r, r) += -3.*beta*D;
			b_(r) += -fermiDirac(kmax_)*(alpha+3.*beta)*D;
			// a_(r, r) += -3.*C;
			// b_(r) += -3.*fermiDirac(kmax_)*C;
		}
		else if (j==nk_-2) {
			a_(r, r-1) += -alpha*D;
			a_(r, r) += -3.*beta*D;
			a_(r, r+1) += (alpha+4.*beta)*D;
			b_(r) += fermiDirac(kmax_)*beta*D;
			// a_(r, r) += -3.*C;
			// a_(r, r+1) += 4.*C;
			// b_(r) += fermiDirac(kmax_)*C;
		}
		else {
			a_(r, r-1) += -alpha*D;
			a_(r, r) += -3.*beta*D;
			a_(r, r+1) += (alpha+4.*beta)*D;
			a_(r, r+2) += -beta*D;
		}
	}
	else if (F > 0) {
		if (j==nk_-1) {
			b_(r) += -fermiDirac(-kmax_)*alpha*D;
			a_(r, r) += 3.*beta*D;
			a_(r, r-1) += -(alpha+4.*beta)*D;
			a_(r, r-2) += beta*D;
			// a_(r, r) += 3.*C;
			// a_(r, r-1) += -4.*C;
			// a_(r, r-2) += C;
		}
		else if (j==0) {
			a_(r, r+1) += alpha*D;
			a_(r, r) += 3.*beta*D;
			b_(r) += fermiDirac(-kmax_)*(alpha+3.*beta)*D;
			// a_(r, r) += 3.*C;
			// b_(r) += 3.*fermiDirac(-kmax_)*C;
		}
		else if (j==1) {
			a_(r, r+1) += alpha*D;
			a_(r, r) += 3.*beta*D;
			a_(r, r-1) += -(alpha+4.*beta)*D;
			b_(r) += -fermiDirac(-kmax_)*beta*D;
			// a_(r, r) += 3.*C;
			// a_(r, r-1) += -4.*C;
			// b_(r) += -fermiDirac(-kmax_)*C;
		}
		else {
			a_(r, r+1) += alpha*D;
			a_(r, r) += 3.*beta*D;
			a_(r, r-1) += -(alpha+4.*beta)*D;
			a_(r, r-2) += beta*D;
		}
	}
}


// #################### Quantum correction term ####################
void WignerFunction::quantumCorrTerm(size_t i, size_t j, double dt){
	// Fills Boltzmann equation matrix with drift terms
	size_t r = i*nk_ + j;
	double C = d3u_(i)/dk_/dk_/dk_/24.;
	if (dt > 0) C *= dt/2.;
	// // UDS2
	// if (d3u_(i)<=0.) {
	// 	if (j==nk_-1) {
	// 		a_(r, r) += -C;
	// 	}
	// 	else if (j==nk_-2) {
	// 		a_(r, r) += -C;
	// 		a_(r, r+1) += 3*C;
	// 	}
	// 	else if (j==nk_-3) {
	// 		a_(r, r) += -C;
	// 		a_(r, r+1) += 3*C;
	// 		a_(r, r+2) += -3*C;
	// 	}
	// 	else {
	// 		a_(r, r) += -C;
	// 		a_(r, r+1) += 3*C;
	// 		a_(r, r+2) += -3*C;
	// 		a_(r, r+3) += C;
	// 	}
	// }
	// if (d3u_(i)>0.) {
	// 	if (j==0) {
	// 		a_(r, r) += C;
	// 	}
	// 	else if (j==1) {
	// 		a_(r, r) += C;
	// 		a_(r, r-1) += -3*C;
	// 	}
	// 	else if (j==2) {
	// 		a_(r, r) += C;
	// 		a_(r, r-1) += -3*C;
	// 		a_(r, r-2) += 3*C;
	// 	}
	// 	else {
	// 		a_(r, r) += C;
	// 		a_(r, r-1) += -3*C;
	// 		a_(r, r-2) += 3*C;
	// 		a_(r, r-3) += -C;
	// 	}
	// }
	// CDS1
	if (j==nk_-1)
	{
		a_(r, r-1) += C;
		a_(r, r-2) += -C/2.;
	}
	else if (j==nk_-2)
	{
		a_(r, r+1) += -C;
		a_(r, r-1) += C;
		a_(r, r-2) += -C/2.;
	}
	else if (j==0)
	{
		a_(r, r+2) += C/2.;
		a_(r, r+1) += -C;
	}
	else if (j==1)
	{
		a_(r, r+2) += C/2.;
		a_(r, r+1) += -C;
		a_(r, r-1) += C;
	}
	else
	{
		a_(r, r+2) += C/2.;
		a_(r, r+1) += -C;
		a_(r, r-1) += C;
		a_(r, r-2) += -C/2.;
	}
}
