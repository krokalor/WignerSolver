#include "lib.hpp"
#include "WignerFunction.hpp"

using namespace wigner;


// ############################## Equilibrium WF ##############################

void WignerFunction::initEq() {
	setBoundCond();
	// setEquilibriumFunction();
	u_ = uB_ + uC_;
	du_ = calcDer(u_, dx_);
	dddu_ = calcThirdDer(u_, dx_);
	a_ = sp_mat(nxk_, nxk_);
	b_ = vec(nxk_, fill::zeros);
}

// ############################## Numerical WF solution ##############################

void WignerFunction::solveWignerEq(){

	initEq();
	// readPotential();

	vec x(nxk_, fill::zeros);
	// setBoundCond();
	// setEquilibriumFunction();

	// Setting up solver options
	superlu_opts opts;
	// opts.symmetric = false;
	// opts.equilibrate = true;
	// opts.refine = superlu_opts::REF_EXTRA;  // 	iterative refinement in extra precision
	// opts.allow_ugly  = true;

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

	spsolve(x, a_, b_, "superlu", opts);  // use SuperLU solver
	// spsolve(x, a, b, "superlu");

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

	vec x(nxk_, fill::zeros);

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

	superlu_opts opts;
	// opts.symmetric = false;
	// opts.equilibrate = true;
	// opts.refine = superlu_opts::REF_NONE;
	// opts.refine = superlu_opts::REF_EXTRA;  // 	iterative refinement in extra precision
    // opts.allow_ugly  = true;
	spsolve(x, a_, b_, "superlu", opts);  // use SuperLU solver
	// spsolve(x, a, b, "superlu");

	for (size_t i=0; i<nx_; ++i)
		for (size_t j=0; j<nk_; ++j)
			f_(i,j) = x(i*nk_+j) - f_(i,j);
}


// #################### diffusionTerm ####################
void WignerFunction::diffusionTerm(size_t i, size_t j, double dt){
	// Fills Boltzmann equation matrix with diffusion term
	size_t r = i*nk_ + j;
	double C = k_(j)/m_/dx_/2.;  // /2. for HDS22 and UDS2!!!
	double bj = bc_(j)*C;
	if (dt > 0) C *= dt/2., bj *= dt;
	// // HDS22
	// double alpha = 2., beta = 1.;
	// double D = C/(alpha+beta);
	// if (k_(j)<0.) {
	// 	if (i==0) {
	// 		a_(r, r) += -3.*C;
	// 		a_(r, (i+1)*nk_ + j) += 4.*C;
	// 		a_(r, (i+2)*nk_ + j) += -C;
	// 	}
	// 	else if (i==nx_-1) {
	// 		a_(r, r) += -3.*C;
	// 		b_(r) += -3.*bj;
	// 	}
	// 	else if (i==nx_-2) {
	// 		a_(r, r) += -3.*C;
	// 		a_(r, (i+1)*nk_ + j) += 4.*C;
	// 		b_(r) += bj;
	// 	}
	// 	else {
	// 		a_(r, (i-1)*nk_ + j) += -alpha*D;
	// 		a_(r, r) += -3.*beta*D;
	// 		a_(r, (i+1)*nk_ + j) += (alpha+4.*beta)*D;
	// 		a_(r, (i+2)*nk_ + j) += -beta*D;
	// 	}
	// }
	// if (k_(j)>0.) {
	// 	if (i==nx_-1) {
	// 		a_(r, r) += 3.*C;
	// 		a_(r, (i-1)*nk_ + j) += -4.*C;
	// 		a_(r, (i-2)*nk_ + j) += C;
	// 	}
	// 	else if (i==0) {
	// 		a_(r, r) += 3.*C;
	// 		b_(r) += 3.*bj;
	// 	}
	// 	else if (i==1) {
	// 		a_(r, r) += 3.*C;
	// 		a_(r, (i-1)*nk_ + j) += -4.*C;
	// 		b_(r) += -bj;
	// 	}
	// 	else {
	// 		a_(r, (i+1)*nk_ + j) += alpha*D;
	// 		a_(r, r) += 3.*beta*D;
	// 		a_(r, (i-1)*nk_ + j) += -(alpha+4.*beta)*D;
	// 		a_(r, (i-2)*nk_ + j) += beta*D;
	// 	}
	// }
	// UDS2
	if (k_(j)<0.) {
		if (i==nx_-1) {
			a_(r, r) += -3*C;
			b_(r) += -3*bj;
		}
		else if (i==nx_-2) {
			a_(r, r) += -3*C;
			a_(r, (i+1)*nk_+j) += 4*C;
			b_(r) += bj;
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
			b_(r) += 3*bj;
		}
		else if (i==1) {
			a_(r, r) += 3*C;
			a_(r, (i-1)*nk_+j) += -4*C;
			b_(r) += -bj;
		}
		else {
			a_(r, r) += 3*C;
			a_(r, (i-1)*nk_+j) += -4*C;
			a_(r, (i-2)*nk_+j) += C;
		}
	}
	// // UDS1
	// if (k_(j)<0.) {
	// 	if (i==nx_-1) {
	// 		a_(r, r) += -C;
	// 		b_(r) += -bc_(j)*C;
	// 	}
	// 	else {
	// 		a_(r, r) += -C;
	// 		a_(r, (i+1)*nk_+j) += C;
	// 	}
	// }
	// if (k_(j)>0.) {
	// 	if (i==0) {
	// 		a_(r, r) += C;
	// 		b_(r) += bc_(j)*C;
	// 	}
	// 	else {
	// 		a_(r, r) += C;
	// 		a_(r, (i-1)*nk_+j) += -C;
	// 	}
	// }
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
		double F = -du_(i);  // klasyczna siła równa -du/dx
		double C = F/dk_/2.;
		if (dt > 0) C *= dt/2.;
		// // HDS22
		// double alpha = 2., beta = 1.;
		// double D = C/(alpha+beta);
		// if (F <= 0.) {
		//   if (j==0) {
		//     a_(r, r) += -3.*C;
		//     a_(r, r+1) += 4.*C;
		//     a_(r, r+2) += -C;
		//   }
		//   else if (j==nk_-1) {
		//     a_(r, r) += -3.*C;
		//     // b_(r) += -3.*bj;
		//   }
		//   else if (j==nk_-2) {
		//     a_(r, r) += -3.*C;
		//     a_(r, r+1) += 4.*C;
		//     // b_(r) += bj;
		//   }
		//   else {
		//     a_(r, r-1) += -alpha*D;
		//     a_(r, r) += -3.*beta*D;
		//     a_(r, r+1) += (alpha+4.*beta)*D;
		//     a_(r, r+2) += -beta*D;
		//   }
		// }
		// if (F > 0.) {
		//   if (j==nk_-1) {
		//     a_(r, r) += 3.*C;
		//     a_(r, r-1) += -4.*C;
		//     a_(r, r-2) += C;
		//   }
		//   else if (j==0) {
		//     a_(r, r) += 3.*C;
		//     // b_(r) += 3.*bj;
		//   }
		//   else if (j==1) {
		//     a_(r, r) += 3.*C;
		//     a_(r, r-1) += -4.*C;
		//     // b_(r) += -bj;
		//   }
		//   else {
		//     a_(r, r+1) += alpha*D;
		//     a_(r, r) += 3.*beta*D;
		//     a_(r, r-1) += -(alpha+4.*beta)*D;
		//     a_(r, r-2) += beta*D;
		//   }
		// }
		// UDS2
		if (F<0.) {
			if (j==nk_-1) {
				a_(r, r) += -3*C;
				b_(r) += -3*fermiDirac(kmax_);
			}
			else if (j==nk_-2) {
				a_(r, r) += -3*C;
				a_(r, r+1) += 4*C;
				b_(r) += fermiDirac(kmax_);
			}
			else {
				a_(r, r) += -3*C;
				a_(r, r+1) += 4*C;
				a_(r, r+2) += -C;
			}
		}
		if (F>=0.) {
			if (j==0) {
				a_(r, r) += 3*C;
				b_(r) += 3*fermiDirac(-kmax_);
			}
			else if (j==1) {
				a_(r, r) += 3*C;
				a_(r, r-1) += -4*C;
				b_(r) += -fermiDirac(-kmax_);
			}
			else {
				a_(r, r) += 3*C;
				a_(r, r-1) += -4*C;
				a_(r, r-2) += C;
			}
		}
		// // UDS1
		// if (F <= 0) {
		// 	if (j == nk_-1) {
		// 		a_(r, r) += -C;
		// 		b_(r) += -fermiDirac(kmax_);
		// 	}
		// 	else{
		// 		a_(r, r) += -C;
		// 		a_(r, r+1) += C;
		// 	}
		// }
		// else if (F > 0) {
		// 	if (j == 0) {
		// 		a_(r, r) += C;
		// 		b_(r) += fermiDirac(-kmax_);
		// 	}
		// 	else {
		// 		a_(r, r) += C;
		// 		a_(r, r-1) += -C;
		// 	}
		// }
		// // CDS1
		// C /= 2;
		// if (j==nk_-1)
		// {
		// 	a_(r, r-1) += -C;
		// }
		// else if (j==0)
		// {
		// 	a_(r, r+1) += C;
		// }
		// else
		// {
		// 	a_(r, r+1) += C;
		// 	a_(r, r-1) += -C;
		// }
	}
}

// #################### scatteringTerm ####################
inline void WignerFunction::scatteringTerm(size_t i, size_t j, double dt){
		size_t r = i*nk_ + j;
		double cR = rR_, cM = rM_, cL = lambda_/dk_/dk_;
		if (dt > 0) cR *= dt/2., cM *= dt/2., cL *= dt/2.;
		// #################### rR term ####################
		a_(r, r) += cR;
		b_(r) += fe_(i, j)*cR;
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
		a_(r,r) += C;
		C *= k_(j)/dk_;
		// // CDS1
		// if (j==0) {
		//   // a_(r, r) += -3.*C;
		//   // a_(r, r+1) += 4.*C;
		//   // a_(r, r+2) += -C;
		//   a_(r, r+1) += C/2.;
		// }
		// else if (j==nk_-1) {
		//   // a_(r, r) += 3.*C;
		//   // a_(r, r-1) += -4.*C;
		//   // a_(r, r-2) += C;
		//   a_(r, r-1) += -C/2.;
		// }
		// else {
		//   a_(r, r-1) += -C/2.;
		//   a_(r, r+1) += C/2.;
		// }
		// double F = -du_(i);  // klasyczna siła równa -du/dx
		// // UDS1
		// if (F > 0){
		// 		if (j == 0)
		// 				a_(r, r) += C;
		// 		else{
		// 				a_(r, r) += C;
		// 				a_(r, r-1) += -C;
		// 		}
		// }
		// else if (F <= 0){
		// 		if (j == nk_-1)
		// 				a_(r, r) += -C;
		// 		else{
		// 				a_(r, r) += -C;
		// 				a_(r, r+1) += C;
		// 		}
		// }
		// HDS1
		double F = -du_(i);  // klasyczna siła równa -du/dx
		double alpha = 2., beta = 1.;
		double D = C/(alpha+beta);
		if (F < 0.) {
		  if (j==0) {
		    a_(r, r) += -3.*C;
		    a_(r, r+1) += 4.*C;
		    a_(r, r+2) += -C;
		  }
		  else if (j==nk_-1) {
		    a_(r, r) += -3.*C;
		    // b_(r) += -3.*bj;
		  }
		  else if (j==nk_-2) {
		    a_(r, r) += -3.*C;
		    a_(r, r+1) += 4.*C;
		    // b_(r) += bj;
		  }
		  else {
		    a_(r, r-1) += -alpha*D;
		    a_(r, r) += -3.*beta*D;
		    a_(r, r+1) += (alpha+4.*beta)*D;
		    a_(r, r+2) += -beta*D;
		  }
		}
		if (F > 0.) {
		  if (j==nk_-1) {
		    a_(r, r) += 3.*C;
		    a_(r, r-1) += -4.*C;
		    a_(r, r-2) += C;
		  }
		  else if (j==0) {
		    a_(r, r) += 3.*C;
		    // b_(r) += 3.*bj;
		  }
		  else if (j==1) {
		    a_(r, r) += 3.*C;
		    a_(r, r-1) += -4.*C;
		    // b_(r) += -bj;
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
	double C = dddu_(i)/dk_/dk_/dk_/24.;
	if (dt > 0) C *= dt/2.;
	// // UDS2
	// if (dddu_(i)<=0.) {
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
	// if (dddu_(i)>0.) {
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
