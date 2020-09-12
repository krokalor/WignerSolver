#include "lib.hpp"
#include "WignerFunction.hpp"

using namespace wigner;


// ############################## Equilibrium WF ##############################

void WignerFunction::initEq() {
	setBoundCond();
	// setEquilibriumFunction();
	du_ = calcDer(u_, dx_);
	a_ = sp_mat(nxk_, nxk_);
	b_ = vec(nxk_, fill::zeros);
}

// ############################## Numerical WF solution ##############################

void WignerFunction::solveWignerEq(){

	initEq();
	// readPotential();

	vec x(nxk_);
	// setBoundCond();
	// setEquilibriumFunction();

	// TODO: OpenMP parallel calculations
	// #pragma omp parallel for collapse(2) default(none)
	for (size_t i=0; i<nx_; ++i) {
		for (size_t j=0; j<nk_; ++j) {
			diffusionTerm(i, j, -1);
			driftTerm(i, j, -1);
			scatteringTerm(i, j, -1);
		}  // end j loop
	}  // end i loop

	superlu_opts settings;
	settings.symmetric = false;
//    settings.refine = superlu_opts::REF_DOUBLE;
	spsolve(x, a_, b_, "superlu", settings);
//     spsolve(x, a, b, "superlu");

	f_.zero();
	for (size_t i=0; i<nx_; ++i)
		for (size_t j=0; j<nk_; ++j)
			f_(i,j) = x(i*nk_+j);

	/*
	std::ofstream file("wyniki/dane/a_mat.out");
	for (i=0; i<nxk_; ++i){
			for (j=0; j<nxk_; ++j){
					if (a(i,j) != 0)
							file<<0<<' ';
					else
							file<<'-'<<' ';
			}
			file<<endl;
	}
	file.close();
	*/

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

	vec x(nxk_);

	// setBoundCond();
	// setEquilibriumFunction();

	// #pragma omp parallel for collapse(2)
	for (size_t i=0; i<nx_; ++i) {
		for (size_t j=0; j<nk_; ++j) {
			a_(i*nk_+j, i*nk_+j) += 1.;
			b_(i*nk_+j) = 2.*f_(i, j);
			diffusionTerm(i, j, dt_);
			driftTerm(i, j, dt_);
			scatteringTerm(i, j, dt_);
		}  // end j loop
	}  // end i loop

	superlu_opts opts;
	opts.symmetric = false;
//     opts.allow_ugly  = true;
	spsolve(x, a_, b_, "superlu", opts);
	// spsolve(x, a, b, "superlu");

	for (size_t i=0; i<nx_; ++i)
			for (size_t j=0; j<nk_; ++j)
					f_(i,j) = x(i*nk_+j) - f_(i,j);
}


// #################### diffusionTerm ####################
void WignerFunction::diffusionTerm(size_t i, size_t j, double dt){
	// Fills Boltzmann equation matrix with diffusion term
	size_t r = i*nk_ + j;
	double C = k_(j)/m_/dx_/2., b_ij = bc_(j)*C;
	if (dt > 0) C *= dt/2., b_ij *= dt;
	double alpha = 2., beta = 1.;
	double D = C/(alpha+beta);
	// HDS22
	if (k_(j)<0.) {
		if (i==0) {
			a_(r, r) += -3.*C;
			a_(r, (i+1)*nk_ + j) += 4.*C;
			a_(r, (i+2)*nk_ + j) += -C;
		}
		else if (i==nx_-1) {
			a_(r, r) += -3.*C;
			b_(r) += -3.*b_ij;
		}
		else if (i==nx_-2) {
			a_(r, r) += -3.*C;
			a_(r, (i+1)*nk_ + j) += 4.*C;
			b_(r) += b_ij;
		}
		else {
			a_(r, (i-1)*nk_ + j) += -alpha*D;
			a_(r, r) += -3.*beta*D;
			a_(r, (i+1)*nk_ + j) += (alpha+4.*beta)*D;
			a_(r, (i+2)*nk_ + j) += -beta*D;
		}
	}
	if (k_(j)>0.) {
		if (i==nx_-1) {
			a_(r, r) += 3.*C;
			a_(r, (i-1)*nk_ + j) += -4.*C;
			a_(r, (i-2)*nk_ + j) += C;
		}
		else if (i==0) {
			a_(r, r) += 3.*C;
			b_(r) += 3.*b_ij;
		}
		else if (i==1) {
			a_(r, r) += 3.*C;
			a_(r, (i-1)*nk_ + j) += -4.*C;
			b_(r) += -b_ij;
		}
		else {
			a_(r, (i+1)*nk_ + j) += alpha*D;
			a_(r, r) += 3.*beta*D;
			a_(r, (i-1)*nk_ + j) += -(alpha+4.*beta)*D;
			a_(r, (i-2)*nk_ + j) += beta*D;
		}
	}
}


// #################### Drift term ####################
void WignerFunction::driftTerm(size_t i, size_t j, double dt){
	// Fills Boltzmann equation matrix with drift terms
	size_t r = i*nk_ + j;
	if (useNLP_) {
		// Non-local potential
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
		double C = F/dk_/1.;
		if (dt > 0) C *= dt/2.;
		// HDS1
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
		    // b_(r) += -3.*b_ij;
		  }
		  else if (j==nk_-2) {
		    a_(r, r) += -3.*C;
		    a_(r, r+1) += 4.*C;
		    // b_(r) += b_ij;
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
		    // b_(r) += 3.*b_ij;
		  }
		  else if (j==1) {
		    a_(r, r) += 3.*C;
		    a_(r, r-1) += -4.*C;
		    // b_(r) += -b_ij;
		  }
		  else {
		    a_(r, r+1) += alpha*D;
		    a_(r, r) += 3.*beta*D;
		    a_(r, r-1) += -(alpha+4.*beta)*D;
		    a_(r, r-2) += beta*D;
		  }
		}
		/*
		// UDS1
		if (F > 0){
				if (j == 0)
						a_(r, r) += C;
				else{
						a_(r, r) += C;
						a_(r, r-1) += -C;
				}
		}
		else if (F <= 0){
				if (j == nk_-1)
						a_(r, r) += -C;
				else{
						a_(r, r) += -C;
						a_(r, r+1) += C;
				}
		}
		*/
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
		C *= k_(j)/dk_/1.;
		if (dt > 0) C *= dt/2.;
		// CDS1
		if (j==0) {
		  // a_(r, r) += -3.*C;
		  // a_(r, r+1) += 4.*C;
		  // a_(r, r+2) += -C;
		  a_(r, r+1) += C/2.;
		}
		else if (j==nk_-1) {
		  // a_(r, r) += 3.*C;
		  // a_(r, r-1) += -4.*C;
		  // a_(r, r-2) += C;
		  a_(r, r-1) += -C/2.;
		}
		else {
		  a_(r, r-1) += -C/2.;
		  a_(r, r+1) += C/2.;
		}
		/*
		double F = -du_(i);  // klasyczna siła równa -du/dx
		// UDS1
		if (F > 0){
				if (j == 0)
						a_(r, r) += C;
				else{
						a_(r, r) += C;
						a_(r, r-1) += -C;
				}
		}
		else if (F <= 0){
				if (j == nk_-1)
						a_(r, r) += -C;
				else{
						a_(r, r) += -C;
						a_(r, r+1) += C;
				}
		}
		*/
		/*
		// HDS1
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
		    // b_(r) += -3.*b_ij;
		  }
		  else if (j==nk_-2) {
		    a_(r, r) += -3.*C;
		    a_(r, r+1) += 4.*C;
		    // b_(r) += b_ij;
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
		    // b_(r) += 3.*b_ij;
		  }
		  else if (j==1) {
		    a_(r, r) += 3.*C;
		    a_(r, r-1) += -4.*C;
		    // b_(r) += -b_ij;
		  }
		  else {
		    a_(r, r+1) += alpha*D;
		    a_(r, r) += 3.*beta*D;
		    a_(r, r-1) += -(alpha+4.*beta)*D;
		    a_(r, r-2) += beta*D;
		  }
		}
		*/
}
