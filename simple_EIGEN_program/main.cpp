#include <iostream>
#include <cmath>
#include "Eigen/Dense"
#include "Eigen/Sparse"

using namespace std;
using namespace Eigen;

typedef Eigen::SparseMatrix<double> SpMat;

#define AU_eV 27.211
#define KB 8.617E-5
#define AU_Acm2 2.364e14    // [A/cm**2]  au_A/au_nm/au_nm
#define AU_cm3 1.48e-25     // [cm**3]  (AU_nm*1e-7)**3

// Supply function (E(k))
inline double sf(double k){
	double uF = 0.00293827, T = 300, m = 0.067;
    double mu = k > 0 ? uF+0.1/AU_eV : uF;
	double c = m/M_PI*KB/AU_eV*T, ex = -(k*k/m/2.-mu)/(KB/AU_eV*T);	 // [au]
	// double c = 4*M_PI*m*KB/AU_eV*temp_, ex = -(k*k/m/2.-mu)/(KB/AU_eV*temp_);	 // [au]
	if (ex < 700)
		return c * log(exp(ex)+1);
	else{
		if (ex > 0)
			return c * ex;
		else return 0;
	}
}

void SolveWigner(int nx, int nk, double l, double k_max) {
    int nxk = nx*nk, nk2 = int(nk/2.), r = 0;
    double dx = l/float(nx), dk = 2*k_max/float(nk);
    double C = 0, k = 0, m = 0.067;
    SpMat a_sp(nxk,nxk);
    MatrixXd a(nxk,nxk), f(nx,nk);
    VectorXd b(nxk), x(nxk);
	for (int i=0; i<nx; ++i) {
		for (int j=0; j<nk; ++j) {
            // Diffusion term
            r = i*nk + j;
            k = dk*(j-(nk-1)*.5);
            C = k/m/dx;
            if (k<0.) {
            	if (i==nx-1) {
            		a(r, r) += -C;
            		b(r) += -sf(k)*C;
            	}
            	else {
            		a(r, r) += -C;
            		a(r, (i+1)*nk+j) += C;
            	}
            }
            if (k>0.) {
            	if (i==0) {
            		a(r, r) += C;
            		b(r) += sf(k)*C;
            	}
            	else {
            		a(r, r) += C;
            		a(r, (i-1)*nk+j) += -C;
            	}
            }
        } // end j loop
    } // end i loop
    for (int i=0; i<nx; ++i)
		for (int j=0; j<nk; ++j)
            if (a(i,j) != 0) a_sp.insert(i,j) = a(i,j);
    // SimplicialCholesky<SpMat> solver(a_sp);  // SimplicialCholesky
    // x = solver.solve(b);
    x =  a.lu().solve(b);
    for (int i=0; i<nx; ++i)
		for (int j=0; j<nk; ++j)
			f(i,j) = x(i*nk+j);
    // Current density
    VectorXd currD(nx);
    for (int i=1; i<nx-2; ++i) {
		for (int j=0; j<nk2; ++j){  // k < 0
            k = dk*(j-(nk-1)*.5);
			currD(i) += k*(3*f(i+1,j)-f(i+2,j))*dk/2./m;  // j(i+1/2)
		}
		for (size_t j=nk2; j<nk; ++j){  // k > 0
            k = dk*(j-(nk-1)*.5);
			currD(i) += k*(3*f(i,j)-f(i-1,j))*dk/2./m;  // j(i+1/2)
		}
		currD(i) /= 2.*M_PI;
	}
	currD(0) = currD(2);
	// currD_(1) = currD_(2);
	currD(nx-1) = currD(nx-3);
	currD(nx-2) = currD(nx-3);
    // Carrier concentration
    VectorXd cdX(nx);
    for (int i=0; i<nx; ++i) {
        for (int j=1; j<nk/2; ++j)
            cdX(i) += (f(i,2*j-2)+4*f(i,2*j-1)+f(i,2*j))*dk/3.;  // simpson
        cdX(i) /= 2.*M_PI;
    }
    cout<<"Current density\n"<<currD*AU_Acm2<<endl;
    cout<<"Carrier concentration\n"<<cdX/AU_cm3<<endl;
}

int main() {
    SolveWigner(50, 50, 756.144, 0.15);
    return 0;
}
