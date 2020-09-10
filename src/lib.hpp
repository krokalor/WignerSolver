#ifndef LIB_H
#define LIB_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <ctime>
#include <map>

#include <armadillo>
// #include<omp.h>

// using namespace std;
using std::cout;
using std::endl;

// ############################## TYPEDEFs ##############################
// TODO: put into a namespace
typedef std::vector<double>::size_type size_t_vec_d;

// ############################## GLOBAL VARIABLES ##############################
// #define ARMA_USE_SUPERLU 1
// #define N_THREADS 2

// ############################## CONSTANTS ##############################

double const PI {3.141592653589793238};
double const KB {8.617E-5};				// Boltzmann constant [eV/K]
double const E0 {1.602E-19};				// Elementary charge [C]
double const HBAR_J {1.05E-34};   			// Planck constant [Js]
double const HBAR_eV {6.55E-16};   			// Planck constant [eVs]
double const M0 {9.11E-31};    			// Electron mass [kg]
double const A_GaAs {0.565};				// GaAs lattice constant [nm]

// ############################## atomic units ##############################

double const AU_nm {0.0529};				// Bohr radius [nm] (0.7896 nm)
double const AU_eV {27.211};     		// Hartree energy [eV]
double const AU_s {2.42E-17};				// Time [ps = 10^-12 s] (HBAR_eV/AU_eV*1e12)
double const AU_A {6.62e-3};        // [A]  _e0/_tau0
double const AU_Acm2 {2.364e14};    // [A/cm**2]  au_A/au_nm/au_nm
double const AU_cm3 {1.48e-25};     // [cm**3]  (AU_nm*1e-7)**3
double const AU_cm2 {2.8e-17};      // [cm**2]  (au_nm*1e-7)**2
double const AU_cm {5.29e-9};      // [cm]  au_nm*1e-7

// #################### klasa array ####################
template <class T>
class array {
	private:
		std::vector<T> a;
		size_t n_row;
	public:
		array() {}
		array(size_t n) : a (std::vector<T>(n)), n_row (n) {}
		array(size_t n, T v) : a (std::vector<T>(n, v)), n_row (n) {}
		~array() { a.clear(); }
		T& operator () (size_t i) { return a.at(i); }
		// Vector addition
		array<T> operator-(array<T> ai) {
			array<T> af(n_row);
			for (size_t i=0; i<n_row; ++i) af.at(i) = a.at(i)-ai.at(i);
			return af;
		}
		// Vector substraction
		array<T> operator+(array<T> ai) {
			array<T> af(n_row);
			for (size_t i=0; i<n_row; ++i) af.at(i) = a.at(i)+ai.at(i);
			return af;
		}
		// Multiplying vector by scalar
		// array<T> operator*(Q s) {
		//   array<T> af(n_row);
		//   for (size_t i=0; i<n_row; ++i) af.at(i) = s*a.at(i);
		//   return af;
		// }
		T& at(size_t i) { return a.at(i); }
		size_t size() { return a.size(); }
		void add(T x) { a.push_back(x); }
		void copy(array<T> ac) { for (size_t i=0; i<n_row; ++i) a.at(i) = ac(i); }
		void zero() { for (size_t i=0; i<n_row; ++i) a.at(i) = 0.; }
};

// #################### klasa matrix ####################
template <class T>
class matrix {
	private:
		array< array<T> > m;
		size_t n_row, n_col;
	public:
		matrix() {}
		matrix(size_t n1, size_t n2) : n_row (n1), n_col (n2)
				{ for (size_t i=0; i<n_row; ++i) m.add( array<T>(n_col) ); }
		~matrix() {}
		array<T>& operator () (size_t i) { return m.at(i); }
		T& operator () (size_t i, size_t j) { return m.at(i).at(j); }
		array<T>& at(size_t i) { return m.at(i); }
		matrix<T> operator-(matrix<T> mi) {
			matrix<T> mf(n_row, n_col);
			for (size_t i=0; i<n_row; ++i)
				for (size_t j=0; j<n_col; ++j)
					mf.at(i) = m.at(i)-mi.at(i);
			return mf;
		}
		matrix<T> operator+(matrix<T> mi) {
			matrix<T> mf(n_row, n_col);
			for (size_t i=0; i<n_row; ++i)
				for (size_t j=0; j<n_col; ++j)
					mf.at(i) = m.at(i)+mi.at(i);
			return mf;
		}
		T& at(size_t i, size_t j) { return m.at(i).at(j); }
		size_t size() { return n_col*n_row; }
		void copy(matrix<T> mc) { for (size_t i=0; i<n_row; ++i) m(i).copy(mc(i)); }
		// void add(T x){ m.push_back(x); }
		void zero() { for (size_t i=0; i<n_row; ++i) m(i).zero(); }
};


// #################### funkcja do liczenia całki ####################
template <class T, class Q>
T calcInt(array<T> f, Q h){
	size_t n = f.size();
	T ig = 0;
	for (size_t i=1; i<n; ++i)
		ig += (f(i-1) + f(i))*h/2.;
	return ig;
}

// #################### funkcja do liczenia pochodnej ####################
template <class T, class Q>
array<T> calcDer(array<T> f, Q h){
	size_t n = f.size();
	array<T> df(n);
	for (size_t i=1; i<n-1; ++i)
		df(i) =(-f(i-1)+f(i+1))/h/2.;
	df(0) = (-3.*f(0)+4.*f(1)-f(2))/h/2.;
	df(n-1) = (3.*f(n-1)-4.*f(n-2)+f(n-3))/h/2.;
	return df;
}


#endif
