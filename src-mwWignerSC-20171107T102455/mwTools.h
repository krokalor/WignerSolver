// DO NOT edit outside SYNC/Sci/src/mwTools/ !

// Copyright 2008-2015 Maciej Woloszyn <woloszyn@agh.edu.pl>


#pragma once

#ifndef MW_TOOLS_H
#define MW_TOOLS_H

// TODO tests for EACH function!!!

#define MW_MAX_DATE_LEN 22
#define MW_CSTR_TMP_LEN 200

#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <complex>

#ifndef DEBUG
  #define DEBUG_PRINT(x)
  #define DEBUG_PRINT2(x)
  #define DEBUG_RUN(x)
#else
  #define DEBUG_PRINT(x) \
    std::cerr <<  "#> " << __FILE__ << ":" << __LINE__ << ": " \
      << #x << "= " << x << std::endl;
  #define DEBUG_PRINT2(x) \
  std::cerr <<  "#> " << mw::time() << " " << x << std::endl;
  #define DEBUG_RUN(x) x
#endif

  

/**
 * A complex datatype for use by the C interfaces to ACML routines.
 * Outside mw namespace to comply with acml.h
 */
#ifndef _ACML_COMPLEX
#define _ACML_COMPLEX
/* from acml.h : */
typedef struct
{
  float real, imag;
} complex;
typedef struct
{
  double real, imag;
} doublecomplex;
#endif /* !defined(_ACML_COMPLEX) */


namespace mw {

// typedefs for vector and matrix data
typedef std::vector<double> v1D; //!< 1D array of doubles
typedef std::vector<v1D> v2D; //!< 2D array


/**
 * @brief operator<< for ACML's doublecomplex
 * @todo TODO rewrite as a doublecomplexToString() function?
 */
//template<class charT, class traits>
//inline std::basic_ostream<charT,traits>&
//operator<< (std::basic_ostream<charT,traits>& ostr, const doublecomplex& rhs) {
template <typename T>
inline std::ostream& operator<< (std::ostream& ostr, const doublecomplex& rhs) {
    ostr << "(" << rhs.real << "," << rhs.imag << ")";
    return ostr;
}

/**
 * @brief operator<< for vector<T> (T must have its own operator<<)
 * @todo TODO add vectorToString(vector v, string sep=", "....) function?
 */
template <typename T>
inline std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
  typename std::vector<T>::const_iterator it;
  for (it=v.begin(); it!=v.end(); ++it) {
    out << *it << " ";
  }
  return out;
}

/**
 * @brief operator<< for vector< vector<T> > (T must have its own operator<<)
 * @todo TODO add vectorvectorToString(v, string sep=", " ...) function?
 */
template <typename T>
inline std::ostream& operator<< (std::ostream& out, const std::vector< std::vector<T> >& v) {
  typename std::vector< std::vector<T> >::const_iterator it;
  for (it=v.begin(); it!=v.end(); ++it) {
    out << *it << std::endl;
  }
  return out;
}



  

// Misc

/**
 * @brief Current date and time.
 * @return string YYYY-MM-DD,HH:MM:SS containing current date and time
 */
std::string time();

/**
 * @brief Hostname.
 * @return the system's host name
 */
std::string host();

/**
 * @brief Command line arguments.
 * @return CL args as a single string.
 */
std::string commandLine(int argc, char *argv[]);

/**
 * @brief Information usually printed at the starting point.
 * @return "## Start (hostname) time"
 */
std::string startTime();

/**
 * @brief Information usually printed at the end.
 * @return "## Stop (hostname) time"
 */
std::string stopTime();


/**
 * @brief Converts @c double to @c std::string. DEPRECATED (c++11): use std::to_string
 * @param x double value.
 * @return string created from @a x.
 */
std::string doubleToString(double x);

/**
 * @brief Returns transposed 2D array.
 * Note: assumes that all inner vectors have the same size.
 */
v2D transpose(const v2D& data);


// Input / Output

/**
 * @brief Loads floating-point data stored in one-column format.
 * All comments, invalid data, extra columns etc. are silently ignored.
 * @return the first read line, if it was commented (by # in the first column);
 *   empty string otherwise
 * @param in input stream
 * @param v values from the 1st column (output)
 *
 * @todo TODO make it work for empty files and other special cases
 */
std::string loadData(std::istream& in, std::vector<double>& v);

/**
 * @brief Loads floating-point data stored in two-column format.
 * All comments, invalid data, extra columns etc. are silently ignored.
 * @return the first read line, if it was commented (by # in the first column);
 *   empty string otherwise
 * @param in input stream
 * @param vx values from the 1st column (output)
 * @param vy values from the 2nd column (output)
 *
 * @todo TODO make it work for empty files and other special cases
 */
std::string loadData(std::istream& in, std::vector<double>& vx,
                       std::vector<double>& vy);
/** @brief (deprecated) alias for backward comp. */
inline std::string load2dData(std::istream& in, std::vector<double>& vx,
                       std::vector<double>& vy) {
  return loadData(in, vx, vy);
}


/**
 * @brief Loads floating-point data stored in multi-column format.
 * All comments, invalid data etc. are silently ignored.
 * @return the first read line, if it was commented (by # in the first column);
 *   empty string otherwise
 * @param in input stream
 * @param t loaded values (output)
 *
 * @todo TODO make it work for empty files and other special cases
 * @todo TODO optimize: pre-alloc vectors
 */
std::string loadData(std::istream& in, std::vector< std::vector<double> >& t);

/**
 * @brief Saves floating-point data in one-column format.
 * @param out output stream
 * @param v values for the 1st column
 * @param info comment, will be saved in the first line as "# info"
 */
void saveData(std::ostream& out, const std::vector<double>& v,
                std::string info = std::string(""));

/**
 * @brief Saves floating-point data in two-column format.
 * Number of saved data lines = min(vx.size,vy.size).
 * @param out output stream
 * @param vx values for the 1st column
 * @param vy values for the 2nd column
 * @param info comment, will be saved in the first line as "# info"
 */
void saveData(std::ostream &out, const std::vector<double>& vx,
                const std::vector<double>& vy, std::string info = "");

/**
 * @brief Saves floating-point data in matrix format.
 * @param out output stream
 * @param v values
 * @param info comment, will be saved in the first line as "# info"
 */
void saveData(std::ostream &out, const std::vector<std::vector<double> >& v,
                std::string info = "");




// Numerics

/**
 * @brief Linear splines interpolation.
 * Assumptions: array xnodes is sorted; xnodes.size()==ynodes.size();
 * If x < xnodes[0] or x > xnodes[size-1]) returns 0.
 * @param x x-value
 * @param xnodes values x_i (interpolation nodes)
 * @param ynodes values y_i (interpolation nodes)
 * @return y(x) - interpolating value at x; for x less than xnodes[0] or more than xnodes[size-1] returns xnodes[0] or xnodes[size-1] respectively
 *
 * @todo TODO separate function for checking assumptions
 * @todo TODO implement a faster algorithm, or use GSL? see mw_lin_spline.{h,c}
 */
double linearSplineValue(double x, const std::vector<double>& xnodes,
                    const std::vector<double>& ynodes);

/**
 * @brief Points obtained using linear splines interpolation.
 * Assumptions: array xnodes is sorted; xnodes.size()==ynodes.size();
 *
 * @param dx discretisation step
 * @param xnodes values x_i (interpolation nodes)
 * @param ynodes values y_i (interpolation nodes)
 * @param vy values of y(x) at points xnodes[0], xnodes[0]+dx, ... ,xnodes[n-1] (n = size of xnodes)
 *
 * @todo TODO implement a faster algorithm, or use GSL? see mw_lin_spline.{h,c}
 * @todo return vy insead of having it as an argument
 */
void discreteLinearSpline(const std::vector<double>& xnodes,
                          const std::vector<double>& ynodes,
                          double dx, std::vector<double>& vy);

/**
 * @brief Least-squares linear fitting.
 * Find coefficients of y=ax+b line best fitted to points x[i],y[i].
 *
 * @param xnodes vector of x-values, length >= 2
 * @param ynodes vector of y-values, y.size() has to equal x.size
 * @return vector [a, b, u(a), u(b)], where u() is an estimate of the error
 * and equals NaN if only two points are given.
 */
std::vector<double> linearFit(const std::vector<double>& xnodes,
                              const std::vector<double>& ynodes);

/**
 * @brief Returns 1st derivative of the data series.
 * 
 * @todo implement order
 */
v1D derivative(const v1D& data, double dx, int order = 1);

// Statistics

/**
 * @brief Average of a vector<double>.
 *
 * @todo throw exc if t.size() == 0
 * @todo rewrite as a template (?)
 */
inline double average(const std::vector<double>& t) {
  double a = 0.;
  for (size_t i = 0; i!=t.size(); ++i) {
    a += t[i];
  }
  return a/(double)t.size();
}

/**
 * @brief Standard deviation of the mean
 * http://pl.wikipedia.org/wiki/Odchylenie_standardowe#Odchylenie_.C5.9Bredniej
 * sigma_<x> = sqrt ( sum((xi - <x>)^2) / (n*(n-1)) )
 * NOTE: 1 pass is not safe!!! --> sigma_<x>  =  sqrt( (<x^2> - <x>^2) / (n-1) )
 * 
 * @todo throw exc if t.size() < 2
 * @todo rewrite as a template
 * @todo see 
 * http://en.wikipedia.org/wiki/Standard_deviation#Rapid_calculation_methods
 * http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
 * 1) http://www.johndcook.com/blog/standard_deviation/
 * 2) Knuth actually has one pass algorithm like this
   --
   for N=1
   _M = x
   _C = 0.0;
   for N>1
   delta = x - _M;
   _M += delta / N;
   _C += delta*(x-_M)
   --
   _M  gets you mean, _C / (N-1) is variance
   --
 */
inline double averageStdDev(const std::vector<double>& t) {
  double x = 0.;
  double x2 = 0.;
  for (size_t i = 0; i!=t.size(); ++i) {
    x += t[i];
    x2 += t[i]*t[i];
  }
  x = x/(double)t.size();
  x2 = x2/(double)t.size();
  return sqrt( (x2-x*x) / ((double)t.size()-1.0) );
}


/**
 * @brief Sample standard deviation
 * 
 * sqrt( sum((x_i-<x>)^2) / (N-1) )
 * http://en.wikipedia.org/wiki/Standard_deviation#Corrected_sample_standard_deviation
 */
inline double stdDevSample(const std::vector<double>& t) {
  double m = average(t);
  double d = 0.;
  for (size_t i = 0; i!=t.size(); ++i) {
    d += (t[i]-m)*(t[i]-m);
  }
  return sqrt( d / ((double)t.size()-1.0) );
}

/**
 * @brief (Population) standard deviation
 * 
 * sqrt( sum((x_i-<x>)^2) / N )
 */
inline double stdDevPopul(const std::vector<double>& t) {
  return stdDevSample(t)*sqrt(((double)t.size()-1.0)/(double)t.size());
}

/**
 * @brief Standard deviation of the mean (standard error of the mean)
 * 
 * sqrt( sum((x_i-<x>)^2) / (N*(N-1)) )
 */
inline double stdDevMean(const std::vector<double>& t) {
  return stdDevSample(t)/sqrt((double)t.size());
}

/**
 * @brief Absolute Correlation Function (aka autocorrelation function)
 * 
 * Formula based on Mathematica's AbsoluteCorrelationFunction
 * R(k) = sum_{i=1}^{N-k} ( x_{i} x_{i+k} )
 * //R(k) = sum_{n=0}^{N-1} ( x_{n} x_{n-k} )
 * see also http://en.wikipedia.org/wiki/Autocorrelation#Signal_processing
 */
v1D autocorrAbs(const v1D& t);

/**
 * @brief Correlation Function (aka autocorrelation function)
 * 
 * Formula based on Mathematica's CorrelationFunction
 * 
 * r(k) = sum_{i=1}^{N-k}((x_i-<x>)(x_{i+k}-<x>)) / sum_{i=1}^{N}((x_i-<x>)^2)
 * see also
 * http://www.itl.nist.gov/div898/handbook/eda/section3/eda35c.htm
 * http://en.wikipedia.org/wiki/Autocorrelation#Estimation
 */
v1D autocorr(const v1D& t);


// Random numbers

/**
 *  pseudo-random double value in the range [min;max)
 */
inline double mwRandDouble(double min, double max) {
  return ( (max-min)*rand()/(RAND_MAX+1.0) + min );
}

/**
 *  pseudo-random int value in the range [0;max)
 */
inline int mwRandInt(int max) {
  return (int) ((double)max*rand()/(RAND_MAX+1.0));
}

/**
 * Faster rand value from [0;max)
 * $ man 3 rand :
 * "The  versions of rand() and srand() in the Linux C Library use the same
 * random number generator as random() and srandom(), so  the  lower-order
 * bits  should  be as random as the higher-order bits."
 */
inline int mwFastRandInt(int max) {
  return rand() % max ;
}

/**
 * Class for generating integer values randomly with given probablilities.
 */
class CMwRandIntChooser {
  protected:
    /**
     * Ranges calculated from weights of real values (lower range
     * boundaries are stored).
     */
    std::vector<int> m_range;

  public:
    /** A constructor.
     * @param weight weights for consecutive integers 0,1,...N-1
     * (normalization is not necessary).
     *
     * @todo check INT_MAX vs RAND_MAX vs types (int, uint, size_t,...)
     */
    CMwRandIntChooser(const std::vector<double> & weight);

    /** Get a random integer.
     * @return value from [0;N)
     */
    unsigned int getRandInt() const;
};


class CMwData {
  protected:
    std::string m_info_pre;
    std::string m_info_columns;
    std::string m_info_post;
    // TODO : abstract class + subclasses for n-column or matrix data
};

// Complex numbers for ACML

#ifdef _ACML_COMPLEX
/**
  @brief Convert C++ complex<double> to ACML's doublecomplex
  */
inline doublecomplex mwConv(const std::complex<double> & w) {
    doublecomplex z;
    z.real = w.real();
    z.imag = w.imag();
    return z;
}

/**
  @brief Convert ACML's doublecomplex to C++ complex<double>
  */
inline std::complex<double> mwConv(const doublecomplex & w) {
    return std::complex<double>(w.real, w.imag);
}

/**
  @brief Create ACML's doublecomplex structure.
  */
inline doublecomplex mwDC(double re, double im) {
    doublecomplex z;
    z.real = re;
    z.imag = im;
    return z;
}


#endif /* defined(_ACML_COMPLEX) */


// Physics

#ifndef ATOMIC_UNITS_H
  // Boltzmann constant
  // $ k_B = 3.16681520371153 \times {{10}^{-6}} {{\hbar }^{2}} {K}^{-1} a_{0}^{-2}\ {m_e}^{-1}$
  #define AU_kB 3.16681520371153e-6
#endif

/**
  @brief The Fermi-Dirac distribution.
  @param E energy [a.u.]
  @param EF = el-chem pot. [au]
  @param T temperature [K]
  @return Fermi-Dirac distr. value
 */
inline double f_FD(double E, double EF, double T) {
    if(T>0.0)
        return 1. / (1. + std::exp((E-EF)/(AU_kB*T)) );
    return ( E > EF ) ? 0.0 : 1.0 ;
}

/**
  @brief Derivative of the Fermi-Dirac distribution.
  @todo TODO analytical formula
  @param E energy [a.u.]
  @param dE delta en. [a.u.]
  @param EF = el-chem pot. [au]
  @param T temperature [K]
  @return simple approx.of the Fermi-Dirac distr. derivative
  */
inline double df_FD(double E, double EF, double T, double dE) {
  double f1 = f_FD(E-dE/2.,EF,T);
  double f2 = f_FD(E+dE/2.,EF,T);
  return (f2-f1)/dE;
}

// Geometry

/**
  @brief Checks if (x,y) is inside a hexagon of radius R, placed to have 
  the diagonal (and 2 vertices) on y-axis.
  */
inline bool isInsideHex(double x, double y, double R = 1.0) {
    const double b = (sqrt(3.)/2.)*R;
    if( ( fabs(x) > b ) || ( fabs(y) + fabs(x)/sqrt(3.) > R ) )
        return false;
    return true;
}

} // namespace mw

#endif // MW_TOOLS_H
