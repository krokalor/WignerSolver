#include "lib.hpp"
#include "WignerFunction.hpp"

using namespace wigner;


// TODO: wskaźniki na funkcję

void WignerFunction::setBoundCond(){
	Gamma_ = rG_*.5;
	// TODO: uL_, uR_ -> uF_
	// double cD = cD_ * 2*lC_/l_;
	uL_ = calcFermiEn(cD_, m_, temp_);  // Fermi levels
  uR_ = calcFermiEn(cD_, m_, temp_);
	// cout<<"# cD = "<<cD/AU_cm3<<", uL = uR = "<<uL_*AU_eV<<endl;
	double k;
	if (bcType_ == 0)
		for (size_t j=0; j<nk_; ++j) {
			k = k_(j);
			bc_(j) = supplyFunction(k);
		}
	else if (bcType_ == -1)
		bc_.zero();
	else {
		if (rG_ > 0) {
			if (bcType_ == 1)
				for (size_t j=0; j<nk_; ++j) {
					k = k_(j);
					bc_(j) = lorentz(k);
				}
			else if (bcType_ == 2)
				for (size_t j=0; j<nk_; ++j) {
					k = k_(j);
					bc_(j) = gauss(k);
				}
			else if (bcType_ == 3)
				for (size_t j=0; j<nk_; ++j) {
					k = k_(j);
					bc_(j) = voigt(k);
				}
			else {
				cout << "# bcType_ = " << bcType_
					<< " IS WRONG BOUNDARY CONDITION TYPE INT, SETTING BC TO SUPPLY FUNCTION" << endl;
				for (size_t j=0; j<nk_; ++j) {
					k = k_(j);
					bc_(j) = supplyFunction(k);
				}
			}
		}
		else {
			cout << "# rG_ = " << rG_<< " SCATTERING RATE SHOULD BE GREATER THAN 0, SETTING BC TO SUPPLY FUNCTION" << endl;
			for (size_t j=0; j<nk_; ++j) {
				k = k_(j);
				bc_(j) = supplyFunction(k);
			}
		}
	}
	double sum = 0;
	std::ofstream file;
	file.open("data/data_files/BC.out", std::ios::out);
	for (size_t j=0; j<nk_; ++j) {
		file<<k_(j)<<' '<<bc_(j)<<'\n';
		sum += bc_(j);
	}
	file<<"# "<<sum;
	file.close();
}


// Supply function (E(k))
double WignerFunction::supplyFunction(double k){
	double mu = k > 0 ? uL_ : uR_;
	double m = m_;
	double c = m/M_PI*KB/AU_eV*temp_, ex = -(k*k/m/2.-mu)/(KB/AU_eV*temp_);	 // [au]
	if (ex < 700)
		return c * log(exp(ex)+1);
	else{
		if (ex > 0)
			return c * ex;
		else return 0;
	}
}


// Supply function (x)
double WignerFunction::sf_x(double mu, double x){
double m = m_;
double c = m/M_PI*KB/AU_eV*temp_, ex = -(x-mu)/(KB/AU_eV*temp_);	 // [au]
	if (ex < 700)
		return c * log(exp(ex)+1);
	else{
		if (ex > 0)
			return c * ex;
		else return 0;
	}
}


// Lorentz * SF convolution
double WignerFunction::lorentz(double k){
	double mu = k > 0 ? uL_ : uR_;
	double m = m_;
	double u = k*k/m/2., g = Gamma_;
	// Lorentz profile
	auto f = [g](double x) { return g/(x*x+g*g)/M_PI; };
	// convolution
	int N = 6E3, i;
	double fb = 0;
	double h = (40/beta_+mu)/float(N);  // 40 from exp(x) -> 0 in SF
	double x0, x1, xm1, xm2, xm3;
	for (i=1; i<N-1; i++){
		x0 = i*h, x1 = (i+1)*h;
		xm2 = (x0+x1)/2., xm1 = (x0+xm2)/2., xm3 = (xm2+x1)/2.;
		fb += 7*f(x0-u) * sf_x(mu, x0) +
					32*f(xm1-u) * sf_x(mu, xm1) +
					12*f(xm2-u) * sf_x(mu, xm2) +
					32*f(xm3-u) * sf_x(mu, xm3) +
					7*f(x1-u) * sf_x(mu, x1);
	}
	return fb*2*h/4./45.;
}


// Gauss * SF convolution
double WignerFunction::gauss(double k){
	double mu = k > 0 ? uL_ : uR_;
	double m = m_;
	double u = k*k/m/2., g = Gamma_;
	// Gauss profile
	auto f = [g](double x) { return 1/g/sqrt(2*M_PI)*exp(-x*x/2./g/g); };
	// convolution
	int N = 6E3, i;
	double fb = 0;
	double h = (40/beta_+mu)/float(N);  // 40 from exp(x) -> 0 in SF
	double x0, x1, xm1, xm2, xm3;
	for (i=1; i<N-1; i++){
		x0 = i*h, x1 = (i+1)*h;
		xm2 = (x0+x1)/2., xm1 = (x0+xm2)/2., xm3 = (xm2+x1)/2.;
		fb += 7*f(x0-u) * sf_x(mu, x0) +
					32*f(xm1-u) * sf_x(mu, xm1) +
					12*f(xm2-u) * sf_x(mu, xm2) +
					32*f(xm3-u) * sf_x(mu, xm3) +
					7*f(x1-u) * sf_x(mu, x1);
	}
	return fb*2*h/4./45.;
}


// Voigt * SF convolution
double WignerFunction::voigt(double k) {
	double mu = k > 0 ? uL_ : uR_;
	double m = m_;
	double u = k*k/m/2., g = Gamma_;
	// Calculating eta - mixing parameter in pseudo-Voigt profile
	double gammaL = 2*Gamma_, gammaG = 2*sqrt(2*log(2))*Gamma_;
	double gamma = pow(pow(gammaG,5) +
		2.69269*pow(gammaG,4)*gammaL +
		2.42843*pow(gammaG,3)*pow(gammaL,2) +
		4.47163*pow(gammaG,2)*pow(gammaL,3) +
		0.07842*gammaG*pow(gammaL,4) + pow(gammaL,5), 0.2);
	double eta = 1.36603*(gammaL/gamma) - 0.47719*pow(gammaL/gamma,2) + 0.11116*pow(gammaL/gamma, 3);
	// (pseudo-)Voigt profile
	auto f = [eta, g](double x)
		{ return eta*( g/(x*x+g*g)/M_PI ) +  // lorentz
			(1-eta)*( 1/g/sqrt(2*M_PI)*exp(-x*x/2./g/g) ); };  // gauss
	// convolution
	int N = 6E3, i;
	double fb = 0;
	double h = (40/beta_+mu)/float(N);  // 40 from exp(x) -> 0 in SF
	double x0, x1, xm1, xm2, xm3;
	for (i=1; i<N-1; i++){
		x0 = i*h, x1 = (i+1)*h;
		xm2 = (x0+x1)/2., xm1 = (x0+xm2)/2., xm3 = (xm2+x1)/2.;
		fb += 7*f(x0-u) * sf_x(mu, x0) +
					32*f(xm1-u) * sf_x(mu, xm1) +
					12*f(xm2-u) * sf_x(mu, xm2) +
					32*f(xm3-u) * sf_x(mu, xm3) +
					7*f(x1-u) * sf_x(mu, x1);
	}
	return fb*2*h/4./45.;
}
