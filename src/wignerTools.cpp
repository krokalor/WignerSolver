#include "lib.hpp"
#include "WignerFunction.hpp"

using namespace wigner;


void WignerFunction::setEquilibriumFunction(std::string input_pot, bool read_pot){
	setPotBias(0.);
	if (read_pot) {
		readPotential(input_pot);
		uC_ = uStart_;
	}
	solveWignerEq();
	fe_ = f_;
	f_.zeros();
}


//
// ############################## Current and current density calculation ##############################
//
double WignerFunction::calcCurr() {
	size_t_vec_d i, j;        // iterators
	double cur1 = 0, cur2 = 0, cur = 0;  // Current density
	//  cur1 - current density in node i-1/2
	//  cur2 - current density in node i+1/2
	std::ofstream file;
	file.open("data/data_files/CD.out", std::ios::out);
	for (i=2; i<nx_-2; ++i) {
		cur1 = 0, cur2 = 0;
		// if (i*dx_>lC_ && i*dx_<l_-lC_) {
		// #################### k<0 ####################
		for (j=0; j<nk2_; ++j){
			// cur1 += k_(j)*f_(i+1, j);  // j(i-1/2)
			// cur2 += k_(j)*f_(i, j);  // j(i+1/2)
			cur1 += k_(j)*(3*f_(i,j)-f_(i+1,j));  // j(i-1/2)
			cur2 += k_(j)*(3*f_(i+1,j)-f_(i+2,j));  // j(i+1/2)
		}
		for (j=nk2_; j<nk_; ++j){
			// cur1 += k_(j)*f_(i, j);  // j(i-1/2)
			// cur2 += k_(j)*f_(i-1, j);  // j(i+1/2)
			cur1 += k_(j)*(3*f_(i-1,j)-f_(i-2,j));
			cur2 += k_(j)*(3*f_(i,j)-f_(i-1,j));
		}
		cur1 *= dk_/(2*m_), cur2 *= dk_/(2*m_);  // /(4*M_PI*m_)
		cur += (cur1+cur2)*dx_/2./l_;
		file<<x_(i)*AU_nm<<' '<<cur1<<' '<<cur2<<' '<<cur<<'\n';
	}
	file<<"# Current: "<<cur<<endl;
	file.close();
	return cur;
}


vec WignerFunction::calcCurrArr() {
	size_t_vec_d i, j;        // iterators
	double cur1 = 0, cur2 = 0;
	vec cur(nx_, fill::zeros);  // Current density  // array<double>
	//  cur1 - current density in node i-1/2
	//  cur2 - current density in node i+1/2
	for (i=2; i<nx_-2; ++i) {
		cur1 = 0, cur2 = 0;
		// #################### k<0 ####################
		for (j=0; j<nk2_; ++j){
			// cur1 += k_(j)*f_(i+1, j);  // j(i-1/2)
			// cur2 += k_(j)*f_(i, j);  // j(i+1/2)
			cur1 += k_(j)*(3*f_(i,j)-f_(i+1,j));  // j(i-1/2)
			cur2 += k_(j)*(3*f_(i+1,j)-f_(i+2,j));  // j(i+1/2)
		}
		for (j=nk2_; j<nk_; ++j){
			// cur1 += k_(j)*f_(i, j);  // j(i-1/2)
			// cur2 += k_(j)*f_(i-1, j);  // j(i+1/2)
			cur1 += k_(j)*(3*f_(i-1,j)-f_(i-2,j));
			cur2 += k_(j)*(3*f_(i,j)-f_(i-1,j));
		}
		cur(i) = (cur1+cur2)*dk_/(2*m_);  // /(4*M_PI*m_)
	}
	cur(0) = cur(2), cur(1) = cur(2);
	cur(nx_-1) = cur(nx_-3), cur(nx_-2) = cur(nx_-3);
	return cur;
}


vec WignerFunction::calcCD_X(){
	// Calculates carrier density in x space
	vec cd(nx_, fill::zeros);  // array<double>
	for (size_t i=0; i<nx_; ++i)
		for (size_t j=1; j<nk_/2; ++j)
			// cd(i) += (f_(i,j-1) + f_(i,j))*dk_/2.;  //  / 2./M_PI  // trapezoid
			cd(i) += (f_(i,2*j-2)+4*f_(i,2*j-1)+f_(i,2*j))*dk_/3.;  // simpson
	return cd;
}


vec WignerFunction::calcCD_K(){
	// Calculates carrier density in k space
	vec cd(nk_, fill::zeros);  // array<double>
	for (size_t j=0; j<nk_; ++j)
		for (size_t i=1; i<nx_/2; ++i)
			// cd(j) += (f_(i-1,j) + f_(i,j))*dx_/2.;  // trapezoid
			cd(j) += (f_(2*i-2,j)+4*f_(2*i-1,j)+f_(2*i,j))*dx_/3.;  // simpson
	return cd;
}


double WignerFunction::calcNorm(){
	// Calculates density function norm (integral over whole space)
	double sum_x = 0, sum_k;
	for (size_t i=0; i<nx_; ++i) {
		sum_k = 0;
		for (size_t j=1; j<nk_; ++j)
			sum_k += ( f_(i,j) + f_(i,j-1) ) * dk_/2.;
		if (i == 0 || i == nx_-1)
			sum_x += sum_k * dx_/2.;
		else
			sum_x += sum_k * dx_;
	}
	return sum_x;  // /2./M_PI
}


double WignerFunction::calcEX(){
	double sum_x = 0, sum_k;
	for (size_t i=0; i<nx_; ++i) {
		sum_k = 0;
		for (size_t j=1; j<nk_; ++j)
			sum_k += ( f_(i,j) + f_(i,j-1) ) * x_(i) * dk_/2.;
		if (i == 0 || i == nx_-1)
			sum_x += sum_k * dx_/2.;
		else
			sum_x += sum_k * dx_;
	}
	return sum_x/2./M_PI / calcNorm();
}


double WignerFunction::calcEK(){
	double sum_x = 0, sum_k;
	for (size_t i=0; i<nx_; ++i) {
		sum_k = 0;
		for (size_t j=1; j<nk_; ++j)
			sum_k += ( f_(i,j) + f_(i,j-1) ) * k_(j) * dk_/2.;
		if (i == 0 || i == nx_-1)
			sum_x += sum_k * dx_/2.;
		else
			sum_x += sum_k * dx_;
	}
	return sum_x/2./M_PI / calcNorm();
}


double WignerFunction::calcSDK(){
	double ev = 0, ev2 = 0, s_ev, s_ev2;
	/*
	double sum_x = 0, sum_k;
	// E[k]
	sum_x = 0;
	for (size_t i=0; i<nx_; ++i) {
		sum_k = 0;
		for (size_t j=1; j<nk_; ++j)
			sum_k += ( f_(i,j) + f_(i,j-1) ) * k_(j) * dk_/2.;
		if (i == 0 || i == nx_-1)
			sum_x += sum_k * dx_/2.;
		else
			sum_x += sum_k * dx_;
	}
	ev = sum_x/2./M_PI;  //  / calcNorm()
	// E[k^2]
	sum_x = 0;
	for (size_t i=0; i<nx_; ++i) {
		sum_k = 0;
		for (size_t j=1; j<nk_; ++j)
			sum_k += ( f_(i,j) + f_(i,j-1) ) * k_(j) * k_(j) * dk_/2.;
		if (i == 0 || i == nx_-1)
			sum_x += sum_k * dx_/2.;
		else
			sum_x += sum_k * dx_;
	}
	ev2 = sum_x/2./M_PI;  //  / calcNorm()
	*/
	for (size_t i=0; i<nx_; ++i) {
		s_ev = 0, s_ev2 = 0;
		for (size_t j=1; j<nk_; ++j) {
			s_ev += ( f_(i,j) + f_(i,j-1) ) * k_(j) * dk_/2.;
			s_ev2 += ( f_(i,j) + f_(i,j-1) ) * k_(j) * k_(j) * dk_/2.;
		}
		if (i == 0 || i == nx_-1) {
			ev += s_ev * dx_/2.;
			ev2 += s_ev2 * dx_/2.;
		}
		else {
			ev += s_ev * dx_;
			ev2 += s_ev2 * dx_;
		}
	}
	ev = ev/2./M_PI / calcNorm();
	ev2 = ev2/2./M_PI / calcNorm();
	return sqrt(ev2-ev*ev);
}


double WignerFunction::calcSDX(){
	double ev = 0, ev2 = 0, s_ev, s_ev2;
	for (size_t i=0; i<nx_; ++i) {
		s_ev = 0, s_ev2 = 0;
		for (size_t j=1; j<nk_; ++j) {
			s_ev += ( f_(i,j) + f_(i,j-1) ) * x_(i) * dk_/2.;
			s_ev2 += ( f_(i,j) + f_(i,j-1) ) * x_(i) * x_(i) * dk_/2.;
		}
		if (i == 0 || i == nx_-1) {
			ev += s_ev * dx_/2.;
			ev2 += s_ev2 * dx_/2.;
		}
		else {
			ev += s_ev * dx_;
			ev2 += s_ev2 * dx_;
		}
	}
	ev = ev/2./M_PI / calcNorm();
	ev2 = ev2/2./M_PI / calcNorm();
	return sqrt(ev2-ev*ev);
}


//
// #################### Set system potential - constant ####################
//
void WignerFunction::setPotBias(double uBias) {
	uBias_ = uBias;
	uC_.zeros();
	double x;
	for (size_t i = 0; i < nx_; ++i) {
		x = x_(i)-lC_;
		if (x<0)
			uC_(i) = 0;
		else if (x>=0 && x<=lD_)
			uC_(i) = -uBias_*x/lD_;
		else if (x>lD_)
			uC_(i) = -uBias_;
	}
}  // End of setPotBias


void WignerFunction::addGaussPot(double v0 = 0.3/AU_eV, double x0 = 1000, double sX = 100) {
	for (size_t i=0; i<nx_; ++i)
		uB_(i) += exp(-(x_(i)-x0)*(x_(i)-x0)/sX/sX)*v0;
}


//
// ############################## Resonant Tunelling Diode ##############################
//
void WignerFunction::setRTD(double w1i, double w2i, double w3i, double w4i,  double V0i) {
	double w1 = w1i/AU_nm, w2 = w2i/AU_nm, w3 = w3i/AU_nm, w4 = w4i/AU_nm, V0 = V0i/AU_eV;
	// double V = uL_ - uR_;
	int p = 10;
	double a = w1, b = w1 + 2*w2 + 2*w3 + w4;
	double sig = w3/2.;
	double b1 = a+w2+sig, b2 = b-w2-sig;
	lD_ = 2*w2 + 2*w3 + w4, lC_ = w1, l_ = lD_+2*lC_ ;
	dx_ = l_/nx_;
	kmax_ = M_PI/2./dx_;
	dk_ = 2*kmax_/float(nk_);
	double e1, e2, x;
	for (size_t i = 0; i < nx_; ++i) {
		x = x_(i);
		e1 = -pow(x-b1, 2*p)/2./pow(sig, 2*p);
		e2 = -pow(x-b2, 2*p)/2./pow(sig, 2*p);
		u_(i) = V0 * ( exp(e1) + exp(e2) );
		if (x < a) u_(i) += 0;
		else if (x >=  a && x <= b) u_(i) += (uL_-uR_)*(x-a)/(a-b);
		else if (x > b) u_(i) += -(uL_-uR_);
	}
	for (size_t i = 0; i < nx_; ++i) {
		x = x_(i);
		e1 = -pow(x-b1, 2*p)/2./pow(sig, 2*p);
		e2 = -pow(x-b2, 2*p)/2./pow(sig, 2*p);
		du_(i) = -V0 * ( 2*p*e1/(x-b1) * exp(e1) + 2*p*e2/(x-b2) * exp(e2) );
		if (x < a || x > b) du_(i) += 0;
		else if (x >=  a && x <= b) du_(i) += (uL_-uR_)/(a-b);
	}
	std::ofstream file;
	file.open("data/data_files/pot.out", std::ios::out);
	for (size_t i = 0; i < nx_; ++i)
		file << x_(i)*AU_nm << ' ' << u_(i) * AU_eV << ' ' << du_(i) * AU_eV/AU_nm<< '\n';
	file.close();
}

//
// ############################## WP initial conditions ##############################
//
void WignerFunction::addWavePacket() {
	/*
	double deviceHW = 1e+6/a0_;
	gwp_tl_ = sqrt(2*M_PI*beta_/m_);
	gwp_dx_ = 0.28 * gwp_tl_;  // 0.28
	gwp_dp_ = sqrt(2*m_/beta_);
	gwp_x0_ = 3*gwp_dx_;
	gwp_A_ = part_num_ /
					(2.*M_PI*gwp_dx_*gwp_dp_*deviceHW*deviceHW*gwp_p0_*gwp_p0_);
			*/
	// gwp_dp_ = 1./(2*gwp_dx_);
	// gwp_p0_ = sqrt(2*m_*uF_);
	// double deviceHW = 1e3*1e3 / AU_nm/AU_nm;  // Pole przekroju (szer. x wysokość)
	// gwp_A_ = part_num_ /
	// 	(2.*M_PI*gwp_dx_*gwp_dp_*deviceHW*gwp_p0_*gwp_p0_);
	gwp_A_ = cD_*lC_ / (gwp_dx_*gwp_dp_*2*M_PI);
	cout<<"# Setting up GWP with parameters:\n"
	<<"# gwp_x0 = "<<gwp_x0_*AU_nm<<" nm, "<<gwp_x0_<<" a.u.\n"
	<<"# gwp_dx = "<<gwp_dx_*AU_nm<<" nm, "<<gwp_dx_<<" a.u.\n"
	<<"# gwp_p0 = "<<gwp_p0_<<" a.u.\n"
	<<"# gwp_dp = "<<gwp_dp_<<" a.u.\n"
	<<"# gwp_A = "<<gwp_A_/AU_cm2<<" cm^-2, "<<gwp_A_<<" a.u.\n"
	<<"# cD_*lC_ = "<<cD_*lC_/AU_cm2<<" cm^-2\n"<<endl;
	// double s2 = 2*gwp_dx_*gwp_dx_;
	double sx = 2*gwp_dx_*gwp_dx_, sp = 2*gwp_dp_*gwp_dp_;
	for (size_t i=0; i<nx_; ++i) {
		for (size_t j=0; j<nk_; ++j)
			f_(i,j) += exp(
				- (k_(j)-gwp_p0_)*(k_(j)-gwp_p0_)/sp
				- (x_(i)-gwp_x0_)*(x_(i)-gwp_x0_)/sx  ) * gwp_A_; // * gwp_A_, / M_PI
	}
}


//
//  WP evolution (no potential)
//
double WignerFunction::WavePacket_TE(double x, double k)
	{ return exp( -2*gwp_dx_*gwp_dx_*(k-gwp_p0_)*(k-gwp_p0_)
								-(x-gwp_x0_)*(x-gwp_x0_)/2./gwp_dx_/gwp_dx_ )/M_PI; }


void WignerFunction::calc_IVchar(double v_min, double v_max, size_t nv){
	double dv = (v_max-v_min)/(nv-1), v = 0, curr;
	double uBias = uBias_;
	iv_v_ = vec(nv, fill::zeros), iv_i_ = vec(nv, fill::zeros);
	// double E;
	std::ofstream ivChar, vpMap;
	ivChar.open("wyniki/dane/ivChar.out", std::ios::out);
	vpMap.open("wyniki/dane/vpMap.out", std::ios::out);
	vpMap<<"# u_B [eV]  p [a.u]  1/4  1/2  3/4\n";
	ivChar<<"# u_B [eV]  J [Acm^{-2}]\n";
	// cout<<"# u_B [eV]  J [Acm^{-2}]\n";
	for (size_t i = 0; i < nv; ++i) {
		// addGaussPot(i*dv);
		v = v_min+i*dv;
		setPotBias(v);
		solveWignerEq();
		// solveWignerPoisson();
		curr = calcCurr();
		// E = i*dv/lD_ * AU_eV/1e3 / AU_cm;
		ivChar<<v*AU_eV<<' '<<curr*AU_Acm2<<endl;
		// cout<<v*AU_eV<<' '<<curr*AU_Acm2<<endl;  // <<' '<<calcNorm()/AU_cm2<<endl;
		iv_v_(i) = v, iv_i_(i) = curr;
		for (size_t j=0; j<nk_; ++j) {
				vpMap<<v*AU_eV<<' '<<k_(j)
					<<' '<<f_(size_t(nx_/4.),j)  // col. 3
					<<' '<<f_(size_t(nx_/2.),j)  // col. 4
					<<' '<<f_(size_t(nx_*3/4.),j)<<'\n';  // col. 5
		}
		vpMap<<'\n';
	}
	ivChar.close();
	vpMap.close();
	uBias_ = uBias;
}


void WignerFunction::calcMobility() {

	vec ne = calcCD_X();
	vec dnedx = calcDer(ne, dx_);
	vec jn = calcCurrArr();
	vec el_f(nx_, fill::zeros);
	vec mob(nx_, fill::zeros);

	for (size_t i = 0; i < nx_; ++i)
		el_f(i) = -du_(i);

	for (size_t i = 0; i < nx_; ++i) {
		double m = ne(i)*el_f(i)+KB/AU_eV*temp_*dnedx(i);
		mob(i) = jn(i)/m * AU_cm2/AU_eV/AU_s;
	}

	/*
	for (size_t i = 0; i < nx_; ++i)
		cout<<x_(i)<<' '<<mob(i)<<' '<<ne(i)<<' '<<dnedx(i) \
			<<' '<<u_(i)*AU_eV<<' '<<el_f(i)<<' '<<jn(i)<<' '<<mob(i)*KB*temp_<<endl;
			*/

	// Average mobility
	double mob_av = 0;
	for (size_t i = 1; i < nx_; ++i)
		mob_av += (mob(i-1) + mob(i))*dx_/2./lD_;
	cout<<m_<<' '<<mob_av<<endl;

}


//
// FERMI INTEGRAL
//


//
// Density of states
//
inline double WignerFunction::nC(double m, double T)
	{ return 2.*pow(m*KB*T/AU_eV/2./M_PI, 3./2.); }  // [au]


//
// Calculates Fermi integral, eta - relative Fermi level - (Ec-Ef)/kb/T
//
inline double WignerFunction::fermiInt(double n, double eta) {
	// Calculate Fermi integral
	auto fi = [n, eta](double x) { return pow(x, n)*exp(eta-x)/(exp(eta-x)+1); } ;
	size_t N = 100, j;
	double h = 1./float(N);
	int i1 = 0, i2 = 0;
	for (size_t i=1; i<N; ++i){
		j = i*h;
		// # i1 += (fi(j)+4*fi((j+1)*h)+fi((j+2)*h))*h/3
		// # i2 += (fi(1/j/h)+4*fi(1/(j+1)/h)+fi(1/(j+2)/h))/(j*h)**2*h/3
		i1 += (fi(j)+fi(j+h))*h/2.;
		i2 += pow( (fi(1./j)+fi(1./(j+h)))/j, 2. )*h/2.;
	}
	double integral = i1+i2;
	// std::tgamma(n+1) <- calculates gamma function for n+1
	return 1./std::tgamma(n+1.)*integral;
}


//
// Calc Fermi integral
//
double WignerFunction::calcFermiEn(double n0, double m, double T) {
	// Calculate Fermi energy from el. density using Fermi integral
	double eta = 0., f = 0.;
	double x = n0/nC(m, T);  // For GaAs, 300 K -> x = 4.6
	//  Looking for eta
	if (x > fermiInt(.5, 5.)) {
		// TODO: Zamienić wyrażenie asymptotyczne na aproksymację
		eta = pow(x*std::tgamma(5./2.), 2./3.);
		/*
		eta = log(x) + 0.25*pow(2., .5)*x \
			+ (3/16.-1/9.*pow(3., .5))*x*x \
			+ (1/8.+5/48.*pow(2., .5)-1/9.*pow(6., .5))*x*x*x \
			+ (1585/6912+5/32.*pow(2., .5)-5/24.*pow(3., .5)-1/25.*pow(5., .5))*x*x*x*x; */
		// eta = log(x) + 3.53553E-1*x - 4.95009E-3*pow(x, 2.) + 1.48386E-4*pow(x, 3.) - 4.42563E-6*pow(x, 4.);
		/*
		while (f < x) {
			eta += 0.01;
			f = pow(eta, 3./2.)/std::tgamma(5/2);
			}
		*/
	}
	else if (x < fermiInt(.5, -2.))
		eta = log(x);  // *pow(2,-3/2.)*x
	else {
		if (x > fermiInt(.5, 0.)) {
			while (f < x) {
				eta += 1e-3;
				f = fermiInt(.5, eta);
			}
		}
		else {
			while (f > x) {
				eta -= 1e-3;
				f = fermiInt(.5, eta);
			}
		}
	}
	eta = eta*KB*T - M_PI*M_PI*KB*T/12./eta;
	// cout<<"eta = "<<eta/(KB*T)<<", x = "<<x<<endl;
	return eta/AU_eV;
}


double WignerFunction::calcFermiEn_MB(double n0, double m, double T) {
	double kBT = KB/AU_eV*T;
	double A = pow(m,3/2.)/sqrt(2)/pow(M_PI,3/2.)*pow(kBT,3/2.);
	return log(n0/A)*kBT;
}


//
// Warunek brzegowy
//


// TODO: wskaźniki na funkcję

void WignerFunction::setBoundCond(){
	Gamma_ = rG_*.5;
	// double cD = cD_ * 2*lC_/l_;
	uL_ = uF_;  // Fermi levels
	uR_ = uF_;
	// cout<<"# cD = "<<cD/AU_cm3<<", uL = uR = "<<uL_*AU_eV<<endl;
	double k;
	if (bcType_ == 0)
		for (size_t j=0; j<nk_; ++j) {
			k = k_(j);
			bc_(j) = supplyFunction(k);
		}
	else if (bcType_ == -1)
		bc_.zeros();
	else if (bcType_ == -2)
		for (size_t j=0; j<nk_; ++j) {
			k = k_(j);
			bc_(j) = maxwell_boltzmann(k);
		}
	else if (bcType_ == -3)
		for (size_t j=0; j<nk_; ++j) {
			k = k_(j);
			bc_(j) = gaussian_bc(k);
		}
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
			// TODO: Armadillo convolution
			// else if (bcType_ == 4)
			// 	vec a = vec(nk_, fill::zeros), b = vec(nk_, fill::zeros), c;
			// 	for (size_t j=0; j<nk_; ++j) {
			// 		a(j) = supplyFunction(k_(j));
			// 		b(j) = lorentz(k_(j));
			// 	}
			// 	c = conv(a, b);
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
	std::ofstream file;
	file.open("data/data_files/BC.out", std::ios::out);
	for (size_t j=0; j<nk_; ++j) {
		file<<k_(j)<<' '<<bc_(j)<<'\n';
	}
	file<<"# "<<calcInt(bc_, dk_)/2./M_PI/AU_cm3;
	file.close();
}


// Supply function (E(k))
double WignerFunction::supplyFunction(double k){
	// double mu = k > 0 ? uL_ : uR_;
	double mu = uF_;
	double m = m_;
	double c = m/M_PI*KB/AU_eV*temp_, ex = -(k*k/m/2.-mu)/(KB/AU_eV*temp_);	 // [au]
	// double c = 4*M_PI*m*KB/AU_eV*temp_, ex = -(k*k/m/2.-mu)/(KB/AU_eV*temp_);	 // [au]
	if (ex < 700)
		return c * log(exp(ex)+1);
	else{
		if (ex > 0)
			return c * ex;
		else return 0;
	}
}


// Maxwell-Boltzmann distribution
double WignerFunction::maxwell_boltzmann(double k){;
	double mu = uF_+uBias_;
	// double mu = k > 0 ? uL_ : uR_
	double m = m_;
	double v = k/m;
	double c = pow(2.*M_PI*m*(KB/AU_eV*temp_), 0.5);  // * cD_
	double ex = exp(-(m*v*v/2.-mu)/(KB/AU_eV*temp_));	 // [au]
	// double c = 4*M_PI*m*KB/AU_eV*temp_, ex = -(k*k/m/2.-mu)/(KB/AU_eV*temp_);	 // [au]
	return c * ex;
}

// Gaussian with sigma parameter
double WignerFunction::gaussian_bc(double k){
	double mu = uF_;
	// double mu = k > 0 ? uL_ : uR_;
	double m = m_;
	double pF = sqrt(uF_*2.*m);
	double kBT = KB/AU_eV*temp_;
	double sigma = 1;
	// double c = 1./sqrt(2.*M_PI*m*(KB/AU_eV*temp_)*sigma*sigma);
	// double ex = exp(-(k*k/m/2.-mu)/(KB/AU_eV*temp_)/sigma/sigma);	 // [au]
	// double c = 4*M_PI*m*KB/AU_eV*temp_, ex = -(k*k/m/2.-mu)/(KB/AU_eV*temp_);	 // [au]
	// double c = 1./sqrt(2.*M_PI*m*kBT)/sigma*cD_;
	double c = pow((2.*M_PI*m*kBT)*sigma*sigma,-1/2.)*cD_;
	double ex = exp(-k*k/2./m/kBT/sigma/sigma);  // (k-pF)*(k-pF)
	return c * ex;
}



// Supply function (x)
inline double WignerFunction::sf_x(double mu, double x){
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


// Gaussian with sigma parameter
double WignerFunction::gaussian_x(double mu, double x){
	// double mu = uF_+uBias_;
	// double mu = k > 0 ? uL_ : uR_;
	double sigma = 1;
	double kBT = KB/AU_eV*temp_;
	double m = m_;
	double c = pow((2.*M_PI*m*kBT)*sigma*sigma,-1/2.)*cD_;
	double ex = exp(-x/kBT/sigma/sigma);	 // [au]
	// double c = 4*M_PI*m*KB/AU_eV*temp_, ex = -(k*k/m/2.-mu)/(KB/AU_eV*temp_);	 // [au]
	return c * ex;
}


// Lorentz * SF convolution
double WignerFunction::lorentz(double k){
	// double mu = k > 0 ? uL_ : uR_;
	double mu = uF_;
	double m = m_;
	double u = k*k/m/2., g = Gamma_;
	// Lorentz profile
	auto f = [g](double x) { return g/(x*x+g*g)/M_PI; };
	// convolution
	size_t N = 1e4, i;
	double h = (40/beta_+mu)/float(N);  // 40 from exp(x) -> 0 in SF
	double fb = 0;
	double x0, x1, xm1, xm2, xm3;
	for (i=1; i<N-1; i++){
		x0 = i*h, x1 = (i+1)*h;
		xm2 = (x0+x1)/2., xm1 = (x0+xm2)/2., xm3 = (xm2+x1)/2.;
		fb += 7*f(x0-u) * gaussian_x(mu, x0) +
					32*f(xm1-u) * gaussian_x(mu, xm1) +
					12*f(xm2-u) * gaussian_x(mu, xm2) +
					32*f(xm3-u) * gaussian_x(mu, xm3) +
					7*f(x1-u) * gaussian_x(mu, x1);
	}
	fb *= 2*h/4./45.;
	// vec a(N, fill::zeros), b(N, fill::zeros), c(N, fill::zeros);
	// for (i=0; i<N; ++i) {
	// 	b(i) = f(i*h-u), a(i) = sf_x(mu, i*h);
	// }
	// c = conv(a, b, "same");
	// fb = sum(c);
	return fb;
}


// Gauss * SF convolution
double WignerFunction::gauss(double k){
	// double mu = k > 0 ? uL_ : uR_;
	double mu = uF_;
	double m = m_;
	double u = k*k/m/2., g = Gamma_;
	// Gauss profile
	auto f = [g](double x) { return 1/g/sqrt(2*M_PI)*exp(-x*x/2./g/g); };
	// convolution
	int N = 1e4, i;
	double fb = 0;
	double h = (40/beta_+mu)/float(N);  // 40 from exp(x) -> 0 in SF
	double x0, x1, xm1, xm2, xm3;
	for (i=1; i<N-1; i++){
		x0 = i*h, x1 = (i+1)*h;
		xm2 = (x0+x1)/2., xm1 = (x0+xm2)/2., xm3 = (xm2+x1)/2.;
		fb += 7*f(x0-u) * gaussian_x(mu, x0) +
					32*f(xm1-u) * gaussian_x(mu, xm1) +
					12*f(xm2-u) * gaussian_x(mu, xm2) +
					32*f(xm3-u) * gaussian_x(mu, xm3) +
					7*f(x1-u) * gaussian_x(mu, x1);
	}
	return fb*2*h/4./45.;
}


// Voigt * SF convolution
double WignerFunction::voigt(double k) {
	// double mu = k > 0 ? uL_ : uR_;
	double mu = uF_;
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
	int N = 1e4, i;
	double fb = 0;
	double h = (40/beta_+mu)/float(N);  // 40 from exp(x) -> 0 in SF
	double x0, x1, xm1, xm2, xm3;
	for (i=1; i<N-1; i++){
		x0 = i*h, x1 = (i+1)*h;
		xm2 = (x0+x1)/2., xm1 = (x0+xm2)/2., xm3 = (xm2+x1)/2.;
		fb += 7*f(x0-u) * gaussian_x(mu, x0) +
			32*f(xm1-u) * gaussian_x(mu, xm1) +
			12*f(xm2-u) * gaussian_x(mu, xm2) +
			32*f(xm3-u) * gaussian_x(mu, xm3) +
			7*f(x1-u) * gaussian_x(mu, x1);
	}
	return fb*2*h/4./45.;
}
