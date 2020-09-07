#include "lib.hpp"
#include "WignerFunction.hpp"

using namespace wigner;


void WignerFunction::setEquilibriumFunction(){

    for (size_t i=0; i<nx_; ++i){
      if (uL_ >= uR_) {
          for (size_t j=0; j<nk2_; ++j){
              fe_(i, j) = bc_(j);
              fe_(i, nk_-j-1) = bc_(j);
          }
      }
      else if (uL_ < uR_) {
          for (size_t j=nk2_; j<nk_; ++j){
              fe_(i, j) = bc_(j);
              fe_(i, nk_-j-1) = bc_(j);
          }
      }
    }

    /*
    for (size_t i=0; i<nx_; ++i) {
        for (size_t j=0; j<nk2_; ++j)
            if (x_(i)>l_-lC_) f_(i,j) = bc_(j);
        for (size_t j=nk2_; j<nk_; ++j)
            if (x_(i)<lC_) f_(i,j) = bc_(j);
    }
    */

    std::ofstream file;
    file.open("data/data_files/FE.out", std::ios::out);
    for (size_t j=0; j<nk_; ++j)
        file<<k_(j)<<' '<<fe_(size_t(nx_/2), j)<<'\n';
    file.close();

}

// ############################## Current and current density calculation ##############################

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
    cur1 *= dk_/(4*M_PI*m_), cur2 *= dk_/(4*M_PI*m_);
    cur += (cur1+cur2)*dx_/2./l_;
    file<<x_(i)*AU_nm<<' '<<cur1<<' '<<cur2<<' '<<cur<<'\n';
  }
  file<<"# Current: "<<cur<<endl;
  file.close();
  return cur;
}

array<double> WignerFunction::calcCurrArr() {
  size_t_vec_d i, j;        // iterators
  double cur1 = 0, cur2 = 0;
  array<double> cur(nx_);  // Current density
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
    cur(i) = (cur1+cur2)*dk_/(4*M_PI*m_);
  }
  cur(0) = cur(2), cur(1) = cur(2);
  cur(nx_-1) = cur(nx_-3), cur(nx_-2) = cur(nx_-3);
  return cur;
}

array<double> WignerFunction::calcCD_X(){
  // Calculates carrier density in x space
  array<double> cd(nx_);
  for (size_t i=0; i<nx_; ++i)
    for (size_t j=1; j<nk_; ++j)
      cd(i) += (f_(i,j-1) + f_(i,j))*dk_/2. / 2./M_PI;
  return cd;
}

array<double> WignerFunction::calcCD_K(){
  // Calculates carrier density in k space
  array<double> cd(nk_);
  for (size_t j=0; j<nk_; ++j)
    for (size_t i=1; i<nx_; ++i)
      cd(j) += (f_(i-1,j) + f_(i,j))*dx_/2.;
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
  return sum_x/2./M_PI;
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


// #################### Set system potential - constant ####################
void WignerFunction::setLinPot(double uB) {
  uB_ = uB;
  u_.zero();
  double x;
  for (size_t i = 0; i < nx_; ++i) {
    x = x_(i)-lC_;
    if (x>=0 && x<=lD_)
      u_(i) = uL_ - uB_*x/lD_;
    else if (x<0)
      u_(i) = uL_;
    else if (x>lD_)
      u_(i) = uR_ - uB_;
  }
  du_ = calcDer(u_, dx_);
}  // End of setLinPot


void WignerFunction::setGaussPot(double v0 = 0.3/AU_eV, double x0 = 1000, double sX = 100) {
  // uB_ = vB;
  // u_.zero();
  for (size_t i=0; i<nx_; ++i)
    u_(i) += exp(-(x_(i)-x0)*(x_(i)-x0)/sX/sX)*v0; // * gwp_A_;
  // double x;
  // for (size_t i = 0; i < nx_; ++i) {
  //   x = x_(i)-lC_;
  //   if (x>=0 && x<=lD_)
  //     u_(i) += uL_ + uB_*x/lD_;
  //   else if (x<0)
  //     u_(i) += uL_;
  //   else if (x>lD_)
  //     u_(i) += uR_ + uB_;
  // }
  du_ = calcDer(u_, dx_);
}



// ############################## Resonant Tunelling Diode ##############################
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

// ############################## WP initial conditions ##############################
void WignerFunction::setWavePacket() {
  /*
  double deviceHW = 1e+6/a0_;
  gwp_tl_ = sqrt(2*M_PI*beta_/m_);
  gwp_dx_ = 0.28 * gwp_tl_;  // 0.28
  gwp_dp_ = sqrt(2*m_/beta_);
  gwp_x0_ = 3*gwp_dx_;
  gwp_p0_ = sqrt(2*m_*uL_);
  gwp_A_ = part_num_ /
          (2.*M_PI*gwp_dx_*gwp_dp_*deviceHW*deviceHW*gwp_p0_*gwp_p0_);
  cout<<"Setting up GWP with parameters:\n"
      <<"gwp_x0 = "<<gwp_x0_*a0_<<" nm, "<<gwp_x0_<<" a.u.\n"
      <<"gwp_dx = "<<gwp_dx_*a0_<<" nm, "<<gwp_dx_<<" a.u.\n"
      <<"gwp_p0 = "<<gwp_p0_*hbar_/a0_*1e9<<" kg*m*s^-1, "<<gwp_p0_<<" a.u.\n"
      <<"gwp_dp = "<<gwp_dp_*hbar_/a0_*1e9<<" kg*m*s^-1, "<<gwp_dp_<<" a.u.\n"
      <<"gwp_A = "<<gwp_A_/hbar_<<" (Js)^-1, "<<gwp_A_<<endl;
      */
  gwp_dp_ = 1./(2*gwp_dx_);
  double s2 = 2*gwp_dx_*gwp_dx_;
  for (size_t i=0; i<nx_; ++i) {
    for (size_t j=0; j<nk_; ++j)
      f_(i,j) += exp(
        - (k_(j)-gwp_p0_)*(k_(j)-gwp_p0_)*s2
        - (x_(i)-gwp_x0_)*(x_(i)-gwp_x0_)/s2  )/M_PI; // * gwp_A_;
  }
}

//  WP evolution (no potential)
double WignerFunction::WavePacket_TE(double x, double k)
  { return exp( -2*gwp_dx_*gwp_dx_*(k-gwp_p0_)*(k-gwp_p0_)
                -(x-gwp_x0_)*(x-gwp_x0_)/2./gwp_dx_/gwp_dx_ )/M_PI; }


void WignerFunction::calc_IVchar(){
  double dv = (v_max_-v_min_)/(nv_-1), v = 0;
  // double E;
  std::ofstream ivChar, vpMap;
  vpMap.open("wyniki/vpMap.out", std::ios::out);
  ivChar.open("wyniki/ivChar.out", std::ios::out);
  vpMap<<"# v_B [V]  p [a.u]  1/4  1/2  3/4\n";
  ivChar<<"# v_B [V]  J [Acm^{-2}]\n";
  cout<<"# v_B [V]  J [Acm^{-2}]\n";
  for (size_t i = 0; i < nv_; ++i) {
    // setGaussPot(i*dv);
    v = v_min_+i*dv;
    setLinPot(v);
    // solveWignerEq();
    solveWignerPoisson();
    // E = i*dv/lD_ * AU_eV/1e3 / AU_cm;
    ivChar<<v*AU_eV<<' '<<calcCurr()*AU_Acm2<<endl;  // <<' '<<calcNorm()/AU_cm2<<endl;
    cout<<v*AU_eV<<' '<<calcCurr()*AU_Acm2<<endl;  // <<' '<<calcNorm()/AU_cm2<<endl;
    for (size_t j=0; j<nk_; ++j) {
        vpMap<<v*AU_eV<<' '<<k_(j)
          <<' '<<f_(size_t(nx_/4.),j)
          <<' '<<f_(size_t(nx_/2.),j)
          <<' '<<f_(size_t(nx_*3/4.),j)<<'\n';
    }
    vpMap<<'\n';
  }
  ivChar.close();
  vpMap.close();
}


void WignerFunction::calcMobility() {

  array<double> ne = calcCD_X();
  array<double> dnedx = calcDer(ne, dx_);
  array<double> jn = calcCurrArr();
  array<double> el_f(nx_);
  array<double> mob(nx_);

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


// Density of states
inline double WignerFunction::nC(double m, double T)
  { return 2.*pow(m*KB*T/AU_eV/2./M_PI, 3./2.); }  // [au]

// Calculates Fermi integral, eta - relative Fermi level - (Ec-Ef)/kb/T
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

// Calc Fermi integral
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
