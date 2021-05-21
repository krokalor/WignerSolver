// File main.cpp - Copyright 2008 Maciej Woloszyn

/**
 * WignerEquation
 * @author Maciej Woloszyn <woloszyn@agh.edu.pl>, Pawel Wojcik
 *
 * Note: values in atomic units [a.u.], in other case variables have
 * information in their name (e.g. x_nm is in [nm], but x in [a.u.]).
 *
 * @todo TODO --input-file option (if not given, read from cin) or simple positional parameter (./exe "file.dat")
 * @todo TODO options for output files (names?)
 * @todo TODO (?) options and parameters in a separate class
 * @todo TODO read parameters also from a file, not only via CLI options
 *
 * @todo TODO see TODO list in the @c WignerEquationSolver class
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <boost/program_options.hpp>
#include <omp.h>

#include "mwAtomicUnits.h"
#include "WignerEquationSolver.h"
#include "Poisson1D.h"

#include "mwTools.h"
#include "ver.h"

using namespace std;
using namespace mw;
namespace po = boost::program_options;


int main(int argc, char *argv[]) {

  //omp_set_num_threads(4);

  string infoProgram;
  {
    stringstream s;
    s << "### mwWignerSC [version "<< VERSION << "]";
    s << " [cmdline " << mw::commandLine(argc,argv) << "]";
    infoProgram = s.str();
  }

  // Print info
  cout << infoProgram;
  cerr << "\n#E> " << infoProgram;

#ifdef DEBUG
  cerr << "\n\n";
  cerr << "#E> !!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  cerr << "#E> !!! DEBUG macro defined !!!\n";
  cerr << "#E> !!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n";
#endif

  // Parameters with default values:
  size_t N_k = 200;
  double T = 4.2; //temperature [K]
  double m_eff = 1.0; // effective mass [a.u.]
  double epsilon = 1.0; // relative permittivity of the material
  double p = 0.004; // N_D region "smoothness"

  double tau = -1.0; // no scattering if (tau <= 0)
  double tauGamma = -1.0; // no scattering at contacts if (tauGamma <= 0)
  double delta_x = 1*AU_1nm; // discr. step in the real space [a.u.]
  double contact_length = 10*AU_1nm; // length of left and right contacts [a.u.]
  size_t contact_N_x = (size_t)(contact_length/delta_x); // number of contact sites

  double E_F = 0.1*AU_1eV; // Fermi en. for the contacts [a.u.]
  double V_min = 0.0; // min. voltage [a.u.]
  double V_max = AU_1V; // max. voltage  [a.u.]
  double delta_V = (V_max-V_min)/10.5; // discr. step for voltage  [a.u.] ( .../10.5 ==> 10 points)
  size_t self_consistent_max_steps = 51;
  double self_consistent_epsilon = 1e-3;
  double self_consistent_omega = 0.05;
  int bc_type = 0; // 0==>SUPPLYFUNCTION; 1==>FERMIDIRAC;
  WignerEquationSolver::BoundaryConditionType bcType = WignerEquationSolver::SUPPLYFUNCTION;
  bool calcExpectValues = false;
  string info; // potential description
  vector<double> potEnergyUB; // pot energy profile of cond band bottom
  //string potentialFileName;

  // parse command line options using Boost lib.
  try {
    // Declare the supported options.
    po::options_description desc("Program reads potential from stdin.\nAllowed options (default values in parentheses)");
    desc.add_options()
        ("help,h", "produce help message")
//        ("input-file,i", po::value(&potentialFileName), "input file with potential values in [a.u.] (OBLIGATORY!)" )
        ("temperature,T", po::value(&T), ("set temperature [K] (" + std::to_string(T) + ")" ).c_str() )
        ("effective-mass,m", po::value(&m_eff), ("effective mass [a.u.] (" + std::to_string(m_eff) + ")").c_str() )
        ("epsilon,e", po::value(&epsilon), ("relative permittivity (" + std::to_string(epsilon) + ")").c_str() )
        ("N-k,N", po::value(&N_k), ( "number of discr. points in the k-space (" + std::to_string(N_k) + ")").c_str() )
        ("delta-x,d", po::value(&delta_x), ("discr. step in the real space [a.u.] (" + std::to_string(delta_x) + ")").c_str() )
        ("contact-length,c", po::value(&contact_length),  ("length of left and right contacts (to allow \"spacers\") [a.u.] (" + std::to_string(contact_length) + ")").c_str() )
//        ("N-x-spacer,S", po::value(&N_x_space),  ("left and right margins (\"spacers\") [sites=dx] (" + std::to_string(N_x_space) + ")").c_str() )
//        ("x-space-right", po::value(&x_space_right_nm), (string("simulated space to the right from the device [nm] (") + std::to_string(x_space_right_nm) + ")").c_str() )
        ("fermi-energy,E", po::value(&E_F),  ("Fermi en. for contacts [a.u.] (" + std::to_string(E_F) + ")").c_str() )
        ("doping-smooth,p", po::value(&p),  ("Fermi en. for contacts [a.u.] (" + std::to_string(E_F) + ")").c_str() )
        ("min-voltage,v", po::value(&V_min), ("min. voltage [a.u.] (" + std::to_string(V_min) + ")").c_str() )
        ("max-voltage,V", po::value(&V_max), ("max. voltage [a.u.] (" + std::to_string(V_max) + ")").c_str() )
        ("delta-voltage,D", po::value(&delta_V), ("discr. step for voltage [a.u.] (" + std::to_string(delta_V) + " => 10 points)").c_str() )
        ("relax-time,t", po::value(&tau), ("enables scattering with given relaxation time (negative=disabled) [a.u.] (" + std::to_string(tau) + ")").c_str() )
        ("relax-time-cont,G", po::value(&tauGamma), (string("enables scattering at contacts with given relaxation time [a.u.] (negative=disabled)")).c_str() )
        ("self-consistent,P", po::value(&self_consistent_max_steps), ("max. no. of iterations in self-consistent procedure (Poisson-Wigner) (" + std::to_string(self_consistent_max_steps) + ")").c_str())
        ("self-cons-epsilon,s", po::value(&self_consistent_epsilon), ("convergence parameter (" + std::to_string(self_consistent_epsilon) + ")").c_str() )
        ("self-cons-omega,w", po::value(&self_consistent_omega), ("self-consistent underrelaxation parameter (" + std::to_string(self_consistent_omega) + ")").c_str() )
        ("bc-type,b", po::value(&bc_type), ( "boundary condition type: 0=suppl.f.,1=FD (" + std::to_string(bc_type) + ")").c_str() )
        ("expect-values,x", "calculate expect.values <x>,<x^2>,<k>,<k^2>,<xk>" )
        ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help")) {
      cout << "\n" << desc << "\n";
      cout << "\nExample:\n  ./WignerEquationBandStorage -d 0.00074 < pot.dat\n";
      return 0;
    }
    if (vm.count("expect-values")) {
      calcExpectValues = true;
    }
    if(bc_type==1) bcType=WignerEquationSolver::FERMIDIRAC; // SUPPLYFUNCTION is default
    contact_N_x = (size_t)(contact_length/delta_x);
    //if (!vm.count("input-file")) {
    //  throw logic_error(string("Input file not specified! Use '-i' option."));
    //}
  }
  catch(exception& e) {
    cerr << "\nError: " << e.what() << "\n\n";
    return 1;
  }

  {
    cerr << "#E> main() : Loading potential data from stdin... \n";
    vector<double> vx, vy;
    info = mw::load2dData(cin, vx, vy);
    mw::discreteLinearSpline(vx, vy, delta_x, potEnergyUB);
  } // vx and vy no longer needed

  // allow only even N_k
  if(N_k%2 == 1) {
    N_k++;
    cerr << "#E> N_k was odd, increasing N_k to " << N_k  << endl;
  }

  // create solver object
  WignerEquationSolver w(N_k, T, m_eff, delta_x, E_F, potEnergyUB,
                         contact_N_x, tau, tauGamma, info, bcType);
  w.p = p;
  w.epsilon = epsilon;
  w.self_consistent_max_steps = self_consistent_max_steps;
  w.self_consistent_epsilon = self_consistent_epsilon;
  w.self_consistent_omega = self_consistent_omega;




  // print info (including start time and hostname)
  cout << w.parameters() << endl;
  cerr << "#E> " << w.parameters() << endl;

  if(V_max<=V_min) w.save_relax_steps = true;

  // loop over bias voltage values
  cout << left; // adjust output to the left
  cout << setw(10) << "# V[a.u.]" <<  " currentDensity[a.u.]    non-class.param.   sc_diff sc_iter";
  if(calcExpectValues) cout << " <x> <x^2> <k> <k^2> <xk> norm";
  cout << endl;
  for( double V_bias = V_min; V_bias <= V_max; V_bias += delta_V ) {
    cerr << "#E> V_bias= " << V_bias << endl;
    //w.setV(V_bias);
    //w.fillMatrices(); ///< @todo TODO call fillMatrices() from setV() (???)
    //w.solve_one_iter(); // solve eq...
    w.solve(V_bias);
    cout << setw(10) << w.getV() << " " << setw(10) << w.currentDensity() ;
    cout << " " << setw(10) << w.nonclParam() ;
    cout << "\t" << w.sc_diff << "\t" << w.sc_iter ;
    if(calcExpectValues) {
        auto ex = w.expectVals();
        cout << "\t" << ex.at(0) << "\t" << ex.at(1) << "\t" << ex.at(2) << "\t" << ex.at(3) << "\t" << ex.at(4) << "\t" << ex.at(5) ;
    }
    cout << endl;
  } // end loop over V_bias

  {
    cerr << "#E> main() : saving misc. data (only for the last value of V)..\n";

    ///< @todo TODO check I/O operations

    ofstream fw("W.out");
    fw << infoProgram << " " << w.parameters() << endl;
    w.printWignerFunction(fw);
    fw.close();

    ofstream fj("J.out");
    fj << infoProgram << " " << w.parameters() << endl;
    w.printCurrentDensity(fj);
    fj.close();

    ofstream fn("N.out");
    fn << infoProgram << " " << w.parameters() << endl;
    w.printChargeDensity(fn);
    fn.close();

    ofstream fu("U.out");
    fu << infoProgram << " " << w.parameters() << endl;
    w.printU(fu);
    fu.close();
  }

/*
  // BEGIN tmp Poisson -------------------
  // see also /home/mw-tmp/DEL/BJS/nu/TBRTD
  vector<double> n = w.chargeDensity();
  vector<double> rho = n;
  for(size_t i=0; i<w.Nx; ++i) {
      rho.at(i) = w.nD(i) - n.at(i);
  }
  Poisson1D p(rho);
  p.epsilonc = 12.9*AU_epsilon0; // 12.9 epsilon_0
  p.phiL = 0.;
  p.phiR =  V_max; //0.15*AU_1V; //0.;
  p.dx = delta_x;
  p.solve();
  ofstream fphi("Phi.out");
  //mw::saveData(fphi,p.phi,"# solution of Poisson Eq.");
  fphi << infoProgram << " " << w.parameters() << endl;
  fphi << "#\n#x phi(x) rho(x) n(x)\n";
  for(size_t i=0; i<w.Nx; ++i) {
      fphi << (int)i*delta_x << "\t" << p.phiAt(i) << "\t" << rho.at(i) << "\t" << n.at(i) << "\n";
    }
  fphi.close();
  // END tmp Poisson -------------------
*/

  // final info
  cout << "#-------------------------------" << endl;
  cout << "## Stop  (" << mw::host() << ") " << mw::time() << endl;
  cerr << "## Stop  (" << mw::host() << ") " << mw::time() << endl;

  return 0;
}
