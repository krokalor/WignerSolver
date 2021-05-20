/**
 * DO NOT EDIT outside SYNC/Sci/mwTools/
 * 
 * M.Woloszyn, 2008-2015
 *
 * Constants and SI units in atomic units: a_0 = m_e = hbar = e = 1
 * References:
 * - H. Shull and G. G. Hall, Atomic Units, Nature, volume 184, no. 4698, page 1559 (Nov. 14, 1959)
 * - http://folk.uio.no/michalj/node72.html
 * - http://en.wikipedia.org/wiki/Atomic_units
 * - http://www.phys.ubbcluj.ro/~vchis/cursuri/cspm/course2.pdf
 * - http://ilan.schnell-web.net/physics/rydberg.pdf
 * 
 * Test results:
 * =============

#include "mwAtomicUnits.h"
#include <iostream>

int main() {
  
 std::cout <<  " 1 a.u = " << 1. / AU_1eV << "eV\n";
 std::cout <<  " 1 a.u = " << 1. / AU_1nm << " nm\n" ;
 std::cout <<  " 1 a.u = " << 1. / AU_1s << " s\n";
 std::cout <<  " 1 a.u = " << 1. / AU_1V_m << " V/m\n";

 return 0;
}

 * prints:
   1 a.u = 27.2114 eV
   1 a.u = 0.0529177 nm
   1 a.u = 2.41888e-17 s
   1 a.u = 5.14221e+11 V/m
 
*/

#pragma once

#ifndef MW_ATOMIC_UNITS_H
#define MW_ATOMIC_UNITS_H 1


// CONSTANTS

// Boltzmann constant
// $ k_B = 3.16681520371153 \times {{10}^{-6}} {{\hbar }^{2}} {K}^{-1} a_{0}^{-2}\ {m_e}^{-1}$
#define AU_kB 3.16681520371153e-6

// speed of light
// $ c = 299792458 m/s$
#define AU_c 137.035999

// vacuum permittivity
// $ epsilon_0 = 1/(4*pi)
#define AU_epsilon0 (1.0/(4.*M_PI))


// UNITS

// electronvolt
// $ 1 eV = 0.03674932587122423 {{\hbar }^{2}} a_{0}^{-2} {m_e}^{-1}$
#define AU_1eV 0.03674932587122423

// angstrom 1 A = 10^{-10} m
// $ 1 A = 1.889726133921252 {a_{0}}$
#define AU_1AO 1.889726133921252

// nanometer 1 nm = 10^{-9} m
#define AU_1nm 18.89726133921252

// volt
// $ 1 V = 0.03674932587122423 {{\hbar }^{2}} e^{-1} a_{0}^{-2} {m_e}^{-1}$
#define AU_1V 0.03674932587122423

// electric field
// $ 1 V/m = 1.944690567144141 \times {{10}^{-12}} {{\hbar }^{2}} e^{-1} a_{0}^{-3} {m_e}^{-1}$
#define AU_1V_m 1.944690567144141e-12

// magnetic field
// $ 1 T = 4.254382547308656 \times {{10}^{-6}} {\hbar } e^{-1} a_{0}^{-2}$
#define AU_1T 4.254382547308656e-6

// ampere
// $ 1 A = 150.9749008100334 e \hbar a_{0}^{-2} {m_e}^{-1}$
#define AU_1A 150.9749008100334

// second
// $ 1 s = 4.134137337414122 \times {{10}^{16}} a_{0}^{2} {m_e}}{\hbar^{-1}$
#define AU_1s 4.134137337414122e16

// cm^{-3} = 1/(10^7 nm)^3
#define AU_1cmr3 1.481847092854201e-25

// A*cm^{-2}
#define AU_1Acmr2 4.22772776964475e-15

#endif /* atomic_units.h */

