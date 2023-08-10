\mainpage MODUL Oriented Dynamics United Library

Framework for developing dynamics.

1. Molecule dynamics and path integral molecule dynamics.
  - MD and MD based techniques 
  - PIMD (done & to be optimized)
  - real-time dynamics: s-PILD (todo & tocheck)
  - MES-PIMD (optimzed but need to check)
2. Nonadiabatic dynamics:
  - Ehrenfest, CMM, wMM (weighted-MM), fMM (focued-MM), pMM (perturbation-MM).
  - Surface hopping.
  - Symmetrical quasi-classical.
3. System-bath dynamics: Hierarchy equation of motion, stochastic Schrödinger equation.
4. Ab Initio calculation (a basic Hartree-Fock is realized).

interfaced to:

1. toy models (harmonic, quadratic, quartic, morse, spin-boson, FMO, etc.)
2. few-bodies forcefield for small molecule in fortran90.
3. interface to gaussian16 and mndo99 for on-the-fly calculation.
