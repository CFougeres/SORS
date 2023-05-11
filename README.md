# Simulation Of Resonant Scattering

User guide of SORS code

Author: C. Fougères (2023) Contact: cfougeres@anl.gov

References:

##############
 Requirements:

- CERN/ROOT version 6 (see https://root.cern.ch/)

##############
 Run:

Execute MC15O4He.C file in root terminal (or MC15N4He.C)

##############
 Inputs:
- user paths of sp_srim/ and scattering/ folder.
- sp_srim/ == stopping power tables from SRIM [J. F. Ziegler et al., SRIM Co., United States of America 6th (2013)] with S.P. in keV.micron
- scattering/ == angular distribution of scattered (beam, alphas) in windows and gas (from SRIM Monte Carlo simulations)
- Pol2 function to determine center-of-mass energies from measured alphas in Si detectors (from code poly_Eameas_Eacm.C)
