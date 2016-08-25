# D3D
This is a distribution of D3D, a software package for computing diffusion-reaction simulations in 3D using either finite-difference or Monte Carlo algorithms. This version of D3D has been specifically formatted to run simulations for the following eLife paper:

"Physical determinants of vesicle mobility and supply at a central synapse"
Rothman JS, Kocsis L, Herzog E, Nusser Z, Silver RA
eLife. 2016 Aug 19.
pii: e15133.
doi: 10.7554/eLife.15133.

Once the java code has been compiled and executed, a GUI will be created with an "eLife" menu at the top. Users can select which simulation to initialize via the following menu items:

1. Init Finite-Difference Frap (figures 2C 3D; Torok iPSF).
2. Init Monte Carlo FRAP (figures 3C,D 4A,B,D; Torok iPSF).
3. Init Monte Carlo FRAP Movie (no figure; default demo; Gaussian iPSF).
4. Init Monte Carlo MSD (figures 4C 5A).
5. Init Monte Carlo AZ (figures 7 8).

Note, simulations with "Torok iPSF" may take a while to initialize.

D3D will run on most machines which run Java. Install Java J2SE 5 or higher. It's better to download the JDK (Java Development Kit) which includes command line tools for Java.

D3D was created by the The Silver Lab -
R. Angus Silver - Department of Neuroscience, Physiology & Pharmacology -
University College London - Gower Street - London WC1E 6BT - UK

http://silverlab.org/

Java code was written by Jason Rothman.

For more information contact j.rothman@ucl.ac.uk

