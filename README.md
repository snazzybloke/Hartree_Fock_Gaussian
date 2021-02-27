This repo shows 10 Fortran files which implement the self-consistent solution of
the Hartre-Fock secular equation (or Roothaan-Halls eqs) in the basis of Gaussian type orbitals.

There is a lib folder with a file libnag.a which only provides the
gamma and incomplete gamma functions (s14aaf and s14baf, respectively).

The PDF file HartreFock.pdf provides the description of the 
To compile the code, customize the variables in the Makefile and type,
place the libnag.a in a suitable directory (and/or change the "-L/... -lnag" line appopriately)
and type:
> make abmol

To run the code, use a file from the inputs directory. e.g.:
> abmol inputs/CH4_631.au
