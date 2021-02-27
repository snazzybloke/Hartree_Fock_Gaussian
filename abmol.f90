! Summary:
! A Hartree-Fock program using Gaussian type orbital (GTO) basis sets
!###################################################
! Author: Ante Bilic                               !
! Since: Mar 30, 2011                              !
! Project: Hartree-Fock-GTO                        !
! Version: N/A                                     !
! Maintainer: Ante Bilic                           !
! Email: ante.bilic.mr@gmail.com                   !
! Status: N/A                                      !
!###################################################


PROGRAM abmol
! This is a driver program with inner subroutines readmol, nuclrep, matsize, atomics
! which read in input file, evaluate the Coulomb repulsive energy between nuclei,
! the size of the matrices (msize), and initial charge density as the sum of the charge
! density of the constituent atoms, respectively.
! Then the external subroutines integral (which evaluates the overlap matrix si,
! the electronic kinetic energy matrix ti, 1-e v1ei and 2-e v2ei electrostaic
! contributions), collect and scf are invoked to find the self consistent solution
! of the Roothan-Hall eqs.

  USE abtype
  USE ktype
  IMPLICIT NONE
  CHARACTER (LEN=80) :: finput
  INTEGER(i1b) :: natoms, charge, msize, istat, nel
  REAL(dp) :: enucl
  CHARACTER (LEN=7) :: bchoice
  TYPE(koord), DIMENSION(:), ALLOCATABLE :: x
  INTEGER(i1b), DIMENSION(:), ALLOCATABLE ::  mat_i

  IF (iargc() /= 1) THEN
    WRITE (6, *) 'Syntax: abmol infile'
    CALL Exit(0)
  ELSE
    CALL getarg(1, finput)
  ENDIF

  CALL readmol
  CALL nuclrep
  WRITE(*,'(1X,A," SCF calculation" //)') bchoice
  CALL matsize
  CALL atomics
  CALL integral(natoms, x, bchoice, msize)
  CALL collect(natoms, msize)
  CALL scf(natoms, nel, msize, enucl)

  Deallocate(x, STAT=istat)

  CONTAINS


  SUBROUTINE readmol
! Reads (from the input file):
! -----
! natoms (int):      the total number of atoms in the file
! bchoice (char):    basis set choice (e.g., STO-3G, 3-21G, 4-31G, 6-31G**)
! charge (int):      the total charge of the molecule (typically 0, neutral)
! x(1D koord array): atom chem symbol, xyz coordinates and atomic number
! 
! Evaluates:
! ----------
! nel (int):      the total electronic charge (i.e., the number of electrons)

    IMPLICIT NONE
    INTEGER(i1b) :: eof, i, j
    CHARACTER (LEN=80) :: line
    CHARACTER (LEN=2) :: celem
    OPEN(1, FILE=finput)
    READ(1, *) natoms  ! Number of atoms
    Allocate(x(natoms), STAT=istat)
    Allocate(mat_i(natoms), STAT=istat)
    READ(1, *) bchoice, charge
    nel = 0
    DO i=1, natoms
      READ(1,*) x(i)%clem, (x(i)%xyz(j), j=1, 3), x(i)%anum
      nel = nel + x(i)%anum
    END DO
    nel = nel - charge ! the total electronic charge
    CLOSE(1)
  END SUBROUTINE readmol


  SUBROUTINE nuclrep
! Calculates the nuclear repulsion
!
! Evaluates:
! ----------
! enucl (real):     the Coulomb energy from positive nuclear charges

    IMPLICIT NONE
    INTEGER(i1b) :: i, j
    enucl = 0.0_dp
    ! loop over all pairs of atoms
    DO i=1, natoms
      DO j=i+1, natoms
        enucl = enucl + x(i)%anum * x(j)%anum / Sqrt(f_dstnc2(x(i), x(j)))
      END DO
    END DO
  END SUBROUTINE nuclrep


  SUBROUTINE matsize
! From the basis set choice bchoice and atomic numbers it determines
! the size of the matrices (msize) spanned by the set of the GTOs.
! For the initial density guess, it does the same job for the individual
! atoms, saving the size of their matrices in the array mat_i.
!
! Evaluates:
! ----------
! msize (int):           the size of the system matrices (overlap, Fock, etc)
! mat_i (1D int array):  the size of the system matrices for each atom

    INTEGER(i1b) :: i
    msize = 0
    DO i=1, natoms
      IF (bchoice == 'STO-3G ') THEN
        IF (x(i)%anum <= 2) THEN   ! H and He
          msize = msize + 1  ! 1s
          mat_i(i) = 1
        ELSE IF (x(i)%anum <= 10) THEN   ! Li to Ne
          msize = msize + 5  ! 1s 2s 2px 2py 2pz
          mat_i(i) = 5
        ELSE IF (x(i)%anum <= 18) THEN  ! Na to Ar
          msize = msize + 9 ! 1s 2s 2px 2py 2pz 3s 3px 3py 3pz
          mat_i(i) = 9
        END IF
      ELSE IF (bchoice == '3-21G ' .OR. bchoice == '4-31G ') THEN
        IF (x(i)%anum <= 2) THEN   ! H and He
          msize = msize + 2 ! 1s 1s'
          mat_i(i) = 2
        ELSE IF (x(i)%anum <= 10) THEN   ! Li to Ne
          msize = msize + 9  ! 1s  2s 2px 2py 2pz 2s' 2px' 2py' 2pz'
          mat_i(i) = 9
        ELSE IF (x(i)%anum <= 18) THEN   ! Na to Ar
          msize = msize + 13 ! 1s 2s 2px 2py 2pz  3s 3px 3py 3pz 3s' 3px' 3py' 3pz'
          mat_i(i) = 13
        END IF
      ELSE IF (bchoice == '6-31G**') THEN
        IF (x(i)%anum <= 2) THEN   ! H and He
          msize = msize + 5  ! 1s 1s' 2px 2py 2pz
          mat_i(i) = 5
        ELSE IF (x(i)%anum <= 10) THEN   ! Li to Ne
          msize = msize + 15  ! 1s  2s 2px 2py 2pz 2s' 2px' 2py' 2pz' + 6d
          mat_i(i) = 15
        ELSE IF (x(i)%anum <= 18) THEN   ! Na to Ar
          msize = msize + 19  ! 1s 2spxyz   3spxyz3s'p'xyz + 6 d orbitals
          mat_i(i) = 19
        END IF
      END IF
    END DO
  END SUBROUTINE matsize


  SUBROUTINE atomics
! Solves the HF problem for individual atoms, to set up initial guess for density matrix mP.
    INTEGER(i1b) :: i
    DO i=1, natoms
      CALL integral(1_i1b, x(i), bchoice, mat_i(i))
      CALL collect(1_i1b, mat_i(i))
      CALL scf(1_i1b, x(i)%anum, mat_i(i), 0.0_dp)
    END DO
  END SUBROUTINE atomics


END PROGRAM abmol
