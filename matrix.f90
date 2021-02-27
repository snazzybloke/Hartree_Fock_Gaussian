!###################################################
! Author: Ante Bilic                               !
! Since: Mar 30, 2011                              !
! Project: Hartree-Fock-GTO                        !
! Version: N/A                                     !
! Maintainer: Ante Bilic                           !
! Email: ante.bilic.mr@gmail.com                   !
! Status: N/A                                      !
!###################################################


MODULE matrixs
! Defines the set of square matrices of essential interest (overlap mS, canonical transformation mX,
! its transpose mXt, "core" Hamiltonian mH, Fock operator mF, 2-e contributions direct + exchange mG,
! orbital coefficients mC, transformed Fock operator mFp, transformed coefficients mCp, density mP
! old denisty mP_old, MO energy eigenvalues mEig) and mVee (same as v2ei).

  USE abtype
  IMPLICIT NONE
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: mS, mX, mXt, mH, mF, mG, mC, &
    mFp, mCp, mP, mP_old, mEig
  REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: mVee

  CONTAINS

  SUBROUTINE setmxs(msize)
    INTEGER(i1b), INTENT(IN) :: msize
    INTEGER(i1b) :: istat
    Allocate(mS(msize, msize), STAT=istat)
    Allocate(mX(msize, msize), STAT=istat)
    Allocate(mXt(msize, msize), STAT=istat)
    Allocate(mH(msize, msize), STAT=istat)
    Allocate(mF(msize, msize), STAT=istat)
    Allocate(mG(msize, msize), STAT=istat)
    Allocate(mC(msize, msize), STAT=istat)
    Allocate(mFp(msize, msize), STAT=istat)
    Allocate(mCp(msize, msize), STAT=istat)
    Allocate(mP(msize, msize), STAT=istat)
    Allocate(mP_old(msize, msize), STAT=istat)
    Allocate(mEig(msize, msize), STAT=istat)
    Allocate(mVee(msize, msize, msize, msize), STAT=istat)
  END SUBROUTINE setmxs

  SUBROUTINE unsetmxs
    INTEGER(i1b) :: istat
    Deallocate(mS, STAT=istat)
    Deallocate(mX, STAT=istat)
    Deallocate(mXt, STAT=istat)
    Deallocate(mH, STAT=istat)
    Deallocate(mF, STAT=istat)
    Deallocate(mG, STAT=istat)
    Deallocate(mC, STAT=istat)
    Deallocate(mFp, STAT=istat)
    Deallocate(mCp, STAT=istat)
    Deallocate(mP, STAT=istat)
    Deallocate(mP_old, STAT=istat)
    Deallocate(mEig, STAT=istat)
    Deallocate(mVee, STAT=istat)
  END SUBROUTINE unsetmxs

END MODULE matrixs
