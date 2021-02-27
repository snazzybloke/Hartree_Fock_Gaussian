!###################################################
! Author: Ante Bilic                               !
! Since: Mar 30, 2011                              !
! Project: Hartree-Fock-GTO                        !
! Version: N/A                                     !
! Maintainer: Ante Bilic                           !
! Email: ante.bilic.mr@gmail.com                   !
! Status: N/A                                      !
!###################################################


MODULE ints
! Defines the matrices that store all the required integrals (overlap si, kinetic
! energy ti, 1-e Coulomb contributions v1ei, and 2-e Coulomb contributions v2ei)

  USE abtype
  IMPLICIT NONE
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: si, ti
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: v1ei
  REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: v2ei

  CONTAINS

  SUBROUTINE initints(na, msize)
! Simply allocates the memory space for the matrices and initializes them with zeros.
!
! Parameters:
! -----------
! na (int):       the number of atoms
! msize (int):    the size of the matrices

    INTEGER(i1b), INTENT(IN) :: na, msize
    INTEGER(i1b) :: istat
    Allocate(si(msize, msize), STAT=istat)
    Allocate(ti(msize, msize), STAT=istat)
    Allocate(v1ei(na, msize, msize), STAT=istat)
    Allocate(v2ei(msize, msize, msize, msize), STAT=istat)
    si = 0.0_dp
    ti = 0.0_dp
    v1ei = 0.0_dp
    v2ei = 0.0_dp
  END SUBROUTINE initints

  SUBROUTINE undefints
    INTEGER(i1b) :: istat
    Deallocate(si, STAT=istat)
    Deallocate(ti, STAT=istat)
    Deallocate(v1ei, STAT=istat)
    Deallocate(v2ei, STAT=istat)
  END SUBROUTINE undefints

END MODULE ints
