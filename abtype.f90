!###################################################
! Author: Ante Bilic                               !
! Since: Mar 30, 2011                              !
! Project: Hartree-Fock-GTO                        !
! Version: N/A                                     !
! Maintainer: Ante Bilic                           !
! Email: ante.bilic.mr@gmail.com                   !
! Status: N/A                                      !
!###################################################


MODULE abtype
! Simply defines the integer and real kinds used by the code
! and some common constants.

  IMPLICIT NONE
! INTEGER, PARAMETER :: i1b = Selected_Int_Kind(2)
  INTEGER, PARAMETER :: i1b = Selected_Int_Kind(4)
  INTEGER, PARAMETER :: i2b = Selected_Int_Kind(4)
  INTEGER, PARAMETER :: i4b = Selected_Int_Kind(9)
  INTEGER, PARAMETER :: sp = Kind(1.0)
  INTEGER, PARAMETER :: dp = Kind(1.0D0)

  REAL(dp), PARAMETER :: pi=3.141592653589793238462643383279502884197D0

! REAL(dp), PARAMETER :: a0_d=0.529177249D0   !Angstroms
! REAL(dp), PARAMETER :: rd2deg=57.29578D0
END MODULE abtype
