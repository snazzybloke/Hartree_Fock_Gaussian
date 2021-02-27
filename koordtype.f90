!###################################################
! Author: Ante Bilic                               !
! Since: Mar 30, 2011                              !
! Project: Hartree-Fock-GTO                        !
! Version: N/A                                     !
! Maintainer: Ante Bilic                           !
! Email: ante.bilic.mr@gmail.com                   !
! Status: N/A                                      !
!###################################################


MODULE ktype
! Defines the koord type as the chemical symbol, xyz coordinates and atomic number.
                                                                                
  USE abtype
  IMPLICIT NONE

  TYPE koord
    CHARACTER (LEN=2) :: clem
    REAL(dp), DIMENSION(3) :: xyz
    INTEGER(i1b) :: anum
  END TYPE koord

  CONTAINS

  FUNCTION f_dstnc2 (a1, a2)
! Calculates squared distance between atoms with koords a1 and a2
!
! Parameters:
! -----------
! a1 (koord):       the coordinates of the atom 1
! a2 (koord):       the coordinates of the atom 2
!
! Returns:
! ----------
! f_dstnc2 (int):    the square of the total distance (in a.u.) between 2 atoms

    IMPLICIT NONE
    TYPE(koord), INTENT(IN) :: a1, a2
    REAL(dp) :: f_dstnc2
    REAL(dp), DIMENSION(3) :: dx
    dx = a1%xyz(1:3) - a2%xyz(1:3)
    f_dstnc2 = Sum(dx * dx)
  END FUNCTION f_dstnc2

END MODULE ktype
