!###################################################
! Author: Ante Bilic                               !
! Since: Mar 30, 2011                              !
! Project: Hartree-Fock-GTO                        !
! Version: N/A                                     !
! Maintainer: Ante Bilic                           !
! Email: ante.bilic.mr@gmail.com                   !
! Status: N/A                                      !
!###################################################


SUBROUTINE collect(na, msize)
! Combines the kinetic energy ti and 1-e v1ei Coulomb energy into the core Hamiltonian matrix mH.
! Forms the orthogonal basis sets of orbitals by diagonalizing the overlap matrix si (mS from here on)
! via canonical orthogonalization.
!
! Parameters:
! ----------
! na (int):     the number of atoms
! msize (int):  the size of the matrices
!
! Evaluates:
! ----------
! mH (2D real array):  the transpose of the transformation matrix
! mS (2D real array):  the orbital overlap matrix (a copy of si)
! mX (2D array):  the transformation matrix
! mXt (2D array): the transpose of the transformation matrix

  USE abtype
  USE ints
  USE matrixs
  USE lapack95       ! from MKL LAPACK
! USE f95_precision  ! from MKL LAPACK
  IMPLICIT NONE
  INTEGER(i1b), INTENT(IN) :: na, msize
! REAL(dp), DIMENSION(msize*(msize+1)/2,msize) :: mTmp
  REAL(dp), DIMENSION(msize,msize) :: mTmp, ms_12
  REAL(dp), DIMENSION(msize) :: vTmp
  INTEGER(i1b) :: nsize, i, j, k, l

  INTERFACE
    FUNCTION v2diag(msize, vect)
      USE abtype
      INTEGER(i1b), INTENT(IN) :: msize
      REAL(dp), DIMENSION(:), INTENT(IN) :: vect
      REAL(dp), DIMENSION(msize, msize) :: v2diag
    END FUNCTION v2diag
  END INTERFACE

  CALL setmxs(msize)
! form the core Hamiltonian matrix mH, A. Szabo & N.S. Ostlund, "Modern Quantum Chemistry", Eq. (3.153):
  mH = ti
  DO i=1, na
    mH = mH + v1ei(i,:,:)
  END DO
! form the overlap matrix
  mS = si
  mTmp = mS
  CALL syev(mTmp, vTmp, 'V')  ! SYEV LAPACK routine
  CALL matprint(mTmp, msize, msize, msize, msize, "U ")
  PRINT *,""
  PRINT *, "Eigenvalues of S: ", vTmp
  vTmp = 1.0_dp/Sqrt(vTmp)
  ms_12 = v2diag(msize, vTmp) ! turn the vector into a diagonal matrix
  CALL matprint(ms_12, msize, msize, msize, msize, "s^-1/2 ")
! use the canonical orthogonalization
  mX = Matmul(mTmp, Matmul(ms_12, Transpose(mTmp)))
  mXt = Transpose(mX)

  mVee = v2ei

  CALL matprint(mS, msize, msize, msize, msize, "S ")
  CALL matprint(mX, msize, msize, msize, msize, "X ")
  CALL matprint(mH, msize, msize, msize, msize, "H ")
  WRITE(*, *)
  CALL undefints
END SUBROUTINE collect
