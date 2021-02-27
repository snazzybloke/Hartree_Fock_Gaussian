!###################################################
! Author: Ante Bilic                               !
! Since: Mar 30, 2011                              !
! Project: Hartree-Fock-GTO                        !
! Version: N/A                                     !
! Maintainer: Ante Bilic                           !
! Email: ante.bilic.mr@gmail.com                   !
! Status: N/A                                      !
!###################################################


SUBROUTINE diag(mF, mC, mE, msize)
! Diagonalizes mF to return the egienvectors in mC and eigenvalues in mE.
!
! Parameters:
! ----------
! mF (2D real array):  the matrix of the Fock operator in the orthogonal basis
! mC (2D real array):  any msize*msize matrix on the input
! mE (2D real array):  any msize*msize matrix on the input
! msize (int):         the size of the matrices
!
! Evaluates:
! ----------
! mC (2D real array):  the orbital coefficients (in the orthogonal basis)
! mE (2D real array):  the MO energies in a diagonal matrix

  USE abtype
  USE lapack95
  IMPLICIT NONE
  INTEGER(i1b), INTENT(IN) :: msize
  REAL(dp), DIMENSION(msize,msize), INTENT(IN) :: mF
  REAL(dp), DIMENSION(msize,msize), INTENT(OUT) :: mC, mE
  REAL(dp), DIMENSION(msize) :: vTmp

  INTERFACE
    FUNCTION v2diag(msize, vect)
      USE abtype
      INTEGER(i1b), INTENT(IN) :: msize
      REAL(dp), DIMENSION(:), INTENT(IN) :: vect
      REAL(dp), DIMENSION(msize,msize) :: v2diag
    END FUNCTION v2diag
  END INTERFACE

  mC = mF
  CALL syev(mC, vTmp, 'V')  ! LAPACK SYEV routine
  mE = v2diag(msize, vTmp)

END SUBROUTINE diag


FUNCTION v2diag(msize, vect)
! Adapted from Martin Counihan "Fortran 95", p 88 (CORRECTED!)
! forming a diagonal matrix from the elements of a vector
!
! Parameters:
! ----------
! msize (int):           the size of the matrices
! vect (1D real array):  a real 1D vector
!
! Evaluates:
! ----------
! v2diag (2D real array): the vector as a diagonal matrix

  USE abtype
  IMPLICIT NONE
  INTEGER(i1b), INTENT(IN) :: msize
  REAL(dp), DIMENSION(:), INTENT(IN) :: vect
  REAL(dp), DIMENSION(msize,msize) :: v2diag
  INTEGER(i1b) :: i
  LOGICAL, DIMENSION(msize) :: v1, v2
  LOGICAL, DIMENSION(msize,msize) :: m1, m2
  REAL(dp), DIMENSION(msize,msize) :: m3
! make a column vector v1, all FALSE:
  v1 = .FALSE.
! shift elements from the bottom by 1, add TRUE at the bottom:
  v2 = Eoshift(ARRAY=v1, DIM=1, SHIFT=1, BOUNDARY=.TRUE.)
! copy the new column vector2 msize times to make a matrix m1:
  m1 = Spread(SOURCE=v2, DIM=2, NCOPIES=msize)
! shift each column by step -i, i=1, msize, so that TRUE sit along the diagonal
  m2 = Cshift(ARRAY=m1, DIM=1, SHIFT=(/ (-i,i=1, msize) /))
! m2 = Cshift(ARRAY=m1, DIM=1, SHIFT=(/ (msize-i,i=1, msize) /)) ! also works
! matrix m3 all zeros:
  m3 = 0.0_dp
! diagonal matrix (where m2=TRUE, otherwise zero):
  v2diag = Unpack(VECTOR=vect, MASK=m2, FIELD=m3)
END FUNCTION v2diag


SUBROUTINE matprint(mA, jm, jn, m, n, label)
! Prints the matrix mA in the desired format.
!
! Parameters:
! ----------
! mA (2D real array):  a real 2D matrix
! jm (int):            the number of rows in the matrix
! jn (int):            the number of columns in the matrix
! m (int):             the number of rows to print
! n (int):             the number of columns to print
! label (char):        the string description of the matrix to print

  USE abtype
  IMPLICIT NONE
  INTEGER(i1b), INTENT(IN) :: jm, jn, m, n
  CHARACTER(LEN=*), INTENT(IN) :: label
  REAL(dp), DIMENSION(jm, jn), INTENT(IN) :: mA
  INTEGER(i1b) :: iHigh, iLow, i, j
  iHigh = 0
  DO
    iLow = iHigh + 1
    iHigh = iHigh + 10
    iHigh = Min(iHigh, n)
    WRITE (*,'(//,3X," The ",A," array",/,15X,10(5X,I3,6X)//)') label, (i, i=iLow, iHigh)
    DO i=1, m
      WRITE (*,'(I6,3X,10(2X,E12.6))') i, (mA(i, j), j=iLow, iHigh)
    END DO
    IF (n <= iHigh) EXIT
  END DO
END SUBROUTINE matprint
