!###################################################
! Author: Ante Bilic                               !
! Since: Mar 30, 2011                              !
! Project: Hartree-Fock-GTO                        !
! Version: N/A                                     !
! Maintainer: Ante Bilic                           !
! Email: ante.bilic.mr@gmail.com                   !
! Status: N/A                                      !
!###################################################


SUBROUTINE scf(na, nel, msize, enucl)
! Carries out the iterative (self-consistent) calculations to solve the Roothaan-Hall eqs.
! The loop stops if the density matrix mP has converged or after a max number of allowed iterations.
! Parameters:
! ----------
! na (int):     the number of atoms
! nel (int):    the total number of electrons
! msize (int):  the size of the matrices
! enucl (real): the Coulomb energy contribution from the nuclear charges
! Evaluates:
! ----------
! mP (2D real array):  the density matrix
! mP_old (2D real array):  the density matrix from the previous iteration,
!                          eventually used to store the Mulliken charges
! mG (2D real array):  the matrix of 2-e electrostatic contributions (direct & exchange)
! mF (2D real array):  the Fock operator matrix (mH + mG)
! mFp (2D real array):  the Fock operator in the transformed orthogonal orbital basis
! mCp (2D real array):  the orbital coefficients (the solution) in the orthogonal orbital basis
! mC (2D real array):  the orbital coefficients (the solution) in the original orbital basis
! mEig (2D real array):  the diagonal matrix of MO energies

  USE abtype
  USE matrixs
  IMPLICIT NONE
  INTEGER(i1b), INTENT(IN) :: na, nel, msize
  REAL(dp), INTENT(IN) :: enucl
  INTEGER(i1b) :: i, j, k, ims, istat, mindx, iter
  INTEGER(i1b), SAVE :: iatomcnt = 1
  REAL(dp) :: en = 0.0_dp, en_tot, deltaP = 0.0_dp
  REAL(dp), PARAMETER :: convcrit = 1.0D-6
  INTEGER(i1b), PARAMETER :: maxiter = 100
  REAL(dp), DIMENSION(:, :), ALLOCATABLE :: mtempDM  ! temporarily store the last atom DM
  REAL(dp), DIMENSION(:, :), ALLOCATABLE, SAVE :: matomsDM  ! merged DM's from atomic calculations

  iter = 0  ! reset for each calculation
  mP = 0.0_dp  ! initial guess zero, pretty bad, used only for individual atoms
  IF (na > 1_i1b) THEN  ! molecular calculation
    PRINT *, "********************* MOLECULAR CALCULATION *********************"
    mP = matomsDM  ! an improved initial guess: taken from individual atomic Density Matrices
  END IF
  CALL matprint(mP, msize, msize, msize, msize, "P")

! start of SCf iterations:
scf_loop:  DO
    iter = iter + 1
    WRITE(*,'(/,4X,"Start of iteration number = ", I3)') iter
! for 2-electron part of the Fock matrix mG:
    CALL formG
    CALL matprint(mG, msize, msize, msize, msize, "G")
    mF = mH + mG
! electronic energy, A. Szabo & N.S. Ostlund, "Modern Quantum Chemistry", Eq. (3.184)
    en = 0.5_dp * Sum(mP * (mH + mF))
    CALL matprint(mF, msize, msize, msize, msize, "F")
    WRITE(*,'(///,4X,"Electronic energy = ", E20.12)') en
! transformed Fock matrix F' = X^+ F X
    mFp = Matmul(mXt, Matmul(mF, mX))
! diagonalize transformed Fock matrix
    CALL diag(mFp, mCp, mEig, msize)
! transformed eigenvectors
    mC = Matmul(mX, mCp)
! form the NEW density matrix, A. Szabo & N.S. Ostlund, "Modern Quantum Chemistry", Eq. (3.145)
! save the old one first
    mP_old = mP
    DO i=1, msize
      DO j=1, msize
        mP(i,j) = 0.0_dp
        DO k=1, nel/2
          mP(i,j) = mP(i,j) + 2.0_dp * mC(i,k) * mC(j,k)
        END DO
      END DO
    END DO
    CALL matprint(mFp, msize, msize, msize, msize, "F' ")
    CALL matprint(mCp, msize, msize, msize, msize, "C' ")
    CALL matprint(mEig, msize, msize, msize, msize, "E ")
    CALL matprint(mC, msize, msize, msize, msize, "C ")
    CALL matprint(mP, msize, msize, msize, msize, "P ")
! calculate deltaP
    deltaP = Sqrt(Sum((mP - mP_old)*(mP - mP_old))/4.0_dp)
    WRITE(*,'(/,4X,"DeltaP (convergence of density matrix) = ",F10.6)') deltaP
    IF (deltaP > convcrit) THEN
! not yet converged
! test for maximum number of iterations
      IF (iter >= maxiter) THEN
        WRITE(*,'(4X,"No convergence in SCF")')
        EXIT scf_loop
      END IF
    ELSE
      EXIT scf_loop
    END IF
  END DO scf_loop
! calculation converged if it reaches here
! add nuclear repulsion
  en_tot = en + enucl
  WRITE(*,'(//,4X,"Calculation converged",//,4X,"Electronic energy = ",E20.12,//,4X,"Total energy = ",E20.12)') en, en_tot
! print out the final results
  CALL matprint(mG, msize, msize, msize, msize, "G ")
  CALL matprint(mF, msize, msize, msize, msize, "F ")
  CALL matprint(mEig, msize, msize, msize, msize, "E ")
  CALL matprint(mC, msize, msize, msize, msize, "C ")
  CALL matprint(mP, msize, msize, msize, msize, "P ")
! Mulliken population, A. Szabo & N.S. Ostlund, "Modern Quantum Chemistry", Eq. (3.195)
  mP_old = Matmul(mP, mS)
  CALL matprint(mP_old, msize, msize, msize, msize, "PS ")
  IF (iatomcnt == 1_i1b) THEN  ! firs atom
    Allocate(matomsDM(msize,msize), STAT=istat)
    matomsDM = mP
  ELSE IF (na == 1_i1b) THEN  ! if this was an individul atom calculations, save the Density Matrix
    ims = Size(matomsDM(:,1))
    Allocate(mtempDM(ims, ims),  STAT=istat)
    mtempDM = matomsDM  ! save the contents in mtempDM
    Deallocate(matomsDM, STAT=istat)
    mindx=ims+msize
    Allocate(matomsDM(mindx,mindx), STAT=istat) ! build a new, bigger, block-diagonal DM:
    matomsDM(1:ims, 1:ims) = mtempDM
    matomsDM(ims+1:mindx, ims+1:mindx) = mP
    matomsDM(ims+1:mindx, 1:ims) = 0.0_dp
    matomsDM(1:ims, ims+1:mindx) = 0.0_dp
    Deallocate(mtempDM, STAT=istat)
  END IF
  iatomcnt = iatomcnt + 1
  CALL unsetmxs

  CONTAINS

  SUBROUTINE formG
! calculates G matrix from the density matrix and 2-electron integrals
! A. Szabo & N.S. Ostlund, "Modern Quantum Chemistry", Eq. (3.154)
    IMPLICIT NONE
    INTEGER(i1b) :: i, j, k, l
    DO i=1, msize
      DO j=1, msize
        mG(i,j) = 0.0_dp
        DO k=1, msize
          DO l=1, msize
            mG(i, j) = mG(i, j) + mP(k, l) * (mVee(i, j, k, l) - 0.5_dp * mVee(i, l, k, j))
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE formG
  
END SUBROUTINE scf
