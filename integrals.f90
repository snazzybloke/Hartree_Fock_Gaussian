!###################################################
! Author: Ante Bilic                               !
! Since: Mar 30, 2011                              !
! Project: Hartree-Fock-GTO                        !
! Version: N/A                                     !
! Maintainer: Ante Bilic                           !
! Email: ante.bilic.mr@gmail.com                   !
! Status: N/A                                      !
!###################################################


SUBROUTINE integral(na, x, bname, msize)
! Calculates all the 1- and 2-electron integrals. This the most difficult part.
! The method of L.E. McMurchie  E.R. Davidson, J.Comp.Phys 26 (1978) is implemented below.
!
! Parameters:
! ----------
! na (int):          the total number of atoms
! x(1D koord array): atom chem symbol, xyz coordinates and atomic number of each atom
! bname (char):      basis set choice (e.g., STO-3G, 3-21G, 4-31G, 6-31G**)
! msize (int):       the size of the system matrices
!
! Evaluates:
! ----------
! si (2D real array):     the matrix of GTO overlaps (S)
! ti (2D real array):     the matrix of electronic kinetic energy operator (T)
! v1ei (3D real array):   the matrices of 1-electron Coulomb integrals (one for each atom, V_ne)
! v2ei (4D real array):   the contributions from the 2-electron Coulomb integrals (V_ee)

  USE abtype
  USE ktype
  USE ints
  USE bsets
  IMPLICIT NONE
  INTEGER(i1b), INTENT(IN) :: na, msize
  TYPE(koord), DIMENSION(na), INTENT(IN) :: x
  CHARACTER (LEN=7), INTENT(IN) :: bname
  INTEGER(i1b) :: i, iM, ii, i1, i2, i3, i3M, j, jM, jj, j1, j2, j3, j3M, k, kM, kk, k1, k2, k3, k3M, l, lM, ll, l1, l2, l3, l3M, istat
  INTEGER(i1b), DIMENSION(:,:), ALLOCATABLE :: inlm, jnlm, knlm, lnlm  ! exponents n, l, m of the cartesians: (x-A_x)^n, (y-A_y)^l, (z-A_z)^m
  REAL(dp) :: rAB2, rPC2, rCD2, rPQ2, alphaP, alphaQ, fint=0.0_dp
  REAL(dp), DIMENSION(3) :: xP, xQ
  TYPE(basis), DIMENSION(na) :: bas
  
! Normalize the basis functions:
  DO i=1, na
    bas(i) = defb(x(i)%anum, bname)
    CALL normalize(bas(i))
  END DO
  CALL initints(na, msize)
! Evaluate the 1-e integrals
! Loop over atoms (i), its basis functions (i1), and their primitives (i3), the rows of the matrices:
  ii = 0
  DO i=1, na  ! loop over atoms i
    iM = Size(bas(i)%nf)  ! the total # of basis functions for atom i
    DO i1=1, iM  ! loop over the GTO basis
      CALL cart_expos(bas(i)%lam(i1), inlm, i3M)
      DO i2=1, i3M  ! loop over px,py,pz (p orbitals); dxy, dyz, dzx, dxx, dyy, dzz (d orbitals)
        ii = ii + 1
        DO i3=1, bas(i)%nf(i1)  ! primitive Gaussians in the i1-th basis fcn on the atom i
          jj = 0
          ! Repeat the 3 loops above indexed with j, j1, and j3, the columns of the matrices:
          DO j=1, na
            jM = Size(bas(j)%nf)
            rAB2 = f_dstnc2(x(i), x(j))
            DO j1=1, jM
              CALL cart_expos(bas(j)%lam(j1), jnlm, j3M)
              DO j2=1, j3M
                jj = jj + 1
                DO j3=1, bas(j)%nf(j1)
! The matrices S, T, V are indexed via (ii, jj)
              ! fundamental overlap integral S=<s,A|s,B>:
                fint = sfun(bas(i)%a(i1)%rp(i3), bas(j)%a(j1)%rp(j3), rAB2) * bas(i)%d(i1)%rp(i3) * bas(j)%d(j1)%rp(j3)
              ! P point along the X(A)-X(B) line:
                xP = f_findP(x(i), x(j), bas(i)%a(i1)%rp(i3), bas(j)%a(j1)%rp(j3))
              ! The Gaussian exponent at the P point:
                alphaP = bas(i)%a(i1)%rp(i3) + bas(j)%a(j1)%rp(j3)
                  si(ii, jj) = si(ii, jj) + fint * def_NLM(inlm(i2, :), jnlm(j2, :), x(i)%xyz, x(j)%xyz, xP, alphaP)
                  ti(ii, jj) = ti(ii, jj) + fint * tkNLM(inlm(i2, :), jnlm(j2, :), bas(i)%a(i1)%rp(i3), bas(j)%a(j1)%rp(j3), x(i)%xyz, x(j)%xyz, xP, alphaP)
                  DO k=1, na  ! loop over atoms k, to calculate the Vne attraction <A|-Zc/R1c|B>
                    v1ei(k, ii, jj) = v1ei(k, ii, jj) + vfun(bas(i)%a(i1)%rp(i3), bas(j)%a(j1)%rp(j3), rAB2, x(k)%anum) * bas(i)%d(i1)%rp(i3) * bas(j)%d(j1)%rp(j3) * sumRNLMdef(inlm(i2, :), jnlm(j2, :), x(i)%xyz, x(j)%xyz, xP, x(k)%xyz, alphaP)
                  END DO
                END DO
              END DO
              Deallocate(jnlm, STAT=istat)
            END DO
          END DO
        END DO
      END DO
      Deallocate(inlm, STAT=istat)
    END DO
  END DO

! evaluate 2-e integrals
! Loop over atoms (i), its basis functions (i1), and their primitives (i3), the 1st dim out of 4:
  ii = 0
  DO i=1, na
   iM = Size(bas(i)%nf)
   DO i1=1, iM
   CALL cart_expos(bas(i)%lam(i1), inlm, i3M)
    DO i2=1, i3M
     ii = ii + 1
     DO i3=1, bas(i)%nf(i1)
      jj = 0
      ! Repeat the 3 loops above indexed with j, j1, and j3, the 2nd dim out of 4:
      DO j=1, na
       jM = Size(bas(j)%nf)
       rAB2 = f_dstnc2(x(i), x(j))
       DO j1=1, jM
        CALL cart_expos(bas(j)%lam(j1), jnlm, j3M)
        DO j2=1, j3M
         jj = jj + 1
         DO j3=1, bas(j)%nf(j1)
         ! The P point along the X(A)-X(B) line and its Gaussian exponent:
          xP = f_findP(x(i), x(j), bas(i)%a(i1)%rp(i3), bas(j)%a(j1)%rp(j3))
          alphaP = bas(i)%a(i1)%rp(i3) + bas(j)%a(j1)%rp(j3)
          kk = 0
          ! Repeat the 3 loops above indexed with k, k1, and k3, the 3nd dim out of 4:
          DO k=1, na
           kM = Size(bas(k)%nf)
           DO k1=1, kM
            CALL cart_expos(bas(k)%lam(k1), knlm, k3M)
            DO k2=1,k3M
             kk = kk + 1
             DO k3=1, bas(k)%nf(k1)
              ll = 0
              ! Repeat the 3 loops above indexed with l, l1, and l3, the 4th dim out of 4:
              DO l=1, na
               lM = Size(bas(l)%nf)
               rCD2 = f_dstnc2(x(k), x(l))
               DO l1=1, lM
                CALL cart_expos(bas(l)%lam(l1), lnlm, l3M)
                DO l2=1, l3M
                 ll = ll + 1
                 DO l3=1, bas(l)%nf(l1)
                 ! The Q point along the X(C)-X(D) line and its Gaussian exponent:
                  xQ = f_findP(x(k), x(l), bas(k)%a(k1)%rp(k3), bas(l)%a(l1)%rp(l3))
                  alphaQ = bas(k)%a(k1)%rp(k3) + bas(l)%a(l1)%rp(l3)
                  rPQ2 = Sum((xP-xQ)*(xP-xQ))
                  v2ei(ii,jj,kk,ll) = v2ei(ii,jj,kk,ll) + twoe(bas(i)%a(i1)%rp(i3), bas(j)%a(j1)%rp(j3), bas(k)%a(k1)%rp(k3), bas(l)%a(l1)%rp(l3), rAB2, rCD2, rPQ2) * bas(i)%d(i1)%rp(i3) * bas(j)%d(j1)%rp(j3) * bas(k)%d(k1)%rp(k3) * bas(l)%d(l1)%rp(l3) * sumRNLMdefx2(inlm(i2,:), jnlm(j2,:), knlm(k2,:), lnlm(l2,:), x(i)%xyz, x(j)%xyz, x(k)%xyz, x(l)%xyz, xP, xQ, alphaP, alphaQ)
                 END DO
                END DO
                Deallocate(lnlm, STAT=istat)
               END DO
              END DO
             END DO
            END DO
            Deallocate(knlm, STAT=istat)
           END DO
          END DO
         END DO
        END DO
        Deallocate(jnlm, STAT=istat)
       END DO
      END DO
     END DO
    END DO
    Deallocate(inlm, STAT=istat)
   END DO
  END DO


  CONTAINS


  SUBROUTINE cart_expos(lam, nlm, idx3M)
! Provides the exponents of the cartesian prefactors in the GTOs based on the orbital angular momentum.
!
! Parameters:
! ----------
! lam (char):          'S' or 'P' or 'D (for orbital angular momentum L=0, 1, 2, respectively)
! nlm (2D int array):  an allocatable 2D integer array
! idx3M (int):         any integer
!
! Returns:
! ----------
! nlm (1D int array):  the n, l, m exponents in (x-A_x)^n, (y-A_y)^l, (z-A_z)^m
! idx3M (int):         the number of the orbitals for the given angular momentum (usually 2L+1, but 6 for L=2)

    IMPLICIT NONE
    CHARACTER(LEN=1), INTENT(IN) :: lam  ! angular momentum l symbol ('S', 'P' or 'D)
    INTEGER(i1b), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: nlm  ! exponents n, l, m of the cartesians: (x-A_x)^n, (y-A_y)^l, (z-A_z)^m
    INTEGER(i1b), INTENT(OUT) :: idx3M  ! i2, j2, k2, l2 loop max index (3 for p orbitals, 6 for d orbitals)
    IF(lam == 'S') THEN   ! s-type orbital, angular momentum L=0
      idx3M = 1  ! (2L+1) projections of ang momentum L. Their x,y,z cartesian exponents:
      Allocate(nlm(1, 3), STAT=istat)
      nlm(1, :) = (/0, 0, 0/)
    ELSE IF(lam == 'P') THEN  ! p-type orbital, angular momentum L=1
      idx3M = 3  ! (2L+1) projections of ang momentum L. Their x,y,z cartesian exponents:
      Allocate(nlm(3, 3), STAT=istat)
      nlm(1, :) = (/1, 0, 0/)
      nlm(2, :) = (/0, 1, 0/)
      nlm(3, :) = (/0, 0, 1/)
    ELSE IF(lam == 'D') THEN  ! d-type orbital, angular momentum L=2
      idx3M = 6  ! 6 d orbitals (i.e., not 2L+1 = 5). Their x,y,z cartesian exponents:
      Allocate(nlm(6, 3), STAT=istat)
      nlm(1, :) = (/1, 1, 0/)
      nlm(2, :) = (/0, 1, 1/)
      nlm(3, :) = (/1, 0, 1/)
      nlm(4, :) = (/2, 0, 0/)
      nlm(5, :) = (/0, 2, 0/)
      nlm(6, :) = (/0, 0, 2/)
    END IF
  END SUBROUTINE cart_expos


  SUBROUTINE normalize(bz)
! Normalizes the contraction coefficients d (prefactors) of the basis set
!
! Parameters:
! ----------
! bz (basis):   GTO basis functions for a given atom 
!
! Returns:
! ----------
! bz (basis):   the same basis normalized, with scaled contraction coefficients d 

    IMPLICIT NONE
    TYPE(basis), INTENT(INOUT) :: bz
    INTEGER(i1b) :: i, j
    DO i=1, Size(bz%d)
      DO j=1, Size(bz%d(i)%rp)
        IF (bz%lam(i) == 'S') THEN          ! s-type Gaussian:
          bz%d(i)%rp(j) = bz%d(i)%rp(j) * (2.0_dp * bz%a(i)%rp(j) / pi)**0.75
        ELSE IF (bz%lam(i) == 'P') THEN     ! p-type Gaussian:
          bz%d(i)%rp(j) = bz%d(i)%rp(j) * (128.0_dp * bz%a(i)%rp(j)**5 / pi**3)**0.25
        ELSE IF (bz%lam(i) == 'D') THEN     ! d-type Gaussian:
          bz%d(i)%rp(j) = bz%d(i)%rp(j) * (2048.0_dp * bz%a(i)%rp(j)**7 / pi**3)**0.25
        END IF
      END DO
    END DO
  END SUBROUTINE normalize


  FUNCTION f_findP (x1, x2, al1, al2)
! The product of two GTOs on two centers can be replaced by a single GTO on a point P along the connecting line,
! see A. Szabo & N.S. Ostlund, "Modern Quantum Chemistry", Appendix A, Eq. (A.4).
!
! Parameters:
! ----------
! x1 (koord):  the coordinates of the 1st center
! x2 (koord):  the coordinates of the 2nd center
! al1 (real):  the alpha coefficent in the Gaussian exponent on the 1st center
! al2 (real):  the alpha coefficent in the Gaussian exponent on the 2nd center
!
! Returns:
! ----------
! f_findP (1D real array):  the xyz coordinates of the point P

    IMPLICIT NONE
    TYPE(koord), INTENT(IN) :: x1, x2
    REAL(dp), INTENT(IN) :: al1, al2
    REAL(dp), DIMENSION(3) :: f_findP
    f_findP = (al1 * x1%xyz + al2 * x2%xyz) / (al1 + al2)
  END FUNCTION f_findP


  FUNCTION sfun(a, b, rab2)
! Overlap between un-normalized Gaussian primitves
! A. Szabo & N.S. Ostlund, "Modern Quantum Chemistry", Appendix A, Eq. (A.9)
!
! Parameters:
! ----------
! a (real):     the exponent alpha from the 1st Gaussian
! b (real):     the exponent alpha from the 2nd Gaussian
! rab2 (real):  the squared distance between the two centres
!
! Returns:
! ----------
! sfun (real):  the overlap integral between two Gaussians

    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: a, b, rab2
    REAL(dp) :: sfun
    sfun = (pi / (a + b))**1.5_dp * Exp(-a * b * rab2 / (a + b))
  END FUNCTION sfun


  FUNCTION vfun(a, b, r2, z)
! Nuclear attraction integrals between un-normalized Gaussian primitves
! A. Szabo & N.S. Ostlund, "Modern Quantum Chemistry", Appendix A, Eq. (A.33)
!
! Parameters:
! ----------
! a (real):     the exponent alpha from the 1st Gaussian
! b (real):     the exponent alpha from the 2nd Gaussian
! r2 (real):  the squared distance between the two centres
! z (int):      the positive nuclear charge on the atom centre
!
! Returns:
! ----------
! vfun (real):  the 

    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: a, b, r2
    INTEGER(i1b), INTENT(IN) :: z
    REAL(dp) :: vfun
    vfun = -z * 2.0_dp * pi / (a + b) * Exp(-a * b * r2 / (a + b))
  END FUNCTION vfun


  FUNCTION ffun(m,x)
! the F_m function
    IMPLICIT NONE
    INTEGER(i1b), INTENT(IN) :: m
    REAL(dp), INTENT(IN) :: x
    REAL(dp) :: ffun, gmlw, gmup, tm
    INTEGER :: ifail

    INTERFACE
      FUNCTION s14aaf (x, ifail)
        USE abtype
        REAL(dp) :: a, x, tol, p, q, s14aaf
        INTEGER :: ifail
      END FUNCTION s14aaf

      SUBROUTINE s14baf (a, x, tol, p, q, ifail)
        USE abtype
        REAL(dp) :: a, x, tol, p, q
        INTEGER :: ifail
      END SUBROUTINE s14baf
    END INTERFACE

    IF (x < 1.0D-6) THEN
    ! lim F(m,x) = 1/2m+1 [1 - (2m+1)*x/(2m+3)]
      ffun = 1.0_dp / (2.0 * m + 1.0) * (1 - x * (2.0 * m + 1.0) / (2.0 * m + 3.0))
    ELSE IF (m == 0) THEN
      ffun = 0.5_dp * Sqrt(pi/x) * Erf(Sqrt(x))
    ELSE
    ! if m>0 one can use recursive formula:
    ! I. Shavitt, in "Methods of Computational Physics", Vol 2 (1963) p1, Eq(22)
    ! Fm(x) = 1/(2*x**(m+0.5)) * gamma(m+0.5, x). However, NOT ADVISABLE!!!
    ! call NAG S14BAF for incomplete gamma function (NORMALIZED, multiply by S14AAF)
      tm = m + 0.5_dp
      CALL s14baf(tm, x, 1.0D-16, gmlw, gmup, ifail)
      ffun = 0.5_dp / x**tm * gmlw * s14aaf(tm, ifail)
    END IF
  END FUNCTION ffun


  FUNCTION twoe(a, b, c, d, r2, p2, q2)
! 2-electron integrals for un-normalized Gaussian primitves
! a, b, c, d are exponents, r2, p2, q2 squared distances
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: a, b, c, d, r2, p2, q2
    REAL(dp) :: twoe
    twoe = 2.0_dp*(pi**2.5_dp)/((a+b)*(c+d)*Sqrt(a+b+c+d)) * Exp(-a*b*r2/(a+b)-c*d*p2/(c+d))
  END FUNCTION twoe


  FUNCTION def_NLM(l1, l2, x1, x2, xpc, alf)
! The last factor in the overlap between two basis functions.          
! McMurchie-Davidson formulas, J.Comp.Phys 26 (1978) 218:
! Eq's (2.26), (2.27)
!
! Parameters:
! ----------
! l1 (1D int array):     the cartesian exponents from the 1st atom centre
! l2 (1D int array):     the cartesian exponents from the 2nd atom centre
! x1 (1D real array):    the cartesian xyz coordinates of the 1st atom centre
! x2 (1D real array):    the cartesian xyz coordinates of the 2nd atom centre
! xpc (1D real array):   the cartesian xyz coordinates of the P point along the connecting line
! alf (real):            the exponent alpha from the Gaussian on the P point
!
! Returns:
! ----------
! def_NLM (real):        the final factor in the overlap integral

    IMPLICIT NONE
    INTEGER(i1b), DIMENSION(3), INTENT(IN) :: l1, l2
    REAL(dp), DIMENSION(3), INTENT(IN) :: x1, x2, xpc
    REAL(dp), INTENT(IN) :: alf
    REAL(dp) :: def_NLM 
    INTEGER(i1b) :: ixyz
    def_NLM = 1.0_dp
    DO ixyz=1, 3
      def_NLM = def_NLM * d_N(l1(ixyz), l2(ixyz), 0, x1, x2, xpc, alf, ixyz)
    END DO
  END FUNCTION def_NLM


  FUNCTION tkNLM(l1, l2, a1, a2, x1, x2, xpc, alf)
! The last factor in the kinetic energy between two basis functions.          
! McMurchie-Davidson formulas, J.Comp.Phys 26 (1978) 218:
! Eq's (2.35), (2.36), (2.37)
!
! Parameters:
! ----------
! l1 (1D int array):     the cartesian exponents from the 1st atom centre
! l2 (1D int array):     the cartesian exponents from the 2nd atom centre
! a1 (real):             the exponent alpha from the Gaussian on the 1st centre
! a2 (real):             the exponent alpha from the Gaussian on the 2nd centre
! x1 (1D real array):    the cartesian xyz coordinates of the 1st atom centre
! x2 (1D real array):    the cartesian xyz coordinates of the 2nd atom centre
! xpc (1D real array):   the cartesian xyz coordinates of the P point along the connecting line
! alf (real):            the exponent alpha from the Gaussian on the P point
!
! Returns:
! ----------
! tk_NLM (real):         the final factor in the kinetic energy integral

    IMPLICIT NONE
    INTEGER(i1b), DIMENSION(3), INTENT(IN) :: l1, l2
    REAL(dp), INTENT(IN) :: a1, a2, alf
    REAL(dp), DIMENSION(3), INTENT(IN) :: x1, x2, xpc
    REAL(dp) :: tkxyz, tkNLM
    REAL(dp), DIMENSION(3) :: def
    INTEGER(i1b) :: n1, n2, ixyz, jxyz, i1=1_i1b
    tkNLM = 0.0_dp
    DO ixyz=1, 3
      def(ixyz) = d_N(l1(ixyz), l2(ixyz), 0, x1, x2, xpc, alf, ixyz)
    END DO
    DO ixyz=1, 3
      n1 = l1(ixyz); n2 = l2(ixyz)
      tkxyz = n1 * n2 * d_N(n1-i1, n2-i1, 0, x1, x2, xpc, alf, ixyz) - 2.0 * n1 * a2 * d_N(n1-i1, n2+i1, 0, x1, x2, xpc, alf, ixyz) - 2.0 * n2 * a1 * d_N(n1+i1, n2-i1, 0, x1, x2, xpc, alf, ixyz) + 4.0 * a1 * a2 * d_N(n1+i1, n2+i1, 0, x1, x2, xpc, alf, ixyz)
      DO jxyz=1, 3
        IF (jxyz /= ixyz) THEN
          tkxyz = tkxyz * def(jxyz)
        END IF
      END DO
      tkNLM = tkNLM + tkxyz
    END DO
    tkNLM = 0.5_dp * tkNLM
  END FUNCTION tkNLM


  FUNCTION sumRNLMdef(l1, l2, x1, x2, xp, xc, alf)
! The last factor in the Vne 1-electron Coulomb energy between two basis functions.          
! McMurchie-Davidson formulas, J.Comp.Phys 26 (1978) 218:
! Eq's (3.9)-(3.14)
! V_nucl = - Z_c * 2 Pi/alphaP E_AB  Sigma_NLM  R_NLM d^nn_N e^ll_L f^mm_M
!
! Parameters:
! ----------
! l1 (1D int array):     the cartesian exponents from the 1st atom centre
! l2 (1D int array):     the cartesian exponents from the 2nd atom centre
! x1 (1D real array):    the cartesian xyz coordinates of the 1st atom centre
! x2 (1D real array):    the cartesian xyz coordinates of the 2nd atom centre
! xp (1D real array):    the cartesian xyz coordinates of the P point along the connecting line
! xc (1D real array):    the cartesian xyz coordinates of the nucleus
! alf (real):            the exponent alpha from the Gaussian on the P point
!
! Returns:
! ----------
! sumRNLMdef (real):     the final factor in the 1-e Coulomb energy integral

    IMPLICIT NONE
    INTEGER(i1b), DIMENSION(3), INTENT(IN) :: l1, l2
    REAL(dp), DIMENSION(3), INTENT(IN) :: x1, x2, xp, xc
    REAL(dp), INTENT(IN) :: alf
    REAL(dp) :: sumRNLMdef
    REAL(dp), DIMENSION(3) :: def
    INTEGER(i1b) :: kN, kL, kM
    sumRNLMdef = 0.0_dp
    DO kN = 0, l1(1) + l2(1)  ! sum over N, 0<=N<=n1+n2
      DO kL = 0, l1(2) + l2(2)  ! sum over L, 0<=L<=l1+l2
        DO kM = 0, l1(3) + l2(3)  ! sum over M, 0<=M<=m1+m2
          sumRNLMdef = sumRNLMdef + r_NLM(kN, kL, kM, xp-xc, alf) * d_N(l1(1), l2(1), kN, x1, x2, xp, alf, 1) * d_N(l1(2), l2(2), kL, x1, x2, xp, alf, 2) * d_N(l1(3), l2(3), kM, x1, x2, xp, alf, 3)
        END DO
      END DO
    END DO
  END FUNCTION sumRNLMdef


  FUNCTION r_NLM(kN, kL, kM, xpq, alf)
! An auxiliary function, called by sumRNLM() and sumRNLMdefx2(),
! which in turns calls r_NLMj(), thus reducing the code complexity.
! 
! Parameters:
! ----------
! kN (int):                the sum of the cartesian exponents of the x prefactor
! kL (int):                the sum of the cartesian exponents of the y prefactor
! kM (int):                the sum of the cartesian exponents of the z prefactor
! xpq (1D real array):     the xyz vector of the difference P-c (for 1-e), P-Q (for 2-e) 
! alf (real):              the exponent alpha from the Gaussian on the P point
!
! Returns:
! ----------
! r_NLM (real)             the value of the recursive function r_NLMj() with j=0

    IMPLICIT NONE
    INTEGER(i1b), INTENT(IN) :: kN, kL, kM
    REAL(dp), DIMENSION(3), INTENT(IN) :: xpq
    REAL(dp), INTENT(IN) :: alf
    REAL(dp) :: r_NLM
    r_NLM = r_NLMj(kN, kL, kM, 0, xpq, alf)
  END FUNCTION r_NLM


  RECURSIVE FUNCTION r_NLMj(kN, kL, kM, j, xpq, alf) RESULT(rN)
! Implements the recursive formulas for the Eq's (4.1) and (4.2) from
! McMurchie-Davidson, J.Comp.Phys 26 (1978) 218,
! for a general case, with the extra j index, Eq's (4.6), (4.7), and (4.8).
! This function is specifically for the Eq (4.8).
! 
! Parameters:
! ----------
! kN (int):                the sum of the cartesian exponents of the x prefactor
! kL (int):                the sum of the cartesian exponents of the y prefactor
! kM (int):                the sum of the cartesian exponents of the z prefactor
! j (int):                 the extra index added after Eq (4.1) for generalisation, here 0
! xpq (1D real array):     the xyz vector of the difference P-c (for 1-e), P-Q (for 2-e) 
! alf (real):              the exponent alpha from the Gaussian on the P point
!
! Returns:
! ----------
! rN (real)               the value of the recursive function r_NLMj() with j=0

    IMPLICIT NONE
    INTEGER(i1b), INTENT(IN) :: kN, kL, kM, j
    REAL(dp), DIMENSION(3), INTENT(IN) :: xpq
    REAL(dp), INTENT(IN) :: alf
    REAL(dp) :: rN
    INTEGER(i1b) :: i1=1_i1b, i2=2_i1b
    IF (kN<0) THEN
      rN = 0.0_dp
    ELSE IF (kN>0) THEN
! McMurchie-Davidson Eq 4.8
      rN = xpq(1) * r_NLMj(kN-i1, kL, kM, j+i1, xpq, alf) + (kN-1) * r_NLMj(kN-i2, kL, kM, j+i1, xpq, alf)
    ELSE IF (kN==0) THEN
      rN = r_0LMj(kL, kM, j, xpq, alf)
    END IF
  END FUNCTION r_NLMj


  RECURSIVE FUNCTION r_0LMj(kL, kM, j, xpq, alf) RESULT(rL)
! Implements the recursive formulas for the Eq's (4.1) and (4.2) from
! McMurchie-Davidson, J.Comp.Phys 26 (1978) 218,
! for a general case, with the extra j index.
! This function is specifically for the Eq (4.7), when N=0.
! 
! Parameters:
! ----------
! kL (int):                the sum of the cartesian exponents of the y prefactor
! kM (int):                the sum of the cartesian exponents of the z prefactor
! j (int):                 the extra index added after Eq (4.1) for generalisation, here 0
! xpq (1D real array):     the xyz vector of the difference P-c (for 1-e), P-Q (for 2-e) 
! alf (real):              the exponent alpha from the Gaussian on the P point
!
! Returns:
! ----------
! rL (real)                the value of the recursive function r_0LMj() with j=0

    IMPLICIT NONE
    INTEGER(i1b), INTENT(IN) :: kL, kM, j
    REAL(dp), DIMENSION(3), INTENT(IN) :: xpq
    REAL(dp), INTENT(IN) :: alf
    REAL(dp) :: rL
    INTEGER(i1b) :: i1=1_i1b, i2=2_i1b
    IF (kL<0) THEN
      rL = 0.0_dp
    ELSE IF (kL>0) THEN
! McMurchie-Davidson Eq 4.7
      rL =  xpq(2) * r_0LMj(kL-i1, kM, j+i1, xpq, alf) + (kL-1) * r_0LMj(kL-i2, kM, j+i1, xpq, alf)
    ELSE IF (kL==0) THEN
      rL =  r_00Mj(kM, j, xpq, alf)
    END IF
  END FUNCTION r_0LMj


  RECURSIVE FUNCTION r_00Mj(kM, j, xpq, alf) RESULT(rM)
! Implements the recursive formulas for the Eq's (4.1) and (4.2) from
! McMurchie-Davidson, J.Comp.Phys 26 (1978) 218,
! for a general case, with the extra j index.
! This function is specifically for the Eq (4.6), when N=L=0.
! 
! Parameters:
! ----------
! kL (int):                the sum of the cartesian exponents of the y prefactor
! kM (int):                the sum of the cartesian exponents of the z prefactor
! j (int):                 the extra index added after Eq (4.1) for generalisation, here 0
! xpq (1D real array):     the xyz vector of the difference P-c (for 1-e), P-Q (for 2-e) 
! alf (real):              the exponent alpha from the Gaussian on the P point
!
! Returns:
! ----------
! rM (real)                the value of the recursive function r_00Mj() with j=0

    IMPLICIT NONE
    INTEGER(i1b), INTENT(IN) :: kM, j
    REAL(dp), DIMENSION(3), INTENT(IN) :: xpq
    REAL(dp), INTENT(IN) :: alf
    REAL(dp) :: rM, rpq2
    INTEGER(i1b) :: i1=1_i1b, i2=2_i1b
    IF (kM<0) THEN
      rM = 0.0_dp
    ELSE IF (kM > 0) THEN
! McMurchie-Davidson Eq 4.6
      rM =  xpq(3) * r_00Mj(kM-i1, j+i1, xpq, alf) + (kM-1) * r_00Mj(kM-i2, j+i1, xpq, alf)
    ELSE IF (kM==0) THEN
      rpq2 = Sum(xpq * xpq)
      rM = (-2.0_dp*alf)**j * ffun(j, alf*rpq2)
    END IF
  END FUNCTION r_00Mj


  RECURSIVE FUNCTION d_N(n1, n2, N, x1, x2, xp, alf, ixyz) RESULT(dd)
! The recursion relations
! McMurchie-Davidson formulas, J.Comp.Phys 26 (1978) 218, Eq's 2.20, 2.21, and 2.22
!
! Parameters:
! ----------
! n1 (int):            the cartesian exponent from either x, y or z factor on the 1st atom centre
! n2 (int):            the cartesian exponent from either x, y or z factor on the 2nd atom centre
! N (int):             the subsript (the order of the corresponding Hermite polynomial, from 0 to the n1+n2)
! x1 (1D real array):  the cartesian xyz coordinates of the 1st atom centre
! x2 (1D real array):  the cartesian xyz coordinates of the 2nd atom centre
! xp (1D real array):  the cartesian xyz coordinates of the P point along the connecting line
! alf (real):          the exponent alpha from the Gaussian on the P point
! ixyz (int):          either 1, 2 or 3 for x, y, z, respectively
!
! Returns:
! ----------
! dd (real):        the result of the final recursion

    IMPLICIT NONE
    INTEGER(i1b), INTENT(IN) :: n1, n2, N, ixyz
    REAL(dp), DIMENSION(3), INTENT(IN) :: x1, x2, xp
    REAL(dp), INTENT(IN) :: alf
! also requires xP and alphaP defined in the calling subroutine initints
    REAL(dp) :: Pab, dd, xx1, xx2
    INTEGER(i1b) :: i1=1_i1b
    IF (n1<0 .OR. n2<0 .OR. N<0 .OR. N>(n1+n2)) THEN
      dd = 0.0_dp
    ELSE IF (n1==0 .AND. n2==0 .AND. N==0) THEN
      dd = 1.0_dp
    ELSE
      xx1 = x1(ixyz); xx2 = x2(ixyz)
      IF (n1 >= n2) THEN
        Pab = xp(ixyz) - xx1
        dd = 0.5_dp/alf * d_N(n1-i1, n2, N-i1, x1, x2, xp, alf, ixyz) + Pab * d_N(n1-i1, n2, N, x1, x2, xp, alf, ixyz) + (N+1) * d_N(n1-i1, n2, N+i1, x1, x2, xp, alf, ixyz)
      ELSE
        Pab = xp(ixyz) - xx2
        dd = 0.5_dp/alf * d_N(n1, n2-i1, N-i1, x1, x2, xp, alf, ixyz) + Pab * d_N(n1, n2-i1, N, x1, x2, xp, alf, ixyz) + (N+1) * d_N(n1, n2-i1, N+i1, x1, x2, xp, alf, ixyz)
      END IF
    END IF
  END FUNCTION d_N


  FUNCTION sumRNLMdefx2(l1, l2, l3, l4, x1, x2, x3, x4, xp, xq, ap, aq)
! The last factor in the Vee 2-electron Coulomb energy between four basis functions.          
! McMurchie-Davidson formulas, J.Comp.Phys 26 (1978) 218:
! Eq's (3.29)-(3.34)
! V_ee ~ Sum_NLMN'L'M'  R_N+N'L+L'M+M' d^nn_N e^ll_L f^mm_M d^nn_N' e^ll_L' f^mm_M'
!
! Parameters:
! ----------
! l1 (1D int array):     the cartesian exponents from the 1st atom centre
! l2 (1D int array):     the cartesian exponents from the 2nd atom centre
! l3 (1D int array):     the cartesian exponents from the 3rd atom centre
! l4 (1D int array):     the cartesian exponents from the 4rd atom centre
! x1 (1D real array):    the cartesian xyz coordinates of the 1st atom centre
! x2 (1D real array):    the cartesian xyz coordinates of the 2nd atom centre
! x3 (1D real array):    the cartesian xyz coordinates of the 3rd atom centre
! x4 (1D real array):    the cartesian xyz coordinates of the 4th atom centre
! xp (1D real array):    the cartesian xyz coordinates of the P point along the connecting line
! xq (1D real array):    the cartesian xyz coordinates of the Q point along the connecting line
! ap (real):             the exponent alpha from the Gaussian on the P point
! aq (real):             the exponent alpha from the Gaussian on the Q point
!
! Returns:
! ----------
! sumRNLMdefx2 (real):   the final factor in the 2-e Coulomb energy integral

    IMPLICIT NONE
    INTEGER(i1b), DIMENSION(3), INTENT(IN) :: l1, l2, l3, l4
    REAL(dp), DIMENSION(3), INTENT(IN) :: x1, x2, x3, x4, xp, xq
    REAL(dp), INTENT(IN) :: ap, aq
    REAL(dp) :: tkxyz, sumRNLMdefx2
    REAL(dp), DIMENSION(3) :: def
    INTEGER(i1b) :: kN, kL, kM, lN, lL, lM
    sumRNLMdefx2 = 0.0_dp
    DO kN = 0, l1(1) + l2(1)  ! sum over N, 0<=N<=n1+n2
     DO kL = 0, l1(2) + l2(2)  ! sum over L, 0<=L<=l1+l2
      DO kM = 0, l1(3) + l2(3)  ! sum over M, 0<=M<=m1+m2
       DO lN = 0, l3(1) + l4(1)  ! sum over N', 0<=N'<=n'1+n'2
        DO lL = 0, l3(2) + l4(2)  ! sum over L', 0<=L'<=l'1+l'2
         DO lM = 0, l3(3) + l4(3)  ! sum over M', 0<=M'<=m'1+m'2
          sumRNLMdefx2 = sumRNLMdefx2 + (-1)**(lN+lL+lM) * r_NLM(kN+lN, kL+lL, kM+lM, xp-xq, ap*aq/(ap+aq)) * d_N(l1(1), l2(1), kN, x1, x2, xp, ap, 1) * d_N(l1(2), l2(2), kL, x1, x2, xp, ap, 2) * d_N(l1(3), l2(3), kM, x1, x2, xp, ap, 3) * d_N(l3(1), l4(1), lN, x3, x4, xq, aq, 1) * d_N(l3(2), l4(2), lL, x3, x4, xq, aq, 2) * d_N(l3(3), l4(3), lM, x3, x4, xq, aq, 3)
         END DO
        END DO
       END DO
      END DO
     END DO
    END DO
  END FUNCTION sumRNLMdefx2

END SUBROUTINE integral
