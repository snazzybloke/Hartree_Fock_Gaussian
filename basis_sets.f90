!###################################################
! Author: Ante Bilic                               !
! Since: Mar 30, 2011                              !
! Project: Hartree-Fock-GTO                        !
! Version: N/A                                     !
! Maintainer: Ante Bilic                           !
! Email: ante.bilic.mr@gmail.com                   !
! Status: N/A                                      !
!###################################################


MODULE bsets
! Provides the definitions of the basis sets
                                                                                
  USE abtype
  IMPLICIT NONE

  TYPE int_pntr
! To contruct an array of pointers.
! It's needed for a 2-D array whose dimensions vary
! varies from basis set to another,
! so array of pointers is better than a 2-D allocatable array.
    INTEGER(i1b), DIMENSION(:), POINTER :: ip
  END TYPE int_pntr

  TYPE real_pntr
    REAL(dp), DIMENSION(:), POINTER :: rp
  END TYPE real_pntr

  TYPE basis
    INTEGER(i1b) :: anum  ! atomic number
    CHARACTER (LEN=7) :: bset  ! name, e.g. STO-3G, 3-21G, 4-31G, 6-31G**
    INTEGER (i1b), DIMENSION(:), ALLOCATABLE :: nf  ! number of basis functions
    TYPE(real_pntr), DIMENSION(:), ALLOCATABLE :: d, a  ! contraction coeff's (d) and exponents (a) of the primitive orbitals
    CHARACTER(LEN=1), DIMENSION(:), ALLOCATABLE :: lam  ! angular momentum l
  END TYPE basis

  CONTAINS

  FUNCTION defb(atn, bname)
! Calculates squared distance between atoms with koords a1 and a2
!
! Parameters:
! -----------
! atn (int):       the atomic number
! bname (char):    the desired GTO basis
!
! Returns:
! ----------
! defb (basis):    the square of the total distance (in a.u.) between 2 atoms

    IMPLICIT NONE
    INTEGER(i1b), INTENT(IN) :: atn
    CHARACTER (LEN=7), INTENT(IN) :: bname
    TYPE(basis) :: defb
    INTEGER(i1b) :: istat
    defb%anum = atn
    defb%bset = bname
    IF (atn <= 2) THEN   ! H and He
      IF (bname == 'STO-3G ') THEN
        Allocate(defb%nf(1), STAT=istat)  ! 1 basis function for He and He
        defb%nf(1) = 3    ! # of primitive gaussians
        Allocate(defb%lam(1),STAT=istat)
        defb%lam(1) = 'S'    ! s-type, l=0
        Allocate(defb%d(1), STAT=istat) ; Allocate(defb%d(1)%rp(3), STAT=istat)
        Allocate(defb%a(1), STAT=istat) ; Allocate(defb%a(1)%rp(3), STAT=istat)
! angular momentum for each Gaussian:
        defb%d(1)%rp = (/ 0.15432897_dp, 0.53532814_dp, 0.44463454_dp /)
        IF (atn == 1) THEN ! Hydrogen
          defb%a(1)%rp = (/ 3.42525091_dp, 0.62391373_dp, 0.16885540_dp /)
        ELSE  ! He
          defb%a(1)%rp = (/ 6.36242139_dp, 1.15892300_dp, 0.31364979_dp /)
        END IF
      ELSE IF (bname == '3-21G ') THEN
        Allocate(defb%nf(2), STAT=istat)  ! 2 basis functions for H and He
        defb%nf(1) = 2; defb%nf(2) = 1 ! (3)-21G
        Allocate(defb%lam(2), STAT=istat) ! angular momentum for each gaussian
        defb%lam(1) = 'S'; defb%lam(2) = 'S' ! s-type, l=0
        Allocate(defb%d(2), STAT=istat)
        Allocate(defb%d(1)%rp(2),STAT=istat); Allocate(defb%d(2)%rp(1),STAT=istat)
        Allocate(defb%a(2), STAT=istat)
        Allocate(defb%a(1)%rp(2), STAT=istat); Allocate(defb%a(2)%rp(1),STAT=istat)
        IF (atn == 1) THEN ! Hydrogen
          defb%d(1)%rp = (/ 0.1562850_dp, 0.9046910_dp /)
          defb%a(1)%rp = (/ 5.4471780_dp, 0.8245470_dp /)
          defb%d(2)%rp = (/ 1.0_dp /)
          defb%a(2)%rp = (/ 0.1831920_dp /)
        ELSE  ! He
          defb%d(1)%rp = (/ 0.1752300_dp, 0.8934830_dp /)
          defb%a(1)%rp = (/ 13.6267000_dp, 1.9993500_dp /)
          defb%d(2)%rp = (/ 1.0_dp /)
          defb%a(2)%rp = (/ 0.3829930_dp /)
        END IF
      ELSE IF (bname == '4-31G ') THEN
        Allocate(defb%nf(2), STAT=istat)  ! 2 basis functions for He and He
        defb%nf(1) = 3; defb%nf(2) = 1  ! (4)-31G
        Allocate(defb%lam(2), STAT=istat) ! angular momentum for each gaussian
        defb%lam(1) = 'S'; defb%lam(2) = 'S' ! s-type, l=0
        Allocate(defb%d(2), STAT=istat)
        Allocate(defb%d(1)%rp(3),STAT=istat); Allocate(defb%d(2)%rp(1),STAT=istat)
        Allocate(defb%a(2), STAT=istat)
        Allocate(defb%a(1)%rp(3), STAT=istat); Allocate(defb%a(2)%rp(1),STAT=istat)
        IF (atn == 1) THEN ! Hydrogen
          defb%d(1)%rp = (/ 0.0334946_dp, 0.2347269_dp, 0.8137573_dp /)
          defb%a(1)%rp = (/ 18.7311370_dp, 2.8253944_dp, 0.6401217_dp /)
          defb%d(2)%rp = (/ 1.0_dp /)
          defb%a(2)%rp = (/ 0.1612778_dp /)
        ELSE  ! He
          defb%d(1)%rp = (/ 0.0237660_dp, 0.1546790_dp, 0.4696300_dp /)
          defb%a(1)%rp = (/ 38.4216340_dp, 5.7780300_dp, 1.2417740_dp /)
          defb%d(2)%rp = (/ 1.0_dp /)
          defb%a(2)%rp = (/ 0.2979640_dp /)
        END IF
      ELSE IF (bname == '6-31G**') THEN
        Allocate(defb%nf(3), STAT=istat)  !3 basis functions (5: 1s,2s,px,py,pz)
        defb%nf(1) = 3; defb%nf(2) = 1; defb%nf(3) = 1 ! (6)-31G, p polarization
        Allocate(defb%lam(3), STAT=istat) ! angular momentum for each gaussian
        defb%lam(1) = 'S'; defb%lam(2) = 'S'; defb%lam(3) = 'P'
        Allocate(defb%d(3), STAT=istat)
        Allocate(defb%d(1)%rp(3),STAT=istat); Allocate(defb%d(2)%rp(1),STAT=istat)
        Allocate(defb%d(3)%rp(1),STAT=istat)
        Allocate(defb%a(3), STAT=istat)
        Allocate(defb%a(1)%rp(3), STAT=istat); Allocate(defb%a(2)%rp(1),STAT=istat)
        Allocate(defb%a(3)%rp(1),STAT=istat)
        IF (atn == 1) THEN ! Hydrogen
          defb%d(1)%rp = (/ 0.03349460_dp, 0.23472695_dp, 0.81375733_dp /)
          defb%a(1)%rp = (/ 18.7311370_dp, 2.8253937_dp, 0.6401217_dp /)
          defb%d(2)%rp = (/ 1.0_dp /)
          defb%a(2)%rp = (/ 0.1612778_dp /)
          defb%d(3)%rp = (/ 1.0_dp /)
          defb%a(3)%rp = (/ 1.1_dp /)
        ELSE  ! He
          defb%d(1)%rp = (/ 0.0237660_dp, 0.1546790_dp, 0.4696300_dp /)
          defb%a(1)%rp = (/ 38.4216340_dp, 5.7780300_dp, 1.2417740_dp /)
          defb%d(2)%rp = (/ 1.0_dp /)
          defb%a(2)%rp = (/ 0.2979640_dp /)
          defb%d(3)%rp = (/ 1.0_dp /)
          defb%a(3)%rp = (/ 1.1_dp /)
        END IF
      END IF
    ELSE IF (atn <= 10) THEN   ! Li to Ne
      IF (bname == 'STO-3G ') THEN
        Allocate(defb%nf(3), STAT=istat)  ! 5 basis functions, 1s,2s,2p(x y z)
        defb%nf(1) = 3; defb%nf(2) = 3; defb%nf(3) = 3 !# of primitive gaussians
        Allocate(defb%lam(3),STAT=istat)
        defb%lam(1) = 'S'; defb%lam(2) = 'S'; defb%lam(3) = 'P' ! s-type, l=0
        Allocate(defb%d(3), STAT=istat)
        Allocate(defb%d(1)%rp(3), STAT=istat)
        Allocate(defb%d(2)%rp(3), STAT=istat)
        Allocate(defb%d(3)%rp(3), STAT=istat)
        Allocate(defb%a(3), STAT=istat)
        Allocate(defb%a(1)%rp(3), STAT=istat)
        Allocate(defb%a(2)%rp(3), STAT=istat)
        Allocate(defb%a(3)%rp(3), STAT=istat)
! angular momentum for each Gaussian:
        defb%d(1)%rp = (/ 0.15432897_dp, 0.53532814_dp, 0.44463454_dp /)
        defb%d(2)%rp = (/ -0.099967237_dp, 0.39951283_dp, 0.70011547_dp /)
        defb%d(3)%rp = (/ 0.15591627_dp, 0.60768372_dp, 0.39195739_dp /)
        IF (atn == 6) THEN ! C
          defb%a(1)%rp = (/ 71.6168370_dp, 13.0450960_dp, 3.5305122_dp /)
          defb%a(2)%rp = (/ 2.9412494_dp, 0.6834831_dp, 0.2222899_dp /)
          defb%a(3)%rp = (/ 2.9412494_dp, 0.6834831_dp, 0.2222899_dp /)
        ELSE IF (atn == 7) THEN ! N
          defb%a(1)%rp = (/ 99.1061690_dp, 18.0523120_dp, 4.8856602_dp /)
          defb%a(2)%rp = (/ 3.7804559_dp, 0.8784966_dp, 0.2857144_dp /)
          defb%a(3)%rp = (/ 3.7804559_dp, 0.8784966_dp, 0.2857144_dp /)
        ELSE IF (atn == 8) THEN ! O
          defb%a(1)%rp = (/ 130.7093200_dp, 23.8088610_dp, 6.4436083_dp /)
          defb%a(2)%rp = (/ 5.0331513_dp, 1.1695961_dp, 0.3803890_dp /)
          defb%a(3)%rp = (/ 5.0331513_dp, 1.1695961_dp, 0.3803890_dp /)
        ELSE IF (atn == 9) THEN ! F
          defb%a(1)%rp = (/ 166.6791300_dp, 30.3608120_dp, 8.2168207_dp /)
          defb%a(2)%rp = (/ 6.4648032_dp, 1.5022812_dp, 0.4885885_dp /)
          defb%a(3)%rp = (/ 6.4648032_dp, 1.5022812_dp, 0.4885885_dp /)
        END IF
      ELSE IF (bname == '3-21G ') THEN
        Allocate(defb%nf(5), STAT=istat) ! 9 bf's: 1s 2s 2p(xyz) 2s' 2p'(xyz)
        defb%nf(1) = 3; defb%nf(2) = 2; defb%nf(3) = 2;
        defb%nf(4) = 1; defb%nf(5) = 1  !  3-21G
        Allocate(defb%lam(5), STAT=istat) ! angular momentum for each gaussian
        defb%lam(1) = 'S'; defb%lam(2) = 'S'; defb%lam(3) = 'P'
        defb%lam(4) = 'S'; defb%lam(5) = 'P'
        Allocate(defb%d(5), STAT=istat)
        Allocate(defb%d(1)%rp(3),STAT=istat)
        Allocate(defb%d(2)%rp(2),STAT=istat)
        Allocate(defb%d(3)%rp(2),STAT=istat)
        Allocate(defb%d(4)%rp(1),STAT=istat)
        Allocate(defb%d(5)%rp(1),STAT=istat)
        Allocate(defb%a(5), STAT=istat)
        Allocate(defb%a(1)%rp(3), STAT=istat)
        Allocate(defb%a(2)%rp(2),STAT=istat)
        Allocate(defb%a(3)%rp(2),STAT=istat)
        Allocate(defb%a(4)%rp(1),STAT=istat)
        Allocate(defb%a(5)%rp(1),STAT=istat)
        IF (atn == 6) THEN ! C
          defb%d(1)%rp = (/ 0.0617669_dp, 0.3587940_dp, 0.7007130_dp /)
          defb%a(1)%rp = (/ 172.2560000_dp, 25.9109000_dp, 5.5333500_dp /)
          defb%d(2)%rp = (/ -0.3958970_dp, 1.2158400_dp /)
          defb%a(2)%rp = (/ 3.6649800_dp, 0.7705450_dp /)
          defb%d(3)%rp = (/ 0.2364600_dp, 0.8606190_dp /)
          defb%a(3)%rp = (/ 3.6649800_dp, 0.7705450_dp /)
          defb%d(4)%rp = (/ 1.0_dp /)
          defb%a(4)%rp = (/ 0.1958570_dp /)
          defb%d(5)%rp = (/ 1.0_dp /)
          defb%a(5)%rp = (/ 0.1958570_dp /)
        ELSE IF (atn == 7) THEN ! N
          defb%d(1)%rp = (/ 0.0598657_dp, 0.3529550_dp, 0.7065130_dp /)
          defb%a(1)%rp = (/ 242.7660000_dp, 36.4851000_dp, 7.8144900_dp /)
          defb%d(2)%rp = (/ -0.4133010_dp, 1.2244200_dp /)
          defb%a(2)%rp = (/ 5.4252200_dp, 1.1491500_dp /)
          defb%d(3)%rp = (/ 0.2379720_dp, 0.8589530_dp /)
          defb%a(3)%rp = (/ 5.4252200_dp, 1.1491500_dp /)
          defb%d(4)%rp = (/ 1.0_dp /)
          defb%a(4)%rp = (/ 0.2832050_dp /)
          defb%d(5)%rp = (/ 1.0_dp /)
          defb%a(5)%rp = (/ 0.2832050_dp /)
        ELSE IF (atn == 8) THEN ! O
          defb%d(1)%rp = (/ 0.0592394_dp, 0.3515000_dp, 0.7076580_dp /)
          defb%a(1)%rp = (/ 322.0370000_dp, 48.4308000_dp, 10.4206000_dp /)
          defb%d(2)%rp = (/ -0.4044530_dp, 1.2215600_dp /)
          defb%a(2)%rp = (/ 7.4029400_dp, 1.5762000_dp /)
          defb%d(3)%rp = (/ 0.2445860_dp, 0.8539550_dp /)
          defb%a(3)%rp = (/ 7.4029400_dp, 1.5762000_dp /)
          defb%d(4)%rp = (/ 1.0_dp /)
          defb%a(4)%rp = (/ 0.3736840_dp /)
          defb%d(5)%rp = (/ 1.0_dp /)
          defb%a(5)%rp = (/ 0.3736840_dp /)
        END IF
      ELSE IF (bname == '4-31G ') THEN
        Allocate(defb%nf(5), STAT=istat) ! 9 bf's: 1s 2s 2p(xyz) 2s' 2p'(xyz)
        defb%nf(1) = 4; defb%nf(2) = 3; defb%nf(3) = 3;
        defb%nf(4) = 1; defb%nf(5) = 1  !  4-31G
        Allocate(defb%lam(5), STAT=istat) ! angular momentum for each gaussian
        defb%lam(1) = 'S'; defb%lam(2) = 'S'; defb%lam(3) = 'P'
        defb%lam(4) = 'S'; defb%lam(5) = 'P'
        Allocate(defb%d(5), STAT=istat)
        Allocate(defb%d(1)%rp(4),STAT=istat)
        Allocate(defb%d(2)%rp(3),STAT=istat)
        Allocate(defb%d(3)%rp(3),STAT=istat)
        Allocate(defb%d(4)%rp(1),STAT=istat)
        Allocate(defb%d(5)%rp(1),STAT=istat)
        Allocate(defb%a(5), STAT=istat)
        Allocate(defb%a(1)%rp(4), STAT=istat)
        Allocate(defb%a(2)%rp(3),STAT=istat)
        Allocate(defb%a(3)%rp(3),STAT=istat)
        Allocate(defb%a(4)%rp(1),STAT=istat)
        Allocate(defb%a(5)%rp(1),STAT=istat)
        IF (atn == 6) THEN ! C
          defb%d(1)%rp = (/ 0.0177258_dp, 0.1234787_dp, 0.4338754_dp, 0.5615042_dp /)
          defb%a(1)%rp = (/ 486.9669300_dp, 73.3710940_dp, 16.4134580_dp, 4.3449836_dp /)
          defb%d(2)%rp = (/ -0.1213837_dp, -0.2273385_dp, 1.1851739_dp /)
          defb%a(2)%rp = (/ 8.6735253_dp, 2.0966193_dp, 0.6046513_dp /)
          defb%d(3)%rp = (/ 0.0635454_dp, 0.2982678_dp, 0.7621032_dp /)
          defb%a(3)%rp = (/ 8.6735253_dp, 2.0966193_dp, 0.6046513_dp /)
          defb%d(4)%rp = (/ 1.0_dp /)
          defb%a(4)%rp = (/ 0.1835578_dp /)
          defb%d(5)%rp = (/ 1.0_dp /)
          defb%a(5)%rp = (/ 0.1835578_dp /)
        ELSE IF (atn == 7) THEN ! N
          defb%d(1)%rp = (/ 0.0175982511_dp, 0.1228462410_dp, 0.4337821410_dp, 0.5614182170_dp /)
          defb%a(1)%rp = (/ 671.2795000_dp, 101.2017000_dp, 22.6999700_dp, 6.0406090_dp /)
          defb%d(2)%rp = (/ -0.1174892990_dp, -0.2139940160_dp, 1.1745021100_dp /)
          defb%a(2)%rp = (/ 12.3935997_dp, 2.9223828_dp, 0.83252808_dp /)
          defb%d(3)%rp = (/ 0.0640203443_dp, 0.3112025550_dp, 0.7527482390_dp /)
          defb%a(3)%rp = (/ 12.3935997_dp, 2.9223828_dp, 0.83252808_dp /)
          defb%d(4)%rp = (/ 1.0_dp /)
          defb%a(4)%rp = (/ 0.2259640_dp /)
          defb%d(5)%rp = (/ 1.0_dp /)
          defb%a(5)%rp = (/ 0.2259640_dp /)
        ELSE IF (atn == 8) THEN ! O
          defb%d(1)%rp = (/ 0.0175506_dp, 0.1228292_dp, 0.4348836_dp, 0.5600108_dp /)
          defb%a(1)%rp = (/ 883.2728600_dp, 133.1292800_dp, 29.9064080_dp, 7.9786772_dp /)
          defb%d(2)%rp = (/ -0.1134010_dp, -0.1772865_dp, 1.1504079_dp /)
          defb%a(2)%rp = (/ 16.1944470_dp, 3.7800860_dp, 1.0709836_dp /)
          defb%d(3)%rp = (/ 0.0685453_dp, 0.3312254_dp, 0.7346079_dp /)
          defb%a(3)%rp = (/ 16.1944470_dp, 3.7800860_dp, 1.0709836_dp /)
          defb%d(4)%rp = (/ 1.0_dp /)
          defb%a(4)%rp = (/ 0.2838798_dp /)
          defb%d(5)%rp = (/ 1.0_dp /)
          defb%a(5)%rp = (/ 0.2838798_dp /)
        END IF
      ELSE IF (bname == '6-31G**') THEN
        Allocate(defb%nf(6), STAT=istat) ! 13 bf's: 1s 2s 2p(xyz) 2s' 2p'(xyz) d(6)
        defb%nf(1) = 6; defb%nf(2) = 3; defb%nf(3) = 3;
        defb%nf(4) = 1; defb%nf(5) = 1; defb%nf(6) = 1  !  6-31G**
        Allocate(defb%lam(6), STAT=istat) ! angular momentum for each gaussian
        defb%lam(1) = 'S'; defb%lam(2) = 'S'; defb%lam(3) = 'P'
        defb%lam(4) = 'S'; defb%lam(5) = 'P'; defb%lam(6) = 'D'
        Allocate(defb%d(6), STAT=istat)
        Allocate(defb%d(1)%rp(6),STAT=istat)
        Allocate(defb%d(2)%rp(3),STAT=istat)
        Allocate(defb%d(3)%rp(3),STAT=istat)
        Allocate(defb%d(4)%rp(1),STAT=istat)
        Allocate(defb%d(5)%rp(1),STAT=istat)
        Allocate(defb%d(6)%rp(1),STAT=istat)
        Allocate(defb%a(6), STAT=istat)
        Allocate(defb%a(1)%rp(6), STAT=istat)
        Allocate(defb%a(2)%rp(3),STAT=istat)
        Allocate(defb%a(3)%rp(3),STAT=istat)
        Allocate(defb%a(4)%rp(1),STAT=istat)
        Allocate(defb%a(5)%rp(1),STAT=istat)
        Allocate(defb%a(6)%rp(1),STAT=istat)
        IF (atn == 6) THEN ! C
          defb%d(1)%rp = (/ 0.0018347_dp, 0.0140373_dp, 0.0688426_dp, 0.2321844_dp, 0.4679413_dp, 0.3623120_dp /)
          defb%a(1)%rp = (/ 3047.5249000_dp, 457.3695100_dp, 103.9486900_dp, 29.2101550_dp, 9.2866630_dp, 3.1639270_dp /)
          defb%d(2)%rp = (/ -0.1193324_dp, -0.1608542_dp, 1.1434564_dp /)
          defb%a(2)%rp = (/ 7.8682724_dp, 1.8812885_dp, 0.5442493_dp /)
          defb%d(3)%rp = (/ 0.0689991_dp, 0.3164240_dp, 0.7443083_dp /)
          defb%a(3)%rp = (/ 7.8682724_dp, 1.8812885_dp, 0.5442493_dp /)
          defb%d(4)%rp = (/ 1.0_dp /)
          defb%a(4)%rp = (/ 0.1687144_dp /)
          defb%d(5)%rp = (/ 1.0_dp /)
          defb%a(5)%rp = (/ 0.1687144_dp /)
          defb%d(6)%rp = (/ 1.0_dp /)
          defb%a(6)%rp = (/ 0.8_dp /)
        ELSE IF (atn == 7) THEN ! N
          defb%d(1)%rp = (/ 0.0018348_dp, 0.0139950_dp, 0.0685870_dp, 0.2322410_dp, 0.4690700_dp, 0.3604550_dp /)
          defb%a(1)%rp = (/ 4173.5110000_dp, 627.4579000_dp, 142.9021000_dp, 40.2343300_dp, 12.8202100_dp, 4.3904370_dp /)
          defb%d(2)%rp = (/ -0.1149610_dp, -0.1691180_dp, 1.1458520_dp /)
          defb%a(2)%rp = (/ 11.6263580_dp, 2.7162800_dp, 0.7722180_dp /)
          defb%d(3)%rp = (/ 0.0675800_dp, 0.3239070_dp, 0.7408950_dp /)
          defb%a(3)%rp = (/ 11.6263580_dp, 2.7162800_dp, 0.7722180_dp /)
          defb%d(4)%rp = (/ 1.0_dp /)
          defb%a(4)%rp = (/ 0.2120313_dp /)
          defb%d(5)%rp = (/ 1.0_dp /)
          defb%a(5)%rp = (/ 0.2120313_dp /)
          defb%d(6)%rp = (/ 1.0_dp /)
          defb%a(6)%rp = (/ 0.8_dp /)
        ELSE IF (atn == 8) THEN ! O
          defb%d(1)%rp = (/ 0.0018311_dp, 0.0139501_dp, 0.0684451_dp, 0.2327143_dp, 0.4701930_dp, 0.3585209_dp /)
          defb%a(1)%rp = (/ 5484.6717000_dp, 825.2349500_dp, 188.0469600_dp, 52.9645000_dp, 16.8975700_dp, 5.7996353_dp /)
          defb%d(2)%rp = (/ -0.1107775_dp, -0.1480263_dp, 1.1307670_dp /)
          defb%a(2)%rp = (/ 15.5396160_dp, 3.5999336_dp, 1.0137618_dp /)
          defb%d(3)%rp = (/ 0.0708743_dp, 0.3397528_dp, 0.7271586_dp /)
          defb%a(3)%rp = (/ 15.5396160_dp, 3.5999336_dp, 1.0137618_dp /)
          defb%d(4)%rp = (/ 1.0_dp /)
          defb%a(4)%rp = (/ 0.2700058_dp /)
          defb%d(5)%rp = (/ 1.0_dp /)
          defb%a(5)%rp = (/ 0.2700058_dp /)
          defb%d(6)%rp = (/ 1.0_dp /)
          defb%a(6)%rp = (/ 0.8_dp /)
        END IF
      END IF
    END IF
  END FUNCTION defb

  SUBROUTINE undefb(defb)
    TYPE(basis) :: defb
    INTEGER(i1b) :: istat
    Deallocate(defb%nf, STAT=istat)
    Deallocate(defb%d, STAT=istat)
    Deallocate(defb%a, STAT=istat)
    Deallocate(defb%lam, STAT=istat)
  END SUBROUTINE undefb

END MODULE bsets
