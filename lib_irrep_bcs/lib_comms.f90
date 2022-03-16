
module lib_comms 

    implicit none 

    integer,     public, parameter :: dp = 8 

    real(dp),    public, parameter :: pi = 3.141592653589793238462643383279_dp

    complex(dp), public, parameter :: cmplx_i = (0.0_dp, 1.0_dp)
  
    complex(dp), public, parameter :: cmplx_0 = (0.0_dp, 0.0_dp)
  
    complex(dp), public, parameter :: cmplx_1 = (1.0_dp, 0.0_dp)

    integer,     public, parameter :: MAXDG = 16

    real(dp),    public, parameter :: TOLDG = 1.E-3

    real(dp),    public, parameter :: TOLREP = 0.05

    integer,     public, parameter :: MAXSYM = 96

    integer,     public, parameter :: MAXIR = 48

    integer,     public, parameter :: MAXIRDG = 4

    integer,     public, parameter :: MAXIRSUM = 4
    
    integer,     public, parameter :: MAXKTYPE = 40

    integer,     public, parameter :: MAXKSAVE = 10

    real(dp),    public, parameter :: epsil = 0.0000001_dp

    integer,                  save :: rot_bilbao(3,3,MAXSYM)
    integer,                  save :: invrot_bilbao(3,3,MAXSYM)
    real(dp),                 save :: tau_bilbao(3,MAXSYM)
    complex(dp),              save :: SU2_bilbao(2,2,MAXSYM)

    
    ! Coordinates of Bilbao k
    real(dp),         save :: samplek(3, MAXKTYPE)
    ! Standard name of Bilbao k
    character(len=2), save :: samplekname(MAXKTYPE)
    ! Number of irreps of little group
    integer,          save :: nirreps(MAXKTYPE)
    ! Number of irreps for single little group
    integer,          save :: sirreps(MAXKTYPE)
    ! The existence of antiunitary symmetrys
    integer,          save :: antisym(MAXKTYPE)
    ! Name of irreps
    character(len=5), save :: Irrepsname(MAXIR, MAXKTYPE)
    ! Reality of irreps
    integer,          save :: Herringrule(MAXIR,MAXKTYPE) 

    ! number of elements in little groups
    integer,          save :: nelelittle(MAXKTYPE)

    ! The first lable can take two values:
    !       0: symmetry element does not belong to the little group
    !       1:                  does belong to
    ! The second lable can take two values:
    !       1: the phase does not contain the u,v,w parameters
    !       2:           does contain
    integer,          save :: labels(2,MAXSYM,MAXIR,MAXKTYPE)
    ! for comparison
    complex(dp),      save :: tableTraces(MAXSYM,MAXIR,MAXKTYPE)
    complex(dp),      save :: chartTraces(MAXSYM,MAXIR,MAXKTYPE)
    real(dp),         save :: coeff_uvw(3,MAXSYM,MAXIR,MAXKTYPE)
    character(len=12),save :: factsTraces(MAXSYM,MAXIR,MAXKTYPE)

    character(len=15),save :: chkpoint(MAXKTYPE)
    character(len=10),save :: spacegroupsymbol
     
    integer, allocatable, save :: reorder(:) 

    real(dp),         save :: difftau_inputbilb(3,MAXSYM)
    complex(dp),      save :: kphase_difftau(MAXSYM)

    integer,          save :: num_doub_sym

    integer,          save :: num_ktype


    ! save data, used for a summary output
    integer,             save :: save_kcount
    integer, allocatable,save :: save_ktype(:)
    integer, allocatable,save :: save_chrct(:,:)
    integer, allocatable,save :: save_numdeg(:,:)
    integer, allocatable,save :: save_numrep(:,:)

    ! 1/2   ->  0.5
    ! 1/3   ->  0.333333
    ! u     ->  0.123  : su
    ! v     ->  0.313  : sv
    ! w     ->  0.427  : sw
    ! 1 - u ->  0.877  : s1_u
    ! -2 u  -> -0.246  : s_2u
    ! -u    -> -0.123  : s_u
    ! 1 + u ->  1.123  : s1u
    ! 1 - v ->  0.687  : s1_v
    !------------codes for u,v,w
    real(dp),    parameter  :: su=0.123_dp,s1_u=0.877_dp,s_2u=-0.246_dp
    real(dp),    parameter  :: sv=0.313_dp,s1_v=0.687_dp
    real(dp),    parameter  :: sw=0.427_dp,s_u=-0.123_dp,s1u = 1.123_dp
    real(dp),    parameter  :: ds1=1._dp/sqrt(3._dp), ds2=0.5_dp/sqrt(3._dp)

    real(dp),    parameter  :: inv3=1._dp/3._dp,dinv3=2._dp/3._dp
    real(dp),    parameter  :: Pabc(3,3) = RESHAPE( &
                (/  1.00000 ,  0.00000 ,  0.00000   &
                  , 0.00000 ,  1.00000 ,  0.00000   &
                  , 0.00000 ,  0.00000 ,  1.00000/) &
                  ,(/3,3/) )
    real(dp),    parameter  :: Cabc(3,3) = RESHAPE( &
                (/  0.50000 , -0.50000 ,  0.00000   &
                  , 0.50000 ,  0.50000 ,  0.00000   &
                  , 0.00000 ,  0.00000 ,  1.00000/) &
                  ,(/3,3/) )
    real(dp),    parameter  :: Cabc68(3,3) = RESHAPE( &
                (/  0.50000 ,  0.50000 ,  0.00000   &
                  ,-0.50000 ,  0.50000 ,  0.00000   &
                  , 0.00000 ,  0.00000 ,  1.00000/) &
                  ,(/3,3/) )
    real(dp),    parameter  :: Babc(3,3) = RESHAPE( &
! SG#5, #8, #9, #12, #15.
                (/  0.50000 ,  0.50000 ,  0.00000   &
                  ,-0.50000 ,  0.50000 ,  0.00000   &
                  , 0.00000 ,  0.00000 ,  1.00000/) &
                  ,(/3,3/) )
    real(dp),    parameter  :: Aabc(3,3) = RESHAPE( &
! SG#38-41.
                (/  1.00000 ,  0.00000 ,  0.00000   &
                  , 0.00000 ,  0.50000 ,  0.50000   &
                  , 0.00000 , -0.50000 ,  0.50000/) &
                  ,(/3,3/) )
   !real(dp),    parameter  :: Rabc(3,3) = RESHAPE( &
   !            (/    ds2   , -0.500_dp,   inv3     &
   !              ,   ds2   ,  0.500_dp,   inv3     &
   !              ,  -ds1   ,  0.000_dp,   inv3  /) &
   !              ,(/3,3/) )
    real(dp),    parameter  :: Rabc(3,3) = RESHAPE( &
                (/  dinv3   ,  inv3    ,   inv3     &
                  , -inv3   ,  inv3    ,   inv3     &
                  , -inv3   ,-dinv3    ,   inv3  /) &
                  ,(/3,3/) )
    real(dp),    parameter  :: Fabc(3,3) = RESHAPE( &
                (/  0.00000 ,  0.50000 ,  0.50000   &
                  , 0.50000 ,  0.00000 ,  0.50000   &
                  , 0.50000 ,  0.50000 ,  0.00000/) &
                  ,(/3,3/) )
    real(dp),    parameter  :: Iabc(3,3) = RESHAPE( &
                (/ -0.50000 ,  0.50000 ,  0.50000   &
                  , 0.50000 , -0.50000 ,  0.50000   &
                  , 0.50000 ,  0.50000 , -0.50000/) &
                  ,(/3,3/) )



    ! used for setup
    real(dp),    allocatable, save :: SO3_save(:,:,:)
    complex(dp), allocatable, save :: SU2_save(:,:,:)
    integer,     allocatable, save :: rot_save(:,:,:)
    integer,     allocatable, save :: invrot_save(:,:,:)

    real(dp),    allocatable, save :: rot_orb(:,:,:,:)
    integer,     allocatable, save :: startorb(:)
    integer,     allocatable, save :: rot_atom(:,:)
    real(dp),    allocatable, save :: rot_phivec(:,:,:)
    
    ! d(5) + f(7)
    integer,      parameter :: mix2l = 12

    complex(dp),  parameter :: pwann(3,3)=RESHAPE(     &
    !inputform just transpose matrix
    !2basis:  pz, px, py
    (/ ( 0.d0,  0.d0) ,( 0.d0,  0.d0) ,( 1.d0,  0.d0)  &
      ,( 1.d0,  0.d0) ,( 0.d0,  0.d0) ,( 0.d0,  0.d0)  &
      ,( 0.d0,  0.d0) ,( 1.d0,  0.d0) ,( 0.d0,  0.d0)/)&
      ,(/3,3/) )
    complex(dp),  parameter :: dwann(5,5)=RESHAPE(     &
    !inputform just transpose matrix
    !basis: 3z2-r2, xz, yz, x2-y2, xy
    (/ ( 0.d0,  0.d0) ,( 0.d0,  0.d0) ,( 0.d0,  0.d0) ,( 0.d0,  0.d0) ,( 1.d0,  0.d0)  &
      ,( 0.d0,  0.d0) ,( 0.d0,  0.d0) ,( 1.d0,  0.d0) ,( 0.d0,  0.d0) ,( 0.d0,  0.d0)  &
      ,( 0.d0,  0.d0) ,( 1.d0,  0.d0) ,( 0.d0,  0.d0) ,( 0.d0,  0.d0) ,( 0.d0,  0.d0)  &
      ,( 0.d0,  0.d0) ,( 0.d0,  0.d0) ,( 0.d0,  0.d0) ,( 1.d0,  0.d0) ,( 0.d0,  0.d0)  &
      ,( 1.d0,  0.d0) ,( 0.d0,  0.d0) ,( 0.d0,  0.d0) ,( 0.d0,  0.d0) ,( 0.d0,  0.d0)/)&
      ,(/5,5/) )
    !added on 12-19-2019 by zjwang
    real(dp), parameter :: sqr38 = dsqrt(0.375d0)
    real(dp), parameter :: sqr58 = dsqrt(0.625d0)
    complex(dp),  parameter :: fwann(7,7)=RESHAPE(     &
    !inputform just transpose matrix
    !basis: 3z2-r2, xz, yz, x2-y2, xy
    (/ ( 0.d0,0.d0) ,( 0.d0,0.d0) ,( 0.d0,0.d0) ,( 1.d0,0.d0) ,( 0.d0,0.d0), ( 0.d0,0.d0 ), ( 0.d0,0.d0 )  &
      ,( 0.d0,0.d0) ,(-sqr38,0.d0) ,( 0.d0,0.d0) ,( 0.d0,0.d0) ,(-sqr58,0.d0), (0.d0,0d0),(0.d0,0.d0) &
      ,( 0.d0,0.d0) ,( 0.d0,0.d0) ,(-sqr38,0.d0) ,( 0.d0,0.d0) ,( 0.d0,0.d0), ( sqr58,0.d0), (0.d0,0.d0)  &
      ,( 0.d0,0.d0) ,( 0.d0,0.d0) ,( 0.d0,0.d0) ,( 0.d0,0.d0) ,( 0.d0,0.d0), ( 0.d0,0.d0), ( 1.d0,0.d0)  &
      ,( 1.d0,0.d0) ,( 0.d0,0.d0) ,( 0.d0,0.d0) ,( 0.d0,0.d0) ,( 0.d0,0.d0), ( 0.d0,0.d0 ), (0.d0,0.d0)  &
      ,( 0.d0,0.d0) ,( sqr58,0.d0) ,(0.d0,0.d0) ,( 0.d0,0.d0) ,(-sqr38,0.d0), ( 0.d0,0.d0), (0.d0,0.d0)  &
      ,( 0.d0,0.d0) ,( 0.d0,0.d0) ,(-sqr58,0.d0) ,( 0.d0,0.d0) ,( 0.d0,0.d0), (-sqr38,0.d0), ( 0.d0,0.d0)/)&
      ,(/7,7/) )
    ! (f1,f2,f3,f4.f5,f6,f7) = (Y3,3  Y3,2  Y3,1  Y3,0  Y3,-1  Y3,-2  Y3,-3) * ss
    ! the default f-orbitals under real sphere harmonic functions
    !       0      -/(5)/4    -i/(5)/4      0      /(3)/4      -i/(3)/4      0
    !  -i//(2)        0           0         0        0             0        1//(2)
    !       0       /(3)/4     -i/(3)/4     0      /(5)/4       i/(5)/4      0
    !       0         0           0         1        0             0         0        = ss1
    !       0      -/(3)/4     -i/(3)/4     0      -/(5)/4      i/(5)/4      0
    !   i//(2)        0           0         0        0             0        1//(2)
    !       0        /(5)/4   -i/(5)/4      0     -/(3)/4      -i/(3)/4      0
    ! the wannier f-orbitals under real sphere harmonic functions
    !       0         0           0         0         0        1//(2)     i//(2)
    !       0         0           0      1//(2)    i//(2)          0         0
    !       0       1//(2)      i//(2)      0         0            0         0
    !       1         0           0         0         0            0         0        = ss2
    !       0      -1//(2)      i//(2)      0         0            0         0
    !       0         0           0      1//(2)   -i//(2)          0         0
    !       0         0           0         0         0        -1//(2)    i//(2)
    ! the transform matrix from wannier orbitals to default orbitals: ss1^{-1}*ss2
    
end module lib_comms
