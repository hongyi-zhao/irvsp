module comms  
    
    implicit none 

    integer,     public, parameter :: dp = 8
    
    real(dp),    public, parameter :: pi = 3.141592653589793238462643383279d0

    integer,     public, parameter :: MAXSYM = 96

    integer,     public, parameter :: MAXIR = 48

    integer,     public, parameter :: MAXIRDG = 4

    integer,     public, parameter :: MAXKTYPE = 40

    integer,     public            :: sgn 

    logical,     public            :: isSymmorphic
    logical,     public            :: isComplexWF
    logical,     public            :: isSpinor
    logical,     public            :: isInv
    logical,     public            :: isSpinPola

    real(dp),    public            :: br2(3,3), br4(3,3)

    integer,     public            :: num_sym
    
    integer,     public            :: num_k

    integer,     public            :: num_bands

    integer,     public            :: bot_band

    integer,     public            :: top_band
    
    integer,     public            :: nspin 

    logical,     public            :: spinor

    integer,     public            :: max_plane

    real(dp),    public            :: lattice(3,3)
    
    real(dp),    save, public ::  det_read(MAXSYM)
    real(dp),    save, public ::  angle_read(MAXSYM)
    real(dp),    save, public ::  axis_read(3,MAXSYM)
    real(dp),    save, public ::  tau_read(3,MAXSYM)
    real(dp),    allocatable, save, public ::  det_input(:)
    real(dp),    allocatable, save, public ::  angle_input(:)
    real(dp),    allocatable, save, public ::  axis_input(:,:)
    real(dp),    allocatable, save, public ::  tau_input(:,:)

    ! from WAVECAR
    integer    , allocatable, save, public ::  igall(:,:) 
    real(dp)   , allocatable, save, public ::  KV(:,:) 

    complex(dp), allocatable, save, public ::  coeffa(:,:),coeffb(:,:)
    real(dp)   , allocatable, save, public ::  EE(:)

    real(dp),                 save, public ::  WK(3)

    integer,                  save, public :: ncnt, nplane
   

end module comms
