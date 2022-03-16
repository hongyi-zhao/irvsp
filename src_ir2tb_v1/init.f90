!=========================================================================!
! project : struct_data
! history : 4/7/2015
! authors : Zhijun Wang  ( zjwang11@hotmail.com )
! purpose : Get input from TBbox.in
! status  : robust
! comment : These programs are distributed in the hope that they will be 
!           useful, but WITHOUT ANY WARRANTY; without even the implied 
!           warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
!=========================================================================!
module  struct_data
implicit none
private

integer   ,    public  ,  parameter ::   dp  =  8
                                         
integer   ,    public  ,  save      ::   NKPTS ! - max. no. of kpoints
integer   ,    public  ,  save      ::   NSPIN ! - spin components 
integer   ,    public  ,  save      ::   nrot  ! - number of op.s 
character(15), public  ,  save      ::   casename
character(20), public  ,  save      ::   hrfile

!--- angle l,atoms,orbits,spin
integer   ,    public  ,  save      ::  orbt  ! orbital order
integer   ,    public  ,  save      ::  ntau  ! number of atoms
integer   ,    public  ,  allocatable, save :: n2obtau(:) ! type of orbitals
integer   ,    public  ,  allocatable, save :: n2ibtau(:) ! index of bands
integer   ,    public  ,  allocatable, save :: rot2tau(:,:) ! index of bands
integer   ,    public  ,  allocatable, save :: rot2vec(:,:,:) ! index of bands
real(dp)  ,    public  ,  allocatable, save :: rot2phi(:,:,:) ! index of bands
real(dp)  ,    public  ,  allocatable, save :: rotjorb(:,:,:,:) !rot orb
real(dp)  ,    public  ,  allocatable, save :: len_k(:),k(:,:)

!----------crystal
real(dp),  save     :: br2(3,3),br4(3,3)
real(dp), parameter :: PI = 3.141592653589793238462643383279d0
integer , parameter :: mix2l  = 12  ! max orbitals of 2*L+1
integer , allocatable, save :: n2ktau(:) ! species of atom
real(dp), allocatable, save :: pos(:,:)  ! index of bands
real(dp), allocatable, save :: kpoints(:,:)

!added on 12-19-2019 by zjwang
!!        s, px, py, pz, xy, yz, zx, x2-y2, 3z2-r2
!         fxyz, f5x3-xr2, f5y3-yr2, f5z3-zr2, fx(y2-z2), fy(z2-x2), fz(x2-y2)
!2wannier!s, pz, px, py, 3z2-r2, xz, yz, x2-y2, xy
!            3   1   2    5      3   2   4      1
!         fz3, fxz2, fyz2, fz(x2-y2), fxyz, fx(x2-3y2), fy(3x2-y2)
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


public   ::  init

CONTAINS

subroutine init( NSYM,FLMAX,FL,IORD,IZ,TAU,DZ2 )
  implicit none
  integer  ,intent(IN)   :: NSYM
  integer  ,intent(IN)   :: FLMAX  ! FLMAX - size of flag (FL) array (should be 4)
  logical  ,intent(OUT)  :: FL(FLMAX)
  integer  ,intent(OUT)  :: IORD
  integer  ,intent(OUT)  :: IZ(3,3,NSYM)
  real(dp) ,intent(OUT)  :: TAU(3,NSYM)
  real(dp) ,intent(OUT)  :: DZ2(3,3,NSYM)  

  integer    :: IIZ(3,3,NSYM)
  integer  , parameter   :: nfst=20
  logical    :: FL1

  integer :: irot,i,j,itau,jtau
  integer :: I1,I2
  integer :: Nk,kmesh

  character*120 :: chaps
  character*10  :: title
  character*5   :: chtp5
  character*15  :: chtp15
  character*35  :: chtp35
  integer       :: itmp,ierr
  real(dp)      :: rtmp(8),rotmt(3,3)
  complex(dp)   :: crotmt(3,3),protmt(3,3),drotmt(5,5),frotmt(7,7)
  complex(dp)   :: cmat3(3,3), cmat5(5,5), cmat7(7,7)
  real(dp):: tmp(3),tmpr(3)


 !**********************************************************
  FL=.FALSE.;FL(1)=.TRUE.;FL1=.FALSE.
  IZ=0;TAU=0.d0;DZ2=0.d0;NSPIN=1
!
  OPEN(unit=nfst,file='tbbox.in',form='formatted',status='old')

       chtp15='case    '; call get_key_para_cht(chtp15,nfst,casename)
      !write(0,'(2X,A10,A5,A15)') chtp15,':  ',               casename
       hrfile=trim(casename)//'_hr.dat'
       title=trim(casename)//'wann_TB'
       IF(casename=='soc') THEN
       NSPIN=2
       FL(2)=.TRUE.
       ENDIF


! set ntau
       chtp15='proj    '; call get_key_para_loc(chtp15,nfst)
       chtp15='orbt    '; call get_key_para_intct(chtp15,nfst, orbt )
       IF (orbt==1) THEN
          PRINT*, "Orbital convertion 1: s, px, py, pz, xy, yz, zx, x2-y2, 3z2-r2 !!!"
       ELSEIF(orbt==2) THEN
          PRINT*, "Orbital convertion 2: s, pz, px, py, 3z2-r2, xz, yz, x2-y2, xy (Wannier90) !!!"
       ELSE
          STOP "Error: orbt !!!!"
       ENDIF
       chtp15='ntau    '; call get_key_para_intct(chtp15,nfst, ntau )
       allocate(n2ktau(ntau),n2obtau(ntau),n2ibtau(ntau))
       allocate(pos(3,ntau))
       n2ibtau=0
       do itau=1,ntau
         read(nfst,*) pos(:,itau),n2ktau(itau),n2obtau(itau)
         IF(itau/=1) n2ibtau(itau)=n2ibtau(itau-1)+n2obtau(itau-1)
        !write(6,"(3F12.6,3I5)") pos(:,itau), n2ktau(itau) &
        !                    , n2obtau(itau), n2ibtau(itau)
       enddo

       allocate(rotjorb(mix2l,mix2l,ntau,NSYM)); rotjorb=0.d0

! set br2,br4
  chtp15='unit_cell'; call get_key_para_loc(chtp15,nfst)
  do i=1,3
    read(nfst,*) br2(:,i),(br4(i,j),j=1,3)
  enddo

!
! set IORD,TAU,DZ2
  do irot=1,NSYM
     read(nfst,*,iostat=ierr)  itmp,rtmp
     if (ierr /= 0) exit
     IORD =irot
     TAU(:,irot)=rtmp(6:8)
     IF( abs( rtmp(1)+1.d0+rtmp(2)   ).lt.1.d-5 ) FL(3)=.TRUE.
     IF( abs( rtmp(6)+rtmp(7)+rtmp(8)).gt.1.d-5 ) FL1  =.TRUE.
     call Dmatrix(rnx=rtmp(3), rny=rtmp(4), rnz=rtmp(5), degree=-rtmp(2), twoja1=3, Dmat=crotmt)
     protmt(:,:)=cmplx(0._dp,0._dp,dp); protmt=crotmt
          IF(orbt==2) THEN
             protmt(:,:)=matmul(protmt,pwann)
             cmat3(:,:)=cmplx(0._dp,0._dp,dp); cmat3=transpose(pwann)
             protmt(:,:)=matmul(cmat3,protmt)
          ENDIF
     call Dmatrix(rnx=rtmp(3), rny=rtmp(4), rnz=rtmp(5), degree=-rtmp(2), twoja1=5, Dmat=drotmt)
          IF(orbt==2) THEN
             drotmt(:,:)=matmul(drotmt,dwann)
             cmat5(:,:)=cmplx(0._dp,0._dp,dp); cmat5=transpose(dwann)
             drotmt(:,:)=matmul(cmat5,drotmt)
          ENDIF
     call Dmatrix(rnx=rtmp(3), rny=rtmp(4), rnz=rtmp(5), degree=-rtmp(2), twoja1=7, Dmat=frotmt)
          IF(orbt==2) THEN
             frotmt(:,:)=matmul(frotmt,fwann)
             cmat7(:,:)=cmplx(0._dp,0._dp,dp); cmat7=transpose(fwann)
             frotmt(:,:)=matmul(cmat7,frotmt)
          ENDIF

!--------set rotjorb
      do itau=1,ntau
       if(n2obtau(itau)==1) then
         rotjorb(1,1,itau,irot)=1.d0!cmplx(1.d0,0.d0,dp)
       elseif(n2obtau(itau)==3) then
         rotjorb(1:3,1:3,itau,irot)=real(protmt,dp)
         IF(rtmp(1).lt. 1.d-2) rotjorb(1:3,1:3,itau,irot)=-real(protmt,dp)
       !----------------
       elseif(n2obtau(itau)==4) then
         IF(irot==1) THEN
           write(0,"(A63,I3,A5)") &
         "WARNING!!! Two kinds of L-orbtals (s and p) are considered for",itau,"-atom"
         ENDIF
         rotjorb(1,1,itau,irot)=1.d0!cmplx(1.d0,0.d0,dp)
         rotjorb(2:4,2:4,itau,irot)= real(protmt,dp)
         IF(rtmp(1).lt. 1.d-2) rotjorb(2:4,2:4,itau,irot)=-real(protmt,dp)
       !-------------
       elseif(n2obtau(itau)==5) then
         rotjorb(1:5,1:5,itau,irot)=real(drotmt,dp)
       !----------------
       elseif(n2obtau(itau)==6) then
         rotjorb(1,1,itau,irot)=1.d0
         rotjorb(2:6,2:6,itau,irot)=real(drotmt,dp)
         IF(irot==1) THEN
            write(0,"(A63,I3,A5)") &
         "WARNING!!! Two kinds of L-orbtals (s and d) are considered for",itau,"-atom"
         ENDIF
         write(0,*) irot
         write(0,"(6F10.6)")  rotjorb(:,:,itau,irot)
         !-------------
       elseif(n2obtau(itau)==7) then
         rotjorb(1:7,1:7,itau,irot)=real(frotmt,dp)
         IF(rtmp(1).lt. 1.d-2) rotjorb(1:7,1:7,itau,irot)=-real(frotmt,dp)
         !-------------
       elseif(n2obtau(itau)==8) then
         rotjorb(1:3,1:3,itau,irot)=real(protmt,dp)
         IF(rtmp(1).lt. 1.d-2) rotjorb(1:3,1:3,itau,irot)=-real(protmt,dp)
         rotjorb(4:8,4:8,itau,irot)=real(drotmt,dp) !  xy, yz, zx, x2-y2, 3z2-r2
       elseif(n2obtau(itau)==9) then
         rotjorb(1,1,itau,irot)=1.d0 
         rotjorb(2:4,2:4,itau,irot)=real(protmt,dp)
         IF(rtmp(1).lt. 1.d-2) rotjorb(2:4,2:4,itau,irot)=-real(protmt,dp)
         rotjorb(5:9,5:9,itau,irot)=real(drotmt,dp) !  xy, yz, zx, x2-y2, 3z2-r2
       elseif(n2obtau(itau)==12) then
         rotjorb(1:5,1:5,itau,irot)=real(drotmt,dp)
         rotjorb(6:12,6:12,itau,irot)=real(frotmt,dp)
         IF(rtmp(1).lt. 1.d-2) rotjorb(6:12,6:12,itau,irot)=-real(frotmt,dp)
       else
         STOP" Please set rotjorb for each atom"
       endif
    
      enddo 

!--------set rotjorb
     rotmt=real(crotmt)
     IF(rtmp(1).lt. 1.d-2) rotmt=-rotmt
     DZ2(:,:,irot)=rotmt(:,:)

  enddo

  nrot=IORD
!
! set IZ
  IZ=0
  do irot=1,IORD
     IZ(:,:,irot)= nint(matmul(matmul(br4, DZ2(:,:,irot) ),br2))
  enddo

! set IIZ
  IIZ=0
  DO irot=1,IORD
  CALL INVMATI(IZ(1,1,irot),IIZ(1,1,irot))
  ENDDO

       allocate(rot2tau(ntau,NSYM),rot2vec(3,ntau,NSYM))
       rot2tau=0;rot2vec=0
       allocate(rot2phi(3,ntau,NSYM));rot2phi=0.d0

  do irot=1,IORD
     do itau=1,ntau
     tmp(:)=matmul(IZ(:,:,irot),pos(:,itau))+tau(:,irot)
     do jtau=1,ntau
       tmpr(:)=tmp(:)-pos(:,jtau)
       IF(abs(tmpr(1)-nint(tmpr(1))) .lt.0.1d-2   &
     .and.abs(tmpr(2)-nint(tmpr(2))) .lt.0.1d-2   &
     .and.abs(tmpr(3)-nint(tmpr(3))) .lt.0.1d-2 ) THEN
       EXIT
       ENDIF
     enddo
      IF(jtau==ntau+1) STOP "rotation error"
      rot2tau(itau,irot)=jtau
      rot2vec(:,itau,irot)=nint(tmpr(:))
      !write(0,'(4I5,3F12.6)') irot,rot2vec(:,itau,irot),tmpr(:)
      !write(0,'(I5,I2,A5,I2,3F12.6)') irot,itau,"->",jtau,tmpr(:)
     !---be careful
      rot2vec(:,itau,irot)=matmul(IIZ(:,:,irot),rot2vec(:,itau,irot))
      !-------------------------------------------------
      !-------PHASE NOTE--------------------------------
      !-----USUALLY R|TAU  GIVEs US A PHASE OF e^-ik{R^-1}tau--------
      !-----DUE TO SYMMORPHIC CONDITION----------------
      !------e^-ik{R^-1}tau === e^-ik.tau----------------
      !-------------------------------------------------
      !-----******************************* ------------
      !-----IN CONVENIENCE WE USE e^-ik.tau ------------
      !-----******************************* ------------
      rot2phi(:,itau,irot)=-TAU(:,irot)+rot2vec(:,itau,irot)
     !---be careful
     enddo
     !write(0,'(5I5)') irot,rot2tau(:,irot)
  enddo

! set kpath
       chtp15='kpoint   '; call get_key_para_loc(chtp15,nfst)
       chtp15='kmesh    '; call get_key_para_intct(chtp15,nfst,kmesh   )
      !write(6,'(2X,A10,A5,I5)') chtp15,':  ',         kmesh
       chtp15='Nk      '; call get_key_para_intct(chtp15,nfst,Nk      )
      !write(6,'(2X,A10,A5,I5)') chtp15,':  ',         Nk
       allocate(kpoints(3,0:Nk)); kpoints=0.d0
       call get_key_para_nvec(nfst,Nk,kpoints)
      !write(6,"(3(F12.6))") kpoints



! set kpoints
       allocate(len_k(0:Nk*kmesh), k(3,0:Nk*kmesh))
  len_k=0.d0;k=0.d0
  k(:,0)=kpoints(:,0)
  do i=1,Nk
    tmp=(kpoints(:,i)-kpoints(:,i-1))/kmesh
    do j=1,kmesh
    k(:,(i-1)*kmesh+j)=kpoints(:,i-1)+tmp*j
    tmpr=matmul(tmp(:),br4(:,:))
    len_k((i-1)*kmesh+j) =  len_k((i-1)*kmesh+j-1) &
                           +dsqrt(tmpr(1)**2+tmpr(2)**2+tmpr(3)**2)
   !write(0,'(I5,4F12.6)') (i-1)*kmesh+j,k(:,(i-1)*kmesh+j),len_k((i-1)*kmesh+j)
    enddo
  enddo
    len_k=2.d0*PI*len_k

  NKPTS=Nk*kmesh

       deallocate(pos,n2ktau)
       deallocate(kpoints)

!
  CLOSE(nfst)
      IF(.NOT.FL(2) .AND. FL(3)) FL(1) = .FALSE.
 !**********************************************************
      WRITE(6,529) TITLE,' P '
!.....write out flags
      IF(     FL1)   WRITE(6,'(A23)',advance='NO') ' Non-symmorphic crystal'
      IF(.NOT.FL1)   WRITE(6,'(A19)',advance='NO') ' Symmorphic crystal'
      IF(     FL(3)) WRITE(6,'(A24)')   ' with inversion symmetry'
      IF(.NOT.FL(3)) WRITE(6,'(A27)')   ' without inversion symmetry'
      IF(     FL(1)) WRITE(6,'(A23)')   ' Complex eigenfunctions'
      IF(.NOT.FL(1)) WRITE(6,'(A20)')   ' Real eigenfunctions'
      IF(     FL(2)) WRITE(6,'(A45)')    &
                   ' Spin-orbit eigenfunctions (->time inversion)'
      IF(.NOT.FL(2)) WRITE(6,'(A29)')    &
                   ' No spin-orbit eigenfunctions'
      IF(     FL(4)) WRITE(6,'(A18)')    &
                   ' Spin-polarization'
      IF(.NOT.FL(4)) WRITE(6,'(A21)')    &
                   ' No spin-polarization'
!
      WRITE(6,590)
      WRITE(6,592)
      WRITE(6,595) ((br2(I1,I2),I2=1,3),I1=1,3)
      WRITE(6,591)
      WRITE(6,595) ((br4(I1,I2),I2=1,3),I1=1,3)
!  
      RETURN
 529  FORMAT(1X,A10,/,1X,A4,' lattice')
 590  FORMAT(//,' Transformations:',/, &
             ' Direct lattice vectors in Cartesian coord. system (BR2)')
 591  FORMAT(' Recipiical lattice vectors in Cartesian coord. system (BR4)')
 592  FORMAT('         pa:             pb:             pc: ')
 595  FORMAT(3(3F16.8,/))

end subroutine init

!!ref    : https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
!!        s, px, py, pz, xy, yz, zx, x2-y2, 3z2-r2
!!        fxyz, f5x3-xr2, f5y3-yr2, f5z3-zr2, fx(y2-z2), fy(z2-x2), fz(x2-y2)
!!TABLE I: http://journals.aps.org/prb/pdf/10.1103/PhysRevB.79.045107
subroutine Dmatrix(ih, ik, il, rnx, rny, rnz, rtheta, rphi, degree, twoja1, Dmat)
  implicit none
  integer, intent(in), optional :: ih, ik, il
  real(dp), intent(in), optional :: rnx, rny, rnz
  real(dp), intent(in), optional :: rtheta, rphi
  real(dp), intent(in) ::  degree
  integer, intent(in) :: twoja1  ! equal to 2*J + 1
  ! Actually, this rotational Dmat should be a real matrix.
  complex(dp),dimension(twoja1,twoja1),intent(out) :: Dmat
  complex(dp),dimension(twoja1,twoja1)             :: Dtmp
  ! the key parameter for Dmatrix
  real(dp) :: nx, ny, nz, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4, nx5,ny5,nz5,nx6,ny6,nz6
  real(dp) ::  omega
  ! complex(dp),dimension(aint(J)*2+1,aint(J)*2+1) :: Jx,Jy,Jz
  ! integer :: m
  complex(dp),parameter:: ci = (0.d0,1.d0)
  real*8     ,parameter:: PI = 3.141592653589793238462643383279d0
  real(dp) :: tmp
  real(dp) :: sqrt3, sqrt15 
  real(dp) :: cosOmega_2, sinOmega_2
  real(dp) :: cos1omega, sin1omega, cos2omega, sin2omega, cos3omega, sin3omega, cos4omega, sin4omega
  real(dp) :: cos5omega, sin5omega, cos6omega, sin6omega
  real(dp) :: cos3omega_2, sin3omega_2, cos5omega_2, sin5omega_2 
  complex(dp) :: exp2omega, exp_2omega, exp4omega, exp_4omega, expomega_2, exp_omega_2
  complex(dp) :: exp3omega, exp_3omega, exp6omega, exp_6omega, expomega, exp_omega, exp5omega, exp_5omega  
  sqrt3 = dsqrt(3.0d0)
  sqrt15 = dsqrt(15.0d0)
  if(present(ih).and.present(ik).and.present(il))then
    tmp = dsqrt(dble(ih)*dble(ih) + dble(ik)*dble(ik) + dble(il)*dble(il))
    nx = dble(ih)/tmp
    ny = dble(ik)/tmp
    nz = dble(il)/tmp
  else if(present(rtheta).and.present(rphi))then
    nx = dsin(rtheta)*dcos(rphi)
    ny = dsin(rtheta)*dsin(rphi)
    nz = dcos(rtheta)
  else if(present(rnx).and.present(rny).and.present(rnz))then
    tmp = dsqrt(rnx*rnx + rny*rny + rnz*rnz )
    if(abs(tmp-1.d0).gt.1.d-6) then
      write(0,'(A,F16.10)')"WARNING!!! The nx, ny, nz input is not normailized: ",tmp
      write(0,*)"The program will normalize them now."
    end if
      nx = rnx/tmp
      ny = rny/tmp
      nz = rnz/tmp
  else
    write(*,*)" Please choose one of the three type of input: "
    write(*,*)" 1. h k l:"
    write(*,*)" 2. nx ny nz:"
    write(*,*)" 3. theta phi:"
  end if
  omega = dble(nint(degree))*PI/180.d0
  IF(twoja1 == 5) omega=omega/2.d0
  cosOmega_2 = dcos(omega/2.d0)
  sinOmega_2 = dsin(omega/2.d0)
  cos3omega_2 = dcos(3.d0*omega/2.d0)
  sin3omega_2 = dsin(3.d0*omega/2.d0)
  cos5omega_2 = dcos(5.d0*omega/2.d0)
  sin5omega_2 = dsin(5.d0*omega/2.d0)
  cos1omega = dcos(omega)
  sin1omega = dsin(omega)
  cos2omega = dcos(2.d0*omega)
  sin2omega = dsin(2.d0*omega)
  cos3omega = dcos(3.d0*omega)
  sin3omega = dsin(3.d0*omega)
  cos4omega = dcos(4.d0*omega)
  sin4omega = dsin(4.d0*omega)
  cos5omega = dcos(5.d0*omega)
  sin5omega = dsin(5.d0*omega)
  cos6omega = dcos(6.d0*omega)
  sin6omega = dsin(6.d0*omega)
  expomega = cmplx(cos1omega, sin1omega)
  exp_omega = cmplx(cos1omega, -sin1omega)
  exp2omega = cmplx(cos2omega,sin2omega)
  exp_2omega = cmplx(cos2omega,-sin2omega)
  exp4omega = cmplx(cos4omega,sin4omega)
  exp_4omega = cmplx(cos4omega,-sin4omega)
  exp3omega = cmplx(cos3omega,sin3omega)
  exp_3omega = cmplx(cos3omega,-sin3omega)
  exp5omega = cmplx(cos5omega, sin5omega)
  exp_5omega = cmplx(cos5omega, -sin5omega)
  exp6omega = cmplx(cos6omega,sin6omega)
  exp_6omega = cmplx(cos6omega,-sin6omega)
  expomega_2 = cmplx(cosOmega_2, sinOmega_2)
  exp_omega_2 = cmplx(cosOmega_2, -sinOmega_2)
  nx2 = nx**2
  ny2 = ny**2
  nz2 = nz**2
  nx3 = nx**3
  ny3 = ny**3
  nz3 = nz**3
  nx4 = nx**4
  ny4 = ny**4
  nz4 = nz**4
  nx5 = nx**5
  ny5 = ny**5
  nz5 = nz**5
  nx6 = nx**6
  ny6 = ny**6
  nz6 = nz**6
  if(twoja1 == 2) then
    Dmat(1,1) = cmplx(cosOmega_2, -nz*sinOmega_2)
    Dmat(2,2) = cmplx(cosOmega_2, nz*sinOmega_2)
    Dmat(1,2) = cmplx(-ny, -nx)*sinOmega_2
    Dmat(2,1) = cmplx(ny, -nx)*sinOmega_2
  else if(twoja1 == 3) then
    Dmat(1,1) = nx**2+(1-nx**2)*cos1omega
    Dmat(2,2) = ny**2+(1-ny**2)*cos1omega
    Dmat(3,3) = nz**2+(1-nz**2)*cos1omega
    Dmat(1,2) = nx*ny*(1-cos1omega) - nz*sin1omega
    Dmat(2,1) = nx*ny*(1-cos1omega) + nz*sin1omega
    Dmat(1,3) = nx*nz*(1-cos1omega) + ny*sin1omega
    Dmat(3,1) = nx*nz*(1-cos1omega) - ny*sin1omega
    Dmat(2,3) = ny*nz*(1-cos1omega) - nx*sin1omega
    Dmat(3,2) = ny*nz*(1-cos1omega) + nx*sin1omega
  else if(twoja1 == 5) then
    Dmat(1,1) = 3.d0*ny**2*nx**2 + &
                0.5d0*(exp_4omega + exp4omega)*(1-ny**2)*(1-nx**2) + &
                0.5d0*(exp_2omega + exp2omega)*(nz**2*(1-nz**2)+(nx**2-ny**2)**2)
    Dmat(2,2) = 3.d0*ny**2*nz**2 + &
                0.5d0*(exp_4omega + exp4omega)*(1-ny**2)*(1-nz**2) + &
                0.5d0*(exp_2omega + exp2omega)*(nx**2*(1-nx**2)+(ny**2-nz**2)**2)
    Dmat(3,3) = 3.d0*nz**2*nx**2 + &
                0.5d0*(exp_4omega + exp4omega)*(1-nz**2)*(1-nx**2) + &
                0.5d0*(exp_2omega + exp2omega)*(ny**2*(1-ny**2)+(nx**2-nz**2)**2)
    Dmat(4,4) = 0.75d0*(nx**2-ny**2)**2 + &
                0.125d0*(exp_4omega + exp4omega)*(4*nz**2+(nx**2-ny**2)**2) + &
                0.5d0*(exp_2omega + exp2omega)*(nz**2*(1-nz**2) + 4*ny**2*nx**2)
    Dmat(5,5) = 0.25d0*(1-3*nz**2)**2 + &
                1.5d0*(exp_2omega + exp2omega)*nz**2*(1-nz**2) + &
                0.375d0*(exp_4omega + exp4omega)*(1-nz**2)**2

    Dmat(1,2) = -2.d0*ci*sin1omega*(ci*ny**3*cos1omega + ci*ny*(1-ny**2)*cos3omega + &
                                      nx*nz*( 3.d0*ci*ny**2*sin1omega+ci*(1-ny**2)*sin3omega))
    Dmat(2,1) =  2.d0*ci*sin1omega*(ci*ny**3*cos1omega + ci*ny*(1-ny**2)*cos3omega + &
                                      nx*nz*(-3.d0*ci*ny**2*sin1omega-ci*(1-ny**2)*sin3omega))
    Dmat(2,3) = -2.d0*ci*sin1omega*(ci*nz**3*cos1omega + ci*nz*(1-nz**2)*cos3omega + &
                                      nx*ny*( 3.d0*ci*nz**2*sin1omega+ci*(1-nz**2)*sin3omega))
    Dmat(3,2) =  2.d0*ci*sin1omega*(ci*nz**3*cos1omega + ci*nz*(1-nz**2)*cos3omega + &
                                      nx*ny*(-3.d0*ci*nz**2*sin1omega-ci*(1-nz**2)*sin3omega))
    Dmat(1,3) =  2.d0*ci*sin1omega*(ci*nx**3*cos1omega + ci*nx*(1-nx**2)*cos3omega + &
                                      ny*nz*(-3.d0*ci*nx**2*sin1omega-ci*(1-nx**2)*sin3omega))
    Dmat(3,1) = -2.d0*ci*sin1omega*(ci*nx**3*cos1omega + ci*nx*(1-nx**2)*cos3omega + &
                                      ny*nz*( 3.d0*ci*nx**2*sin1omega+ci*(1-nx**2)*sin3omega))

    Dmat(1,4) = -ci*sin1omega*(ci*nz*(3.d0-nz**2)*cos1omega + ci*nz*(1+nz**2)*cos3omega - &
                                 4.d0*ci*nx*ny*(ny**2-nx**2)*sin1omega**3)
    Dmat(4,1) =  ci*sin1omega*(ci*nz*(3.d0-nz**2)*cos1omega + ci*nz*(1+nz**2)*cos3omega + &
                                 4.d0*ci*nx*ny*(ny**2-nx**2)*sin1omega**3)
    Dmat(1,5) = -2.d0*sqrt3*sin1omega**2*(nx*ny*(1.d0-nz**2) + nx*ny*(1+nz**2)*cos2omega - &
                                               nz*(ny**2-nx**2)*sin2omega)
    Dmat(5,1) = -2.d0*sqrt3*sin1omega**2*(nx*ny*(1.d0-nz**2) + nx*ny*(1+nz**2)*cos2omega + &
                                               nz*(ny**2-nx**2)*sin2omega)

    Dmat(2,4) = 1.5d0*(nx**2 - ny**2)*ny*nz + ny*nz*(2.d0*ny**2 - 2.d0*nx**2 - 1.d0)*cos2omega + &
                0.5d0*ny*nz*(2.d0*nx**2 + nz**2 + 1.d0)*cos4omega - &
                nx*( 2.d0*ny**2 - nz**2 + (1.d0 - 2.d0*ny**2 + nz**2)*cos2omega)*sin2omega
    Dmat(3,4) = 1.5d0*(nx**2 - ny**2)*nx*nz + nx*nz*(4.d0*ny**2 + 2.d0*nz**2 - 1.d0)*cos2omega - &
                0.5d0*nx*nz*(2.d0*ny**2 + nz**2 + 1.d0)*cos4omega - &
                ny*( 2.d0*nx**2 - nz**2 + (-1.d0 + 2.d0*ny**2 + 3.d0*nz**2)*cos2omega)*sin2omega
    Dmat(4,2) = 1.5d0*(nx**2 - ny**2)*ny*nz + ny*nz*(2.d0*ny**2 - 2.d0*nx**2 - 1.d0)*cos2omega + &
                0.5d0*ny*nz*(2.d0*nx**2 + nz**2 + 1.d0)*cos4omega + &
                nx*( 2.d0*ny**2 - nz**2 + (1.d0 - 2.d0*ny**2 + nz**2)*cos2omega)*sin2omega
    Dmat(4,3) = 1.5d0*(nx**2 - ny**2)*nx*nz + nx*nz*(4.d0*ny**2 + 2.d0*nz**2 - 1.d0)*cos2omega - &
                0.5d0*nx*nz*(2.d0*ny**2 + nz**2 + 1.d0)*cos4omega + &
                ny*( 2.d0*nx**2 - nz**2 + (-1.d0 + 2.d0*ny**2 + 3.d0*nz**2)*cos2omega)*sin2omega

    Dmat(2,5) = 2.d0*ci*sqrt3*(nz**2 + (1 - nz**2)*cos2omega)*sin1omega*ci*(nx*cos1omega - ny*nz*sin1omega)
    Dmat(3,5) =-2.d0*ci*sqrt3*(nz**2 + (1 - nz**2)*cos2omega)*sin1omega*ci*(ny*cos1omega + nx*nz*sin1omega)
    Dmat(5,3) = 2.d0*ci*sqrt3*(nz**2 + (1 - nz**2)*cos2omega)*sin1omega*ci*(ny*cos1omega - nx*nz*sin1omega)
    Dmat(5,2) =-2.d0*ci*sqrt3*(nz**2 + (1 - nz**2)*cos2omega)*sin1omega*ci*(nx*cos1omega + ny*nz*sin1omega)

    Dmat(4,5) = -sqrt3*sin1omega**2*(nx**4 - ny**4 + (nx**2 - ny**2)*(1.d0 + nz**2)*cos2omega - &
                                          4.d0*nx*ny*nz*sin2omega)
    Dmat(5,4) = -sqrt3*sin1omega**2*(nx**4 - ny**4 + (nx**2 - ny**2)*(1.d0 + nz**2)*cos2omega + &
                                          4.d0*nx*ny*nz*sin2omega)
  else if (twoja1 == 7) then 
    Dmat(1,1) = 0.25d0*exp_3omega*(60d0*exp3omega*nx2*ny2*nz2 + 3d0*(1d0-nx2)*(1d0-ny2)*(1d0-nz2) + &
                3d0*exp6omega*(1d0-nx2)*(1d0-ny2)*(1d0-nz2) + &
                5d0*exp2omega*(nx4*(1d0-nx2)+(1d0-nx2)*ny2*nz2 + nx2*(ny4-6d0*ny2*nz2+nz4)) + &
                5d0*exp4omega*(nx4*(1d0-nx2)+(1d0-nx2)*ny2*nz2 + nx2*(ny4-6d0*ny2*nz2+nz4)) + &
                2d0*expomega *(nx6-nx4*(1d0-nx2)+(1d0-nx2)*(ny2-nz2)**2 - nx2*(ny4-3d0*ny2*nz2+nz4)) + &
                2d0*exp5omega*(nx6-nx4*(1d0-nx2)+(1d0-nx2)*(ny2-nz2)**2 - nx2*(ny4-3d0*ny2*nz2+nz4))   )
    Dmat(2,1) = -2d0*sqrt15*(nx2+(1d0-nx2)*cos1omega)*sinOmega_2**2*((1d0-nx2)*ny*nz+(1d0+nx2)*ny*nz*cos1omega+nx*(nz2-ny2)*sin1omega)
    Dmat(3,1) = -2d0*sqrt15*(ny2+(1d0-ny2)*cos1omega)*sinOmega_2**2*((1d0-ny2)*nx*nz+(1d0+ny2)*nx*nz*cos1omega+ny*(nx2-nz2)*sin1omega)
    Dmat(4,1) = -2d0*sqrt15*(nz2+(1d0-nz2)*cos1omega)*sinOmega_2**2*((1d0-nz2)*nx*ny+(1d0+nz2)*nx*ny*cos1omega+nz*(ny2-nx2)*sin1omega)
    Dmat(5,1) = 0.25d0*ci*exp_omega_2*(expomega-1d0)*(2d0*nx3*(2d0*nx2+5d0*(1d0-nx2))*cosOmega_2 + &
                                                      nx*(4d0*nx4+5d0*(1d0-nx2)**2)*cos3omega_2 + &
                                                      6d0*nx3*ny2*cos5omega_2 + 3d0*nx*ny4*cos5omega_2 + 6d0*nx3*nz2*cos5omega_2 + &
                                                      6d0*nx*ny2*nz2*cos5omega_2 + 3d0*nx*nz4*cos5omega_2 + &
                                                      30d0*nx2*ny*nz*(nz2-ny2)*sinOmega_2 + &
                                                      5d0*(1d0-3d0*nx2)*ny*nz*(nz2-ny2)*sin3omega_2 + &
                                                      3d0*ny*nz*(ny4-nz4)*sin5omega_2 )
    Dmat(6,1) = -0.5d0*sinOmega_2*(2d0*ny3*(5d0*nx2+2d0*ny2+5d0*nz2)*cosOmega_2 + &
                                            ny*(4d0*ny4+5d0*(1d0-ny2)**2)*cos3omega_2 + &
                                            6d0*ny3*nx2*cos5omega_2 + 3d0*ny*nx4*cos5omega_2 + 6d0*nx2*ny*nz2*cos5omega_2 + &
                                            6d0*ny3*nz2*cos5omega_2 + 3d0*ny*nz4*cos5omega_2 + &
                                            30d0*ny2*nx*nz*(nx2-nz2)*sinOmega_2 + &
                                            5d0*(1d0-3d0*ny2)*nx*nz*(nx2-nz2)*sin3omega_2 + &
                                            3d0*nx*nz*(nz4-nx4)*sin5omega_2 )
    Dmat(7,1) = -0.5d0*sinOmega_2*(2d0*nz3*(2d0*nz2+5d0*(1d0-nz2))*cosOmega_2 + &
                                   nz*(4d0*nz4+5d0*(1d0-nz2)**2)*cos3omega_2 + &
                                   6d0*nz3*nx2*cos5omega_2 + 3d0*nz*nx4*cos5omega_2 + 6d0*nx2*nz*ny2*cos5omega_2 + &
                                   6d0*nz3*ny2*cos5omega_2 + 3d0*nz*ny4*cos5omega_2 + &
                                   30d0*nz2*nx*ny*(nx2-ny2)*sinOmega_2 + &
                                   5d0*(1d0-3d0*nz2)*nx*ny*(ny2-nx2)*sin3omega_2 + &
                                   3d0*nx*ny*(nx4-ny4)*sin5omega_2) 
    Dmat(1,2) = -2d0*sqrt15*(nx2+(1d0-nx2)*cos1omega)*sinOmega_2**2*((1d0-nx2)*ny*nz+(1d0+nx2)*ny*nz*cos1omega+nx*(ny2-nz2)*sin1omega)
    Dmat(2,2) = 1d0/16d0*exp_3omega*( 3d0*exp2omega*(1d0-5d0*nx2)**2*(1d0-nx2) + 3d0*exp4omega*(1d0-5d0*nx2)**2*(1d0-nx**2) + &
                                       30d0*expomega*nx2*(1d0-nx2)**2 + 30d0*exp5omega*nx2*(1d0-nx2)**2 + &
                                       5d0*(1d0-nx2)**3 + 5d0*exp6omega*(1d0-nx2)**3 + 4d0*exp3omega*nx2*(2d0*nx2-3d0*(1d0-nx2))**2 )
    Dmat(3,2) = -0.25d0*sinOmega_2*(2d0*nz*(6d0*nx4-3d0*nx2*ny2+6d0*ny4+nz4+7d0*nz2*(1d0-nz2))*cosOmega_2 + &
                                    5d0*nz*(9d0*nx2*ny2+nz4+nz2*(1d0-nz2))*cos3omega_2 - 15d0*nx2*ny2*nz*cos5omega_2 + &
                                    5d0*nx2*nz3*cos5omega_2 + 5d0*ny2*nz3*cos5omega_2 + 5d0*nz5*cos5omega_2 + &
                                    2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*ny*(3d0*nx2-2d0*ny2+3d0*nz2)*sinOmega_2 + &
                                    5d0*nx*ny*(5d0*nx2*ny2-3d0*nz4-3d0*nz2*(1d0-nz2))*sin3omega_2 + &
                                    5d0*nx*ny*(-nx2*ny2+3d0*nz4+3d0*nz2*(1-nz2))*sin5omega_2 )
    Dmat(4,2) = 0.25d0*sinOmega_2*(2d0*ny*(6d0*nx4+nx2*(7d0*ny2-3d0*nz2)+(1d0-nx2)*(ny2+6d0*nz2))*cosOmega_2 + &
                                   5d0*ny*((1-nx2)*ny2+nx2*(ny2+9d0*nz2))*cos3omega_2 - 15d0*nx2*ny*nz2*cos5omega_2 + &
                                   5d0*nx2*ny3*cos5omega_2 + 5d0*ny5*cos5omega_2 + 5d0*ny3*nz2*cos5omega_2 + &
                                   (-2d0)*nx*(2d0*nx2-3d0*(1d0-nx2))*nz*(-2d0*nz2+3d0*(1d0-nz2))*sinOmega_2 + &
                                   5d0*nx*nz*(3d0*(1d0-nx2)*ny2+nx2*(3d0*ny2-5d0*nz2))*sin3omega_2 + &
                                   5d0*nx*nz*(-3d0*(1d0-nx2)*ny2 + nx2*(-3d0*ny2+nz2))*sin5omega_2 )
    Dmat(5,2) = -sqrt15*(nx2+(1d0-nx2)*cos1omega)*sinOmega_2**2*(ny4-nz4+(1d0+nx2)*(ny2-nz2)*cos1omega-4d0*nx*ny*nz*sin1omega)
    Dmat(6,2) = 0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*(2d0*nz*(2d0*nx4+3d0*nx2*ny2+nz4+nz2*(1d0-nz2))*cosOmega_2 + &
                                                              nz*(-nx2*(ny2-7d0*nz2)+(1d0-nx2)*(2d0*ny2+nz2))*cos3omega_2 + &
                                                              3d0*nx2*ny2*nz*cos5omega_2 + 2d0*ny4*nz*cos5omega_2 + &
                                                            (-1d0)*nx2*nz3*cos5omega_2 + 3d0*ny2*nz3*cos5omega_2 + nz5*cos5omega_2 + &
                                                             2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*ny*(nx2-nz2)*sinOmega_2 + &
                                                             nx*ny*(-2d0*ny4-ny2*nz2+nz4+nx2*(3d0*ny2+11d0*nz2))*sin3omega_2 + &
                                                             nx*ny*(2d0*ny4+ny2*nz2-nz4+nx2*(ny2-3d0*nz2))*sin5omega_2 )
    Dmat(7,2) = 0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*(2d0*ny*(2d0*nx4+nx2*ny2+ny4+(3d0*nx2+ny2)*nz2)*cosOmega_2 + &
                                                              ny*(-nx2*(nz2-7d0*ny2)+(1d0-nx2)*(2d0*nz2+ny2))*cos3omega_2 + &
                                                              3d0*nx2*ny*nz2*cos5omega_2 + 2d0*ny*nz4*cos5omega_2 + &
                                                            (-1d0)*nx2*ny3*cos5omega_2 + 3d0*ny3*nz2*cos5omega_2 + ny5*cos5omega_2 + &
                                                             2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*nz*(ny2-nx2)*sinOmega_2 + &
                                                            (-1d0)*nx*nz*((1d0-nx2)*(ny2-2d0*nz2)+nx2*(11d0*ny2+3d0*nz2))*sin3omega_2 + &
                                                             nx*nz*(3d0*nx2*ny2+ny4-2d0*nz4-nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(1,3) = -2d0*sqrt15*(ny2+(1d0-ny2)*cos1omega)*sinOmega_2**2*(nx*(1d0-ny2)*nz+nx*(1d0+ny2)*nz*cos1omega+ny*(nz2-nx2)*sin1omega)
    Dmat(2,3) = 0.25d0*sinOmega_2*(2d0*nz*(6d0*nx4-3d0*nx2*ny2+6d0*ny4+nz4+7d0*nz2*(1d0-nz2))*cosOmega_2 + &
                                   5d0*nz*(9d0*nx2*ny2+nz4+nz2*(1d0-nz2))*cos3omega_2 - 15d0*nx2*ny2*nz*cos5omega_2 + &
                                   5d0*nz3*cos5omega_2 - & 
                                   2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*ny*(3d0*nx2-2d0*ny2+3d0*nz2)*sinOmega_2 + &
                                   5d0*nx*ny*(-5d0*nx2*ny2+3d0*nz4+3d0*nz2*(1d0-nz2))*sin3omega_2 + &
                                   5d0*nx*ny*(nx2*ny2-3d0*nz4-3d0*nz2*(1-nz2))*sin5omega_2 )
    Dmat(3,3) = 1d0/16d0*exp_3omega*(3d0*exp2omega*(1d0-5d0*ny2)**2*(1d0-ny2) + 3d0*exp4omega*(1d0-5d0*ny2)**2*(1d0-ny2) + &
                                     30d0*expomega*ny2*(1d0-ny2)**2 + 30d0*exp5omega*ny2*(1d0-ny2)**2 + 5d0*(1d0-ny2)**3 + &
                                     5d0*exp6omega*(1d0-ny2)**3 + 4d0*exp3omega*ny2*(3d0*nx2-2d0*ny2+3d0*nz2)**2 )
    Dmat(4,3) = 0.125d0*(4d0*ny*nz*(3d0*(3d0*nx4+nx2*(1d0-nx2)-2d0*ny4 + ny2*nz2 - 2d0*nz4) -20d0*ny2*nz2*cos1omega - &
                                    5d0*(3d0*nx4+3d0*nx2*(1d0-nx2)-ny2*nz2)*cos2omega)*sinOmega_2**2 - &
                                    2d0*nx*(nx4+7d0*nx2*(1d0-nx2)+6d0*ny4-33d0*ny2*nz2+6d0*nz4+60d0*ny2*nz2*cos1omega + &
                                    5d0*(nx4+nx2*(1d0-nx2)-3d0*ny2*nz2)*cos2omega) * sin1omega )
    Dmat(5,3) = 1d0/16d0*sqrt15*exp_3omega*(-2d0*(expomega-1d0)**2*(1d0+exp4omega)*nx5*ny + &
                                            2d0*ci*nx4*nz*(-1d0+exp4omega*(2d0*cos2omega-1d0)) + &
                                        8d0*exp3omega*nx3*ny*(-ny2+3d0*nz2+4d0*ny2*cos1omega+(1d0-nx2)*cos2omega)*sinOmega_2**2 - &
                 8d0*exp3omega*nx*ny*(-2d0*ny4+ny2*nz2-3d0*nz4+nz2*(-8d0*ny2*cos1omega + (3d0*ny2+nz2)*cos2omega))*sinOmega_2**2- &
                 4d0*exp3omega*nx2*nz*(5d0*ny2+nz2-4d0*ny2*cos1omega+3d0*(1d0-nx2)*cos2omega)*sin1omega - &
                 4d0*exp3omega*nz*(8d0*ny2*nz2*cos1omega + (ny2-nz2)*(2d0*ny2-nz2-nz2*cos2omega))*sin1omega )
    Dmat(6,3) = sqrt15*(ny2+(1d0-ny2)*cos1omega)*sinOmega_2**2 * (nx4-nz4+(1+ny2)*(nx2-nz2)*cos1omega+4d0*nx*ny*nz*sin1omega)
    Dmat(7,3) = 0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)* &
                ( 2d0*nx*cosOmega_2 * ( nx4 - 3d0*nx2*ny2 + 2d0*ny4 + nz2*(1d0+4d0*ny2-nz2) + 4d0*ny2*(2d0*nx2-nz2)*cos1omega + &
                                      (nx2*(nx2-ny2)+2d0*nz4+3d0*nz2*(1d0-nz2))*cos2omega ) + &
                  ny*nz * (2d0*(nx2-ny2)*(3d0*nx2-2d0*ny2+3d0*nz2)*sinOmega_2 + &
                          (nx4 + 3d0*ny2*nz2 - 2d0*nz4 + nx2*(11d0*ny2-nz2))*sin3omega_2 + &
                          (2d0*nz4 + nz2*(1d0-nz2) - nx2*(1d0+2d0*ny2-nz2))*sin5omega_2 ))
    Dmat(1,4) = -2d0*sqrt15*(nz2+(1d0-nz2)*cos1omega)*sinOmega_2**2*(nx*ny*(1d0-nz2)+nx*ny*(1d0+nz2)*cos1omega+(nx2-ny2)*nz*sin1omega)
    Dmat(2,4) = -0.25d0*sinOmega_2*(2d0*ny*(6d0*nx4+nx2*(7d0*ny2-3d0*nz2) + (1d0-nx2)*(ny2+6d0*nz2))*cosOmega_2 + &
                                    5d0*ny*((1d0-nx2)*ny2+nx2*(ny2+9d0*nz2))*cos3omega_2 + 5d0*nx2*ny3*cos5omega_2 + &
                                    5d0*ny5*cos5omega_2 -15d0*nx2*ny*nz2*cos5omega_2 + 5d0*ny3*nz2*cos5omega_2 + &
                                    2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*nz*(-2d0*nz2+3d0*(1d0-nz2))*sinOmega_2 - &
                                    5d0*nx*nz*(3d0*(1d0-nx2)*ny2+nx2*(3d0*ny2-5d0*nz2))*sin3omega_2 + &
                                    5d0*nx*nz*(3d0*(1d0-nx2)*ny2+nx2*(3d0*ny2-nz2))*sin5omega_2 )
    Dmat(3,4) = 0.125d0*(4d0*ny*nz*(3d0*(3d0*nx4+nx2*(1d0-nx2)-2d0*ny4+ny2*nz2-2d0*nz4) - 20d0*ny2*nz2*cos1omega - &
                                     5d0*(3d0*nx4+3d0*nx2*(1d0-nx2)-ny2*nz2)*cos2omega)*sinOmega_2**2 + &
                          2d0*nx*(nx4+7d0*nx2*(1d0-nx2)+6d0*ny4-33d0*ny2*nz2+6d0*nz4+60*ny2*nz2*cos1omega + &
                          5d0*(nx4+nx2*(1d0-nx2)-3d0*ny2*nz2)*cos2omega)*sin1omega )
    Dmat(4,4) = 1d0/16d0*exp_3omega*(3d0*exp2omega*(1d0-5d0*nz2)**2*(1d0-nz2)+3d0*exp4omega*(1d0-5d0*nz2)**2*(1d0-nz2) + &
                                     30d0*expomega*nz2*(1d0-nz2)**2+30d0*exp5omega*nz2*(1d0-nz2)**2 + &
                                     5d0*(1d0+exp6omega)*(1d0-nz2)**3 + 4d0*exp3omega*(-2d0*nz3+3d0*nz*(1d0-nz2))**2) 
    Dmat(5,4) = 0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*(2d0*ny*cosOmega_2*(ny4-3d0*ny2*nz2+2d0*nz4+nx2*(ny2+5d0*nz2) - &
                                        4d0*nz2*(1d0-3d0*ny2-nz2)*cos1omega + (2d0*nx4+3d0*nx2*(1d0-nx2)+ny2*(ny2-nz2))*cos2omega)+&
                                        nx*nz*(2d0*(ny2-nz2)*(3d0-5d0*nz2)*sinOmega_2 + &
                                        (-2d0*nx4+ny4+11d0*ny2*nz2-nx2*(ny2-3d0*nz2))*sin3omega_2 + &
                                        (2d0*nx4+nx2*(1d0-nx2)-ny2*(ny2+3d0*nz2))*sin5omega_2))
    Dmat(6,4) = 0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*(2d0*nx*cosOmega_2*(nx4-3d0*nx2*nz2+2d0*nz4+ny2*(nx2+5d0*nz2) + &
                                        4d0*nz2*(2d0*nx2-ny2)*cos1omega + (nx4+2d0*ny4+3d0*ny2*nz2+nx2*(3d0*ny2-nz2))*cos2omega)+&
                                        ny*nz*(-2d0*(nx2-nz2)*(3d0-5d0*nz2)*sinOmega_2 - &
                                        (nx4-nx2*ny2-2d0*ny4+(11d0*nx2+3d0*ny2)*nz2)*sin3omega_2 + &
                                        (nx4-nx2*ny2-2d0*ny4+nz2*(3d0*nx2-ny2))*sin5omega_2))
    Dmat(7,4) = -sqrt15*(nz2+(1d0-nz2)*cos1omega)*sinOmega_2**2*(nx4-ny4+(nx2-ny2)*(1d0+nz2)*cos1omega-4d0*nx*ny*nz*sin1omega)
    Dmat(1,5) = -0.25d0*ci*exp_omega_2*(expomega-1d0)*(2d0*nx3*(2d0*nx2+5d0*(1d0-nx2))*cosOmega_2 + &
                                                      nx*(4d0*nx4+5d0*(1d0-nx2)**2)*cos3omega_2 + &
                                                      6d0*nx3*ny2*cos5omega_2 + 3d0*nx*ny4*cos5omega_2 + 6d0*nx3*nz2*cos5omega_2 + &
                                                      6d0*nx*ny2*nz2*cos5omega_2 + 3d0*nx*nz4*cos5omega_2 + &
                                                      40d0*nx2*ny*nz*(ny2-nz2)*sinOmega_2**3 + &
                                                      8d0*nz*(ny5-ny*nz4)*(2d0+3d0*cos1omega)*sinOmega_2**3)
    Dmat(2,5) = -sqrt15*(nx2+(1d0-nx2)*cos1omega)*sinOmega_2**2*(ny4-nz4+(1d0+nx2)*(ny2-nz2)*cos1omega+4d0*nx*ny*nz*sin1omega)
    Dmat(3,5) = 0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*&
        (-2d0*nz*cosOmega_2*(2d0*ny4-3d0*ny2*nz2+nz4+nx2*(5d0*ny2+nz2)-4d0*ny2*(1d0-ny2-3d0*nz2)*cos1omega + &
         (2d0*nx4+3d0*nx2*(1d0-nx2)-ny2*nz2+nz4)*cos2omega) + &
         nx*ny*(2d0*(ny2-nz2)*(3d0*nx2-2d0*ny2+3d0*nz2)*sinOmega_2 + (2d0*nx4+nx2*(-3d0*ny2+nz2)-nz2*(11d0*ny2+nz2))*sin3omega_2 + &
         (-2d0*nx4-nx2*(1d0-nx2)+3d0*ny2*nz2+nz4)*sin5omega_2 ))
    Dmat(4,5) = 1d0/16d0*sqrt15*exp_3omega*( 2d0*(expomega-1d0)**2*(1d0+exp4omega)*nx5*nz - &
                                              2d0*ci*nx4*ny*(exp4omega*(2d0*cos2omega-1d0)-1d0) - &
                                    8d0*exp3omega*nx3*nz*(3d0*ny2-nz2+4d0*nz2*cos1omega + (1d0-nx2)*cos2omega)*sinOmega_2**2 + &
                8d0*exp3omega*nx*nz*(-3d0*ny4+ny2*nz2-2d0*nz4+ny2*(-8d0*nz2*cos1omega+(ny2+3d0*nz2)*cos2omega))*sinOmega_2**2 + &
                4d0*exp3omega*nx2*ny*(ny2+5d0*nz2-4d0*nz2*cos1omega+3d0*(1d0-nx2)*cos2omega)*sin1omega + &
                4d0*exp3omega*ny*(8d0*ny2*nz2*cos1omega+(ny2-nz2)*(ny2-2d0*nz2+ny2*cos2omega))*sin1omega )
    Dmat(5,5) = 1d0/16d0*exp_3omega*(60d0*exp3omega*nx2*(ny2-nz2)**2+3d0*(1d0-nx2)*(4d0*nx4+4d0*nx2*(1d0-nx2)+(ny2-nz2)**2) + &
                                    3d0*exp6omega*(1d0-nx2)*(4d0*nx4+4d0*nx2*(1d0-nx2)+(ny2-nz2)**2) + &
                            2d0*expomega*((-2d0*nx3+nx*ny2)**2-2d0*(2d0*nx4+9d0*nx2*ny2-8d0*ny4)*nz2 + nz4*(1d0+15d0*ny2-nz2)) + &
                            2d0*exp5omega*((-2d0*nx3+nx*ny2)**2-2d0*(2d0*nx4+9d0*nx2*ny2-8d0*ny4)*nz2+ nz4*(1d0+15d0*ny2-nz2)) + &
                            5d0*exp2omega*(4d0*nx4*(1d0-nx2)+(1d0-nx2)*(ny2-nz2)**2 - 4d0*nx2*(ny4-6d0*ny2*nz2+nz4)) + &
                            5d0*exp4omega*(4d0*nx4*(1d0-nx2)+(1d0-nx2)*(ny2-nz2)**2 - 4d0*nx2*(ny4-6d0*ny2*nz2+nz4)))
    Dmat(6,5) = -0.25d0*sinOmega_2*( 2*nz*(-15d0*nx2*ny2+nz4-5d0*nz2*(1d0-nz2))*cosOmega_2 + &
                                    nz*(-5d0*(2d0*nx4-5d0*nx2*ny2+2d0*ny4)-3d0*nz4+5d0*nz2*(1d0-nz2))*cos3omega_2 + &
                                    6d0*nx4*nz*cos5omega_2 -3d0*nx2*ny2*nz*cos5omega_2 + 6d0*ny4*nz*cos5omega_2 - &
                                    3d0*nx2*nz3*cos5omega_2 - &
                                    3d0*ny2*nz3*cos5omega_2 -3d0*nz5*cos5omega_2 + 30d0*nx*ny*(nx2-nz2)*(ny2-nz2)*sinOmega_2 + &
                                    5d0*nx*ny*(2d0*nx4+nx2*ny2+2*ny4-3d0*nz4+5d0*nz2*(1d0-nz2))*sin3omega_2 + &
                                    3d0*nx*ny*(2d0*nx4+5d0*nx2*ny2+2d0*ny4+5d0*nz4+5d0*nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(7,5) = -0.25d0*sinOmega_2*( 2d0*ny*(-ny4+5d0*ny2*nz2+5d0*nx2*(ny2+3d0*nz2))*cosOmega_2 + &
                                    ny*(10d0*nx4-5d0*nx2*ny2+3d0*ny4-5d0*(5d0*nx2+ny2)*nz2+10d0*nz4)*cos3omega_2 - &
                                    6d0*nx4*ny*cos5omega_2 + 3d0*nx2*ny*nz2*cos5omega_2 - 6d0*ny*nz4*cos5omega_2 + &
                                    3d0*nx2*ny3*cos5omega_2 + &
                                    3d0*ny3*nz2*cos5omega_2 + 3d0*ny5*cos5omega_2 - 30d0*nx*nz*(nx2-ny2)*(ny2-nz2)*sinOmega_2 + &
                                    5d0*nx*nz*(2d0*nx4+5d0*nx2*ny2+2*nz4-3d0*ny4+nz2*(1d0+4d0*ny2-nz2))*sin3omega_2 + &
                                    3d0*nx*nz*(2d0*nx4+5d0*nx2*ny2+2d0*nz4+5d0*ny4+5d0*nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(1,6) = 0.5d0*sinOmega_2*(2d0*ny3*(5d0*nx2+2d0*ny2+5d0*nz2)*cosOmega_2 + &
                                            ny*(4d0*ny4+5d0*(1d0-ny2)**2)*cos3omega_2 + &
                                            6d0*ny3*nx2*cos5omega_2 + 3d0*ny*nx4*cos5omega_2 + 6d0*nx2*ny*nz2*cos5omega_2 + &
                                            6d0*ny3*nz2*cos5omega_2 + 3d0*ny*nz4*cos5omega_2 + &
                                            30d0*ny2*nx*nz*(nz2-nx2)*sinOmega_2 + &
                                            5d0*(1d0-3d0*ny2)*nx*nz*(nz2-nx2)*sin3omega_2 + &
                                            3d0*nx*nz*(nx4-nz4)*sin5omega_2 )
    Dmat(2,6) = 0.25d0*sqrt15*sinOmega_2*( 2d0*nz*(2d0*nx4+3d0*nx2*ny2+nz4+nz2*(1d0-nz2))*cosOmega_2 + &
                                           nz*(-nx2*(ny2-7d0*nz2)+(1d0-nx2)*(2d0*ny2+nz2))*cos3omega_2 + &
                                           3d0*nx2*ny2*nz*cos5omega_2 + 2d0*ny4*nz*cos5omega_2 - nx2*nz3*cos5omega_2 + &
                                           3d0*ny2*nz3*cos5omega_2 + nz5*cos5omega_2 - &
                                           2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*ny*(nx2-nz2)*sinOmega_2 - &
                                           nx*ny*(-2d0*ny4-ny2*nz2+nz4+nx2*(3d0*ny2+11d0*nz2))*sin3omega_2 + &
                                           nx*ny*(-2d0*ny4-ny2*nz2+nz4-nx2*(ny2-3d0*nz2))*sin5omega_2 )
    Dmat(3,6) = sqrt15*(ny2+(1d0-ny2)*cos1omega)*sinOmega_2**2*(nx4-nz4+(1d0+ny2)*(nx2-nz2)*cos1omega-4d0*nx*ny*nz*sin1omega)
    Dmat(4,6) =-0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*(2d0*nx*cosOmega_2*(nx4-3d0*nx2*nz2+2d0*nz4+ny2*(nx2+5d0*nz2) + &
                                        4d0*nz2*(2d0*nx2-ny2)*cos1omega + (nx4+2d0*ny4+3d0*ny2*nz2+nx2*(3d0*ny2-nz2))*cos2omega)+&
                                        ny*nz*(2d0*(nx2-nz2)*(3d0-5d0*nz2)*sinOmega_2 + &
                                        (nx4-nx2*ny2-2d0*ny4+(11d0*nx2+3d0*ny2)*nz2)*sin3omega_2 - &
                                        (nx4-nx2*ny2-2d0*ny4+nz2*(3d0*nx2-ny2))*sin5omega_2))
    Dmat(5,6) = -0.25d0*sinOmega_2*(-2*nz*(-15d0*nx2*ny2+nz4-5d0*nz2*(1d0-nz2))*cosOmega_2 + &
                                    nz*(5d0*(2d0*nx4-5d0*nx2*ny2+2d0*ny4)+3d0*nz4-5d0*nz2*(1d0-nz2))*cos3omega_2 - &
                                    6d0*nx4*nz*cos5omega_2 +3d0*nx2*ny2*nz*cos5omega_2 - 6d0*ny4*nz*cos5omega_2 + &
                                    3d0*nx2*nz3*cos5omega_2 + &
                                    3d0*ny2*nz3*cos5omega_2 +3d0*nz5*cos5omega_2 + 30d0*nx*ny*(nx2-nz2)*(ny2-nz2)*sinOmega_2 + &
                                    5d0*nx*ny*(2d0*nx4+nx2*ny2+2*ny4-3d0*nz4+5d0*nz4*(1d0-nz2))*sin3omega_2 + &
                                    3d0*nx*ny*(2d0*nx4+5d0*nx2*ny2+2d0*ny4+5d0*nz4+5d0*nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(6,6) = 1d0/16d0*exp_3omega*(60d0*exp3omega*ny2*(1d0-ny2-2d0*nz2)**2 + &
                                     3d0*(1d0-ny2)*(nx4+nx2*(4d0*ny2-2d0*nz2)+(2d0*ny2+nz2)**2) + &
                                     3d0*exp6omega*(1d0-ny2)*(nx4+nx2*(4d0*ny2-2d0*nz2)+(2d0*ny2+nz2)**2) + &
                                2d0*expomega*(nx4*(ny2+16d0*nz2)+(-2d0*ny3+ny*nz2)**2-2d0*nx2*(2d0*ny4+9d0*ny2*nz2-8d0*nz4)) + &
                                2d0*exp5omega*(nx4*(ny2+16d0*nz2)+(-2d0*ny3+ny*nz2)**2-2d0*nx2*(2d0*ny4+9d0*ny2*nz2-8d0*nz4)) + &
                    5d0*(exp2omega+exp4omega)*(nx6-nx4*(4d0*ny2+nz2)+(-2d0*ny2*nz+nz3)**2+nx2*(4d0*ny4+24d0*ny2*nz2-nz4)) )
    Dmat(7,6) = -0.25d0*sinOmega_2*(2d0*nx*(nx4-5d0*nx2*(1d0-nx2)-15d0*ny2*nz2)*cosOmega_2 + &
                                    nx*(-3d0*nx4+5d0*nx2*(1d0-nx2)-5d0*(2d0*ny4-5d0*ny2*nz2+2d0*nz4))*cos3omega_2 + &
                                    cos5omega_2*(-3d0*nx5-3d0*nx3*ny2+6d0*nx*ny4-3d0*nx3*nz2-3d0*nx*ny2*nz2+6d0*nx*nz4) + &
                                    30d0*ny*nz*(nx2-ny2)*(nx2-nz2)*sinOmega_2 + &
                                    5d0*ny*nz*(-3d0*nx4+5d0*nx2*(1d0-nx2)+2d0*ny4+ny2*nz2+2d0*nz4)*sin3omega_2 + &
                                    3d0*ny*nz*(5d0*nx4+5d0*nx2*ny2+2d0*ny4+2d0*nz4+5d0*nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(1,7) =  0.5d0*sinOmega_2*(2d0*nz3*(2d0*nz2+5d0*(1d0-nz2))*cosOmega_2 + &
                                   nz*(4d0*nz4+5d0*(1d0-nz2)**2)*cos3omega_2 + &
                                   6d0*nz3*nx2*cos5omega_2 + 3d0*nz*nx4*cos5omega_2 + 6d0*nx2*nz*ny2*cos5omega_2 + &
                                   6d0*nz3*ny2*cos5omega_2 + 3d0*nz*ny4*cos5omega_2 + &
                                   30d0*nz2*nx*ny*(nx2-ny2)*sinOmega_2 + &
                                   5d0*(1d0-3d0*nz2)*nx*ny*(nx2-ny2)*sin3omega_2 + &
                                   3d0*nx*ny*(ny4-nx4)*sin5omega_2) 
    Dmat(2,7) =-0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*(2d0*ny*(2d0*nx4+nx2*ny2+ny4+(3d0*nx2+ny2)*nz2)*cosOmega_2 + &
                                                              ny*(-nx2*(nz2-7d0*ny2)+(1d0-nx2)*(2d0*nz2+ny2))*cos3omega_2 + &
                                                              3d0*nx2*ny*nz2*cos5omega_2 + 2d0*ny*nz4*cos5omega_2 + &
                                                            (-1d0)*nx2*ny3*cos5omega_2 + 3d0*ny3*nz2*cos5omega_2 + ny5*cos5omega_2 + &
                                                             2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*nz*(nx2-ny2)*sinOmega_2 + &
                                                            nx*nz*((1d0-nx2)*(ny2-2d0*nz2)+nx2*(11d0*ny2+3d0*nz2))*sin3omega_2 + &
                                                             nx*nz*(-3d0*nx2*ny2-ny4+2d0*nz4+nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(3,7) = 1d0/16d0*sqrt15*exp_3omega*(-ci*(exp2omega-1d0)*(exp2omega+1d0)**2*nx5 - &
                                            (expomega-1d0)**2*(1d0-6d0*exp2omega+exp4omega)*nx4*ny*nz + &
                                8d0*exp3omega*nx2*ny*nz*(ny2-3d0*nz2-8d0*ny2*cos1omega+(3d0*ny2-nz2)*cos2omega)*sinOmega_2**2 - &
                            8d0*exp3omega*ny*nz*(ny2*(2d0*ny2-nz2+4d0*nz2*cos1omega)+nz2*(ny2+2d0*nz2)*cos2omega)*sinOmega_2**2 + &
                            4d0*exp3omega*nx3*(-3d0*ny2+nz2+8d0*ny2*cos1omega-(ny2-3d0*nz2)*cos2omega)*sin1omega + &
                            4d0*exp3omega*nx*(ny2*(2d0*ny2+5d0*nz2-4d0*nz2*cos1omega)+nz2*(3d0*ny2+2d0*nz2)*cos2omega)*sin1omega )
    Dmat(4,7) = -sqrt15*(nz2+(1d0-nz2)*cos1omega)*sinOmega_2**2*(nx4-ny4+(nx2-ny2)*(1d0+nz2)*cos1omega+4d0*nx*ny*nz*sin1omega)
    Dmat(5,7) = -0.25d0*sinOmega_2*(-2d0*ny*(-ny4+5d0*ny2*nz2+5d0*nx2*(ny2+3d0*nz2))*cosOmega_2 - &
                                    ny*(10d0*nx4-5d0*nx2*ny2+3d0*ny4-5d0*(5d0*nx2+ny2)*nz2+10d0*nz4)*cos3omega_2 + &
                                    6d0*nx4*ny*cos5omega_2 - 3d0*nx2*ny*nz2*cos5omega_2 + 6d0*ny*nz4*cos5omega_2 - &
                                    3d0*nx2*ny3*cos5omega_2 - &
                                    3d0*ny3*nz2*cos5omega_2 - 3d0*ny5*cos5omega_2 - 30d0*nx*nz*(nx2-ny2)*(ny2-nz2)*sinOmega_2 + &
                                    5d0*nx*nz*(2d0*nx4+5d0*nx2*ny2+2*nz4-3d0*ny4+nz2*(1d0+4d0*ny2-nz2))*sin3omega_2 + &
                                    3d0*nx*nz*(2d0*nx4+5d0*nx2*ny2+2d0*nz4+5d0*ny4+5d0*nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(6,7) = -0.25d0*sinOmega_2*(-2d0*nx*(nx4-5d0*nx2*(1d0-nx2)-15d0*ny2*nz2)*cosOmega_2 + &
                                    nx*(3d0*nx4-5d0*nx2*(1d0-nx2)+5d0*(2d0*ny4-5d0*ny2*nz2+2d0*nz4))*cos3omega_2 + &
                                    cos5omega_2*(3d0*nx5+3d0*nx3*ny2-6d0*nx*ny4+3d0*nx3*nz2+3d0*nx*ny2*nz2-6d0*nx*nz4) + &
                                    30d0*ny*nz*(nx2-ny2)*(nx2-nz2)*sinOmega_2 + &
                                    5d0*ny*nz*(-3d0*nx4+5d0*nx2*(1d0-nx2)+2d0*ny4+ny2*nz2+2d0*nz4)*sin3omega_2 + &
                                    3d0*ny*nz*(5d0*nx4+5d0*nx2*ny2+2d0*ny4+2d0*nz4+5d0*nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(7,7) = 1d0/16d0*exp_3omega*(60d0*exp3omega*nz2*(1d0-2d0*ny2-nz2)**2 + &
                                     2d0*expomega*((nx4-18d0*nx2*ny2+ny4)*nz2+4d0*nz6+16d0*nx2*ny2*(1d0-nz2)-4d0*nz4*(1d0-nz2)) + &
                                     2d0*exp5omega*((nx4-18d0*nx2*ny2+ny4)*nz2+4d0*nz6+16d0*nx2*ny2*(1d0-nz2)-4d0*nz4*(1d0-nz2)) + &
                                     3d0*(1d0+exp6omega)*(1d0-nz2)*(4d0*nz4+4d0*nz2*(1d0-nz2)+(1d0-2d0*ny2-nz2)**2) + &
                        5d0*(exp2omega+exp4omega)*(-4d0*(nx4-6d0*nx2*ny2+ny4)*nz2+4d0*nz4*(1d0-nz2)+(1d0-nz2)*(1d0-2d0*ny2-nz2)**2) )

  else
    write(*,*) " This subroutine can only deal with J = 1/2, 1, 2, 3 !!! "
  end if
end subroutine

end module
