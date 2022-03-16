!=========================================================================!
! project : struct_data
! history : 06/22/2016
! authors : Zhijun Wang  ( zjwang11@hotmail.com )
! purpose : Get input from OUTCAR
! status  : good  
! comment : These programs are distributed in the hope that they will be 
!           useful, but WITHOUT ANY WARRANTY; without even the implied 
!           warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
!=========================================================================!
module  struct_data
implicit none
private

integer   ,  public , parameter :: dp  =  8
integer   ,  public , save :: NMAT  ! - max. no. of plane waves
integer   ,  public , save :: NUME  ! - max. no. of eigenvalues(bands)
integer   ,  public , save :: NKPTS ! - max. no. of kpoints
integer   ,  public , save :: NSPIN ! - spin poralorization

public   ::  init

CONTAINS

subroutine init( NSYM,FLMAX,FL,IORD,IZ,TAU,SU2 )
  implicit none
  integer  ,intent(IN)   :: NSYM
  integer  ,intent(IN)   :: FLMAX  ! FLMAX - size of flag (FL) array (should be 4)
  logical  ,intent(OUT)  :: FL(FLMAX)
  integer  ,intent(OUT)  :: IORD
  integer  ,intent(OUT)  :: IZ(3,3,NSYM)
  real(dp) ,intent(OUT)  :: TAU(3,NSYM)
  complex(dp),intent(OUT):: SU2(2,2,NSYM)        

  integer  , parameter   :: nfst=20
  logical    :: FL1
  real(dp)   :: br2(3,3),br4(3,3)
  real(dp)   :: DZ2(3,3,NSYM)  

  integer :: irot,i,j
  integer :: I1,I2

  character*120 :: chaps
  character*10  :: title
  character*5   :: chtp5
  character*15  :: chtp15
  character*35  :: chtp35
  integer       :: itmp,ierr
  real(dp)      :: rtmp(8),rotmt(3,3)
  complex(dp)   :: crotmt(3,3),srotmt(2,2)

 !**********************************************************
  FL=.FALSE.;FL(1)=.TRUE.;FL1=.FALSE.
  IZ=0;TAU=0.d0;DZ2=0.d0
!
  OPEN(unit=nfst,file='OUTCAR',form='formatted',status='old')
!
! set IORD,TAU,DZ2
  do 
   read(nfst,"(A90)") chaps
   chtp5=chaps(1:5)
   if(chtp5 =='Space' ) exit
  enddo
   read(nfst,*) 
     srotmt=0.d0
  do irot=1,NSYM
     read(nfst,*,iostat=ierr)  itmp,rtmp
     if (ierr /= 0) exit
     IORD =irot
     TAU(:,irot)=rtmp(6:8)
     IF( abs( rtmp(1)+1.d0+rtmp(2)   ).lt.1.d-5 ) FL(3)=.TRUE.
     IF( abs( rtmp(6)+rtmp(7)+rtmp(8)).gt.1.d-5 ) FL1  =.TRUE.
     call Dmatrix(rnx=rtmp(3), rny=rtmp(4), rnz=rtmp(5), degree=-rtmp(2), twoja1=3, Dmat=crotmt)
     call Dmatrix(rnx=rtmp(3), rny=rtmp(4), rnz=rtmp(5), degree=-rtmp(2), twoja1=2, Dmat=srotmt)
     rotmt=real(crotmt)
     IF(rtmp(1).lt. 1.d-2) rotmt=-rotmt
     DZ2(:,:,irot)= rotmt(:,:)
     SU2(:,:,irot)=srotmt(:,:)
  enddo
!
! read NKPTS NUME
  do 
   read(nfst,"(A120)") chaps
   chtp15=chaps(4:11)
   if(chtp15 =='k-points       ' ) exit
  enddo
   chtp15=chaps(32:38)  ; read(chtp15,*) NKPTS
   chtp15=chaps(105:115); read(chtp15,*) NUME
!
! read title
  do 
   read(nfst,"(A90)") chaps
   chtp15=chaps(2:7)
   if(chtp15 =='SYSTEM         ' ) exit
  enddo
   read(chaps,*) title,title,title
!
! set FL(:),FL1
  do 
   read(nfst,"(A90)") chaps
   chtp5=chaps(4:8)
   if(chtp5 =='ISPIN' ) exit
  enddo
   read(chaps,"(A15I5)") chtp15,NSPIN
   IF(NSPIN==2) FL(4)=.TRUE.
   IF(NSPIN==2) THEN
    WRITE(0,*) 'ERROR!!! ***spin-polarized calculation***'
    WRITE(0,*) 'ERROR!!! Please download the spin-polarized version of vasp2trace at'
    WRITE(0,*) 'ERROR!!! https://github.com/zjwang11/irvsp/blob/master/src_trace_v2.tar.gz!!!'
    STOP ' '
   ENDIF
   read(nfst,"(A90)") chaps
   read(nfst,"(A15L5)") chtp15,FL(2)
!
! set IZ
  do 
     read(nfst,"(A90)") chaps
     chtp5=chaps(7:11)
     if(chtp5 =='direc' ) exit
  enddo
     write(614,"(A40)") '################---input---#############'
     write(614,"(A90)") chaps
  do i=1,3
    !read(nfst,*) br2(:,i),(br4(i,j),j=1,3)
     read(nfst,"(2(3X,3F13.9))") br2(:,i),(br4(i,j),j=1,3)
    !write( 6 ,"(2(3X,3F13.9))") br2(:,i),(br4(i,j),j=1,3)
     write(614,"(2(3X,3F13.9))") br2(:,i),(br4(i,j),j=1,3)
  enddo
  IZ=0
  do irot=1,IORD
     IZ(:,:,irot)= nint(matmul(matmul(br4, DZ2(:,:,irot) ),br2))
  enddo
!
! set NMAT
  do 
   read(nfst,"(A90)") chaps
   chtp15=chaps(1:15)
   if(chtp15 ==' maximum number' ) exit
  enddo
   read(chaps,"(A35I10)") chtp35,NMAT
!
  CLOSE(nfst)
      IF(.NOT.FL(2) .AND. FL(3)) FL(1) = .FALSE.
 !**********************************************************
      IF(     FL(2)) WRITE(614,'(A50)')    &
                   '  1 : Spin-orbit eigenfunctions (->time inversion)'
      IF(.NOT.FL(2)) WRITE(614,'(A34)')    &
                   '  0 : No spin-orbit eigenfunctions'
      IF(     FL(2)) WRITE(624,'(I3)')  1
      IF(.NOT.FL(2)) WRITE(624,'(I3)')  0
      IF(     FL(4)) WRITE(614,'(A23)')    &
                   '  1 : Spin-polarization'
      IF(.NOT.FL(4)) WRITE(614,'(A26)')    &
                   '  0 : No spin-polarization'
!
      RETURN
end subroutine init

subroutine Dmatrix(ih, ik, il, rnx, rny, rnz, rtheta, rphi, degree, twoja1, Dmat)
  implicit none
  integer, intent(in), optional :: ih, ik, il
  real(dp), intent(in), optional :: rnx, rny, rnz
  real(dp), intent(in), optional :: rtheta, rphi
  real(dp), intent(in) ::  degree
  integer, intent(in) :: twoja1  ! equal to 2*J + 1
  complex(dp),dimension(twoja1,twoja1),intent(out) :: Dmat
  ! the key parameter for Dmatrix
  real(dp) :: nx, ny, nz
  real(dp) ::  omega
  ! complex(dp),dimension(aint(J)*2+1,aint(J)*2+1) :: Jx,Jy,Jz
  ! integer :: m
  complex(dp),parameter:: ci = (0.d0,1.d0)
  real*8     ,parameter:: PI = 3.141592653589793238462643383279d0
  real(dp) :: tmp
  real(dp) :: sqrt3
  real(dp) :: cosOmega_2, sinOmega_2
  real(dp) :: cos1omega, sin1omega, cos2omega, sin2omega, cos3omega, sin3omega, cos4omega, sin4omega
  complex(dp) :: exp2omega, exp_2omega, exp4omega, exp_4omega
  sqrt3 = dsqrt(3.0d0)
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
    if(abs(tmp-1.d0).gt.0.5d-6) then
      write(0,'(A,F16.10)')"WARMING!!! The nx, ny, nz input is not normailized: ",tmp
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
  cosOmega_2 = dcos(omega/2.d0)
  sinOmega_2 = dsin(omega/2.d0)
  cos1omega = dcos(omega)
  sin1omega = dsin(omega)
  cos2omega = dcos(2.d0*omega)
  sin2omega = dsin(2.d0*omega)
  cos3omega = dcos(3.d0*omega)
  sin3omega = dsin(3.d0*omega)
  cos4omega = dcos(4.d0*omega)
  sin4omega = dsin(4.d0*omega)
  exp2omega = cmplx(cos2omega,sin2omega)
  exp_2omega = cmplx(cos2omega,-sin2omega)
  exp4omega = cmplx(cos4omega,sin4omega)
  exp_4omega = cmplx(cos4omega,-sin4omega)
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
  else
    write(*,*) " This subroutine can only deal with J = 1/2, 1, 2 !!! "
  end if
end subroutine

end module
