!=========================================================================!
! project : wave_data
! history : 07/18/2014
! authors : Zhijun Wang  ( zjwang11@hotmail.com )
! purpose : Get input from WAVECAR
! status  : robust
! comment : These programs are distributed in the hope that they will be 
!           useful, but WITHOUT ANY WARRANTY; without even the implied 
!           warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
!=========================================================================!
module wave_data
use symm,only : NSYM,ksym
use struct_data
implicit none
private
!

real (dp)   , parameter ::  pi    = 3.141592653589793238462643383279_dp 
complex(dp) , parameter ::  ci    = dcmplx(0.0_dp, 1.0_dp)   
complex(dp) , parameter ::  czero = dcmplx(0.0_dp, 0.0_dp)
!
integer     ,  save     ::  NMAT  ! - max. no. of orbitals
integer     ,  save     ::  NUME  ! - max. no. of eigenvalues(bands)
integer     ,  save     ::  NRPTS

complex(dp) , allocatable , save ::  LMAT(:,:,:) 
                          
complex(dp) , allocatable , save ::  coeffa(:,:),coeffb(:,:) !ncnt,nband
complex(dp) , allocatable , save ::  XM(:,:)         !NSYM,nband
real(dp)    , allocatable , save ::  EE(:)           !nband
                          
integer     , allocatable , save ::  irvec(:,:)      !nband
complex(dp) , allocatable , save ::  HmnR(:,:,:)     !NSYM,nband

public   ::   KPTIN
public   ::   setarray
public   ::   downarray

contains

      subroutine setarray()
       call readhr()  !set NMAT,NUME
       allocate (EE(NUME))         ; EE=0.d0
       allocate (LMAT(NMAT,NMAT,NSYM)); LMAT=0.d0
       allocate (coeffa(NMAT,NUME)); coeffa=czero
       allocate (coeffb(NMAT,NUME)); coeffb=czero
       allocate (XM(NSYM,NUME))    ; XM=czero

      end subroutine

      subroutine downarray()
       if(allocated(EE    ))  deallocate (EE    )
       if(allocated(LMAT  ))  deallocate (LMAT  )
       if(allocated(coeffa))  deallocate (coeffa)
       if(allocated(coeffb))  deallocate (coeffb)
       if(allocated(XM    ))  deallocate (XM    )

       if(allocated(irvec ))  deallocate (irvec )
       if(allocated(HmnR  ))  deallocate (HmnR  )

       if(allocated(n2obtau))  deallocate (n2obtau)
       if(allocated(n2ibtau))  deallocate (n2ibtau)
       if(allocated(rot2tau))  deallocate (rot2tau)
       if(allocated(rot2vec))  deallocate (rot2vec)
       if(allocated(rot2phi))  deallocate (rot2phi)
       if(allocated(rotjorb))  deallocate (rotjorb)
       if(allocated(len_k  ))  deallocate (len_k  )
       if(allocated(k      ))  deallocate (k      )

      end subroutine



!======================================================================
subroutine KPTIN(KKK)
integer,intent(in)::KKK

integer :: i,info
real(dp)  :: wk(3)
complex(dp)  :: Hmnk(NUME,NUME), evec(NUME,NUME)

integer    :: L(NSYM,NMAT)

     info=0
!-----
      wk(:)=k(:,KKK)
      do i=1,3
        IF(abs(3.d0*wk(i)-1.d0) .lt. 0.0002d0)  wk(i)=1.d0/3.d0
      enddo

!.....output
      WRITE(6,500)
      WRITE(6,509) KKK
      WRITE(6,510) (WK(I),I=1,3)
 500  FORMAT(/,80('*'))
 509  FORMAT(/,/,'knum =',I3,4X,'kname= ',A10)
 510  FORMAT('k =',3F9.6,/)


!--------------------------
!
     Hmnk=czero
     call get_Hk(NUME,wk,Hmnk)
     evec=czero 
     call myzheev(NUME,NUME,Hmnk,EE,evec,info)
     write(102,'(1000F12.6)') len_k(KKK),EE(:)
     coeffa=czero; coeffb=czero;
     coeffa(1:NMAT,1:NUME)=evec(1:NMAT,1:NUME)
     IF(NSPIN==2) coeffb(1:NMAT,1:NUME)=evec(NMAT+1:2*NMAT,1:NUME)


     call setL(wk)

     XM=czero
!     SUB  KSYM(KKK,WK,NMAT,NV  ,KV,   A  ,   B  ,NUME, NE ,EE,L,PH,XM)
      !CALL KSYM(KKK,WK,NMAT,ncnt,KV,coeffa,coeffb,NUME,NUME,EE,L,PH,XM)
      CALL KSYM(KKK,WK,NMAT,NMAT,coeffa,coeffb,NUME,NUME,EE,LMAT,XM)
!
!
      return
      end subroutine  KPTIN

!!$*   routine for computing vector cross-product
subroutine vcross(a,b,c)
  implicit real*8(a-h,o-z)
  dimension a(3),b(3),c(3)
  
  a(1)=b(2)*c(3)-b(3)*c(2)
  a(2)=b(3)*c(1)-b(1)*c(3)
  a(3)=b(1)*c(2)-b(2)*c(1)
  return
end subroutine vcross      

!=========+=========+=========+=========+=========+=========+=========+
!  get the data
!=========+=========+=========+=========+=========+=========+=========+
subroutine readhr ( )
  implicit none

  integer :: nbs, NPS
  integer :: ftemp=1000,i,j,iR,Itmp,Jtmp
  integer ,allocatable :: NDEGEN(:)
  real(dp):: r1tmp,r2tmp

!!!--------------------------------------------------------------------
  !write(*,"(A20,2X,A15)")  'read hr-data from', hrfile

  open(ftemp,file=hrfile,status='old')
  read(ftemp,*) 
  read(ftemp,*)  nbs
  read(ftemp,*)  NPS

  NMAT=nbs/NSPIN
  NUME=nbs
  NRPTS=NPS

  allocate (NDEGEN(NRPTS))   ;NDEGEN=0

  allocate (irvec(3,NRPTS))   ; irvec=0
  allocate (HmnR(NUME  ,NUME  ,NRPTS));HmnR=czero

  read(ftemp,"(15I5)") NDEGEN(:)
  do iR=1,NPS
     do j=1,Nbs
     do i=1,Nbs
     !read(ftemp,"(5I5,2F12.6)") irvec(:,iR),Itmp,Jtmp,HmnR(i,j,iR)
      read(ftemp,*) irvec(:,iR),Itmp,Jtmp,r1tmp,r2tmp
      HmnR(i,j,iR)=cmplx(r1tmp,r2tmp,dp)
     enddo
     enddo
     IF(NDEGEN(iR)==1) CYCLE
     HmnR(:,:,IR)=HmnR(:,:,IR)/NDEGEN(iR)
     NDEGEN(iR)=1
  enddo
  close(ftemp)

  deallocate(NDEGEN)
!!!--------------------------------------------------------------------
  return
end subroutine readhr
!=========+=========+=========+=========+=========+=========+=========+


!=========+=========+=========+=========+=========+=========+=========+
!  get the data
!=========+=========+=========+=========+=========+=========+=========+
subroutine get_Hk(nbands,kk,Hk)
implicit none
integer,intent(in) :: nbands
real(dp), intent(in) :: kk(3)
complex(dp), intent(out) :: Hk(nbands,nbands)

integer :: ir
real(dp):: kdotr

Hk=czero
do ir = 1,Nrpts
   kdotr=dot_product(kk(:),irvec(:,ir))
   Hk(1:nbands,1:nbands) =   Hk(1:nbands,1:nbands)   &
                         + HmnR(1:nbands,1:nbands,ir)*exp(CI*kdotr*2.d0*PI)
enddo

return
end subroutine get_Hk

!=========+=========+=========+=========+=========+=========+=========+

  subroutine myzheev(lda,N,Hmn,eval,evec,info)
     implicit none
     integer    , intent(in) :: lda
     integer    , intent(in) :: N
     complex(dp), intent(in) :: Hmn(lda,N)
     real(dp)   , intent(out):: eval(N)
     complex(dp), intent(out):: evec(lda,N)
     integer    , intent(out):: info

   ! dummy integer variables for lapack call
   INTEGER    :: LWORK
   REAL(DP)   :: RWORK(3*N) 
   COMPLEX(DP):: WORK(2*N)
   complex(dp) :: A(N,N)

   LWORK=2*N

   A=czero;A(:,:)=Hmn(1:N,1:N)
   CALL ZHEEV('V', 'U', N, A, N, eval, WORK, LWORK, RWORK, INFO )
   evec(1:N,1:N)=A(:,:)

     return
  end subroutine myzheev


    subroutine setL(wk)
    implicit none
    real(dp),intent(in) :: wk(3)
    integer:: itau,irot,jtau,jm,in,bs,dt
    real(dp):: phk

    LMAT=czero ! allocate (LMAT(NMAT,NMAT,NSYM)); L=0
    do irot=1,nrot
    do itau =1,ntau
       jtau=rot2tau(itau,irot)
       jm=n2obtau(jtau);bs=n2ibtau(jtau)
       in=n2obtau(itau);dt=n2ibtau(itau)
       IF(jm/=in) STOP "error symm in subroutine setL"
      !phk= dot_product(wk(:),rot2vec(:,itau,irot))
       phk= dot_product(wk(:),rot2phi(:,itau,irot))
       LMAT((/1:jm/)+bs,(/1:in/)+dt,irot)=rotjorb((/1:jm/),(/1:in/),itau,irot)*EXP(-CI*phk*2.d0*pi)
    enddo
    !write(*,'(8F12.6)') real(LMAT(:,:,irot))
    !write(*,*)
    enddo

    end subroutine setL

end module wave_data
