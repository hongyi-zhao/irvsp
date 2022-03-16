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
private
!
integer    , allocatable, save ::  igall(:,:) 
real(dp)   , allocatable, save ::  KV(:,:) 
integer    , allocatable, save ::  L(:,:) 
complex(dp), allocatable, save ::  PH(:,:)

complex(dp), allocatable, save ::  coeffa(:,:),coeffb(:,:) !ncnt,nband
complex(dp), allocatable, save ::  XM(:,:)         !NSYM,nband
real(dp)   , allocatable, save ::  EE(:)           !nband

public ::  kptin
public ::  setarray
public ::  downarray

contains


      subroutine setarray()
       allocate (EE(NUME))         ; EE=0.d0
       allocate (igall(3,NMAT))    ; igall=0
       allocate (KV(3,NMAT))       ; KV=0.d0
       allocate (L(NSYM,NMAT))     ; L=0
       allocate (PH(NSYM,NMAT))    ; PH=0.d0
       allocate (coeffa(NMAT,NUME)); coeffa=0.d0
       allocate (coeffb(NMAT,NUME)); coeffb=0.d0
       allocate (XM(NSYM,NUME))    ; XM=0.d0
      end subroutine

      subroutine downarray()
       if(allocated(EE    ))  deallocate (EE    )
       if(allocated(igall ))  deallocate (igall )
       if(allocated(KV    ))  deallocate (KV    )
       if(allocated(L     ))  deallocate (L     )
       if(allocated(PH    ))  deallocate (PH    )
       if(allocated(coeffa))  deallocate (coeffa)
       if(allocated(coeffb))  deallocate (coeffb)
       if(allocated(XM    ))  deallocate (XM    )
      end subroutine



!======================================================================
! project : coefficients read from WAVECAR
! authors : Zhijun Wang  ( zjwang11@hotmail.com )
! compile : ifort -assume byterecl wave_data.f90  
! status  : good
! comment : These programs are distributed in the hope that they will be 
!           useful, but WITHOUT ANY WARRANTY; without even the implied 
!           warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
!======================================================================
subroutine kptin(KKK)
implicit real*8 (a-h, o-z)
integer,intent(in)::KKK

complex*8 , allocatable :: coeff(:)
complex*16, allocatable :: cener(:)
real*8    , allocatable :: occ  (:)

real*8   :: omeg, n_x,n_y,n_z

dimension  a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),a2xa3(3),sumkg(3),vtmp(3)
dimension  wk(3)

character*75  filename

PARAMETER (PI=3.141592653589793d0)      


!-----read WAVECAR HEAD------------
     filename="WAVECAR"
     data c/0.26246582250210965422d0/ 
     data TOLK/1E-4/ 
     data TOLPH/1E-10/ 


!----- set nrecl ------------------
      nrecl=24

!-END-read WAVECAR HEAD------------
      open(unit=10,file=filename,access='direct',recl=nrecl, &
           iostat=iost,status='old')
      if (iost.ne.0) write(6,*) 'open error - iostat =',iost            
      read(unit=10,rec=1) xnrecl,xispin,xnprec
      nrecl=nint(xnrecl)
      ispin=nint(xispin)
      nprec=nint(xnprec)
      if(nprec.eq.45210) then
         stop  '*** error - WAVECAR_double requires complex*16'
      endif
      if(ispin.eq.2) then
         !write(0,*) '*** NOTE - COMPLETE for FM  case, ISPIN =',ispin
      endif
      close(unit=10)
!-END-read WAVECAR HEAD------------

!-----open FILE--------
      open(unit=10,file=filename,access='direct',recl=nrecl, &
         iostat=iost,status='old')
      if (iost.ne.0) write(6,*) 'open error - iostat =',iost
      read(unit=10,rec=2) xnwk,xnband,ecut,(a1(j),j=1,3),(a2(j),j=1,3), &
           (a3(j),j=1,3)
      nwk=nint(xnwk)
      nband=nint(xnband)
      if ( NSPIN .ne. ispin) then
         write(0,*)   '*** error - selected k=',NSPIN,' > max k=',ispin
         stop
      endif
      if ( NKPTS .ne. nwk) then
         write(0,*)   '*** error - selected k=',NKPTS,' > max k=',nwk
         stop
      endif
      if ( NUME .ne. nband ) then
         write(0,*)   '*** error - selected band=',NUME,' > max band=',nband
         stop
      endif
      


      allocate(cener(nband));cener=0.d0
      allocate(occ  (nband)); occ=0.d0
      !------Find desired wavefunction----------
      irec=3+(KKK-1)*(nband+1)
      read(unit=10,rec=irec) xnplane,(wk(i),i=1,3), &
           (cener(iband),occ(iband),iband=1,nband)
      nplane=nint(xnplane)
      EE=0.d0;EE(1:nband)=real(cener(1:nband))

!-----
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


      
!----------------  FROM ECUT  -----------------
!-----compute reciprocal properties------------
      call vcross(a2xa3,a2,a3)
      Vcell=dot_product(a1,a2xa3)
      a3mag=dsqrt(dot_product(a3,a3))
      call vcross(b1,a2,a3)
      call vcross(b2,a3,a1)
      call vcross(b3,a1,a2)
         b1=2.*pi*b1/Vcell
         b2=2.*pi*b2/Vcell
         b3=2.*pi*b3/Vcell
      b1mag=dsqrt(b1(1)**2+b1(2)**2+b1(3)**2)
      b2mag=dsqrt(b2(1)**2+b2(2)**2+b2(3)**2)
      b3mag=dsqrt(b3(1)**2+b3(2)**2+b3(3)**2)

      
      phi12=acos((b1(1)*b2(1)+b1(2)*b2(2)+b1(3)*b2(3))/(b1mag*b2mag))
      call vcross(vtmp,b1,b2)
      vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
      sinphi123=(b3(1)*vtmp(1)+b3(2)*vtmp(2)+b3(3)*vtmp(3))/(vmag*b3mag)
      nb1maxA=(dsqrt(ecut*c)/(b1mag*abs(sin(phi12))))+1
      nb2maxA=(dsqrt(ecut*c)/(b2mag*abs(sin(phi12))))+1
      nb3maxA=(dsqrt(ecut*c)/(b3mag*abs(sinphi123)))+1
      npmaxA=nint(4.*pi*nb1maxA*nb2maxA*nb3maxA/3.)
            
      phi13=acos((b1(1)*b3(1)+b1(2)*b3(2)+b1(3)*b3(3))/(b1mag*b3mag))
      call vcross(vtmp,b1,b3)
      vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
      sinphi123=(b2(1)*vtmp(1)+b2(2)*vtmp(2)+b2(3)*vtmp(3))/(vmag*b2mag)
      phi123=abs(asin(sinphi123))
      nb1maxB=(dsqrt(ecut*c)/(b1mag*abs(sin(phi13))))+1
      nb2maxB=(dsqrt(ecut*c)/(b2mag*abs(sinphi123)))+1
      nb3maxB=(dsqrt(ecut*c)/(b3mag*abs(sin(phi13))))+1
      npmaxB=nint(4.*pi*nb1maxB*nb2maxB*nb3maxB/3.)
            
      phi23=acos((b2(1)*b3(1)+b2(2)*b3(2)+b2(3)*b3(3))/(b2mag*b3mag))
      call vcross(vtmp,b2,b3)
      vmag=dsqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
      sinphi123=(b1(1)*vtmp(1)+b1(2)*vtmp(2)+b1(3)*vtmp(3))/(vmag*b1mag)
      phi123=abs(asin(sinphi123))
      nb1maxC=(dsqrt(ecut*c)/(b1mag*abs(sinphi123)))+1
      nb2maxC=(dsqrt(ecut*c)/(b2mag*abs(sin(phi23))))+1
      nb3maxC=(dsqrt(ecut*c)/(b3mag*abs(sin(phi23))))+1 
      npmaxC=nint(4.*pi*nb1maxC*nb2maxC*nb3maxC/3.)
      
      nb1max=max0(nb1maxA,nb1maxB,nb1maxC)
      nb2max=max0(nb2maxA,nb2maxB,nb2maxC)
      nb3max=max0(nb3maxA,nb3maxB,nb3maxC)
      npmax=min0(npmaxA,npmaxB,npmaxC)

      npmax=2*(1+2*nb1max)*(1+2*nb2max)*(1+2*nb3max)

      igall=0
      KV=0.d0
 
      igall = 0
!-----Calculate plane waves---------------------
!-----FOR a special K point -------------------
      ncnt=0
      do ig3=0,2*nb3max
         ig3p=ig3
         if (ig3.gt.nb3max) ig3p=ig3-2*nb3max-1
         do ig2=0,2*nb2max
            ig2p=ig2
            if (ig2.gt.nb2max) ig2p=ig2-2*nb2max-1
            do ig1=0,2*nb1max
               ig1p=ig1
               if (ig1.gt.nb1max) ig1p=ig1-2*nb1max-1
               do j=1,3
                  sumkg(j)=(wk(1)+ig1p)*b1(j)+ &
                       (wk(2)+ig2p)*b2(j)+(wk(3)+ig3p)*b3(j)
               enddo
               gtot=sqrt(dot_product(sumkg,sumkg))
               etot=gtot**2/c
               if (etot.lt.ecut) then
                  ncnt=ncnt+1
                  igall(1,ncnt)=ig1p
                  igall(2,ncnt)=ig2p
                  igall(3,ncnt)=ig3p
                  KV(:,ncnt)=DBLE(igall(:,ncnt))+WK(:)
               end if
            enddo
         enddo
      enddo

!--------judge SOC----------
      !! 2* to handle two component spinors
      if (2*ncnt.eq.nplane) then
        !write(6,*) 'Spin-orbit Wavefunctions(INCOMPLETE): just for parity'
      elseif ( ncnt .eq. nplane) then
        !write(6,*) 'No spin-orbit Wavefunctions'
      else
         write(0,*) '*** error - computed 2*ncnt=',2*ncnt, &
              ' != input nplane=',nplane
         stop
      endif
!---------------- ----------
!
      allocate (coeff(2*ncnt )); coeff=0.0
      coeffa=0.d0
      coeffb=0.d0
      L=0
      PH=0.d0
      XM=0.d0
      DO iband=1,nband
         !------Read desired wavefunction----------
         coeff=0.0
         irec=3+( KKK -1)*(nband+1)+iband
         read(unit=10,rec=irec) (coeff(iplane), iplane=1,nplane)
!! wzj
!     WRITE(547,500)
!     WRITE(547,507) KKK,iband
!     WRITE(547,510) (WK(I),I=1,3)
!     do iplane=1,nplane
!     write(547,511) iplane,igall(:,iplane),coeff(iplane),abs(coeff(iplane))!coeff(iplane+ncnt)
!     enddo
!507  FORMAT(/,/,'knum =',I3,4X,'nband= ',I3)
!511  FORMAT(4I5,2X,' (',2F10.6,' )',2X,2F12.6)
!! wzj
         coeffa(1:ncnt,iband)=coeff(1:ncnt)
         coeffb(1:ncnt,iband)=coeff(ncnt+1:2*ncnt)
      ENDDO !!!end for band
!     SUB  KSYM(KKK,WK,NMAT,NV  ,KV,   A  ,   B  ,NUME, NE ,EE,L,PH,XM)
      CALL KSYM(KKK,WK,NMAT,ncnt,KV,coeffa,coeffb,NUME,NUME,EE,L,PH,XM)
!
      deallocate(coeff)
      deallocate(cener)
      deallocate(occ)
      close(10)
!
      return
      end subroutine  kptin

!!$*   routine for computing vector cross-product
subroutine vcross(a,b,c)
  implicit real*8(a-h,o-z)
  dimension a(3),b(3),c(3)
  
  a(1)=b(2)*c(3)-b(3)*c(2)
  a(2)=b(3)*c(1)-b(1)*c(3)
  a(3)=b(1)*c(2)-b(2)*c(1)
  return
end subroutine vcross      

end module wave_data
