SUBROUTINE CHRCT(NMAT,NUME,A,B,SU2,EE,L,PH,XM,FL,NE,NV,LKG,IKG,IORD,nele)
use symm
implicit none
integer    , intent(IN)   :: NMAT,NUME,L(NSYM,NMAT),NE,NV,LKG(NSYM),IKG,IORD,nele
complex(dp), intent(inout):: A(NMAT,NUME),B(NMAT,NUME)
complex(dp), intent(in)   :: SU2(2,2,NSYM),PH(NSYM,NMAT)
real(dp)   , intent(in)   :: EE(NUME)
logical    , intent(in)   :: FL(FLMAX)
COMPLEX(dp), INTENT(OUT)  :: XM(NSYM,NUME)

complex(dp),allocatable :: a2(:,:),b2(:,:)
integer     :: i,j,I1,ND,IE,IG,i2,j2,no,np,in
real(dp)    :: sum
complex(dp) :: csum, GAM(MAXDG,MAXDG)

complex(dp), parameter :: czero=cmplx(0._dp,0._dp,dp)

!**********************************************************************
!
      if(allocated(a2))  deallocate (a2,b2)
      allocate (a2(nv,maxdg),b2(nv,maxdg))

      DO I=1,IORD
      DO J=1,NE
         XM(I,J)=CMPLX(1.2345678_dp,1.2345678_dp,dp)
      ENDDO
      ENDDO
!
!.....normalization
      DO J=1,NE
        SUM=0._DP
        DO I=1,NV
           SUM=SUM +ABS(A(I,J))**2 + ABS(B(I,J))**2
        ENDDO
        DO I=1,NV
           A(I,J)=A(I,J)/DSQRT(SUM)
           B(I,J)=B(I,J)/DSQRT(SUM)
        ENDDO
      ENDDO
!
!.....characters for all energies and symm. ops.
      IE=1
      DO WHILE(IE.LE.NE)
      ND=1
      DO 
       IF ((IE+ND).GT.NE ) EXIT
       IF ((EE(IE+ND)-EE(IE)).GE.TOLDG) EXIT
        ND=ND+1
      END DO
      IF(ND.GT.MAXDG) THEN
         IE=IE+ND
         CYCLE
      ENDIF
!
 
!https://doi.org/10.1016/j.cpc.2021.107993
!To make vasp2trace output trace data for all bands, two tiny changes should be made to the source code of vasp2trace :
!changing the nele in the 30th line of wrtir.f90 to ne and deleting the 55th line of chrct.f90 , i.e. IF(IE>nele) EXIT .
!      IF(IE>nele) EXIT

!.....for all classes
      DO IG=1,IKG
!
!.....rotation of spin components; include phase from real space rotation 

      
      IF(FL(2)) THEN 
        J2=0
        DO J=IE,IE+ND-1
        J2=J2+1
        DO I=1,NV
          A2(I,J2)=(SU2(1,1,LKG(IG))*A(I,J)+ &
                    SU2(1,2,LKG(IG))*B(I,J))*PH(IG,I)
          B2(I,J2)=(SU2(2,1,LKG(IG))*A(I,J)+ &
                    SU2(2,2,LKG(IG))*B(I,J))*PH(IG,I)
        ENDDO
        ENDDO
      ELSE
        J2=0
        DO J=IE,IE+ND-1
        J2=J2+1
        DO I=1,NV
          A2(I,J2)=A(I,J)*PH(IG,I)
          B2(I,J2)=CZERO
        ENDDO
        ENDDO
      ENDIF
!
!.....character
      XM(LKG(IG),IE)=CZERO
      DO NO=0,ND-1
      DO NP=0,ND-1 
        CSUM=CZERO
        DO I=1,NV
           CSUM=CSUM+DCONJG(A(L(IG,I),IE+NP))*A2(I,1+NO) &
                    +DCONJG(B(L(IG,I),IE+NP))*B2(I,1+NO)
        enddo
        GAM(NP+1,NO+1)=CSUM
        IF(NO.EQ.NP) XM(LKG(IG),IE)=XM(LKG(IG),IE)+CSUM
      enddo
      enddo
!
!.....for all degenerate states
      DO I1=IE+1,IE+ND-1
         XM(LKG(IG),I1)=XM(LKG(IG),IE)
      enddo
!
      ENDDO
      IE=IE+ND
      END DO
!
!.....output (characters for all bands and all symm ops.)
!       DO 300 IE=1,NE
!       WRITE(5,510) KKK,IE,EE(IE)
!       WRITE(5,530) (DREAL(XM(IN,IE)),DIMAG(XM(IN,IE)),IN=1,IORD)
!300  CONTINUE
!
      RETURN
 510  FORMAT(I5,X,I5,F11.7,$)
 530  FORMAT(16(/,3(2F11.7,2X)))
      END
