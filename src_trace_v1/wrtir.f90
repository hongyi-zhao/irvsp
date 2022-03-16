SUBROUTINE WRTIR(NUME,EE,XM,NE,LKG,IKG,IK,nele)
      use symm
      INTEGER ,INTENT(IN) :: NUME,IKG,IK,nele
      INTEGER ,INTENT(IN) :: NE,LKG(NSYM)
      REAL(DP),INTENT(IN) :: EE(NUME)
      COMPLEX(DP),INTENT(IN)::XM(NSYM,NUME)

      real(DP),    parameter  :: epsil= 0.1D-6
      LOGICAL  :: FG1
      REAL(DP) :: DSUM
      COMPLEX(DP) :: X,ZIR
      INTEGER :: IE,ND,IG,I1,J1,I2,J2,IS
      INTEGER :: NIR,NIR1,NIR2
     
!**********************************************************************
!
!.....output
      IE=1
      DO WHILE(IE.LE.NE)
      ND=1
      DO WHILE((IE+ND).LE.NE)
      IF((EE(IE+ND)-EE(IE)).LT.TOLDG) THEN
        ND=ND+1
      ELSE
        EXIT 
      ENDIF
      END DO
      IF(ND.GT.MAXDG) GOTO 100

     !https://doi.org/10.1016/j.cpc.2021.107993
     !To make vasp2trace output trace data for all bands, two tiny changes should be made to the source code of vasp2trace :
! changing the nele in the 30th line of wrtir.f90 to ne and deleting the 55th line of chrct.f90 , i.e. IF(IE>nele) EXIT .
!     if(IE<=nele) then
     if(IE<=ne) then
      WRITE(614,'(I3,I3,F12.6,$)') IE,ND,EE(IE)
      WRITE(624,'(I3,I3,F12.6,$)') IE,ND,EE(IE)
      DO IG=1,IKG
         X=XM(LKG(IG),IE)
         WRITE(614,590) DREAL(X)+epsil,DIMAG(X)+epsil
         WRITE(624,'(2F12.6,$)') DREAL(X)+epsil,DIMAG(X)+epsil
      ENDDO
         WRITE(624,*)
         WRITE(614,*)
     endif
!
      DSUM=0.D0
!
 100  IE=IE+ND
      END DO

      RETURN
 590  FORMAT(2X,2F12.6,$)
END SUBROUTINE WRTIR
