      SUBROUTINE SU2OP(JIJ,IORD,DZ2,DET,SZ,ANG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NSYM=     96)
!
!       see for instance ref.[2]:
!....................................................................................
!
      COMPLEX*16   SZ(2,2),IMAG,CZERO,SS(2,2),AA,BB
      DIMENSION    DIX(3,3),DZ2(3,3)
      DIMENSION    ANG(3),JIJ(NSYM,NSYM)
      DATA         IMAG/(0.D0,1.D0)/,CZERO/(0.D0,0.D0)/
!**********************************************************************
!
      IS=1
      IF(DET.LT.0.D0) IS=-1 
      DO 10 I=1,3
      DO 10 J=1,3
 10      DIX(I,J)=IS*DZ2(I,J)
!
!   ..determine eulers angles
      CALL EULANG(DIX,ANG)
!
!   ..u(R) for SU(2)
      SIP=DSIN(ANG(1)/2.D0)
      COP=DCOS(ANG(1)/2.D0)
      EPP=(ANG(2)+ANG(3))/2.D0
      ENN=(ANG(2)-ANG(3))/2.D0
!
      SZ(1,1)= COP*CMPLX(DCOS( EPP),DSIN( EPP))
      SZ(1,2)= SIP*CMPLX(DCOS( ENN),DSIN( ENN))
      SZ(2,1)=-SIP*CMPLX(DCOS(-ENN),DSIN(-ENN))
      SZ(2,2)= COP*CMPLX(DCOS(-EPP),DSIN(-EPP))
      RETURN 
!
!-----rotated by zjwang on 11.6.2015
      PI=DACOS(-1.0D0)
      THE=0.D0;PHI=0.D0         !001 spin polarization
      THE=PI/2.D0;PHI=PI/4.D0   !110 spin polarization
      XL=0.D0;XM=0.D0;XN=1.D0   !001 spin polarization
     !XL=1.D0;XM=1.D0;XN=0.D0   !110 spin polarization
     !XL=1.D0;XM=1.D0;XN=1.D0   !111 spin polarization
      XP=DSQRT(XL**2+XM**2)
      THE=DACOS(XN/DSQRT(XP**2+XN**2))
      IF(XP.LT.1D-5) THEN
         PHI=0.D0
      ELSE 
         PHI=ASIN(XM/XP)
      ENDIF

      AA=CMPLX(1.D0+DCOS(THE),DSIN(THE))/2.D0
      BB=CMPLX(DCOS((THE+PHI)/2.D0),-DSIN((THE+PHI)/2.D0))
      SS(1,1)=AA*BB
      SS(1,2)=DSIN(THE/2.D0)*CMPLX(-DCOS(PHI/2.D0), DSIN(PHI/2.D0))
      SS(2,1)=DSIN(THE/2.D0)*CMPLX( DCOS(PHI/2.D0), DSIN(PHI/2.D0))
      SS(2,2)=DCOS(THE/2.D0)*CMPLX( DCOS(PHI/2.D0), DSIN(PHI/2.D0))
      SZ=MATMUL(CONJG(TRANSPOSE(SS)),MATMUL(SZ,SS))
!--- -rotated by zjwang on 11.6.2015
!
      RETURN
      END



