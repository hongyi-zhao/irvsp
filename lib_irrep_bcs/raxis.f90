      SUBROUTINE RAXIS(IP,DET,R,RAN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      DIMENSION        R(3,3),RAN(4),DANG(5)
      DATA             DANG/0.33333333,0.57735027,0.66666667, &
                            0.70710678,0.8660254/
      DATA             TOL/5E-5/
      DATA             PI/3.141592653589793238462643383279/
!****************************************************************************
!
      DO 2 I1=1,4
 2       RAN(I1)=0.D0
!
      do 15 ix=-2,2
      do 15 iy=-2,2
      x=ix*0.5d0
      y=iy*0.5d0
!      DO 15 X=-1.D0,1.D0,0.5
!      DO 15 Y=-1.D0,1.D0,0.5
         IF(ABS(X*X+Y*Y).LE.1.D0) THEN
         Z  = DSQRT(1.0-X*X-Y*Y)
         RX = DET*R(1,1)*X + DET*R(1,2)*Y + DET*R(1,3)*Z
         RY = DET*R(2,1)*X + DET*R(2,2)*Y + DET*R(2,3)*Z
         RZ = DET*R(3,1)*X + DET*R(3,2)*Y + DET*R(3,3)*Z
         AN = DACOS(RX*X + RY*Y + RZ*Z)
         IF(DABS(AN).LT.TOL) THEN
             RAN(1) = X
             RAN(2) = Y
             RAN(3) = Z
             GOTO 399
         ENDIF
         ENDIF
 15   CONTINUE
!
      DO 25 IX=-5,5
      DO 25 IY=-5,5
        IF(IY.NE.0.AND.IX.NE.0) THEN
          X  = (DBLE(IX)/ABS(IX))*DANG(ABS(IX))
          Y  = (DBLE(IY)/ABS(IY))*DANG(ABS(IY))
          IF(ABS(X*X+Y*Y).LE.1) THEN
          Z  = DSQRT(1.0-X*X-Y*Y)
          RX = DET*R(1,1)*X + DET*R(1,2)*Y + DET*R(1,3)*Z
          RY = DET*R(2,1)*X + DET*R(2,2)*Y + DET*R(2,3)*Z
          RZ = DET*R(3,1)*X + DET*R(3,2)*Y + DET*R(3,3)*Z
          AN = DACOS(RX*X + RY*Y + RZ*Z)
          IF(DABS(AN).LT.TOL) THEN
             RAN(1) = X
             RAN(2) = Y
             RAN(3) = Z
             GOTO 399
          ENDIF
          ENDIF
        ENDIF
 25   CONTINUE
!
      DMIN = PI
      DO 35 I1 = 0,200
        TE = (PI/2.D0)*DBLE(I1)/200.D0
        DO 35 I2 = 0,800
          PH = 2.D0*PI*DBLE(I2)/800.D0
          X  = DSIN(TE)*DCOS(PH)
          Y  = DSIN(TE)*DSIN(PH)
          Z  = DCOS(TE)
          RX = DET*R(1,1)*X + DET*R(1,2)*Y + DET*R(1,3)*Z
          RY = DET*R(2,1)*X + DET*R(2,2)*Y + DET*R(2,3)*Z
          RZ = DET*R(3,1)*X + DET*R(3,2)*Y + DET*R(3,3)*Z
          DELT=1E-8
          IF((RX*X + RY*Y + RZ*Z).LT.0.D0) DELT=-DELT
          AN = DACOS(RX*X + RY*Y + RZ*Z - DELT)
          IF(DABS(AN).LT.DMIN) THEN
             RAN(1) = X
             RAN(2) = Y
             RAN(3) = Z
             DMIN   = AN
          ENDIF
 35   CONTINUE
!
 399  CONTINUE
!
!.....positive angles
      IF(RAN(3).LT.TOL) THEN
        IF(RAN(3).LE.-TOL) STOP 'raxis: Z < 0'
        RAN(3)=0.D0
        IF(RAN(1).LE.-TOL) THEN
          DO 401 I1=1,2
 401         RAN(I1) = -RAN(I1)
        ENDIF
        IF(DABS(RAN(1)).LT.TOL.AND.RAN(2).LE.-TOL) THEN
          DO 402 I1=1,2
 402         RAN(I1) = -RAN(I1)
        ENDIF        
      ENDIF
!
!.....direction of the rotation
      IF(ABS(IP).GE.3) THEN
        IF(DABS(RAN(3)).GE.TOL) THEN
          X= 1.D0
          Y= 0.D0
          Z= 0.D0
        ELSE
          X= 0.D0
          Y= 0.D0
          Z= 1.D0
        ENDIF
        RX = DET*R(1,1)*X + DET*R(1,2)*Y + DET*R(1,3)*Z
        RY = DET*R(2,1)*X + DET*R(2,2)*Y + DET*R(2,3)*Z
        RZ = DET*R(3,1)*X + DET*R(3,2)*Y + DET*R(3,3)*Z
!    
        AN = RAN(1)*(Y*RZ-Z*RY) &
           + RAN(2)*(Z*RX-X*RZ) &
           + RAN(3)*(X*RY-Y*RX)
!
        RAN(4)=1.D0
        IF(DABS(AN).LT.TOL) STOP 'raxis: volume=0'
        IF(AN.GT.0) RAN(4)=-1.D0
      ENDIF
!
      RETURN
      END


