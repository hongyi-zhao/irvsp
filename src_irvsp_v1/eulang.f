      SUBROUTINE EULANG(RM,A,IS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      DIMENSION A(3),RM(3,3),RM2(3,3)
      DATA TOL/5.0E-5/,TOLP/5.0E-6/
!**********************************************************************
!
      DO 10 I=1,3
      DO 10 J=1,3
 10      RM2(I,J)=RM(I,J)
!
      PI=DACOS(-1.D0)
      A(1)=DACOS((RM(3,3)))
      DO WHILE(A(1).GT.PI) 
        A(1)=A(1)-PI
      END DO
      DO WHILE(A(1).LT.0.D0) 
        A(1)=A(1)+PI
      END DO
!
      IF(DABS(DABS(RM(3,3))-1.D0).LT.TOL ) THEN
        A(2)=DACOS(RM(3,3)*RM(1,1))
        IF(RM(1,2).LT.0.D0) A(2)=2*PI-A(2)
        A(3)=0.D0
      ELSE
        A(2)=DACOS(-RM(1,3)/DSIN(A(1)))
        IF((RM(2,3)/DSIN(A(1))).LT.0.D0) A(2)=2*PI-A(2)
        A(3)=DACOS(RM(3,1)/DSIN(A(1)))
        IF((RM(3,2)/DSIN(A(1))).LT.0.D0) A(3)=2*PI-A(3)
      ENDIF
!
!.....0<= a(1) <=pi,   0<= a(2), a(3) <=2*pi


      DO 30 I1=1,3
        DO 32 I2=-4,4
 32        IF(ABS(A(I1)-I2*PI).LT.1E-6) A(I1)=I2*PI

        AMAX=PI
        IF(I1.GE.2) AMAX=2*PI
        DO WHILE(A(I1).LT.0.D0)
             A(I1)=A(I1)+AMAX
        END DO
        DO WHILE(A(I1).GT.AMAX) 
            A(I1)=A(I1)-AMAX
        END DO
        IF(DABS(A(I1)).LT.TOLP)      A(I1)=1E-12
        IF(DABS(A(I1)-AMAX).LT.TOLP) A(I1)=AMAX
 30   CONTINUE
!
!   ..check
      RM(1,1)=-DSIN(A(3))*DSIN(A(2))+DCOS(A(1))*DCOS(A(3))*DCOS(A(2))    
      RM(1,2)= DCOS(A(3))*DSIN(A(2))+DCOS(A(1))*DSIN(A(3))*DCOS(A(2)) 
      RM(1,3)=-DSIN(A(1))*DCOS(A(2))
      RM(2,1)=-DSIN(A(3))*DCOS(A(2))-DCOS(A(1))*DCOS(A(3))*DSIN(A(2))
      RM(2,2)= DCOS(A(3))*DCOS(A(2))-DCOS(A(1))*DSIN(A(3))*DSIN(A(2))
      RM(2,3)= DSIN(A(1))*DSIN(A(2))
      RM(3,1)= DSIN(A(1))*DCOS(A(3))
      RM(3,2)= DSIN(A(1))*DSIN(A(3))
      RM(3,3)= DCOS(A(1))
!
      SUMP=0.D0
      DO 20 I=1,3
      DO 20 J=1,3
       SUMP=SUMP+DABS(RM(I,J)-RM2(I,J))
 20   CONTINUE
      IF(SUMP.GT.9*TOL) THEN
         WRITE(6,*) (A(I1),I1=1,3)
         WRITE(6,*) (RM(1,I1),I1=1,3)
         WRITE(6,*) (RM(2,I1),I1=1,3)
         WRITE(6,*) (RM(3,I1),I1=1,3)
         STOP 'Error in eulang: Wrong angles'
      ENDIF
! 
      RETURN
      END


