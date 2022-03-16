      SUBROUTINE ROTKV(NMAT,PH,L,ZKV,LKG,IKG,NV,WK,IZ,IIZ,TAU) 
      use symm,only:NSYM,PI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      INTEGER         NMAT
      COMPLEX*16      PH(NSYM,NMAT)
      DIMENSION       IZ(3,3,NSYM),IIZ(3,3,NSYM),TAU(3,NSYM)
      DIMENSION       WK(3)
      DIMENSION       ZKVR(3),LKG(NSYM) &
                     ,ZKV(3,NMAT),L(NSYM,NMAT)
!**********************************************************************
!
      L=0;PH=(0.D0,0.D0) 
!.....find K' such that (k+K)inv(Ri)=k+K'
      DO 10 IG=1,IKG
      DO 10 N=1,NV
        DO 30 J=1,3
 30      ZKVR(J)=ZKV(1,N)*DBLE(IIZ(1,J,LKG(IG))) &
                +ZKV(2,N)*DBLE(IIZ(2,J,LKG(IG))) &
                +ZKV(3,N)*DBLE(IIZ(3,J,LKG(IG)))
!
        IN=1
        DO WHILE(( (ABS(ZKVR(1)-ZKV(1,IN)) &
                   +ABS(ZKVR(2)-ZKV(2,IN)) &
                   +ABS(ZKVR(3)-ZKV(3,IN))).GT.1.D-2).AND.(IN.LE.NV))
!modified by wzj   +ABS(ZKVR(3)-ZKV(3,IN))).GT.1.D-3).AND.(IN.LE.NV))
          IN=IN+1
        END DO
!
        IF(IN.EQ.NV+1)then
           WRITE(6,*) 'rotkv: cannot find (k+K)inv(Ri)'
           WRITE(6,*) IG,N, NV,(ZKV(J,N),J=1,3)
           WRITE(6,*) IG,IN,NV,(ZKVR(J), J=1,3)
           WRITE(6,*) 'All k-vectors:'
           DO 888 II=1,NV
 888          WRITE(6,*)'k+K=',(ZKV(J,II),J=1,3)
           STOP 'rotkv: cannot find (k+K)inv(Ri)'
        endif
!
        L(IG,N)=IN
        AG=0.D0
        DO 40 J=1,3
!       DO 40 J=1,2 
 40       AG=AG-2*PI*DBLE(ZKV(J,IN)-WK(J))*TAU(J,LKG(IG)) 
!40       AG=AG-2*PI*DBLE(ZKV(J,IN))*TAU(J,LKG(IG)) !by zjwang 4.8.2016
          PH(IG,N)=CMPLX(DCOS(AG),DSIN(AG))
 10   CONTINUE
!
      RETURN
      END


