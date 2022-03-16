      SUBROUTINE RMPROP(IZ,IIZ,TAU,IORD,DZ2,SU2,JIJ,ILC,RAN,RNAM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NSYM=     48)
!
      LOGICAL          FG1,FG2
      CHARACTER*6      RNAM(NSYM)
      DIMENSION        IZ(3,3,NSYM),IIZ(3,3,NSYM),TAU(3,NSYM)
      DIMENSION        RAN(4,NSYM),ANG(3,NSYM)
      COMPLEX*16       SU2(2,2,NSYM),SDET
      COMPLEX*16       SX,SY,SZ(2,2),SZT(2,2)
      DIMENSION        IZ2(3,3),IZ3(3,3)
!     DIMENSION        BR1(3,3),DR1(3,3)
      DIMENSION        JIJ(NSYM,NSYM),IUNI(3,3),ILC(NSYM)
      DIMENSION        DTMP(3,3),DZ2(3,3,NSYM)
      DATA             TOLM/1E-3/
!****************************************************************************
!
      DO 2 I1=1,3
        IUNI(I1,I1)=1
        DO 2 I2=1,3
 2        IF(I1.NE.I2) IUNI(I1,I2)=0
      DO 4 I1=1,NSYM
      DO 3 J1=1,4
 3       RAN(J1,I1)=0
      DO 4 I2=1,NSYM
 4       JIJ(I1,I2)=0
!
!.....calculate JIJ(i,j)=n, where Rn=Rj*Ri*inv(Rj); 
!.....for all Ri and Rj:
      DO 300 I=1,IORD
      DO 200 J=1,IORD
        JIJ(I,J)=0
!
!.......calculate -inv(Rj)*tj
!        DO 250 I1=1,3
! 250       TAU1(I1)=-IZ(1,I1,J)*TAU(1,J)
!     &              -IZ(2,I1,J)*TAU(2,J)
!     &              -IZ(3,I1,J)*TAU(3,J)
!
!.........calculate Pi*inv(Pj)
          DO 202 I1=1,3
          DO 202 I2=1,3
 202         IZ2(I1,I2)=IZ(I1,1,I)*IIZ(1,I2,J) &
                       +IZ(I1,2,I)*IIZ(2,I2,J) &
                       +IZ(I1,3,I)*IIZ(3,I2,J)
!          DO 252 I1=1,3
! 252         TAU2(I1)=IZ(I1,1,I)*TAU1(1)
!     &               +IZ(I1,2,I)*TAU1(2)
!     &               +IZ(I1,3,I)*TAU1(3)
!     &               +TAU(I1,I)
!
!.........calculate Pj*Pi*inv(Pj)
          DO 204 I1=1,3
          DO 204 I2=1,3
 204         IZ3(I1,I2)=IZ(I1,1,J)*IZ2(1,I2) &
                       +IZ(I1,2,J)*IZ2(2,I2) &
                       +IZ(I1,3,J)*IZ2(3,I2)
!
!          DO 254 I1=1,3
! 254         TAU3(I1)=IZ(I1,1,J)*TAU2(1)
!     &               +IZ(I1,2,J)*TAU2(2)
!     &               +IZ(I1,3,J)*TAU2(3)
!     &               +TAU(I1,J)
!
!.........find the smallest non-primitiv translation of tau3
!          DO 260 I1=1,3
!          DO 260 I2=4,1,-1
!            IF(TAU3(I1).GE. I2) TAU3(I1)=TAU3(I1)-I2
!            IF(TAU3(I1).LE.-I2) TAU3(I1)=TAU3(I1)+I2
! 260      CONTINUE
!
!.........find Rn
          FG1=.TRUE.
          K=0
          DO WHILE((FG1).AND.(K.LT.IORD))
            K=K+1
            IDFF=0
            DO 206 I1=1,3
            DO 206 I2=1,3
 206           IDFF=IDFF+ABS(IZ(I1,I2,K)-IZ3(I1,I2))
            IF(IDFF.EQ.0) THEN
              FG1=.FALSE.
              JIJ(I,J)=K
            ENDIF
          END DO
          IF(FG1) STOP 'rmprop: cannot find Rj*Ri*inv(Rj)'
 200  CONTINUE
 300  CONTINUE
!
!
      WRITE(6,*)
      DO 100 I=1,IORD
!
!.....type of rotations, E,C2,C3,C4,C6
      IDET=IZ(1,1,I)*IZ(2,2,I)*IZ(3,3,I) &
          +IZ(1,2,I)*IZ(2,3,I)*IZ(3,1,I)   &
          +IZ(1,3,I)*IZ(2,1,I)*IZ(3,2,I) &
          -IZ(1,1,I)*IZ(2,3,I)*IZ(3,2,I)        &
          -IZ(1,2,I)*IZ(2,1,I)*IZ(3,3,I) &
          -IZ(1,3,I)*IZ(2,2,I)*IZ(3,1,I)
      IF(IDET.NE.1.AND.IDET.NE.-1) STOP 'classes: |det| not 1'
      IDFF=0
      DO 23 I1=1,3
      DO 23 I2=1,3
         IZ3(I1,I2)=IDET*IZ(I1,I2,I)
 23      IDFF=IDFF+ABS(IZ3(I1,I2)-IUNI(I1,I2))
      IP=1
      DO WHILE(IP.LE.6.AND.IDFF.NE.0) 
        IP=IP+1
        IDFF=0
        DO 24 I1=1,3
        DO 24 I2=1,3
           IZ2(I1,I2)=IDET*IZ(I1,1,I)*IZ3(1,I2) &
                     +IDET*IZ(I1,2,I)*IZ3(2,I2) &
                     +IDET*IZ(I1,3,I)*IZ3(3,I2)
 24        IDFF=IDFF+ABS(IZ2(I1,I2)-IUNI(I1,I2))
        DO 25 I1=1,3
        DO 25 I2=1,3
 25       IZ3(I1,I2)=IZ2(I1,I2)
      END DO 
      IF(IDFF.NE.0) STOP 'classes: cannot find type of rot (1)' 
      IF(IP.EQ.5)   STOP 'classes: cannot find type of rot (2)' 
!
      IP=IP*IDET
      IF(IP.EQ.1) WRITE(RNAM(I)(1:6),'(A6)')       '   E  '
      IF(IP.GT.1) WRITE(RNAM(I)(1:6),'(A3,I1,A2)') '  C',IP,'   '
      IF(IP.EQ.-1)WRITE(RNAM(I)(1:6),'(A6)')       '   I  '
      IF(IP.LT.-1)WRITE(RNAM(I)(1:6),'(A3,I1,A2)') ' IC',-IP,'  '
      ILC(I)=IP
!
!.....rotation in the Cartesian coordinate system
!.....IZ_c = (inv(BR1)~)*IZ*(BR1~)
!     DO 320 I1=1,3
!     DO 320 I2=1,3
!320     DTMP(I1,I2)=DFLOAT(IZ(I1,1,I))*BR1(I2,1) &
!                   +DFLOAT(IZ(I1,2,I))*BR1(I2,2) &
!                   +DFLOAT(IZ(I1,3,I))*BR1(I2,3)
!     DO 321 I1=1,3
!     DO 321 I2=1,3
!321     DZ2(I1,I2,I)=DR1(1,I1)*DTMP(1,I2) &
!                    +DR1(2,I1)*DTMP(2,I2) &
!                    +DR1(3,I1)*DTMP(3,I2)
!
      DET= DZ2(1,1,I)*DZ2(2,2,I)*DZ2(3,3,I) &
          +DZ2(1,2,I)*DZ2(2,3,I)*DZ2(3,1,I) &
          +DZ2(1,3,I)*DZ2(2,1,I)*DZ2(3,2,I) &
          -DZ2(1,1,I)*DZ2(2,3,I)*DZ2(3,2,I) &
          -DZ2(1,2,I)*DZ2(2,1,I)*DZ2(3,3,I) &
          -DZ2(1,3,I)*DZ2(2,2,I)*DZ2(3,1,I)
!
      IF(DABS(DET-DBLE(IDET)).GT.TOLM) &
                  STOP 'rmprop: wrong DETERM'
!
!.....rotation axis 
      CALL RAXIS(IP,DET,DZ2(1,1,I),RAN(1,I))
!
!.....determine the operations acting on the spin components
      CALL SU2OP(JIJ,IORD,DZ2(1,1,I),DET,SU2(1,1,I),ANG(1,I)) 
!
 100  CONTINUE
!
!.....use trace(u(Pi)) >= 0
      DO 109 I=1,IORD
      SX=SU2(1,1,I)+SU2(2,2,I)
      IF(DREAL(SX).LT.0.D0) THEN
        DO 90 I1=1,2
        DO 90 J1=1,2
 90        SU2(I1,J1,I)=-SU2(I1,J1,I)
        SX=-SX
      ENDIF
      IF(DREAL(SX).LT.TOLM.AND.DREAL(SU2(1,1,I)).LT.TOLM )THEN
        YY1=ABS(DREAL(SU2(1,1,I)))+ABS(DIMAG(SU2(1,1,I)))
        YY2=ABS(DREAL(SU2(1,2,I)))+YY1
        IF(DREAL(SU2(1,1,I)).LT.-TOLM ) THEN
           FG2=.TRUE.
        ELSEIF(DIMAG(SU2(1,1,I)).LT.-TOLM ) THEN
           FG2=.TRUE.
        ELSEIF(YY1.LT.TOLM.AND.DREAL(SU2(1,2,I)).LT.-TOLM ) THEN
           FG2=.TRUE.
        ELSEIF(YY2.LT.TOLM.AND.DIMAG(SU2(1,2,I)).LT.-TOLM ) THEN
           FG2=.TRUE.   
        ELSE
           FG2=.FALSE.   
        ENDIF 
        IF(FG2) THEN
          DO 91 I1=1,2
          DO 91 J1=1,2
 91          SU2(I1,J1,I)=-SU2(I1,J1,I)
          SX=-SX
        ENDIF
      ENDIF
109  CONTINUE
!
!.....use only unbarred operations [u(Pi)=-u(Pi_bar)]
      DO 110 I=1,IORD
      SX=SU2(1,1,I)+SU2(2,2,I)
      SY=SU2(1,2,I)+SU2(2,1,I)
      SDET=SU2(1,1,I)*SU2(2,2,I)-SU2(2,1,I)*SU2(1,2,I)
      IF(DABS(DIMAG(SX))  .GT.TOLM) STOP 'imagina trace(u)'
      IF(DABS(DIMAG(SDET)).GT.TOLM) STOP 'imagina det(u)  '
      DO 120 J=1,IORD         
         K=JIJ(I,J)
         IF(K.LT.I) THEN
!..........inv(Uj)
           SDET =SU2(1,1,J)*SU2(2,2,J) &
                -SU2(2,1,J)*SU2(1,2,J) 
           SZ(1,1)= SU2(2,2,J)/SDET
           SZ(1,2)=-SU2(1,2,J)/SDET
           SZ(2,1)=-SU2(2,1,J)/SDET
           SZ(2,2)= SU2(1,1,J)/SDET
!las           do 888 i1=1,2
!           do 888 i2=1,2
! 888          ss(i1,i2)=sz(i1,i2)
!..........Ui*inv(Uj)
           SZT(1,1)=SU2(1,1,I)*SZ(1,1)+SU2(1,2,I)*SZ(2,1)  
           SZT(1,2)=SU2(1,1,I)*SZ(1,2)+SU2(1,2,I)*SZ(2,2)  
           SZT(2,1)=SU2(2,1,I)*SZ(1,1)+SU2(2,2,I)*SZ(2,1)  
           SZT(2,2)=SU2(2,1,I)*SZ(1,2)+SU2(2,2,I)*SZ(2,2)  
!..........Uj*Ui*inv(Uj)
           SZ(1,1)=SU2(1,1,J)*SZT(1,1)+SU2(1,2,J)*SZT(2,1)  
           SZ(1,2)=SU2(1,1,J)*SZT(1,2)+SU2(1,2,J)*SZT(2,2)  
           SZ(2,1)=SU2(2,1,J)*SZT(1,1)+SU2(2,2,J)*SZT(2,1)  
           SZ(2,2)=SU2(2,1,J)*SZT(1,2)+SU2(2,2,J)*SZT(2,2)  
!
           DFP=0D0
           DFM=0D0
           DO 121 I1=1,2
           DO 121 I2=1,2
              DFP=DFP+ABS(SZ(I1,I2)+SU2(I1,I2,K))
              DFM=DFM+ABS(SZ(I1,I2)-SU2(I1,I2,K))
 121       CONTINUE
           IF(ABS(DFP).GT.TOLM.AND.ABS(DFM).GT.TOLM)  &
               STOP 'rmprop: not correct un=uj*ui*inv(uj)'
           IF(ABS(DFP).LT.TOLM) THEN
              DO 122 I1=1,2
              DO 122 I2=1,2
 122             SU2(I1,I2,I)=-SU2(I1,I2,I)
           ENDIF
!..........check determinates
           SDETK=SU2(1,1,K)*SU2(2,2,K) &
                -SU2(2,1,K)*SU2(1,2,K)
           SDET =SU2(1,1,I)*SU2(2,2,I) &
                -SU2(2,1,I)*SU2(1,2,I)
           IF(ABS(SDETK-SDET).GT.TOLM)  &
                STOP 'rmprop: |u(i)| not |u(k)| despite same class'
         ENDIF 
 120  CONTINUE
 110  CONTINUE
!
!.....output
       WRITE(6,500)
       DO 150 I=1,IORD
       WRITE(6,516)(IZ( 1,I2,I),I2=1,3),TAU(1,I), &
                   (IIZ(1,I2,I),I2=1,3), &
                   (DZ2(1,I2,I),I2=1,3),ANG(1,I)
       DO 410 I1=2,3                                                 
  410     WRITE(6,515)(IZ(I1,I2,I),I2=1,3),TAU(I1,I), &
            (IIZ(I1,I2,I),I2=1,3), &
            (DZ2(I1,I2,I),I2=1,3),ANG(I1,I),(SU2(I1-1,I2,I),I2=1,2)
       WRITE(6,517) I
       IF(ILC(I).GE.1) THEN
         IF(ILC(I).EQ. 1) WRITE(6,518) '(E)','unity op.   '
         IF(ILC(I).EQ. 2) WRITE(6,519)  ILC(I),180
         IF(ILC(I).EQ. 3) WRITE(6,519)  ILC(I),120
         IF(ILC(I).EQ. 4) WRITE(6,521)  ILC(I), 90
         IF(ILC(I).EQ. 6) WRITE(6,521)  ILC(I), 60
       ELSE
         IF(ILC(I).EQ.-1) WRITE(6,518) '(I)','inverse op. '
         IF(ILC(I).EQ.-2) WRITE(6,525) -ILC(I),180
         IF(ILC(I).EQ.-3) WRITE(6,525) -ILC(I),120
         IF(ILC(I).EQ.-4) WRITE(6,527) -ILC(I), 90
         IF(ILC(I).EQ.-6) WRITE(6,527) -ILC(I), 60
       ENDIF
       IF(RAN(4,I).LT.-0.5) WRITE(6,530) 'clockwise '
       IF(RAN(4,I).GT.0.5)  WRITE(6,531) 'counterclockwise '
       WRITE(6,532) (RAN(I1,I),I1=1,3)
       WRITE(6,*)
  150  CONTINUE
!
!      WRITE(6,570) &
!        'Note:                                                                         ',&
!        '(1) Labelling of IRs depends on atom positions.                               ',&
!        '(2) Labelling of IRs in the character table of Koster et al. as well as of    ',&
!        'Altmann et al. (see below) is defined by the chosen symmetry axes and for     ',&
!        'each point group. Hence, in order to compare the labelling of the present     ',&
!        'IRs with the literature one has to determine the relation, for each k-point,  ',&
!        'between their symmetry axes and the present ones.                             ',&
!        'That is, please check the group elements of the different classes.            '
!
!
      RETURN
 500  FORMAT(/,'SYMMETRY OPERATIONS Pi={Ri|taui+tm}',/, &
            '  Ri     taui ', &
            '  inv(Ri) Ri(Cartesian coord)', &
            '  Eulers angles  spin transf.'/) 
 515  FORMAT(3I2,F8.3,2X,3I2,2X,3F6.3,2X,F6.3,2X,2('(',2F6.3,')'))
 516  FORMAT(3I2,F8.3,2X,3I2,2X,3F6.3,2X,F6.3)
 517  FORMAT(I4,$)
 518  FORMAT(1X,A3,/,A12) 
 519  FORMAT(1X,'(C',I1,')', /,I3,'-degree rotation') 
 521  FORMAT(1X,'(C',I1,')', /,I2,'-degree rotation')
 525  FORMAT(1X,'(IC',I1,')',/,I3,'-degree rotation times inversion')
 527  FORMAT(1X,'(IC',I1,')',/,I2,'-degree rotation times inversion')

 530  FORMAT(A10,$)
 531  FORMAT(A17,$)
 532  FORMAT('rotation through (',F6.3,2(',',F6.3),')')

 555  FORMAT(/)
! 570  FORMAT(/,8(A78,/))
!
      END  
