      SUBROUTINE CLASSE(IZ,LKG,IKG,ILC,NROT,JIJ,RAN,RNAM,CNAM,GRPNAM,IAXC2)
      use symm,only : NSYM,MAXIR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      LOGICAL          FG1,FG2,FG3,FG4
      CHARACTER*6      RNAM(NSYM),CNAM(NSYM)
      CHARACTER*3      GRPNAM, PDX
      DIMENSION        AXXD(3)
      DIMENSION        LKG(NSYM),IAXC2(3)
      DIMENSION        JIJ(NSYM,NSYM),RAN(4,NSYM),ILC(NSYM)
      DIMENSION        IZ(3,3,NSYM),LINC(NSYM),NROT(10)
      DIMENSION        LCL(NSYM),LFST(NSYM)
!
      CHARACTER*80     CTIR(MAXIR),TTIR
      DIMENSION        NTAB(4)
      COMPLEX*16       ZTIR(MAXIR,MAXIR)
      COMMON /CTAB/    NTAB,CTIR,TTIR,ZTIR
!*******************************************************************
!
      PDX = ' `"'
      DO 2 I=1,3
 2       IAXC2(I) = 0
      DO 4 I=1,IKG
 4       CNAM(LKG(I))(1:6)=RNAM(LKG(I))(1:6)
      DO 6 I=1,NSYM
 6       LFST(I)=0
!
!.....the number of elements in the classes
      ILF=0
      FG3=.FALSE.
      DO 12 I=1,IKG
        IF(ILC(LKG(I)).EQ.-1) FG3=.TRUE.
        INC=0
        FG2=.TRUE.
        DO 11 J1=1,NSYM
 11        LINC(J1)=-1
        DO 13 J1=1,IKG  
          ITT=JIJ(LKG(I),LKG(J1))
!
          FG1=.TRUE.
          K=0
          DO WHILE(K.LT.IKG.AND.FG1)
            K=K+1
            IF(ITT.EQ.LKG(K)) FG1=.FALSE.
          END DO
          IF(ITT.EQ.0) STOP 'classes: incorrect (1)'
          IF(FG1)      STOP 'classes: incorrect (2)'
          IF(LKG(K).LT.LKG(I)) FG2=.FALSE.
!
          FG1=.TRUE. 
          J2=0   
          DO WHILE(J2.LT.INC+1.AND.FG1)
            J2=J2+1
            IF(LKG(K).EQ.LINC(J2)) FG1=.FALSE.
          END DO
          IF(FG1) THEN
            INC=INC+1
            LINC(INC)=LKG(K)
          ENDIF
 13     CONTINUE
!        
!.......first element in each class
        IF(FG2) THEN
           ILF=ILF+1
           LFST(ILF)=LKG(I)
        ENDIF
!
!.......write the number of elements in the class name
        IF(INC.GE.2) THEN
          IF(ILC(LKG(I)).GT.0) THEN
            WRITE(CNAM(LKG(I))(2:2),'(I1)') INC
          ELSE
            WRITE(CNAM(LKG(I))(1:1),'(I1)') INC
          ENDIF
        ENDIF
        LCL(I)=INC
 12   CONTINUE
!
!.....main C2-axis (closest to z-axis)
      DFFZ  = 909.0
      DO I=1,IKG
         IF( CNAM(LKG(I))(2:4).EQ.' C2'.OR.  &
           ( CNAM(LKG(I))(1:4).EQ.' IC2'.AND.GRPNAM.EQ.'C3h').OR. &
           ( CNAM(LKG(I))(1:4).EQ.' IC2'.AND.GRPNAM.EQ.'D3h').OR. &
           ( CNAM(LKG(I))(1:4).EQ.' IC2'.AND.GRPNAM.EQ.'Cs ') ) THEN
           DLENA = SQRT(RAN(1,LKG(I))**2 + RAN(2,LKG(I))**2 + RAN(3,LKG(I))**2)
           TFFZ = ACOS( RAN(3,LKG(I))/DLENA )
           IF( TFFZ.LE.DFFZ ) THEN
             IAXC2(1) = LKG(I)
             DFFZ     = TFFZ
           ENDIF
         ENDIF
      ENDDO
!
!.....axes in D2d,D3h
      IF(GRPNAM.EQ.'D2d'.OR.GRPNAM.EQ.'D3h' ) THEN
      DO I=1,IKG 
        IF(CNAM(LKG(I))(1:4).EQ.' 2C2'.OR.CNAM(LKG(I))(1:4).EQ.' 3C2') IAXC2(2)= LKG(I)
        IF(CNAM(LKG(I))(1:4).EQ.'2IC2'.OR.CNAM(LKG(I))(1:4).EQ.'3IC2') IAXC2(3)= LKG(I)
      ENDDO
      ENDIF
!
!.....axes in D2,C2v,D2h (neg. volume) D4,C4v,D4h (closet to any axes) 
!.....and D6,C6v,D6h (closest to y-axis)
      IF(GRPNAM.EQ.'D2 '.OR.GRPNAM.EQ.'C2v'.OR.GRPNAM.EQ.'D2h'.OR. &
         GRPNAM.EQ.'D4 '.OR.GRPNAM.EQ.'C4v'.OR.GRPNAM.EQ.'D4h'.OR. &
         GRPNAM.EQ.'D6 '.OR.GRPNAM.EQ.'C6v'.OR.GRPNAM.EQ.'D6h' ) THEN
      DFFZ  = 909.0      
      DO 21 I=1,IKG      
      IF(.NOT.FG3.OR.CNAM(LKG(I))(2:2).NE.'I' ) THEN
      IF(LKG(I).NE.IAXC2(1).AND.CNAM(LKG(I))(3:4).EQ.'C2') THEN
        DO 22 J=1,IKG      
        IF(.NOT.FG3.OR.CNAM(LKG(J))(2:2).NE.'I' ) THEN
        IF(LKG(J).NE.IAXC2(1).AND.CNAM(LKG(J))(3:4).EQ.'C2') THEN
!
        FG4=.TRUE.
        DO IJ=1,IKG
          IF(JIJ(LKG(J),LKG(IJ)).EQ.IAXC2(1).OR. &
             JIJ(LKG(J),LKG(IJ)).EQ.LKG(I) ) FG4=.FALSE.
        ENDDO
        IF(FG4) THEN
          DLENA  = SQRT(RAN(1,LKG(I))**2 + RAN(2,LKG(I))**2 + RAN(3,LKG(I))**2)
          DO I1=1,3
            AXXD(I1) = ACOS( ABS(RAN(I1,LKG(I)))/DLENA )
          ENDDO
!
          IF(    GRPNAM.EQ.'D2 '.OR.GRPNAM.EQ.'C2v'.OR.GRPNAM.EQ.'D2h') THEN
            TFFZ = RAN(1,IAXC2(1))*( RAN(2,LKG(I))*RAN(3,LKG(J))    &
                                    -RAN(3,LKG(I))*RAN(2,LKG(J)) )  &
                 + RAN(2,IAXC2(1))*( RAN(3,LKG(I))*RAN(1,LKG(J))    &
                                    -RAN(1,LKG(I))*RAN(3,LKG(J)) )  &
                 + RAN(3,IAXC2(1))*( RAN(1,LKG(I))*RAN(2,LKG(J))    &
                                    -RAN(2,LKG(I))*RAN(1,LKG(J)) )  
          ELSEIF(GRPNAM.EQ.'D4 '.OR.GRPNAM.EQ.'C4v'.OR.GRPNAM.EQ.'D4h') THEN
            TFFZ = AXXD(1)
            IF( AXXD(2).LT.TFFZ) TFFZ = AXXD(2)
            IF( AXXD(3).LT.TFFZ) TFFZ = AXXD(3)
          ELSEIF(GRPNAM.EQ.'D6 '.OR.GRPNAM.EQ.'C6v'.OR.GRPNAM.EQ.'D6h') THEN
            TFFZ = AXXD(2)
          ENDIF
!
          IF(TFFZ.LE.DFFZ) THEN
            IAXC2(2)= LKG(I)
            IAXC2(3)= LKG(J)
            DFFZ    = TFFZ
          ENDIF
!
        ENDIF
        ENDIF
        ENDIF
 22     CONTINUE
      ENDIF
      ENDIF
 21   CONTINUE
      ENDIF   
!
!.....prefix for C2 
      IF(IAXC2(2).NE.0) THEN
        DO 30 I=2,3
        DO 30 J=1,IKG  
           ITT=JIJ(IAXC2(I),LKG(J))
           CNAM(ITT)(5:5)=PDX(I:I)
 30     CONTINUE
      ENDIF
!      do i=1,ikg 
!       write(6,*) i,LKG(i),CNAM(LKG(I))
!      enddo
!
!.....prefix for C2-elements in O,Td,Oh (not having C2-symmetry axies)
      IF(GRPNAM.EQ.'O  '.OR.GRPNAM.EQ.'Td '.OR.GRPNAM.EQ.'Oh ' ) THEN
      DO I=1,IKG      
        IF(CNAM(LKG(I))(2:4).EQ.'6C2'.OR.CNAM(LKG(I))(1:4).EQ.'6IC2') CNAM(LKG(I))(5:5)='`'
      ENDDO
      ENDIF
!
!.....prefix for C3, C4, C6-elemnts (positiv and negativ rot.)
      DO I=1,IKG
         IF(CNAM(LKG(I))(2:4).EQ.' C3'.OR.CNAM(LKG(I))(1:4).EQ.' IC3'.OR. &
            CNAM(LKG(I))(2:4).EQ.' C4'.OR.CNAM(LKG(I))(1:4).EQ.' IC4'.OR. &
            CNAM(LKG(I))(2:4).EQ.' C6'.OR.CNAM(LKG(I))(1:4).EQ.' IC6') THEN          
               CNAM(LKG(I))(5:5)=' '
               IF(RAN(4,LKG(I)).LE.-0.5) CNAM(LKG(I))(5:5)='-'
         ENDIF      
      ENDDO
      IF(GRPNAM.EQ.'T  '.OR.GRPNAM.EQ.'Th ') THEN
      DO I=1,IKG
         ITT=LKG(I)!-----added by wzj
         IF(CNAM(LKG(I))(2:4).EQ.'4C3'.OR.CNAM(LKG(I))(1:4).EQ.'4IC3') THEN          
               CNAM(ITT)(5:5)=' '
               IF(RAN(4,LKG(I)).LE.-0.5) CNAM(ITT)(5:5)='-'
         ENDIF
      ENDDO
      ENDIF
!
!.....set -Cx' = ICx' for crystals with inversion symmetry
      IF(FG3) THEN
      DO 34 I=1,IKG
      IF(CNAM(LKG(I))(2:2).EQ.'I') THEN
      DO 35 J=1,IKG
      IF(CNAM(LKG(J))(2:2).NE.'I') THEN   
          IDFF=0
          DO 36 I1=1,3
          DO 36 I2=1,3
 36          IDFF=IDFF+ABS(IZ(I1,I2,LKG(I))+IZ(I1,I2,LKG(J)))
          IF(IDFF.EQ.0) CNAM(LKG(I))(5:5)=CNAM(LKG(J))(5:5)
      ENDIF
 35   CONTINUE
      ENDIF
 34   CONTINUE
      ENDIF
!
!.....make sure group element have same fourfold axes.
      IF(GRPNAM.EQ.'C4 '.OR.GRPNAM.EQ.'S4 '.OR.GRPNAM.EQ.'C4h'.OR. &
         GRPNAM.EQ.'D4 '.OR.GRPNAM.EQ.'C4v'.OR.GRPNAM.EQ.'D2d'.OR. &
         GRPNAM.EQ.'D4h' ) THEN
      DO 45 J=1,IKG
      IF(CNAM(LKG(J))(3:4).EQ.'C4') THEN
          DFF1= ABS(RAN(1,IAXC2(1))-RAN(1,LKG(J))) &
              + ABS(RAN(2,IAXC2(1))-RAN(2,LKG(J))) &
              + ABS(RAN(3,IAXC2(1))-RAN(3,LKG(J))) 
          IF(DFF1.GE.0.001) THEN
            WRITE(6,*) 'WARNING: C2 and C4 have not same axis'
            STOP       'WARNING: C2 and C4 have not same axis'
          ENDIF
      ENDIF
 45   CONTINUE
      ENDIF
!
!.....make sure group element have same sixfold axes.
      IF(GRPNAM.EQ.'C6 '.OR.GRPNAM.EQ.'C3h'.OR.GRPNAM.EQ.'C6h'.OR. &
         GRPNAM.EQ.'D6 '.OR.GRPNAM.EQ.'C6v'.OR.GRPNAM.EQ.'D3h'.OR. &
         GRPNAM.EQ.'D6h' ) THEN
      DO 65 J=1,IKG
      IF(CNAM(LKG(J))(3:4).EQ.'C3'.OR.CNAM(LKG(J))(3:4).EQ.'C6') THEN
          DFF1= ABS(RAN(1,IAXC2(1))-RAN(1,LKG(J))) &
              + ABS(RAN(2,IAXC2(1))-RAN(2,LKG(J))) &
              + ABS(RAN(3,IAXC2(1))-RAN(3,LKG(J)))
          IF(DFF1.GE.0.001) THEN
            WRITE(6,*) 'WARNING: C2, C3, and C6 have not same axis'
            STOP       'WARNING: C2, C3, and C6 have not same axis'
          ENDIF
      ENDIF
 65   CONTINUE
      ENDIF
!
!.....test that all members in a group has same class name
      DO 80 I=1,IKG
      DO 80 J=1,IKG
        ITT=JIJ(LKG(I),LKG(J))
       !IF(CNAM(ITT)(1:6).NE.CNAM(LKG(I))(1:6) ) THEN
        IF(CNAM(ITT)(1:4).NE.CNAM(LKG(I))(1:4) ) THEN !wzj
            WRITE(6,*) 'ERROR: incorrect classes'
            STOP       'ERROR: incorrect classes'
        ENDIF 
 80   CONTINUE
!
      RETURN
      END
