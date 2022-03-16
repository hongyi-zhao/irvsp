      SUBROUTINE WRTIR(NUME,FL,FGT,EE,XM,NE,LKG,IKG,CNAM)
      use symm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      LOGICAL          FL(FLMAX),FG1,FGT 
      CHARACTER*6      CNAM(NSYM)
      CHARACTER*3      GRPNAM 
      COMPLEX*16       X,XOLD,XR(MAXIR),ZIR,XM(NSYM,NUME)
      COMPLEX*16       CZERO
      DIMENSION        LKG(NSYM),EE(NUME)
      DIMENSION        NJ1I(MAXIRDG),SK(3)
!
      CHARACTER*80     CTIR(MAXIR),TTIR
      DIMENSION        NTAB(4)
      COMPLEX*16       ZTIR(MAXIR,MAXIR)
      COMMON /CTAB/    NTAB,CTIR,TTIR,ZTIR
!
      DATA             CZERO/(0.0D0,0.0D0)/
      DATA             TOL/5.0E-2/
!**********************************************************************
!
!.....calc. the characters for all energies and tranformations
     !IE=1
      IE=nmin
      DO WHILE(IE.LE.NE)
      ND=1
      DO WHILE((IE+ND).LE.NE.AND.(EE(IE+ND)-EE(IE)).LT.TOLDG)
        ND=ND+1
      END DO
     !IF(FGT) WRITE(6,570) IE,ND,EE(IE)
              WRITE(6,570) IE,ND,EE(IE)
      IF(ND.GT.MAXDG) GOTO 100
!.....if 'Cornwell condition' is satisfied
      IF(FGT) THEN
!
!.....for all classes
      IF(.NOT.FL(2)) THEN
        NIR1=1
        NIR2=NTAB(4)
      ELSE
        NIR1=NTAB(4)+1
        NIR2=NTAB(3)        
      ENDIF
      NIR=NIR2-NIR1+1
!
      DO 90 IR=1,NTAB(4)
      IWRT=0
      DO 90 IG=1,IKG
      IF(CNAM(LKG(IG)).EQ.TTIR(6*(IR-1)+9:6*IR+8)) THEN
      IWRT=IWRT+1    
      X=XM(LKG(IG),IE)
      
!
!.....test: all elements in a class shall give the same character
      IF(IWRT.EQ.1) THEN
         XR(IR)=X
         WRITE(6,580) DREAL(X),DIMAG(X) 
      ELSE
        IF(ABS(X-XOLD).GT.TOL) THEN
          WRITE(6,581) X,XOLD
          WRITE(6,*) LIGO,LKG(IG),'   ',CNAM(LKG(IG))
          WRITE(6,*) 'X not equal for all elements in the class '
!Clas          STOP 'wrtir: X not equal for all elements in the class '
        ELSE
          XOLD=X
        ENDIF
      ENDIF
!
      XOLD=X
      LIGO=LKG(IG)
      ENDIF
  90  CONTINUE
!
!.....identify the irreducible representation
      FG1=.TRUE.
      DO 400 I1=1,MAXIRDG
      DO 401 J1=1,MAXIRDG
 401     NJ1I(J1)=0
!
         DO 405 I2=0,NIR**I1-1
         NJ1I(1)=MOD(I2,NIR)+1
         DO 410 I3=2,I1
           IF(MOD(I2,NIR**(I3-1)).EQ.0) THEN
             NJ1I(I3)=NJ1I(I3)+1
             IF(NJ1I(I3).GT.NIR) NJ1I(I3)=1
           ENDIF
 410     CONTINUE
!
         DSUM=0.D0
         DO 430 IR=1,NTAB(4)
           ZIR=CZERO
           DO 440 I3=1,I1
 440          ZIR=ZIR+ZTIR(NIR1-1+NJ1I(I3),IR)
           DSUM=DSUM+ABS(XR(IR)-ZIR)
 430     CONTINUE
         IF(DSUM.LT.(TOL*DBLE(NIR))) GOTO 490
 405  CONTINUE
 400  CONTINUE
      WRITE(6,585) '??'
      FG1=.FALSE.
!
!.....output
 490  CONTINUE
      IF(FG1) THEN
      DO 493 I3=I1,1,-1
       IF(I3.EQ.I1) THEN
          WRITE(6,590)
       ELSE
          WRITE(6,591)
       ENDIF

       IF(CTIR(NIR1-1+NJ1I(I3))(3:3).EQ.' ') THEN
          WRITE(6,595) CTIR(NIR1-1+NJ1I(I3))(1:2)
       ELSE
          WRITE(6,596) CTIR(NIR1-1+NJ1I(I3))(1:3)
       ENDIF
 493  CONTINUE
      WRITE(6,*) 
      ENDIF
!
!.....end 'Cornwell condition'
      ELSE
      DO IG=1,IKG
         X=XM(LKG(IG),IE)
         WRITE(6,580) DREAL(X),DIMAG(X) 
         XR(IG)=X
      ENDDO

         DO I1=1,MAXIRDG
         DO J1=1,MAXIRDG
           ZIR=CZERO
          !ZIR=IR(I1)+IR(J1)
         ENDDO
         ENDDO
      
         print*
      ENDIF
!
!.....output for band plotting (case.irrep)
!
      IF(IE.EQ.1) WRITE(5,506) NE
      DO 301 IIJ=1,ND
        WRITE(5,508) IE+IIJ-1,ND,EE(IE+IIJ-1)
        IF(.NOT.FG1) THEN
          WRITE(5,518) 0,0,0,0,0,0,0,0
        ELSEIF( .NOT.FGT )  THEN
          WRITE(5,518) 0,0,0,0,0,0,0,0
        ELSE
          DO 303 I3=I1,1,-1
            IF(CTIR(NIR1-1+NJ1I(I3))(3:3).EQ.' ') THEN
              READ(CTIR(NIR1-1+NJ1I(I3))(2:2),515) IIIR
            ELSEIF(CTIR(NIR1-1+NJ1I(I3))(3:3).EQ.'-') THEN
              READ(CTIR(NIR1-1+NJ1I(I3))(2:2),515) IIIR
              IIIR=-IIIR
            ELSEIF(CTIR(NIR1-1+NJ1I(I3))(3:3).EQ.'+') THEN
              READ(CTIR(NIR1-1+NJ1I(I3))(2:2),515) IIIR
            ELSEIF(CTIR(NIR1-1+NJ1I(I3))(4:4).EQ.' ') THEN
              READ(CTIR(NIR1-1+NJ1I(I3))(2:3),516) IIIR
            ELSEIF(CTIR(NIR1-1+NJ1I(I3))(4:4).EQ.'-') THEN
              IIIR=-IIIR
              READ(CTIR(NIR1-1+NJ1I(I3))(2:3),516) IIIR
            ELSEIF(CTIR(NIR1-1+NJ1I(I3))(4:4).EQ.'+') THEN
              READ(CTIR(NIR1-1+NJ1I(I3))(2:3),516) IIIR
            ELSE
              STOP 'wrtir: cannot find no of irrep'
            ENDIF      
            NODEG = NINT( DREAL( ZTIR(abs(IIIR),1) ) )
            IF( ABS(NODEG-ABS(ZTIR(abs(IIIR),1))).GT.0.1 ) THEN
              STOP 'wrtir: cannot find no of degen'
            ENDIF
            WRITE(5,517)  IIIR, NODEG
 303      CONTINUE          
          DO 304 I3=4,I1+1,-1
 304        WRITE(5,517) 0,0
          WRITE(5,*)
        ENDIF
 301  CONTINUE
!

 100  IE=IE+ND  
      END DO
!
      RETURN
 506  FORMAT(1I10)
 507  FORMAT(5X,A5)
 508  FORMAT(2I10,F10.6,2X,$)
 510  FORMAT(I5,1X,I5,F11.7,$)
 515  FORMAT(I1)
 516  FORMAT(I2)
 517  FORMAT( (2X,2I3),$)
 518  FORMAT(4(2X,2I3))
 530  FORMAT(16(/,3(2F11.7,2X)))
 570  FORMAT(I3,I3,F10.6,$)
 580  FORMAT(F5.2,SP,F5.2,'i ',$)
 581  FORMAT(' STOP: X=(',2F8.5,') Xold=(',2F8.5,')')
 585  FORMAT(A2)
 590  FORMAT('=',$)
 591  FORMAT(' + ',$)
 595  FORMAT(A2,$)
 596  FORMAT(A3,$)
      END



