      SUBROUTINE TRSYMA(FL,LKG,IKG,TAU,CNAM,WK,IORD, &
                        GRPNAM,RAN,IAXC2,KKK)
      use symm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     PARAMETER (PI=3.141592653589793d0)      
!
      LOGICAL          FL(FLMAX),FG1,FG2
      CHARACTER*6      CNAM(NSYM)
      CHARACTER*3      GRPNAM
      COMPLEX*16       CZERO,IMAG,PH(NSYM),PHOLD
      DIMENSION        WK(3),TAU(3,NSYM)
      DIMENSION        LKG(NSYM),IAXC2(3),RAN(4,NSYM)
      CHARACTER*80     CTIR(MAXIR),TTIR
      DIMENSION        NTAB(4)
      COMPLEX*16       ZTIR(MAXIR,MAXIR)
      COMMON /CTAB/    NTAB,CTIR,TTIR,ZTIR
      DATA             TOL/5.0E-8/
! default changed because of accident degeneracy
!      DATA             TOL/5.0E-6/
      DATA             CZERO/(0.0D0,0.0D0)/,IMAG/(0.0D0,1.0D0)/
!**********************************************************************
!
!     Time-reversal symmetry can cause IRp and (IRp)* to belong to the 
!     same eigenvaue, although they span different eigenspaces. Ref [6] p21
!
!     (a) If IRp  ~ (IRp)*  ~ IRq_real:  IRp is potentially real 
!     (b) If IRp !~ (IRp)*            :  IRp is essentially complex
!     (c) If IRp  ~ (IRp)* !~ IRq_real:  IRP is pseudoreal
!
!     single groups: cases (b) and (c) give extra degeneracy
!     double groups: cases (a) and (b) give extra degeneracy
!
      DO 2 I=1,NSYM
 2       PH(I)=CZERO
!
      IF(.NOT.FL(2)) THEN
        NIR1=1
        NIR2=NTAB(4)
      ELSE
        NIR1=NTAB(4)+1
        NIR2=NTAB(3)        
      ENDIF
      NIR=NIR2-NIR1+1
!
!.....for all pairs of IRi and IRj
      IWRE=0
      IDGN=0
      DO 100 I=NIR1,NIR2
      FG1=.FALSE.    
      DO 100 J=I+1,NIR2
!
!.....for all symm.ops. Pg={Rg|tg}
      SUMT=0.D0
      SUMI=0.D0
      IWRE=1
      DO 90 IR=1,NTAB(4)
      IWRT=0
      DO 92 IG=1,IKG
      IF(CNAM(LKG(IG)).EQ.TTIR(6*(IR-1)+9:6*IR+8)) THEN
!
!........sum over |X(i,Pg)-conj(X(j,Pg))|
         IWRT=IWRT+1
         ARG=-2*PI*    ( WK(1)*TAU(1,LKG(IG)) &
                        +WK(2)*TAU(2,LKG(IG)) &
                        +WK(3)*TAU(3,LKG(IG)) )
         PH(LKG(IG))=CMPLX(DCOS(ARG),DSIN(ARG))
         SUMT=SUMT +ABS(DCONJG(PH(LKG(IG))*ZTIR(I,IR)) &
                              -PH(LKG(IG))*ZTIR(J,IR) )
         SUMI=SUMI +ABS( DIMAG(PH(LKG(IG))*ZTIR(J,IR)))
      ENDIF
 92   CONTINUE
      IF(IWRT.GT.IWRE) IWRE=IWRT
 90   CONTINUE
!
!.....check if case (a), (b), or (c)
      FG2=.FALSE.
      TOLD=TOL*IKG
      IF(.NOT.FL(2)) THEN
        IF(ABS(SUMT).LT.TOLD.AND.SUMI.GT.TOLD) FG2=.TRUE.
      ELSE
        IF(ABS(SUMT).LT.TOLD.AND.SUMI.LT.TOLD) FG2=.TRUE.
      ENDIF
!
!.....output
      IF(FG2) THEN
        IDGN=IDGN+1
        IF(IDGN.EQ.1) WRITE(6,500)
        IUI=1
        DO WHILE(CTIR(I)(IUI+1:IUI+1).NE.' '.AND.IUI.LE.7)
          IUI=IUI+1
        ENDDO
        IUJ=1
        DO WHILE(CTIR(J)(IUJ+1:IUJ+1).NE.' '.AND.IUJ.LE.7)
          IUJ=IUJ+1
        ENDDO
        WRITE(6,510) CTIR(I)(1:IUI),CTIR(J)(1:IUJ)
      ENDIF
 100  CONTINUE
!
!.....check the IR of the translation part
      WRITE(5,502) KKK,WK(1),WK(2),WK(3)
      WRITE(5,503) GRPNAM, IKG, NTAB(4)
      IF(IAXC2(1).EQ.0) THEN
        WRITE(6,515)
      ELSE
        WRITE(6,516)       
      ENDIF
!      call flush(5)
!      call flush(6)

      
      DO 200 IR=1,NTAB(4)
      IWRT=0
      INAX=0
      DO 220 IG=1,IKG
      IF(CNAM(LKG(IG)).EQ.TTIR(6*(IR-1)+9:6*IR+8)) THEN
        IWRT=IWRT+1
        IF(IWRT.EQ.1) THEN
          WRITE(5,508) CNAM(LKG(IG)), LKG(IG)
          WRITE(6,520) CNAM(LKG(IG)), LKG(IG)
          PHOLD=PH(LKG(IG))
        ELSE
          WRITE(5,521) LKG(IG)
          WRITE(6,521) LKG(IG)
          IF(ABS(PH(LKG(IG))-PHOLD).GT.TOL)THEN
             WRITE(6,581) PH(LKG(IG)),PHOLD
             WRITE(6,*) IG,LKG(IG),CNAM(LKG(IG))
             WRITE(6,*) 'Phase not equal for all elements in the class '
          ENDIF       
        ENDIF
        DO 221 I1=1,3
 221       IF(LKG(IG).EQ.IAXC2(I1)) INAX=I1
      ENDIF
 220  CONTINUE
      IWR2=14
      IF(3*IWRE.GT.IWR2) IWR2=3*IWRE
      DO 443 I1=3*IWRT,IWR2
 443     WRITE(6,'(A,$)') ' ' 
      WRITE(6,525) DREAL(PHOLD),DIMAG(PHOLD)
      IF(INAX.NE.0) WRITE(6,531) (RAN(I2,IAXC2(INAX)), I2=1,3)
 200  CONTINUE
      WRITE(5,*)
!
      WRITE(6,538)
      DO 202 IR=2,NTAB(4)
      I=1
      DO WHILE(CNAM(LKG(I)).NE.TTIR(6*(IR-1)+9:6*IR+8)) 
        I=I+1
        IF(I.GT.IKG) THEN
          WRITE(6,*) 'trsym:cannot find class', CNAM(LKG(I))
          STOP 'trsym:cannot find class'
        ENDIF
      END DO
      WRITE(6,539) CNAM(LKG(I))
 202  CONTINUE
      WRITE(6,*)
!
!
!
      RETURN
 500  FORMAT(/,7X,'Due to time-reversal symmetry, the following pairs', &
             /,7X,'of irreducible representations are degenerate:', &
             /,7X,'',$)
 502  FORMAT(I6,3F10.6,       ' # k number, k-vector')
 503  FORMAT(3X,A3,I3,I10,17X,  ' # pnt-grp, no symm.ops, no classes')
 507  FORMAT(/,1X,A6,1X,I2,$)
 508  FORMAT(/,1X,A6,/,I3,$)
 510  FORMAT('(',A,' and ',A,'); ',$)
 515  FORMAT(//,'class, symmetry ops, exp(-i*k*taui)',$)
 516  FORMAT(//,'class, symmetry ops, exp(-i*k*taui),  main axes',$)
 520  FORMAT(/,A6,1X,I2,$)
 521  FORMAT(I3,$)
 525  FORMAT(' (',SP,F5.2,F5.2,'i)',$)
 531  FORMAT('  (',F6.3,2(',',F6.3),')',$)
 538  FORMAT(//,'bnd ndg  eigval     E  ',$)
 539  FORMAT(6X,A6,$)
 543  FORMAT(A1,$)
 544  FORMAT(' ',$)
 581  FORMAT(' WARNING: PH=(',2F8.5,') PHold=(',2F8.5,')')
      END


