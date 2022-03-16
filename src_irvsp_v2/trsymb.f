      SUBROUTINE TRSYMB(FL,IZ,LKG,IKG,TAU,SK,IORD,GRPNAM,ILC,RAN,KKK,KPH)
      use symm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      LOGICAL          FL(FLMAX),FG1,FG2
      CHARACTER*3      GRPNAM
      COMPLEX*16       CZERO,IMAG,KPH(NSYM),PHOLD
      DIMENSION        SK(3),TAU(3,NSYM),RT(3)
      DIMENSION        LKG(NSYM),RAN(4,NSYM)
      CHARACTER*80     CTIR(MAXIR),TTIR
      DIMENSION        NTAB(4),IZ(3,3,NSYM),ILC(NSYM)
      COMPLEX*16       ZTIR(MAXIR,MAXIR)
      COMMON /CTAB/    NTAB,CTIR,TTIR,ZTIR
      DATA             TOL/5.0E-6/
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

!
!.....output for spaghetti
!
!.....check the IR of the translation part
      WRITE(5,502) KKK,SK(1),SK(2),SK(3)
      WRITE(5,503) 'Nop', 1, 1
      DO 220 I=1,IKG
        IDFF = ABS(1-IZ(1,1,LKG(I))) + ABS(0-IZ(1,2,LKG(I))) + ABS(0-IZ(1,3,LKG(I))) &
             + ABS(0-IZ(2,1,LKG(I))) + ABS(1-IZ(2,2,LKG(I))) + ABS(0-IZ(2,3,LKG(I))) &
             + ABS(0-IZ(3,1,LKG(I))) + ABS(0-IZ(3,2,LKG(I))) + ABS(1-IZ(3,3,LKG(I))) 
        IF(IDFF.EQ.0) WRITE(5,508) '  E   ', LKG(I)
 220  CONTINUE
      WRITE(5,*) 
!-----wzj
      WRITE(6,516)
      DO I=1,NSYM
         KPH(I)=CZERO
      ENDDO

      DO IG=1,IKG 
        !WRITE(6,520) IG,LKG(IG)
         I=LKG(IG)
       IF(ILC(I).GE.1) THEN
         IF(ILC(I).EQ. 1) WRITE(6,545) 'E'    ,LKG(IG)
         IF(ILC(I).EQ. 2) WRITE(6,546)  ILC(I),LKG(IG)
         IF(ILC(I).EQ. 3) WRITE(6,546)  ILC(I),LKG(IG)
         IF(ILC(I).EQ. 4) WRITE(6,546)  ILC(I),LKG(IG)
         IF(ILC(I).EQ. 6) WRITE(6,546)  ILC(I),LKG(IG)
       ELSE
         IF(ILC(I).EQ.-1) WRITE(6,545) 'I'    ,LKG(IG)
         IF(ILC(I).EQ.-2) WRITE(6,547) -ILC(I),LKG(IG)
         IF(ILC(I).EQ.-3) WRITE(6,547) -ILC(I),LKG(IG)
         IF(ILC(I).EQ.-4) WRITE(6,547) -ILC(I),LKG(IG)
         IF(ILC(I).EQ.-6) WRITE(6,547) -ILC(I),LKG(IG)
       ENDIF

      DO 443 I1=8,19
 443     WRITE(6,'(A,$)') ' '
         ARG=-2*PI* ( SK(1)*TAU(1,LKG(IG)) & 
                     +SK(2)*TAU(2,LKG(IG)) &
                     +SK(3)*TAU(3,LKG(IG)) )
         KPH(LKG(IG))=CMPLX(DCOS(ARG),DSIN(ARG))
         PHOLD=KPH(LKG(IG))
         WRITE(6,525) DREAL(PHOLD),DIMAG(PHOLD)
         WRITE(6,531) (RAN(I2,LKG(IG)), I2=1,3)
      ENDDO
!-----wzj
!
      WRITE(6,538)
!     DO 202 IR=2,NTAB(4)
!     I=1
!     DO WHILE(CNAM(LKG(I)).NE.TTIR(6*(IR-1)+9:6*IR+8))
!       I=I+1
!       IF(I.GT.IKG) THEN
!         WRITE(6,*) 'trsym:cannot find class', CNAM(LKG(I))
!         STOP 'trsym:cannot find class'
!       ENDIF
!     END DO
!     WRITE(6,539) CNAM(LKG(I))
!202  CONTINUE

      do IG=2,IKG
      WRITE(6,579) LKG(IG)
      enddo
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
 515  FORMAT(//,'elemt ,symmetry ops, exp(-i*k*taui)'$)!by zjwang 11.9.2016
 516  FORMAT(//,'elemt ,symmetry ops, exp(-i*k*taui),  main axes',$)
!520  FORMAT(/,A6,1X,I2,$)
 520  FORMAT(/,I3,4X,I2,$)
 521  FORMAT(I3,$)
 525  FORMAT(' (',SP,F5.2,F5.2,'i)',$)
 531  FORMAT('  (',F6.3,2(',',F6.3),')',$)
 538  FORMAT(//,'bnd ndg  eigval     E  ',$)
 539  FORMAT(6X,A6,$)
 543  FORMAT(A1,$)
 544  FORMAT(' ',$)
 579  FORMAT(10X,I2,$)
 581  FORMAT(' WARNING: PH=(',2F8.5,') PHold=(',2F8.5,')')

 545  FORMAT(/,A3     ,4X,I2,$)
 546  FORMAT(/,' C',I1,4X,I2,$)
 547  FORMAT(/,'IC',I1,4X,I2,$)
      END


