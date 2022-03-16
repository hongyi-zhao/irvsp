      SUBROUTINE TRSYMB(FL,IZ,LKG,IKG,TAU,SK,IORD,GRPNAM,RAN,KKK)
      use symm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      LOGICAL          FL(FLMAX),FG1,FG2
      CHARACTER*3      GRPNAM
      COMPLEX*16       CZERO,IMAG,PH(NSYM),PHOLD
      DIMENSION        SK(3),TAU(3,NSYM)
      DIMENSION        LKG(NSYM),RAN(4,NSYM)
      CHARACTER*80     CTIR(MAXIR),TTIR
      DIMENSION        NTAB(4),IZ(3,3,NSYM)
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
 515  FORMAT(//,'class, symmetry ops, exp(-i*k*taui)'$)
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


