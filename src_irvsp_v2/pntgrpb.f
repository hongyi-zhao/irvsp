      SUBROUTINE PNTGRPB(FL,LKG,IKG,ILC,JIJ,RNAM,GRPNAM,SGN,NROT)
!     SUBROUTINE PNTGRPB(FL,LKG,IKG,ILC,    RNAM,       NCC     )
      use symm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      LOGICAL          FL(FLMAX),FG1,FG2,FGD,FGE
      COMPLEX*16       CZERO,IMAG
      CHARACTER        LIN
      CHARACTER*11     BTE
      CHARACTER*3      GRPNAM
      CHARACTER*6      RNAM(NSYM),irreplabs(0:maxir)
      DIMENSION        JIJ(NSYM,NSYM),JTAB(2)
      DIMENSION        NROT(10),LKG(NSYM),LCL(NSYM),ILC(NSYM)
      CHARACTER*80     CTIR(MAXIR),TTIR,IRNOTE
      DIMENSION        NTAB(4)
      COMPLEX*16       ZTIR(MAXIR,MAXIR),irs(NSYM,MAXIR)
      COMMON /CTAB/    NTAB,CTIR,TTIR,ZTIR
      DATA             CZERO/(0.0D0,0.0D0)/,IMAG/(0.0D0,1.0D0)/
!*******************************************************************
!
      LIN='-'
      DO 2 I1=1,MAXIR
      DO 2 I2=1,80
         IRNOTE(I2:I2)  =' '
         TTIR(I2:I2)    =' '
 2       CTIR(I1)(I2:I2)=' '
!
!.....the number of classes in G(k)
      NCC=0
      DO 12 I=1,IKG
        FG1=.FALSE.
        J1=0
        DO WHILE((J1.LT.IKG).AND.(.NOT.FG1))
          J1=J1+1
          ITT=LKG(I)!JIJ(LKG(I),LKG(J1))
          FG2=.FALSE.
          DO 13 K=1,IKG
 13          IF(ITT.EQ.LKG(K)) FG2=.TRUE.
          IF(ITT.EQ.0) STOP 'pntgrp: (1) incorrect k-group'
          IF(.NOT.FG2) STOP 'pntgrp: (2) incorrect k-group'
          IF(ITT.LT.LKG(I)) FG1=.TRUE.
        END DO
        IF(.NOT.FG1)  NCC=NCC+1 
        LCL(I)=NCC
 12   CONTINUE

!.....output
      WRITE(6,500) GRPNAM
      IF(IKG.LT.10) THEN
        WRITE(6,530) IKG,SGN
      ELSE
        WRITE(6,531) IKG,SGN
      ENDIF
      WRITE(6,533) !NTAB(1),NTAB(2)
     !WRITE(6,532) !NTAB(1),NTAB(2)
      WRITE(6,534) !JTAB(1),JTAB(2)

     !!WRITE(6,540) TTIR
     !irreplabs(0)='      '
     !irreplabs(1)='  G1  '
     !irreplabs(2)='  G2  '
     !irreplabs(3)='  G3  '
     !irreplabs(4)='  G4  '
     !irreplabs(5)='  G5  '
 
     !irs(:,:)=CZERO

     !WRITE(6,577) irreplabs(0)
     !DO I=1,IKG
     !WRITE(6,579) I !NTAB(1),NTAB(2)
     !ENDDO
     !do ir=1,5
     !   WRITE(6,578) irreplabs(ir)
     !DO I=1,IKG
     !   WRITE(6,580) irs(LKG(i),ir)
     !ENDDO
     !enddo

     !IF(.NOT.FL(2)) THEN
     !  J1=1
     !  J2=NTAB( 4)
     !ELSE
     !  J1=NTAB( 4)+1
     !  J2=NTAB( 3)
     !ENDIF
!LAS !    DO 460 J=J1,J2  
     !DO 460 J=1,  NTAB( 3)
     !  IF(J.EQ.NTAB(4)+1) WRITE(6,558) (LIN(1:1),I=1,8+NTAB(4)*6)
     !  WRITE(6,559) CTIR(J)(1:4),CTIR(J)(5:14)
     !  WRITE(6,560) (CTIR(J)(8+I*6+1:8+I*6+7),I=1,NTAB( 4)-1)
!460  CONTINUE
     !IF(FGD) WRITE(6,537)
     !IF(FGE) WRITE(6,538)
     !IF(IRNOTE(1:30).NE.' ') WRITE(6,539) IRNOTE(1:30)    
!
!.....blanc characters in order to use only 80 character for ctir
     !DO 470 J=1, NTAB( 3)  
!470     CTIR(J)(9:10)='  ' 
!
     !print*, 'wzj';stop! wzj 
      RETURN
 500  FORMAT(7X,'The k-point name is ',A3)
 530  FORMAT(7X,I1,' symmetry operations (module lattice translations) in space group ', I3)
 531  FORMAT(7X,I2,' symmetry operations (module lattice translations) in space group ', I3)
 577  FORMAT(/,1X,A6,$)
 578  FORMAT(/,7X,A6,6X,$)
 579  FORMAT(10X,I2,$)
 580  FORMAT(F5.2,SP,F5.2,'i ',$)
 532  FORMAT(7X,'The irreps for non-symmorphic groups are implemented &
             by Zhijun Wang',/,7X,'Email: zjwang11@hotmail.com')
 533  FORMAT(7X,'We do NOT classify the elements into classes.')
 534  FORMAT(7X,'Tables can be found on website: http://www.cryst.ehu.es/.',$)
 537  FORMAT(7X,'d=exp(2pi*i/8)')
 538  FORMAT(7X,'e=exp(2pi*i/3)')
 539  FORMAT(/,7X,'labeling of IRs can change due to choice of ',/,7X,'symmetry axes: ',A30)
 540  FORMAT(/,8X,A80)
 558  FORMAT(7X,81A)
 559  FORMAT(7X,A4,1X,A10,$)
 560  FORMAT(11A6)
!
      END
