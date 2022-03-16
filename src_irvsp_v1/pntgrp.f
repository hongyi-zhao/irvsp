      SUBROUTINE PNTGRP(FL,LKG,IKG,ILC,JIJ,RNAM,GRPNAM,NCC,NROT)
      use symm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      LOGICAL          FL(FLMAX),FG1,FG2,FGD,FGE
      COMPLEX*16       CZERO,IMAG
      CHARACTER        LIN
      CHARACTER*11     BTE
      CHARACTER*3      GRPNAM
      CHARACTER*6      RNAM(NSYM)
      DIMENSION        JIJ(NSYM,NSYM),JTAB(2)
      DIMENSION        NROT(10),LKG(NSYM),LCL(NSYM),ILC(NSYM)
      CHARACTER*80     CTIR(MAXIR),TTIR,IRNOTE
      DIMENSION        NTAB(4)
      COMPLEX*16       ZTIR(MAXIR,MAXIR)
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
          ITT=JIJ(LKG(I),LKG(J1))
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
!
!.....number of different rotations (E,C2,C3,C4,...)
      DO 3 I=1,10
 3       NROT(I)=0
      DO 14 I=1,IKG
        IPD=ILC(LKG(I))
        IF(ABS(IPD).EQ.1.AND.IPD.GT.0) NROT( 1)=NROT( 1)+1
        IF(ABS(IPD).EQ.2.AND.IPD.GT.0) NROT( 2)=NROT( 2)+1
        IF(ABS(IPD).EQ.3.AND.IPD.GT.0) NROT( 3)=NROT( 3)+1
        IF(ABS(IPD).EQ.4.AND.IPD.GT.0) NROT( 4)=NROT( 4)+1
        IF(ABS(IPD).EQ.6.AND.IPD.GT.0) NROT( 5)=NROT( 5)+1

        IF(ABS(IPD).EQ.1.AND.IPD.LT.0) NROT( 6)=NROT( 6)+1
        IF(ABS(IPD).EQ.2.AND.IPD.LT.0) NROT( 7)=NROT( 7)+1
        IF(ABS(IPD).EQ.3.AND.IPD.LT.0) NROT( 8)=NROT( 8)+1
        IF(ABS(IPD).EQ.4.AND.IPD.LT.0) NROT( 9)=NROT( 9)+1
        IF(ABS(IPD).EQ.6.AND.IPD.LT.0) NROT(10)=NROT(10)+1
 14   CONTINUE
!
!.....tests
      NSUM=0
      DO 16 I=1,10
 16      NSUM=NSUM+NROT(I)
      IF(NSUM.NE.IKG)  STOP 'pntgrp: incorrect NROT'
      IF(NROT(1).NE.1) STOP 'pntgrp: wrong identity op.'
      IF(NROT(6).GT.1) STOP 'pntgrp: wrong inverse op.'
      IF(NROT(6).EQ.1) THEN
      DO 17 I=2,5
 17      IF(NROT(I).NE.NROT(I+5)) STOP 'classes: Cn not eq ICn'
      ENDIF
!
!.....identify the point group G(k)     
      WRITE(BTE,'(5I1,1X,5I1)') NROT( 1),NROT( 2),NROT( 3),NROT( 4), &
                               NROT( 5),NROT( 6),NROT( 7),NROT( 8), &
                               NROT( 9),NROT(10)
!  
      FG1=.FALSE. 
      IF(IKG.EQ. 1.AND.NCC.EQ. 1.AND.BTE.EQ.'10000 00000')  GOTO  10
      IF(IKG.EQ. 2.AND.NCC.EQ. 2.AND.BTE.EQ.'10000 10000')  GOTO  10
      IF(IKG.EQ. 2.AND.NCC.EQ. 2.AND.BTE.EQ.'11000 00000')  GOTO  30
      IF(IKG.EQ. 2.AND.NCC.EQ. 2.AND.BTE.EQ.'10000 01000')  GOTO  30
      IF(IKG.EQ. 4.AND.NCC.EQ. 4.AND.BTE.EQ.'11000 11000')  GOTO  30
      IF(IKG.EQ. 4.AND.NCC.EQ. 4.AND.BTE.EQ.'13000 00000')  GOTO  60
      IF(IKG.EQ. 4.AND.NCC.EQ. 4.AND.BTE.EQ.'11000 02000')  GOTO  60
      IF(IKG.EQ. 8.AND.NCC.EQ. 8.AND.BTE.EQ.'13000 13000')  GOTO  60
      IF(IKG.EQ. 4.AND.NCC.EQ. 4.AND.BTE.EQ.'11020 00000')  GOTO  90
      IF(IKG.EQ. 4.AND.NCC.EQ. 4.AND.BTE.EQ.'11000 00020')  GOTO  90
      IF(IKG.EQ. 8.AND.NCC.EQ. 8.AND.BTE.EQ.'11020 11020')  GOTO  90
      IF(IKG.EQ. 8.AND.NCC.EQ. 5.AND.BTE.EQ.'15020 00000')  GOTO 120
      IF(IKG.EQ. 8.AND.NCC.EQ. 5.AND.BTE.EQ.'11020 04000')  GOTO 120
      IF(IKG.EQ. 8.AND.NCC.EQ. 5.AND.BTE.EQ.'13000 02020')  GOTO 120
      IF(IKG.EQ.16.AND.NCC.EQ.10.AND.BTE.EQ.'15020 15020')  GOTO 120
      IF(IKG.EQ. 3.AND.NCC.EQ. 3.AND.BTE.EQ.'10200 00000')  GOTO 160
      IF(IKG.EQ. 6.AND.NCC.EQ. 6.AND.BTE.EQ.'10200 10200')  GOTO 160
      IF(IKG.EQ. 6.AND.NCC.EQ. 3.AND.BTE.EQ.'13200 00000')  GOTO 180
      IF(IKG.EQ. 6.AND.NCC.EQ. 3.AND.BTE.EQ.'10200 03000')  GOTO 180
      IF(IKG.EQ.12.AND.NCC.EQ. 6.AND.BTE.EQ.'13200 13200')  GOTO 180
      IF(IKG.EQ. 6.AND.NCC.EQ. 6.AND.BTE.EQ.'11202 00000')  GOTO 210
      IF(IKG.EQ. 6.AND.NCC.EQ. 6.AND.BTE.EQ.'10200 01002')  GOTO 210
      IF(IKG.EQ.12.AND.NCC.EQ.12.AND.BTE.EQ.'11202 11202')  GOTO 210
      IF(IKG.EQ.12.AND.NCC.EQ. 6.AND.BTE.EQ.'17202 00000')  GOTO 240
      IF(IKG.EQ.12.AND.NCC.EQ. 6.AND.BTE.EQ.'11202 06000')  GOTO 240
      IF(IKG.EQ.12.AND.NCC.EQ. 6.AND.BTE.EQ.'13200 04002')  GOTO 240
      IF(IKG.EQ.24.AND.NCC.EQ.12.AND.BTE.EQ.'17202 17202')  GOTO 240
      IF(IKG.EQ.12.AND.NCC.EQ. 4.AND.BTE.EQ.'13800 00000')  GOTO 280
      IF(IKG.EQ.24.AND.NCC.EQ. 8.AND.BTE.EQ.'13800 13800')  GOTO 280
      IF(IKG.EQ.24.AND.NCC.EQ. 5.AND.BTE.EQ.'19860 00000')  GOTO 300
      IF(IKG.EQ.24.AND.NCC.EQ. 5.AND.BTE.EQ.'13800 06060')  GOTO 300
      IF(IKG.EQ.48.AND.NCC.EQ.10.AND.BTE.EQ.'19860 19860')  GOTO 300
      WRITE(6,*) IKG,NCC,'___',BTE
      DO 5 I=1,IKG
 5      write(6,*) i,lkg(i),ilc(lkg((i))),rnam(lkg(i))
      STOP 'classes: Wrong k-group'
!
   10 GRPNAM= 'C1 '
      IGRP= 1
      NTAB( 1) =   3
      NTAB( 2) =  31
      NTAB( 3) =   2
      NTAB( 4) =   1
      JTAB( 1) =   1
      JTAB( 2) = 108
      TTIR     ='           E  ' 
      CTIR( 1) ='G1   A1    1  '
      CTIR( 2) ='G2   A1/2  1  '
      IF(BTE(7:7).EQ.'1') GOTO 20
      GOTO 400
   20 GRPNAM= 'Ci '
      IGRP= 2
      NTAB( 1) =   6
      NTAB( 2) =  32
      JTAB( 1) =  11
      JTAB( 2) = 138
      GOTO 399
   30 GRPNAM= 'C2 '
      IGRP= 3
      NTAB( 1) =   9
      NTAB( 2) =  33
      NTAB( 3) =   4
      NTAB( 4) =   2
      JTAB( 1) =   2
      JTAB( 2) = 110
      TTIR     ='           E    C2  ' 
      CTIR( 1) ='G1   A     1     1  '
      CTIR( 2) ='G2   B     1    -1  '
      CTIR( 3) ='G3  1E1/2  1     i  '
      CTIR( 4) ='G4  2E1/2  1    -i  '
      IF(BTE.EQ.'10000 01000')  GOTO  40
      IF(BTE(7:7).EQ.'1')       GOTO  50
      GOTO 400
   40 GRPNAM= 'Cs '
      IGRP= 4
      JTAB( 1) =  12
      JTAB( 2) = 140
      TTIR     ='           E   IC2  ' 
      CTIR( 1) ='G1   A`    1     1  '
      CTIR( 2) ='G2   A"    1    -1  '
      CTIR( 3) ='G3  1E1/2  1     i  '
      CTIR( 4) ='G4  2E1/2  1    -i  '
      GOTO 400
   50 GRPNAM= 'C2h'
      IGRP= 5
      NTAB( 1) =  15
      NTAB( 2) =  35
      JTAB( 1) =  60
      JTAB( 2) = 532
      GOTO 399
   60 GRPNAM= 'D2 '
      IGRP= 6
      NTAB( 1) =  17
      NTAB( 2) =  36
      NTAB( 3) =   5
      NTAB( 4) =   4
      JTAB( 1) =  22
      JTAB( 2) = 194
      TTIR     ='           E    C2    C2`   C2" '
      CTIR( 1) ='G1   A1    1     1     1     1  '
      CTIR( 2) ='G2   B3    1    -1     1    -1  '
      CTIR( 3) ='G3   B1    1     1    -1    -1  '
      CTIR( 4) ='G4   B2    1    -1    -1     1  '
      CTIR( 5) ='G5   E1/2  2     0     0     0  '
      IRNOTE(1:30)= 'G2 <-> G3 <-> G4              '
      IF(BTE.EQ.'11000 02000')  GOTO  70
      IF(BTE(7:7).EQ.'1')       GOTO  80
      GOTO 400
   70 GRPNAM= 'C2v'
      IGRP= 7
      JTAB( 1) =  50
      JTAB( 2) = 482
      TTIR     ='           E    C2   IC2`  IC2" '
      CTIR( 1) ='G1   A1    1     1     1     1  '
      CTIR( 2) ='G2   B2    1    -1     1    -1  '
      CTIR( 3) ='G3   A2    1     1    -1    -1  '
      CTIR( 4) ='G4   B1    1    -1    -1     1  '
      CTIR( 5) ='G5   E1/2  2     0     0     0  '
      IRNOTE(1:30)= 'G2 <-> G4                     '
      GOTO 400
   80 GRPNAM= 'D2h'
      IGRP= 8
      NTAB( 1) =  23
      NTAB( 2) =  39
      JTAB( 1) =  31
      JTAB( 2) = 247
      GOTO 399
   90 GRPNAM= 'C4 '
      IGRP= 9
      NTAB( 1) =  25
      NTAB( 2) =  40
      NTAB( 3) =   8
      NTAB( 4) =   4
      JTAB( 1) =   4
      JTAB( 2) = 114
      TTIR     ='           E    C4    C2    C4- '
      CTIR( 1) ='G1   A     1     1     1     1  '
      CTIR( 2) ='G2   B     1    -1     1    -1  '
      CTIR( 3) ='G3  2E     1     i    -1    -i  '
      CTIR( 4) ='G4  1E     1    -i    -1     i  '
      CTIR( 5) ='G5  2E1/2  1     d     i     d* '
      CTIR( 6) ='G6  1E1/2  1     d*   -i     d  '
      CTIR( 7) ='G7  2E3/2  1    -d     i    -d* '
      CTIR( 8) ='G8  1E3/2  1    -d*   -i    -d  '
      IF(BTE.EQ.'11000 00020')  GOTO 100
      IF(BTE(7:7).EQ.'1')       GOTO 110
      GOTO 400
  100 GRPNAM= 'S4 '
      IGRP=10
      JTAB( 1) =  13
      JTAB( 2) = 144
      TTIR     ='           E   IC4    C2   IC4- '
      GOTO 400
  110 GRPNAM= 'C4h'
      IGRP=11
      NTAB( 1) =  31
      NTAB( 2) =  43
      JTAB( 1) =  62
      JTAB( 2) = 538
      GOTO 399
  120 GRPNAM= 'D4 '
      IGRP=12
      NTAB( 1) =  33
      NTAB( 2) =  45
      NTAB( 3) =   7
      NTAB( 4) =   5
      JTAB( 1) =  24
      JTAB( 2) = 199
      TTIR     ='           E   2C4    C2   2C2`  2C2" '
      CTIR( 1) ='G1   A1    1     1     1     1     1  '
      CTIR( 2) ='G2   A2    1     1     1    -1    -1  '
      CTIR( 3) ='G3   B1    1    -1     1     1    -1  '
      CTIR( 4) ='G4   B2    1    -1     1    -1     1  '
      CTIR( 5) ='G5   E     2     0    -2     0     0  '
      CTIR( 6) ='G6   E1/2  2     /2    0     0     0  '
      CTIR( 7) ='G7   E3/2  2    -/2    0     0     0  '
      IRNOTE(1:30)= 'G3 <-> G4                     '
      IF(BTE.EQ.'11020 04000')  GOTO 130
      IF(BTE.EQ.'13000 02020')  GOTO 140
      IF(BTE(7:7).EQ.'1')       GOTO 150
      GOTO 400
  130 GRPNAM= 'C4v'
      IGRP=13
      JTAB( 1) =  52
      JTAB( 2) = 489
      TTIR     ='           E   2C4    C2  2IC2` 2IC2" '
      GOTO 400
  140 GRPNAM= 'D2d'
      IGRP=14
      JTAB( 1) =  41
      JTAB( 2) = 366
      TTIR     ='           E  2IC4    C2   2C2` 2IC2" '
      IRNOTE(1:30)= '                              '
      GOTO 400
  150 GRPNAM= 'D4h'
      IGRP=15
      NTAB( 1) =  40
      NTAB( 2) =  50
      JTAB( 1) =  33
      JTAB( 2) = 258
      GOTO 399
  160 GRPNAM= 'C3 '
      IGRP=16
      NTAB( 1) =  42
      NTAB( 2) =  51
      NTAB( 3) =   6
      NTAB( 4) =   3
      JTAB( 1) =   3
      JTAB( 2) = 112
      TTIR     ='           E    C3    C3- '
      CTIR( 1) ='G1   A     1     1     1  '
      CTIR( 2) ='G2  2E     1     e     e* '
      CTIR( 3) ='G3  1E     1     e*    e  '
      CTIR( 4) ='G4  1E1/2  1    -e*   -e  '
      CTIR( 5) ='G5  2E1/2  1    -e    -e* '
      CTIR( 6) ='G6   A3/2  1    -1    -1  '
      IRNOTE(1:30)= 'G2 <-> G3 and G4 <-> G5       '
      IF(BTE(7:7).EQ.'1') GOTO 170
      GOTO 400
  170 GRPNAM= 'C3i'
      IGRP=17
      NTAB( 1) =  47
      NTAB( 2) =  54
      JTAB( 1) =  14
      JTAB( 2) = 146
      GOTO 399
  180 GRPNAM= 'D3 '
      IGRP=18
      NTAB( 1) =  49
      NTAB( 2) =  55
      NTAB( 3) =   6
      NTAB( 4) =   3
      JTAB( 1) =  23
      JTAB( 2) = 196
      TTIR     ='           E   2C3   3C2  '
      CTIR( 1) ='G1   A1    1     1     1  '
      CTIR( 2) ='G2   A2    1     1    -1  '
      CTIR( 3) ='G3   E     2    -1     0  '
      CTIR( 4) ='G4   E1/2  2     1     0  '
      CTIR( 5) ='G5  1E3/2  1    -1     i  '
      CTIR( 6) ='G6  2E3/2  1    -1    -i  '
      IF(BTE.EQ.'10200 03000')  GOTO 190
      IF(BTE(7:7).EQ.'1')       GOTO 200
      GOTO 400
  190 GRPNAM= 'C3v'
      IGRP=19
      JTAB( 1) =  51
      JTAB( 2) = 485
      TTIR     ='           E   2C3  3IC2  '
      GOTO 400
  200 GRPNAM= 'D3d'
      IGRP=20
      NTAB( 1) =  55
      NTAB( 2) =  58
      JTAB( 1) =  42
      JTAB( 2) = 371
      GOTO 399
  210 GRPNAM= 'C6 '
      IGRP=21
      NTAB( 1) =  57
      NTAB( 2) =  59
      NTAB( 3) =  12
      NTAB( 4) =   6
      JTAB( 1) =   6
      JTAB( 2) = 119
      TTIR     ='           E    C6    C3    C2    C3-   C6- '
      CTIR( 1) ='G1   A     1     1     1     1     1     1  '
      CTIR( 2) ='G2  2E2    1     e*    e     1     e*    e  '
      CTIR( 3) ='G3  1E2    1     e     e*    1     e     e* '
      CTIR( 4) ='G4   B     1    -1     1    -1     1    -1  '
      CTIR( 5) ='G5  2E1    1    -e*    e    -1     e*   -e  '
      CTIR( 6) ='G5  1E1    1    -e     e*   -1     e    -e* '
      CTIR( 7) ='G7  1E1/2  1   -ie    -e*    i    -e    ie* '
      CTIR( 8) ='G8  2E1/2  1    ie*   -e    -i    -e*  -ie  '
      CTIR( 9) ='G9  2E5/2  1    ie    -e*   -i    -e   -ie* '
      CTIR(10) ='G10 1E5/2  1   -ie*   -e     i    -e*   ie  '
      CTIR(11) ='G11 1E3/2  1    -i    -1     i    -1     i  '
      CTIR(12) ='G12 2E3/2  1     i    -1    -i    -1    -i  '
      IF(BTE.EQ.'10200 01002')  GOTO 220
      IF(BTE(7:7).EQ.'1')       GOTO 230
      GOTO 400
  220 GRPNAM= 'C3h'
      IGRP=22
      JTAB( 1) =  61
      JTAB( 2) = 534
      TTIR     ='           E   IC6    C3   IC2    C3-  IC6- '
      CTIR( 1) ='G1   A`    1     1     1     1     1     1  '
      CTIR( 2) ='G2  2E`    1     e*    e     1     e*    e  '
      CTIR( 3) ='G3  1E`    1     e     e*    1     e     e* '
      CTIR( 4) ='G4   A"    1    -1     1    -1     1    -1  '
      CTIR( 5) ='G5  2E"    1    -e*    e    -1     e*   -e  '
      CTIR( 6) ='G5  1E"    1    -e     e*   -1     e    -e* '
      CTIR( 7) ='G7  1E1/2  1   -ie    -e*    i    -e    ie* '
      CTIR( 8) ='G8  2E1/2  1    ie*   -e    -i    -e*  -ie  '
      CTIR( 9) ='G9  2E5/2  1    ie    -e*   -i    -e   -ie* '
      CTIR(10) ='G10 1E5/2  1   -ie*   -e     i    -e*   ie  '
      CTIR(11) ='G11 1E3/2  1    -i    -1     i    -1     i  '
      CTIR(12) ='G12 2E3/2  1     i    -1    -i    -1    -i  '
      GOTO 400
  230 GRPNAM= 'C6h'
      IGRP=23
      NTAB( 1) =  63
      NTAB( 2) =  64
      JTAB( 1) =  64
      JTAB( 2) = 546
      GOTO 399
  240 GRPNAM= 'D6 '
      IGRP=24
      NTAB( 1) =  65
      NTAB( 2) =  67
      NTAB( 3) =   9
      NTAB( 4) =   6
      JTAB( 1) =  26
      JTAB( 2) = 207
      TTIR     ='           E    C2   2C3   2C6   3C2`  3C2" '
      CTIR( 1) ='G1   A1    1     1     1     1     1     1  '
      CTIR( 2) ='G2   A2    1     1     1     1    -1    -1  '
      CTIR( 3) ='G3   B1    1    -1     1    -1     1    -1  '
      CTIR( 4) ='G4   B2    1    -1     1    -1    -1     1  '
      CTIR( 5) ='G5   E1    2    -2    -1     1     0     0  '
      CTIR( 6) ='G6   E2    2     2    -1    -1     0     0  '
      CTIR( 7) ='G7   E1/2  2     0     1     /3    0     0  '
      CTIR( 8) ='G8   E5/2  2     0     1    -/3    0     0  '
      CTIR( 9) ='G9   E3/2  2     0    -2     0     0     0  '
      IRNOTE(1:30)= 'G3 <-> G4                     '
      IF(BTE.EQ.'11202 06000')  GOTO 250
      IF(BTE.EQ.'13200 04002')  GOTO 260
      IF(BTE(7:7).EQ.'1')       GOTO 270
      GOTO 400
  250 GRPNAM= 'C6v'
      IGRP=25
      JTAB( 1) =  54
      JTAB( 2) = 497
      TTIR     ='           E    C2   2C3   2C6  3IC2` 3IC2" '
      GOTO 400
  260 GRPNAM= 'D3h'
      IGRP=26
      JTAB( 1) =  32
      JTAB( 2) = 250
      TTIR     ='           E   IC2   2C3  2IC6   3C2` 3IC2" '
      CTIR( 1) ='G1   A1`   1     1     1     1     1     1  '
      CTIR( 2) ='G2   A2`   1     1     1     1    -1    -1  '
      CTIR( 3) ='G3   A1"   1    -1     1    -1     1    -1  '
      CTIR( 4) ='G4   A2"   1    -1     1    -1    -1     1  '
      CTIR( 5) ='G5   E"    2    -2    -1     1     0     0  '
      CTIR( 6) ='G6   E`    2     2    -1    -1     0     0  '
      CTIR( 7) ='G7   E1/2  2     0     1     /3    0     0  '
      CTIR( 8) ='G8   E5/2  2     0     1    -/3    0     0  '
      CTIR( 9) ='G9   E3/2  2     0    -2     0     0     0  '
      IRNOTE(1:30)= '                              '
      GOTO 400
  270 GRPNAM= 'D6h'
      IGRP=27
      NTAB( 1) =  72
      NTAB( 2) =  76
      JTAB( 1) =  35
      JTAB( 2) = 277
      GOTO 399
  280 GRPNAM= 'T  '
      IGRP=28
      NTAB( 1) =  74
      NTAB( 2) =  79
      NTAB( 3) =   7
      NTAB( 4) =   4
      JTAB( 1) =  70
      JTAB( 2) = 590
      TTIR     ='           E   3C2   4C3   4C3- ' 
      CTIR( 1) ='G1   A     1     1     1     1  '
      CTIR( 2) ='G2  1E     1     1     e     e* '
      CTIR( 3) ='G3  2E     1     1     e*    e  '
      CTIR( 4) ='G4   T     3    -1     0     0  '
      CTIR( 5) ='G5   E1/2  2     0     1     1  '
      CTIR( 6) ='G6  1F3/2  2     0     e     e* '
      CTIR( 7) ='G7  2F3/2  2     0     e*    e  '
      IRNOTE(1:30)= 'G2 <-> G3 and G6 <-> G7       '
      IF(BTE(7:7).EQ.'1') GOTO 290
      GOTO 400
  290 GRPNAM= 'Th '
      IGRP=29
      NTAB( 1) =  79
      NTAB( 2) =  86
      JTAB( 1) =  72
      JTAB( 2) = 633
      GOTO 399
  300 GRPNAM= 'O  '
      IGRP=30
      NTAB( 1) =  81
      NTAB( 2) =  88
      NTAB( 3) =   8
      NTAB( 4) =   5
      JTAB( 1) =  72
      JTAB( 2) = 633
      TTIR     ='           E   8C3   3C2   6C4   6C2` ' 
      CTIR( 1) ='G1   A1    1     1     1     1     1  '
      CTIR( 2) ='G2   A2    1     1     1    -1    -1  '
      CTIR( 3) ='G3   E     2    -1     2     0     0  '
      CTIR( 4) ='G4   T1    3     0    -1     1    -1  '
      CTIR( 5) ='G5   T2    3     0    -1    -1     1  '
      CTIR( 6) ='G6   E1/2  2     1     0     /2    0  '
      CTIR( 7) ='G7   E5/2  2     1     0    -/2    0  '
      CTIR( 8) ='G8   F3/2  4    -1     0     0     0  '
      IF(BTE.EQ.'13800 06060')  GOTO 310
      IF(BTE(7:7).EQ.'1')       GOTO 320
      GOTO 400
  310 GRPNAM= 'Td '
      IGRP=31 
      JTAB( 1) =  73
      JTAB( 2) = 637
      TTIR     ='           E   8C3   3C2  6IC4  6IC2` ' 
      GOTO 400 
  320 GRPNAM= 'Oh '
      IGRP=32
      NTAB( 1) =  87
      NTAB( 2) = 103
      JTAB( 1) =  71
      JTAB( 2) = 641
      GOTO 399
! 
 399  CALL DPRODCI(CTIR(1),TTIR,NTAB(1))
 400  CONTINUE
!
!.....evaluate ctir to complex numbers
      FGD=.FALSE.
      FGE=.FALSE.
      DO 410 I1=1,NTAB(3)
      J1=0
      DO 410 I2=8+1,8+NTAB(4)*6,6
        J1=J1+1 
        ZTIR(I1,J1)=CZERO
        IF(    CTIR(I1)(I2+3:I2+3).EQ.'1') THEN
          ZTIR(I1,J1)=1.D0
        ELSEIF(CTIR(I1)(I2+3:I2+3).EQ.'2') THEN
          ZTIR(I1,J1)=2.D0
        ELSEIF(CTIR(I1)(I2+3:I2+3).EQ.'3') THEN
          ZTIR(I1,J1)=3.D0
        ELSEIF(CTIR(I1)(I2+3:I2+3).EQ.'4') THEN
          ZTIR(I1,J1)=4.D0
        ELSEIF(CTIR(I1)(I2+3:I2+3).EQ.'i') THEN
          ZTIR(I1,J1)=IMAG
        ELSEIF(CTIR(I1)(I2+3:I2+3).EQ.'/') THEN
          READ(CTIR(I1)(I2+4:I2+4),'(I1)') IDI
          ZTIR(I1,J1)=DSQRT(DBLE(IDI))
        ELSEIF(CTIR(I1)(I2+3:I2+3).EQ.'e') THEN
          FGE=.TRUE.
          DPH=2.D0*PI/3.D0
          IF(CTIR(I1)(I2+4:I2+4).EQ.'*') THEN
             ZTIR(I1,J1)=CMPLX(DCOS(DPH),-DSIN(DPH))
          ELSE
             ZTIR(I1,J1)=CMPLX(DCOS(DPH), DSIN(DPH))
          ENDIF
          IF(CTIR(I1)(I2+2:I2+2).EQ.'i')  &
                      ZTIR(I1,J1)= ZTIR(I1,J1)*IMAG
        ELSEIF(CTIR(I1)(I2+3:I2+3).EQ.'d') THEN
          FGD=.TRUE.
          DPH=2.D0*PI/8.D0
          IF(CTIR(I1)(I2+4:I2+4).EQ.'*') THEN
             ZTIR(I1,J1)=CMPLX(DCOS(DPH),-DSIN(DPH))
          ELSE
             ZTIR(I1,J1)=CMPLX(DCOS(DPH), DSIN(DPH))
          ENDIF
        ENDIF
        IF(CTIR(I1)(I2+2:I2+2).EQ.'-') ZTIR(I1,J1)=-ZTIR(I1,J1)
        IF(CTIR(I1)(I2+1:I2+1).EQ.'-') ZTIR(I1,J1)=-ZTIR(I1,J1)
 410  CONTINUE
!     
!.....output
      WRITE(6,500) GRPNAM
      IF(IKG.LT.10) THEN
        WRITE(6,530) IKG,NCC
      ELSE
        WRITE(6,531) IKG,NCC
      ENDIF
      WRITE(6,533) NTAB(1),NTAB(2)
      WRITE(6,534) JTAB(1),JTAB(2)

      WRITE(6,540) TTIR
      IF(.NOT.FL(2)) THEN
        J1=1
        J2=NTAB( 4)
      ELSE
        J1=NTAB( 4)+1
        J2=NTAB( 3)
      ENDIF
!LAS      DO 460 J=J1,J2  
      DO 460 J=1,  NTAB( 3)
        IF(J.EQ.NTAB(4)+1) WRITE(6,558) (LIN(1:1),I=1,8+NTAB(4)*6)
        WRITE(6,559) CTIR(J)(1:4),CTIR(J)(5:14)
        WRITE(6,560) (CTIR(J)(8+I*6+1:8+I*6+7),I=1,NTAB( 4)-1)
 460  CONTINUE
      IF(FGD) WRITE(6,537)
      IF(FGE) WRITE(6,538)
      IF(IRNOTE(1:30).NE.' ') WRITE(6,539) IRNOTE(1:30)    
!
!.....blanc characters in order to use only 80 character for ctir
      DO 470 J=1, NTAB( 3)  
 470     CTIR(J)(9:10)='  ' 
!
      RETURN
 500  FORMAT(7X,'The point group is ',A3)
 530  FORMAT(7X,I1,' symmetry operations in ',I2,' classes')
 531  FORMAT(7X,I2,' symmetry operations in ',I2,' classes')
 533  FORMAT(7X,'Table ',I2,'   on page ',I3,' in Koster  et al [7]')
 534  FORMAT(7X,'Table ',I2,'.4 on page ',I3,' in Altmann et al [8]')
 537  FORMAT(7X,'d=exp(2pi*i/8)')
 538  FORMAT(7X,'e=exp(2pi*i/3)')
 539  FORMAT(/,7X,'labeling of IRs can change due to choice of ',/,7X,'symmetry axes: ',A30)
 540  FORMAT(/,8X,A80)
 558  FORMAT(7X,81A)
 559  FORMAT(7X,A4,1X,A10,$)
 560  FORMAT(11A6)
!
      END
