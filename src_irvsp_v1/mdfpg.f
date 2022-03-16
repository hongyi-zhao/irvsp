      SUBROUTINE MDFPG(FL,FGT)
      use symm,only:FLMAX,MAXIR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      LOGICAL          FL(FLMAX),FGT
      COMPLEX*16       CZERO,IMAG
      CHARACTER*80     CTIR(MAXIR),TTIR
      DIMENSION        NTAB(4)
      COMPLEX*16       ZTIR(MAXIR,MAXIR)
      COMMON /CTAB/    NTAB,CTIR,TTIR,ZTIR
      DATA             CZERO/(0.0D0,0.0D0)/,IMAG/(0.0D0,1.0D0)/
!**********************************************************************
!
      WRITE(6,510)
!
!.....output
      WRITE(6,*)
!
 510  FORMAT(/,7X,'WILL BE IMPLEMENTED')
      RETURN
      END


