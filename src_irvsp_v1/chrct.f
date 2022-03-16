      SUBROUTINE CHRCT(NMAT,NUME,A,B,SU2,EE,L,PH,XM,FL,NE,NV,LKG,IKG,IORD,KKK)
      use symm,only:FLMAX,MAXIR,NSYM,MAXDG,TOLDG
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      LOGICAL          FL(FLMAX)
      COMPLEX*16       SU2(2,2,NSYM)
      COMPLEX*16       GAM(MAXDG,MAXDG)
      COMPLEX*16       PH(NSYM,NMAT),XM(NSYM,NUME)
      COMPLEX*16       A(NMAT,NUME),B(NMAT,NUME)
      COMPLEX*16, allocatable ::  A2(:,:),B2(:,:)  !NMAT,MAXDG
      COMPLEX*16       CSUM,CZERO
      DIMENSION        LKG(NSYM)
      DIMENSION        EE(NUME),L(NSYM,NMAT)
!
      DATA             CZERO/(0.0D0,0.0D0)/
!     DATA             TOL/5.0E-3/
!**********************************************************************
!
      if(allocated(a2))  deallocate (a2,b2)
      allocate (a2(nv,maxdg),b2(nv,maxdg))
      DO 4 I=1,IORD
      DO 4 J=1,NE
 4       XM(I,J)=CMPLX(1.2345678,1.2345678)

!
!.....normalization
      DO 10 J=1,NE
        SUM=0.D0
        DO 12 I=1,NV
 12        SUM=SUM +ABS(A(I,J))**2 + ABS(B(I,J))**2
        DO 14 I=1,NV
           A(I,J)=A(I,J)/SQRT(SUM)
           B(I,J)=B(I,J)/SQRT(SUM)
 14   CONTINUE
 10   CONTINUE
!
!.....characters for all energies and symm. ops.
      IE=1
      DO WHILE(IE.LE.NE)
      ND=1
      DO WHILE((IE+ND).LE.NE.AND.(EE(IE+ND)-EE(IE)).LT.TOLDG)
        ND=ND+1
      END DO
      IF(ND.GT.MAXDG) GOTO 100
!
!.....for all classes
      DO 90 IG=1,IKG
!las1
      if(ie.eq.1.and.kkk.eq.1111)then
         write(*,*)ig,ie,nd
         write(*,555) LKG(IG),((   &
           DREAL(SU2(i1,i2,LKG(IG))),  &
           DIMAG(SU2(i1,i2,LKG(IG))),i2=1,2),i1=1,2)
         write(*,*) 'end'
 555     format(i3,/,2(2(F9.6,' I*',F9.6,'     '),/))
      endif
!las2
!
!.....rotation of spin components; include phase from real space rotation 
      
      IF(FL(2)) THEN 
        J2=0
        DO 350 J=IE,IE+ND-1
        J2=J2+1
        DO 350 I=1,NV
          A2(I,J2)=(SU2(1,1,LKG(IG))*A(I,J)+ &
                    SU2(1,2,LKG(IG))*B(I,J))*PH(IG,I)
          B2(I,J2)=(SU2(2,1,LKG(IG))*A(I,J)+ &
                    SU2(2,2,LKG(IG))*B(I,J))*PH(IG,I)
 350    CONTINUE
      ELSE
        J2=0
        DO 360 J=IE,IE+ND-1
        J2=J2+1
        DO 360 I=1,NV
          A2(I,J2)=A(I,J)*PH(IG,I)
          B2(I,J2)=CZERO
 360    CONTINUE
      ENDIF
!
!.....character
      XM(LKG(IG),IE)=CZERO
      DO 370 NO=0,ND-1
      DO 370 NP=0,ND-1 
        CSUM=CZERO
        DO 380 I=1,NV
 380       CSUM=CSUM+DCONJG(A(L(IG,I),IE+NP))*A2(I,1+NO) &
                    +DCONJG(B(L(IG,I),IE+NP))*B2(I,1+NO)
        GAM(NP+1,NO+1)=CSUM
        IF(NO.EQ.NP) XM(LKG(IG),IE)=XM(LKG(IG),IE)+CSUM
 370  CONTINUE
!
!.....for all degenerate states
      DO 375 I1=IE+1,IE+ND-1
 375     XM(LKG(IG),I1)=XM(LKG(IG),IE)
!
  90  CONTINUE
 100  IE=IE+ND
      END DO
!
!.....output (characters for all bands and all symm ops.)
        DO 300 IE=1,NE
        WRITE(5,510) KKK,IE,EE(IE)
        WRITE(5,530) (DREAL(XM(IN,IE)),DIMAG(XM(IN,IE)),IN=1,IORD)
 300  CONTINUE
!
      RETURN
 510  FORMAT(I5,X,I5,F11.7,$)
 530  FORMAT(16(/,3(2F11.7,2X)))
      END



