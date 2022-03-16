      SUBROUTINE KGROUP( IIZ,IORD,WK,LKG,IKG )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NSYM=     48)
!
      DIMENSION        IIZ(3,3,NSYM)
      DIMENSION        WK(3),WKR(3)
      DIMENSION        LKG(NSYM),DB1(3,3),DKR(3)
      DATA             TOLK/1E-4/
!******************************************************************************
!
!.....determine the group G(k) of the allowed wave vector k,
!.....i.e., finding those Ri having the property k*inv(Ri)=k+Km
!.....db1(3,3) transforms Km to units of primitiv vectors.
      LKG(:)=0
      IKG=0
      DO I=1,IORD
         DO I1=1,3
          WKR(I1)=dot_product(WK(:),IIZ(:,I1,I))-WK(I1)
         ENDDO
         DTEST=DABS(NINT(WKR(1))-WKR(1))   &
              +DABS(NINT(WKR(2))-WKR(2))   &
              +DABS(NINT(WKR(3))-WKR(3))
         IF(DTEST.LE.TOLK) THEN
           IKG=IKG+1
           LKG(IKG)=I
         ENDIF
      ENDDO
      RETURN
      END  
