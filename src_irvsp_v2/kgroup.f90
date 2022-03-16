subroutine kgroup( IIZ,IORD,WK,LKG,IKG )
use symm,only: NSYM
implicit none

integer, intent(in)  :: IIZ(3,3,NSYM),IORD
real(8), intent(in)  :: WK(3)
integer, intent(out) :: LKG(NSYM)
integer, intent(out) :: IKG


real(8) :: WKR(3),DB1(3,3),DKR(3),DTEST
real(8), parameter :: TOLK=0.1E-4
integer :: I,J
!******************************************************************************
!
!.....determine the group G(k) of the allowed wave vector k,
!.....i.e., finding those Ri having the property k*inv(Ri)=k+Km
!.....db1(3,3) transforms Km to units of primitiv vectors.

      LKG(:)=0
      IKG=0
      DO I=1,IORD
         DO J=1,3
          WKR(J)=dot_product(WK(:),IIZ(:,J,I))-WK(J)
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

end subroutine kgroup
