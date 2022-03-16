subroutine kgroup(IIZ,IORD,WK,LKG,IKG)
use symm
implicit none
integer, intent(in)  :: IIZ(3,3,NSYM),IORD
real(DP),intent(in)  :: WK(3)
integer, intent(out) :: LKG(NSYM),IKG

real(DP):: dtest,wkr(3)
integer :: i

!******************************************************************************
!
!.....determine the little group G(k) of the give wave vector k,
!.....i.e., finding those Ri having the property k*inv(Ri)=k+G
!.....Kvectors are given in units of primitive vectors.

LKG(:)=0
IKG=0
do i=1,iord
   wkr(:)=  matmul(wk(:),IIZ(:,:,I))-wk(:)
   dtest= dabs(wkr(1)-nint(wkr(1))) &
         +dabs(wkr(2)-nint(wkr(2))) &
         +dabs(wkr(3)-nint(wkr(3))) 
   if (dtest .le. 1.d-4) then
       IKG=IKG+1
       LKG(IKG)=I
   endif
enddo
   
end subroutine kgroup
