subroutine rot_KGvec(NMAT,PH,L,KGvec,LKG,IKG,NV,WK,IZ,IIZ,TAU)
use symm,only: DP,NSYM,PI
implicit none
integer  ,  intent(in)  :: NMAT,NV,IZ(3,3,NSYM),IIZ(3,3,NSYM)
real(dp) ,  intent(in)  :: WK(3),KGvec(3,NMAT),TAU(3,NSYM)
integer  ,  intent(in)  :: LKG(NSYM),IKG
integer  ,  intent(out) :: L(NSYM,NMAT)
complex(dp),intent(out) :: PH(NSYM,NMAT)

real(dp) :: KGvecR(3),AG
integer  :: IG,N,IN,J,II
!**********************************************************************
!
!.....find G' such that (k+G)inv(Ri)=k+G'
!
      L=0;PH=cmplx(0._dp,0._dp,dp)
      DO IG=1,IKG
      DO N=1,NV
        KGvecR(:)=matmul(KGvec(:,N),IIZ(:,:,LKG(IG)))
!
        IN=1
        DO WHILE(( (ABS(KGvecR(1)-KGvec(1,IN)) &
                   +ABS(KGvecR(2)-KGvec(2,IN)) &
                   +ABS(KGvecR(3)-KGvec(3,IN))).GT.1.D-3).AND.(IN.LE.NV))
          IN=IN+1
        END DO
!
        IF(IN.EQ.NV+1)then
           WRITE(6,*) 'rotkv: cannot find (k+G)inv(Ri)'
           WRITE(6,*) IG,N, NV,(KGvec(J,N),J=1,3)
           WRITE(6,*) IG,IN,NV,(KGvecR(J), J=1,3)
           WRITE(6,*) 'All k-vectors:'
           do  II=1,NV
              WRITE(6,*)'k+G=',(KGvec(J,II),J=1,3)
           enddo
           STOP 'rotkv: cannot find (k+G)inv(Ri)'
        endif
!
        L(IG,N)=IN
        AG=0._dp
        DO J=1,3
           AG=AG-2*PI* KGvec(J,IN)       *TAU(J,LKG(IG)) !for luis zjwang 12.6.2017
          !AG=AG-2*PI*(KGvec(J,IN)-WK(J))*TAU(J,LKG(IG)) !by zjwang 12.6.2017
          !AG=AG+2*PI*(KGvec(J,IN)-WK(J))*TAU(J,LKG(IG)) !by zjwang 11.9.2016
           PH(IG,N)=CMPLX(DCOS(AG),DSIN(AG),DP)
        ENDDO
      enddo
      enddo
!
end subroutine rot_KGvec
