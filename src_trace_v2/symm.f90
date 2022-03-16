!=========================================================================!
! project : symmetries of space groups
! history : 06/22/2016
! authors : Zhijun Wang  ( zjwang.phy@mail.com )
! purpose : get the irreps for G(k)/T(k) ~ point group
! status  : good  
! comment : These programs are distributed in the hope that they will be 
!           useful, but WITHOUT ANY WARRANTY; without even the implied 
!           warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
!=========================================================================!
module symm
use struct_data,only:init
implicit none
private

integer ,  public ,  parameter :: DP      = 8

integer ,  public ,  parameter :: MAXDG   = 16 !! MAXDG   =  6
integer ,  public ,  parameter :: FLMAX   =  4 !! FLMAX   =  4
integer ,  public ,  parameter :: NSYM    = 96 !! NSYM    = 96
                                          
real(dp),  public ,  parameter :: TOLDG   = 2.E-3
real(dp),  public ,  parameter :: PI      = 3.141592653589793238462643383279D0

!.....notation: 
!     Pi={Ri|ti}    -  crystallographic symmetry operations {rotation|transl};
!                      ti=taui+tm; 0<=taui<1 and tm is a lattice vectors.
!                      The product Pi*Pj = {Ri*Rj|Ri*tj+ti} 
!     inv(Pi)       -  inverse of Pi; inv(Pi) = {inv(Ri)|-inv(Ri)*ti}
!     Ri~           -  transpose of Ri; Ri~ not always inv(Ri) 
!     G(k)          -  space group of the allowed wave vector k 
!                      (little space group); those {Ri|ti} such that  
!                      k*inv(Ri)=k+K, where K is a reciprocal lattice vector.
!     Go            -  crystallographic point group of G, consisting of the 
!                      elements {Ri|0} 
!     Go(k)         -  point group of the allowed wave vector k
!     T(k)          -  translation group of allowed wave vector k; 
!                      those {1|tm} such  that exp(-i*k*tm)=1; 
!                      T(k) is an invariant (normal) subgroup to G(k)
!
!.....basic theory, Ref. [5] p46-48:
!     <psi_l(k,r)|Pi*psi_l'(k,r)> = exp(-ik*ti)*IRp({Ri|0})_l,l' for k-points 
!                      obeying the Cornwell conditions (Ref [1] p239)
!     psi(k,r)      -  nf*SUM{ c(k+G)*fi(k+G,r) };  fi(k,r)=exp(ikr); 
!                      nf is the normalization factor
!     Pi*psi(k,r)   -  nf*SUM{ c(k+G)*exp( i(k+G)inv(Pi)r]) }
!                     =nf*SUM{ c(k+G)*exp( i(k+G)inv(Ri)(r-ti)) }
!                     =nf*SUM{ c(k+G)*exp(-i(k+G)inv(Ri)ti)*fi((k+G)inv(Ri),r)}
!                     =nf*SUM{ c((k+G')Ri) *exp(-i(k+G')ti) * fi(k+G',r) },
!                      where (k+G)inv(Ri)=k+G' => k+G=(k+G')Ri 
!                      In Ref. [6] p47 orthogonal matrices are presumed,
!                      i.e, Ri~ =inv(Ri)
!
!.....input:
!     FL(4)         -  flags:
!                      (1) true if complex eigenfunctions
!                      (2) true if spinors (i.e., with spin-orbit coupling)
!                      (3) true if crystallographic inversion symmetry  
!                      (4) true if spin-polarized
!     IORD          -  the total number of symm. ops.
!     IZ,TAU,IIZ    -  symm. ops. {Ri|ti} read from OUTCAR. IIZ=inv(IZ)
!     A(*,*)        -  spin-up   part of eigenfunctions 
!     B(*,*)        -  spin-down part of eigenfunctions
!     WK,Kvec       -  k-point in the first Brillouin zone
!     EE(NE)        -  eigenvalues in eV  
!     KGvec(3,NV)   -  (k+G), where G is reciprocal lattice vector
! 
!.....calculated quantities:
!     SU2(2,2,*)    -  spin rotation of double groups
!     PH(*,*)       -  phase factor  
!     LKG(IKG)      -  list of symm.op. in the little k-group G(k)
!
!.....output:
!     XM(IKG,NE)    -  Characters for each eigenstate and each symm. ops.
!                      XM(*,*)=trace(GAM)
!     
!*******************************************************************
!
      integer    ,    save   ::    IORD
      logical    ,    save   ::    FL(FLMAX)
      integer    ,    save   ::    IZ(3,3,NSYM), IIZ(3,3,NSYM)
      real(dp)   ,    save   ::    TAU(3,NSYM)
      real(dp)   ,    save   ::    DZ2(3,3,NSYM)
!
!
      complex(dp),    save   ::    SU2(2,2,NSYM)        
      complex(dp),    save   ::    zerop=cmplx(0.1D-8,0.1D-8,dp)
!
!
      integer    ,    save   ::    LKG(NSYM)   ,  IKG
!      
!
      public :: ssym
      public :: ksym
!
CONTAINS

      SUBROUTINE SSYM()
      INTEGER  :: IK,J
      integer  :: IZtmp(3,3,NSYM)
      real(dp) :: TAUtmp(3,NSYM)

      !----init----- 
      call init(NSYM,FLMAX,FL,IORD,IZ,TAU,SU2)

      IIZ=0
      DO IK=1,IORD
         CALL inv_int33(IZ(1,1,IK),IIZ(1,1,IK))
      ENDDO
 
      write(614,'(I3,A32)') IORD,' : the total number of operators'
      write(624,'(I3)'    ) IORD
      write(625,'(I3)'    ) IORD
      DO IK=1,IORD
         write(624,'(3I3,$)') (IZ(1,J,IK),J=1,3)
         write(624,'(3I3,$)') (IZ(2,J,IK),J=1,3)
         write(624,'(3I3,$)') (IZ(3,J,IK),J=1,3)
         write(624,'(3F12.6,$)') (TAU(J,IK),J=1,3)
         write(624,'(4F12.6,$)') (SU2(1,J,IK)+zerop,J=1,2)
         write(624,'(4F12.6)'  ) (SU2(2,J,IK)+zerop,J=1,2)
         write(625,'(3I3,$)') (IZ(1,J,IK),J=1,3)
         write(625,'(3I3,$)') (IZ(2,J,IK),J=1,3)
         write(625,'(3I3,$)') (IZ(3,J,IK),J=1,3)
         write(625,'(3F12.6,$)') (TAU(J,IK),J=1,3)
         write(625,'(4F12.6,$)') (SU2(1,J,IK)+zerop,J=1,2)
         write(625,'(4F12.6)'  ) (SU2(2,J,IK)+zerop,J=1,2)
         
         write(614,'(2X,3I3,F12.8)') IK
         write(614,'(4X,3I3,F12.8)') (IZ(1,J,IK),J=1,3),TAU(1,IK)
         write(614,"(4X,3I3,F12.8,2(2X,'(',2F12.8,')'))") (IZ(2,J,IK),J=1,3),TAU(2,IK),(SU2(1,J,IK)+zerop,J=1,2)
         write(614,"(4X,3I3,F12.8,2(2X,'(',2F12.8,')'))") (IZ(3,J,IK),J=1,3),TAU(3,IK),(SU2(2,J,IK)+zerop,J=1,2)
      ENDDO
 
      END SUBROUTINE SSYM

      SUBROUTINE KSYM(KKK,WK,NMAT,NV,KV,A,B,NUME,NE,EE,L,PH,XM,nele)
      INTEGER , INTENT(IN) :: KKK,nele
      REAL(DP), INTENT(IN) :: WK(3)
      INTEGER , INTENT(IN) :: NMAT
      INTEGER , INTENT(IN) :: NUME
      INTEGER , INTENT(IN) :: NV
      INTEGER , INTENT(IN) :: NE
      REAL(DP), INTENT(IN) :: KV(3,NMAT)
      REAL(DP), INTENT(IN) :: EE(NUME)
      COMPLEX(DP), INTENT(IN) :: A(NMAT,NUME)
      COMPLEX(DP), INTENT(IN) :: B(NMAT,NUME)

      INTEGER    , INTENT(OUT) ::  L(NSYM,NMAT)
      COMPLEX(DP), INTENT(OUT) :: PH(NSYM,NMAT)
      COMPLEX(DP), INTENT(OUT) :: XM(NSYM,NUME)

      INTEGER    ::  IK,IG
!.....find elements of G(k): Refs [1] p235; [2] p89; [3] p79
      CALL KGROUP(IIZ,IORD,WK,LKG,IKG)
!
!.......rotate reciprocal lattice vector.
      CALL rot_KGvec(NMAT,PH,L,KV,LKG,IKG,NV,WK,IZ,IIZ,TAU) 
!
      IK=KKK
      WRITE(614,*)
      WRITE(614,"(I3,A)") IK," : ik"
      WRITE(614,"(I3,A51)") IKG," : the total number of operations in k little group"
      WRITE(624,"(I3,A51)") IKG
      do IG=1,IKG
         WRITE(614,"(2X,I3)")  LKG(IG)
         WRITE(624,"(I5,$)")   LKG(IG)
      enddo
      WRITE(624,*)

      WRITE(614,548)
      do IG=2,IKG
      WRITE(614,579) LKG(IG)
      enddo
      WRITE(614,*)
 548  FORMAT('bnd ndg  eigval     E  ',$)
 579  FORMAT(10X,I2,$)

!.......calculate the characters: Ref [6] p46
        call CHRCT(NMAT,NUME,A,B,SU2,EE,L,PH,XM,FL,NE,NV,LKG,IKG,IORD,nele)

        CALL WRTIR(NUME,EE,XM,NE,LKG,IKG,IK,nele)
!
      RETURN
      END SUBROUTINE KSYM
end module symm
