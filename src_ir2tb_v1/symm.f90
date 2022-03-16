!=========================================================================!
! project : wave_data
! history : 07/18/2014
! authors : Zhijun Wang  ( zjwang11@hotmail.com )
! purpose : find symmetry for a specific k point
! status  : robust
! comment : These programs are distributed in the hope that they will be 
!           useful, but WITHOUT ANY WARRANTY; without even the implied 
!           warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
!=========================================================================!
module symm
use struct_data,only:init
implicit none
private

integer ,  public ,  parameter  :: MAXDG   =  6
integer ,  public ,  parameter  :: MAXIRDG =  4
integer ,  public ,  parameter  :: FLMAX   =  4
integer ,  public ,  parameter  :: MAXIR   = 48
integer ,  public ,  parameter  :: NSYM    = 48

integer           ,   parameter ::   dp   = 8

real(dp),  public ,  parameter  :: TOLDG  = 0.5E-3
real(dp),  public ,  parameter  :: PI     = 3.141592653589793238462643383279D0

!.....notation: 
!     Pi={Ri|ti}    -  crystallographic symmetry operations {rotation|transl};
!                      ti=taui+tm; 0<=taui<1 and tm is a lattice vectors.
!                      The product Pi*Pj = {Ri*Rj|Ri*tj+ti} 
!     inv(Pi)       -  inverse of Pi; inv(Pi) = {inv(Ri)|-inv(Ri)*ti}
!     Ri~           -  transpose of Ri; Ri~ not always inv(Ri) 
!     G             -  space group, consisting of the elements {Ri|ti} 
!     G(k)          -  space group of the allowed wave vector k 
!                      (little space group); those {Ri|ti} such that  
!                      k*inv(Ri)=k+K, where K is a reciprocal lattice vector.
!     Go            -  crystallographic point group of G, consisting of the 
!                      elements {Ri|0} 
!     Go(k)         -  point group of the allowed wave vector k
!     T(k)          -  translation group of allowed wave vector k; 
!                      those {1|tm} such  that exp(-i*k*tm)=1; 
!                      T(k) is an invariant (normal) subgroup to G(k)
!     F(k)          -  factor group G(k)/T(k), used for special k-points at 
!                      BZ surface in non-symmorphic crystals.
!     Fo(k)         -  point group of F(k)
!     IRp           -  p:th 'relevant' irreducible representation of 
!                      Go(k) or Fo(k). IRp is real if all representation 
!                      matrices are real. IRp can be equivalent to a real 
!                      representation (IRp ~ IRq_real), or equivalent to its 
!                      complex conjugate (IRp ~ (IRp)*), or essentially 
!                      complex (IRp !~ (IRp)* !~ IRq_real)
!                      
!
!.....basic theory, Ref. [5] p46-48:
!     <psi_l(k,r)|Pi*psi_l'(k,r)> = exp(-ik*ti)*IRp({Ri|0})_l,l' for k-points 
!                      obeying the Cornwell conditions (Ref [1] p239)
!     psi(k,r)      -  nf*SUM{ c(k+K)*fi(k+K,r) };  fi(k,r)=exp(ikr); 
!                      nf is the normalization factor
!     Pi*psi(k,r)   -  nf*SUM{ c(k+K)*exp( i(k+K)inv(Pi)r]) }
!                     =nf*SUM{ c(k+K)*exp( i(k+K)inv(Ri)(r-ti)) }
!                     =nf*SUM{ c(k+K)*exp(-i(k+K)inv(Ri)ti)*fi((k+K)inv(Ri),r)}
!                     =nf*SUM{ c((k+K')Ri) *exp(-i(k+K')ti) * fi(k+K',r) },
!                      where (k+K)inv(Ri)=k+K' => k+K=(k+G')Ri 
!                      In Ref. [6] p47 orthogonal matrices are presumed,
!                      i.e, Ri~ =inv(Ri)
!
!.....input:
!     FL(4)         -  flags:
!                      (1) true if complex eigenfunctions
!                      (2) true if spinors (i.e., with spin-orbit coupling)
!                      (3) true if crystallographic inversion symmetry  
!                      (4) true if spin-polarized
!     IZ,TAU,IIZ    -  symm. ops. {Ri|ti} read from case.struct. IIZ=inv(IZ)
!     A(*,*)        -  spin-up   part of eigenfunctions 
!     B(*,*)        -  spin-down part of eigenfunctions
!     WK=ISK/ISKDEN -  k-point; ISK and ISKDEN are integers
!     EE(NE)        -  eigenvalues in Rydberg  
!     KV(3,NV)      -  (k+K)*ISKDEN, where K is reciprocal lattice vector
! 
!...............................................................................
!.....transformations:
!     BR1(3,3)      -  lapw1 k-list coord. into Cartesian coord (recipr. space)
!     BR2(3,3)      -  reciprocal coord.   into Cartesian coord (recipr. space)
!     DR1,DR2       -  inv(BR1), inv(BR2)
!     DB1           -  DR2*BR1; from lapw1 k-list coord. into primitiv coord.
!
!.....calculated quantities:
!     RNAME,CNAME   -  name of rotaions and classes
!     RAN(4,*)      -  rotation direction of Ri
!     IAXC2(3)      -  main C2 axes
!     SU2(2,2,*)    -  spin rotation of double groups
!     JIJ(*,*)      -  equals n in Rn=Rj*Ri*inv(Rj)
!     PH(*,*)       -  phase factor  
!     LKG(IKG)      -  list of symm.op. in the little k-group G(k)
!     FGT           -  true if IR of G(k) cannot simple be associated to an IR
!                      of Go(k), which can happen for non-symmorphic crystals 
!                      if k is a special k-point at BZ surface. In this case 
!                      the local relevant IR is presented.
!     GAM(I,J)      -  IR matrix element (I,J) of Go(k) for a certain energy 
!                      and symm.op. I,J=1,...ND, where ND is the no of 
!                      degenerate states
!
!.....output:
!     XM(IKG,NE)    -  Characters for each eigenstate and each symm. ops.
!                      XM(*,*)=trace(GAM)
!
!.....common block:
!     COMMON /CTAB/    NTAB,CTIR,TTIR,ZTIR
!     CTIR(i,j)     -  Character tables of the 32 point group
!                      String format for each irreducible representation
!     TTIR          -  the title of the classes
!     ZTIR          -  same as CTIR, but complex numbers instead of string
!     NTAB( 1)      -  Table number in G.F. Koster, et al., Ref. 7
!     NTAB( 2)      -  Page number for the table in Ref. 7
!     NTAB( 3)      -  Total no. of IRs 
!     NTAB( 4)      -  No. of IRs in the single group
!     JTAB( 1)      -  Table number in S.L. Altmann, et al., Ref. 8
!     JTAB( 2)      -  Page number for the table in Ref. 8
!     
!*******************************************************************
!
      character*10, parameter::    symmlog='irtb.symm'
      logical    ,    save   ::    FL(FLMAX)
      integer    ,    save   ::    IZ(3,3,NSYM),     IIZ(3,3,NSYM)
      integer    ,    save   ::    IORD
      real(dp)   ,    save   ::    TAU(3,NSYM)
!      real(dp)   ,    save   ::    BR1(3,3)
!      real(dp)   ,    save   ::    DR1(3,3)
      real(dp)   ,    save   ::    DZ2(3,3,NSYM)
!
!
      complex(dp),    save   ::    SU2(2,2,NSYM)        
!
!
      integer    ,    save   ::    JIJ(NSYM,NSYM),  ILC(NSYM)
      real(dp)   ,    save   ::    RAN(4,NSYM)
      character*6,    save   ::    RNAM(NSYM)
      integer    ,    save   ::    LKG(NSYM)     ,  IKG
      character*3,    save   ::    GRPNAM
      integer    ,    save   ::    NCC           ,  NROT(10)
      logical    ,    save   ::    FGT
!      
      character*6,    save   ::    CNAM(NSYM)
      integer    ,    save   ::    IAXC2(3)
!      
!
      public :: ssym
      public :: ksym
!
CONTAINS

      SUBROUTINE SSYM()

      INTEGER  :: IK
      CALL WRTDATE(6)

!----!init----- 
      call init( NSYM,FLMAX,FL,IORD,IZ,TAU,DZ2 )
!
      IIZ=0
      DO IK=1,IORD
      CALL INVMATI(IZ(1,1,IK),IIZ(1,1,IK))
      ENDDO
!
!.....properties of Ri: Refs [1] p25,55; [2] p10,197; [3] p16,169
      SU2=(0.D0,0.D0);JIJ=0;ILC=0;RAN=0.D0
      CALL RMPROP(IZ,IIZ,TAU,IORD,DZ2,SU2,JIJ,ILC,RAN,RNAM)
     !print*,RNAM

      open (unit=5,file=symmlog,status='unknown')
      CALL WRTDATE(5)
      close(unit=5)

!
      END SUBROUTINE SSYM

      SUBROUTINE KSYM(KKK,WK,NMAT,NV,A,B,NUME,NE,EE,LMAT,XM)
      INTEGER , INTENT(IN) :: KKK
      REAL(DP), INTENT(IN) :: WK(3)
      INTEGER , INTENT(IN) :: NMAT
      INTEGER , INTENT(IN) :: NUME
      INTEGER , INTENT(IN) :: NV
      INTEGER , INTENT(IN) :: NE
!     REAL(DP), INTENT(IN) :: KV(3,NV)
      REAL(DP), INTENT(IN) :: EE(NUME)
      COMPLEX(DP), INTENT(IN) :: A(NMAT,NUME)
      COMPLEX(DP), INTENT(IN) :: B(NMAT,NUME)
      COMPLEX(DP), INTENT(IN) :: LMAT(NMAT,NMAT,NSYM)

!     COMPLEX(DP), INTENT(OUT) :: PH(NSYM,NMAT)
      COMPLEX(DP), INTENT(OUT) :: XM(NSYM,NUME)

!.....find elements of G(k): Refs [1] p235; [2] p89; [3] p79
      CALL KGROUP( IIZ,IORD,WK,LKG,IKG )
!     
!.....identify the point group Go(k)
       !only use FL(2)
      CALL PNTGRP(FL,LKG,IKG,ILC,JIJ,RNAM,GRPNAM,NCC,NROT)
!
!.....check the 'Cornwell condition': Ref [1] p239
      CALL CRWCND(WK,IZ,TAU,LKG,IKG,FGT)
!
      open(unit=5,file=symmlog,access='append')
!.....characterize according to Fo(k): Ref [1] p241; [2] p168
      IF(.NOT.FGT) THEN
!       THIS WILL BE IMPLEMENTED SOON     
!.......find the factor group F(k)=G(k)/T(k)
        CALL MDFPG(FL,FGT)
!.......character table of Fo(k)
!       CALL MPNTGRP()
!.......time-reversal symm: Ref [1] p158; [2] p149,210; [3] p177
        CALL TRSYMB(FL,IZ,LKG,IKG,TAU,WK,IORD,GRPNAM,RAN,KKK)             
!
      ELSE
!.......determine classes
!        call flush(6)
        CALL CLASSE(IZ,LKG,IKG,ILC,NROT,JIJ,RAN,RNAM,CNAM,GRPNAM,IAXC2)
     !print*,CNAM
     !print*,IAXC2
!
!.......time-reversal symm: Ref [1] p158; [2] p149,210; [3] p177
      !open(unit=5,file="tmp.top",status='unknown')
       !only use FL(2)
        CALL TRSYMA(FL,LKG,IKG,TAU,CNAM,WK,IORD,GRPNAM,RAN,IAXC2,KKK)
!
!.......rotate reciprocal lattice vector.
!       CALL ROTKV(NMAT,PH,L,KV,LKG,IKG,NV,ISK,ISKDEN,IZ,IIZ,TAU) 
       !CALL ROTKV(NMAT,PH,L,KV,LKG,IKG,NV,WK,IZ,IIZ,TAU) 
!      
!.......calculate the characters: Ref [6] p46
       !only use FL(2)
!       CALL CHRCT(                             FL,NE,NV,LKG,IKG,IORD,KKK)
!       call CHRCT(NMAT,NUME,A,B,SU2,EE,LMAT,PH,XM,FL,NE,NV,LKG,IKG,IORD,KKK)
        call CHRCT(NMAT,NUME,A,B,SU2,EE,LMAT,XM,FL,NE,NV,LKG,IKG,IORD,KKK)
!       call CHRCT(npmax,nband,coeffa,coeffb,FL,NE,NV,LKG,IKG,IORD,KKK)
      ENDIF
!
!.....identify the irreducible representations
       !only use FL(2)
      CALL WRTIR(NUME,FL,FGT,EE,XM,NE,LKG,IKG,CNAM,KKK)
      close(5)

      RETURN
      END SUBROUTINE KSYM

end module symm



