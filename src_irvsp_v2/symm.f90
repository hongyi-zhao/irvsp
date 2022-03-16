!=========================================================================!
! project : symmorphic space groups
! history : 07/18/2014
! authors : Zhijun Wang  ( zjwang11@hotmail.com )
! purpose : get the irreps for G(k)/T(k) ~ point group
! status  : robust
! comment : These programs are distributed in the hope that they will be 
!           useful, but WITHOUT ANY WARRANTY; without even the implied 
!           warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
!=========================================================================!
module symm
use struct_data,only:init
use nonsymm
implicit none
private

integer ,  public ,  parameter  ::  MAXDG   = 16 !! MAXDG   =  6
integer ,  public ,  parameter  ::  MAXIRDG =  4 !! MAXIRDG =  4
integer ,  public ,  parameter  ::  FLMAX   =  4 !! FLMAX   =  4
integer ,  public ,  parameter  ::  MAXIR   = 48 !! MAXIR   = 48
integer ,  public ,  parameter  ::  NSYM    = 96 !! NSYM    = 48
                                    
integer ,            parameter  ::  DP      = 8
                                    
real(dp),  public ,  parameter  ::  TOLDG  = 2.E-3
real(dp),  public ,  parameter  ::  PI     = 3.141592653589793238462643383279_DP

integer ,  public ,    save     ::  nmin   = 1  !! the minimum number of bands

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
      character*10, parameter::    symmlog='vasp.symm'
      integer    ,    save   ::    IORD
      logical    ,    save   ::    FL(FLMAX)
      integer    ,    save   ::    IZ(3,3,NSYM),     IIZ(3,3,NSYM)
      real(dp)   ,    save   ::    TAU(3,NSYM)
!      real(dp)   ,    save   ::    BR1(3,3)
!      real(dp)   ,    save   ::    DR1(3,3)
      real(dp)   ,    save   ::    DZ2(3,3,NSYM)
!
!
      complex(dp),    save   ::    SU2(2,2,NSYM)
      complex(dp),    save   ::    SU2v1(2,2,NSYM)
      complex(dp),    save   ::    zerop=cmplx(0.1D-8,0.1D-8,dp)
!
!
      integer    ,    save   ::    JIJ(NSYM,NSYM),  ILC(NSYM)
      real(dp)   ,    save   ::    RAN(4,NSYM)
      character*6,    save   ::    RNAM(NSYM)
      integer    ,    save   ::    LKG(NSYM)     ,  IKG
      character*3,    save   ::    GRPNAM
      integer    ,    save   ::    NCC           ,  NROT(10)
      logical    ,    save   ::    FGT,FGT2
!      
      character*6,    save   ::    CNAM(NSYM)
      integer    ,    save   ::    IAXC2(3)
!      
      integer    ,    save   ::    ver_n   = 0  !! version control
      integer    ,    save   ::    SGN     = 0  !! space group number

!
      public :: ssym
      public :: ksym
!
CONTAINS

      SUBROUTINE SSYM(sgn_,vn_,nmax)
      use struct_data,only:NUME,NKPTS
      integer  , intent(in) ::  sgn_, vn_
      integer  , intent(out)::  nmax
      INTEGER  :: IK,J
      integer  :: IZtmp(3,3,NSYM)
      real(dp) :: TAUtmp(3,NSYM)
      CALL WRTDATE(6)

      SGN=sgn_;ver_n=vn_
!----!init----- 
      SU2=(0._DP,0._DP) 
      call init(NSYM,FLMAX,FL,IORD,IZ,TAU,DZ2,SU2)
!
      IIZ=0
      DO IK=1,IORD
      CALL INVMATI(IZ(1,1,IK),IIZ(1,1,IK))
      ENDDO

!-----added by zjwang for nonsymmorphic space groups
      CALL MPNTGRP(IORD,IZ,IIZ,TAU,DZ2,SU2,SGN)

      IF(nmax==0) nmax=NUME !wzj 2019.12.30


!注意 irvsp 只有当 KPOINTS 中k点个数少于10个时才会输出trace.txt 文件（可以通过简单修改源码解除该限制）。

!werner@X10DAi-00:~/Public/repo/github.com/zjwang11/irvsp.git$ rg '::[ ]*[IJKLMN][^ ]+' | grep -i point | grep -Po '::[ ]*\K[^ ]+'| sort -u |xargs -I{} rg "{}[ ]*(<|\.lt\.)[ ]*10"
!src_irvsp_v2_release/irrep.f90
!89:      IF(NKPTS<10) THEN
!103:      IF(NKPTS<10) THEN

!src_irvsp_v2_release/wave_data.f90
!129:      if(NKPTS<10.and.kkk==1) then

!src_irvsp_v2_release/symm.f90
!184:      IF(NKPTS<10) THEN


!      IF(NKPTS<10) THEN
      OPEN(unit=624,file='trace.txt',form='formatted',status='unknown')
      IF(nmax/=0) WRITE(624,'(I3)') nmax-nmin+1
      IF(     FL(2)) WRITE(624,'(I3)')  1
      IF(.NOT.FL(2)) WRITE(624,'(I3)')  0
      write(624,'(I3)'    ) IORD
      DO IK=1,IORD
         write(624,'(3I3,$)') (IZ(1,J,IK),J=1,3)
         write(624,'(3I3,$)') (IZ(2,J,IK),J=1,3)
         write(624,'(3I3,$)') (IZ(3,J,IK),J=1,3)
         write(624,'(3F12.6,$)') (TAU(J,IK),J=1,3)
         write(624,'(4F12.6,$)') (SU2(1,J,IK)+zerop,J=1,2)
         write(624,'(4F12.6)'  ) (SU2(2,J,IK)+zerop,J=1,2)
      ENDDO
!      ENDIF

!
!.....properties of Ri: Refs [1] p25,55; [2] p10,197; [3] p16,169
     !SU2=(0.D0,0.D0);JIJ=0;ILC=0;RAN=0.D0
                      JIJ=0;ILC=0;RAN=0.D0
      SU2v1=(0._DP,0._DP) 
      CALL RMPROP(IZ,IIZ,TAU,IORD,DZ2,SU2v1,SU2,JIJ,ILC,RAN,RNAM)
     !print*,RNAM;stop

      open (unit=5,file=symmlog,status='unknown')
      CALL WRTDATE(5)
      close(unit=5)
!
      END SUBROUTINE SSYM

      SUBROUTINE KSYM(KKK,WK,NMAT,NV,KV,A,B,NUME,NE,EE,L,PH,XM)
      INTEGER , INTENT(IN) :: KKK
      REAL(DP), INTENT(IN) :: WK(3)
      INTEGER , INTENT(IN) :: NMAT
      INTEGER , INTENT(IN) :: NUME
      INTEGER , INTENT(IN) :: NV
      INTEGER , INTENT(IN) :: NE
     !REAL(DP), INTENT(IN) :: KV(3,NMAT)
      REAL(DP), INTENT(INOUT) :: KV(3,NMAT)
      REAL(DP), INTENT(IN) :: EE(NUME)
      COMPLEX(DP), INTENT(IN)  ::  A(NMAT,NUME)
      COMPLEX(DP), INTENT(IN)  ::  B(NMAT,NUME)

      INTEGER    , INTENT(OUT) ::  L(NSYM,NMAT)
      COMPLEX(DP), INTENT(OUT) :: PH(NSYM,NMAT)
      COMPLEX(DP), INTENT(OUT) :: XM(NSYM,NUME)

      INTEGER     ::  IK,IG,IR
      REAL(DP)    :: WKt(3)
      COMPLEX(DP) :: KPH(NSYM)

!.....find elements of G(k): Refs [1] p235; [2] p89; [3] p79
      CALL KGROUP( IIZ,IORD,WK,LKG,IKG )
!
!.....check the 'Cornwell condition': Ref [1] p239
      CALL CRWCND(WK,IZ,TAU,LKG,IKG,FGT)
      IF    (ver_n==1) THEN; FGT2=.TRUE.    !  - version I  (v1)
      ELSEIF(ver_n==2) THEN; FGT2=.FALSE.   !  - version II (v2)
      ELSEIF(ver_n==3) THEN; FGT2=FGT       !  - version III(v3)
      ELSE   ; FGT2=.FALSE.; FGT=FGT2       !  - version IV (v4)
      ENDIF
!
!.....characterize according to Fo(k): Ref [1] p241; [2] p168
      IF(.NOT.FGT) THEN
        IF(FGT2 .NEQV. FGT) THEN
           WRITE(6,'(/,7X,A,/)') 'IT IS IMPLEMENTED IN VERSION II(v2) !!!'
           RETURN
        ENDIF
        open(unit=5,file=symmlog,access='append')
!
!--------2019.10.18
!.......identify the nonsymmorphic point group
        WKt(:)=WK(:)
        CALL getkid( FL(3),WKt,IK,GRPNAM,IIZ,NMAT,KV,IR)
        CALL KGROUP( IIZ,IORD,WKt,LKG,IKG )
!--------2019.10.18

        WRITE(624,"(I3,A51)") IKG
        do IG=1,IKG
           WRITE(624,"(I5,$)")   LKG(IG)
        enddo
        WRITE(624,*)
!
!.......rotate reciprocal lattice vector.
!       CALL ROTKV(NMAT,PH,L,KV,LKG,IKG,NV,ISK,ISKDEN,IZ,IIZ,TAU) 
        CALL ROTKV(NMAT,PH,L,KV,LKG,IKG,NV,WKt,IZ,IIZ,TAU) 
!
     !  only use FL(2)
        CALL PNTGRPB(FL,LKG,IKG,ILC,JIJ,RNAM,GRPNAM,SGN,NROT)
     !  CALL TRSYMA(FL,LKG,IKG,TAU,CNAM,WK,IORD,GRPNAM,RAN,IAXC2,KKK)
        call dumptableofIrs(IK,6)
!
!       THIS WILL BE IMPLEMENTED SOON     
!.......find the factor group F(k)=G(k)/T(k)
!.......time-reversal symm: Ref [1] p158; [2] p149,210; [3] p177
        CALL TRSYMB(FL,IZ,LKG,IKG,TAU,WKt,IORD,GRPNAM,ILC,RAN,KKK,KPH)             
     !  CALL TRSYMA(FL,LKG,IKG,TAU,CNAM,WK,IORD,GRPNAM,RAN,IAXC2,KKK)
!.......character table of Fo(k)
!.......calculate the characters: Ref [6] p46
!       CALL CHRCT(                             FL,NE,NV,LKG,IKG,IORD,KKK)
        call CHRCT(NMAT,NUME,A,B,SU2,EE,L,PH,XM,FL,NE,NV,LKG,IKG,IORD,KKK,IR)

       !only use FL(2)
        CALL WRTIRB(NUME,FL,FGT,EE,XM,NE,LKG,IKG,IK,KPH,nmin,IR)
!
      ELSE
        IF(FGT2 .NEQV. FGT) THEN
           WRITE(6,'(/,7X,A,/)') 'THIS OPTION IS CLOSED IN VERSION II(v2) !!!'
           RETURN
        ENDIF
        open(unit=5,file=symmlog,access='append')
!
!.......identify the point group Go(k)
     !  only use FL(2)
        CALL PNTGRP(FL,LKG,IKG,ILC,JIJ,RNAM,GRPNAM,NCC,NROT)
!
!.......determine classes
!       call flush(6)
        CALL CLASSE(IZ,LKG,IKG,ILC,NROT,JIJ,RAN,RNAM,CNAM,GRPNAM,IAXC2)
     !  print*,CNAM
     !  print*,IAXC2
!
!.......rotate reciprocal lattice vector.
!       CALL ROTKV(NMAT,PH,L,KV,LKG,IKG,NV,ISK,ISKDEN,IZ,IIZ,TAU) 
        CALL ROTKV(NMAT,PH,L,KV,LKG,IKG,NV,WK,IZ,IIZ,TAU) 
!
!.......time-reversal symm: Ref [1] p158; [2] p149,210; [3] p177
     !  only use FL(2)
        CALL TRSYMA(FL,LKG,IKG,TAU,CNAM,WK,IORD,GRPNAM,RAN,IAXC2,KKK)
!      
!.......calculate the characters: Ref [6] p46
     !  only use FL(2)
!       CALL CHRCT(                             FL,NE,NV,LKG,IKG,IORD,KKK)
        call CHRCT(NMAT,NUME,A,B,SU2v1,EE,L,PH,XM,FL,NE,NV,LKG,IKG,IORD,KKK,0)
!       call CHRCT(npmax,nband,coeffa,coeffb,FL,NE,NV,LKG,IKG,IORD,KKK)
!.......identify the irreducible representations
     !  only use FL(2)
        CALL WRTIR(NUME,FL,FGT,EE,XM,NE,LKG,IKG,CNAM)
!
      ENDIF
      close(5)
      RETURN
      END SUBROUTINE KSYM

end module symm



