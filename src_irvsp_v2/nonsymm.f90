!=========================================================================!
! project : character tables for nonsymmorphic space groups
! history : 07/18/2014
! authors : Zhijun Wang  ( zjwang11@hotmail.com )
! purpose : find symmetry for a certain k point
! status  : robust
! comment : These programs are distributed in the hope that they will be 
!           useful, but WITHOUT ANY WARRANTY; without even the implied 
!           warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
!=========================================================================!
module nonsymm
implicit none
private
      integer ,    parameter  :: MAXDG  = 16
      integer ,    parameter  :: MAXIRDG=  4
      integer ,    parameter  :: FLMAX  =  4
      integer ,    parameter  :: MAXIR  = 48
      integer ,    parameter  :: NSYM   = 96
      
      integer ,    parameter  ::     dp = 8
      
      real(dp),    parameter  :: TOLDG  = 2.E-3
      real(dp),    parameter  :: PI     = 3.141592653589793238462643383279_DP
      
      real(dp),    parameter  :: epsil  = 0.0000001_dp
      
      complex(dp), parameter  :: czero=(0._dp,0._dp),cone=(1._dp,0._dp)
      
      integer     :: DoubNum            ! the number of elements in the double space group
      integer     :: SymElemR(3,3,NSYM) ! the R part of the symmetry operators {R|t}
      real(dp)    :: SymElemt(3,NSYM)   ! the t part of the symmetry operators {R|t}
      complex(dp) :: SymElemS(2,2,NSYM) ! the SU2 symmetry operators {S}
      
      !---------irreps of little groups
      integer      , parameter  ::                          ktypes=40
      real(dp)     , save       ::                samplek(3,ktypes) ! Coordinates 
      character*2  , save       ::            samplekname  (ktypes) ! stardand name
      integer      , save       ::                  nirreps(ktypes) ! Number of Irreps for little groups
      integer      , save       ::                  sirreps(ktypes) ! Number of Irreps for single little groups
      integer      , save       ::                  antisym(ktypes) ! the existance of antiunitary symmetrys
      character*5  , save       ::        Irrepsname (MAXIR,ktypes) ! Names of Irreps
      integer      , save       ::        Herringrule(MAXIR,ktypes) ! reality of Irreps
          !  \delta(A^2)= 1,     equivalent, real    ! trival for      integer spin
          !  \delta(A^2)=-1,     equivalent, complex ! trival for half-integer spin
          !  \delta(A^2)= 0, not equivalent, pair with itself.
      integer      , save       ::           nelelittle    (ktypes) ! Number of elements in litte groups
      integer      , save       ::      labels(2,NSYM,MAXIR,ktypes)
          !  the first label can take two values: 0 whern the symmetry element does
          !          NOT belong to the little group or 1, when it does belong to it.
          !  the second label can take two values: 1 if the phase does not contain
          !  the u,v,w parameters or 2 if the phase does contain at least one of them.
      complex(dp)  , save       ::   tableTraces(NSYM,MAXIR,ktypes) ! for comparison
      complex(dp)  , save       ::   chartTraces(NSYM,MAXIR,ktypes)
      real(dp)     , save       ::   coeff_uvw(3,NSYM,MAXIR,ktypes)
      character*12 , save       ::   factsTraces(NSYM,MAXIR,ktypes)

      integer      , save       ::   nktp ! the number of k types in the files
      character*15 , save       ::   chkpoint(ktypes) 
      character*10 , save       ::   spacegroupsymbol
      
      ! 1/2   ->  0.5
      ! 1/3   ->  0.333333
      ! u     ->  0.123  : su
      ! v     ->  0.313  : sv
      ! w     ->  0.427  : sw
      ! 1 - u ->  0.877  : s1_u
      ! -2 u  -> -0.246  : s_2u
      ! -u    -> -0.123  : s_u
      ! 1 + u ->  1.123  : s1u
      ! 1 - v ->  0.687  : s1_v
      !------------codes for u,v,w
      real(dp),    parameter  :: su=0.123_dp,s1_u=0.877_dp,s_2u=-0.246_dp
      real(dp),    parameter  :: sv=0.313_dp,s1_v=0.687_dp
      real(dp),    parameter  :: sw=0.427_dp,s_u=-0.123_dp,s1u = 1.123_dp
      
      real(dp),     save      :: Kc2p(3,3),p2cR(3,3)
      
      real(dp),    parameter  :: ds1=1._dp/sqrt(3._dp), ds2=0.5_dp/sqrt(3._dp)
      real(dp),    parameter  :: inv3=1._dp/3._dp,dinv3=2._dp/3._dp
      real(dp),    parameter  :: Pabc(3,3) = RESHAPE( &
                  (/  1.00000 ,  0.00000 ,  0.00000   &
                    , 0.00000 ,  1.00000 ,  0.00000   &
                    , 0.00000 ,  0.00000 ,  1.00000/) &
                    ,(/3,3/) )
      real(dp),    parameter  :: Cabc(3,3) = RESHAPE( &
                  (/  0.50000 , -0.50000 ,  0.00000   &
                    , 0.50000 ,  0.50000 ,  0.00000   &
                    , 0.00000 ,  0.00000 ,  1.00000/) &
                    ,(/3,3/) )
      real(dp),    parameter  :: Cabc68(3,3) = RESHAPE( &
                  (/  0.50000 ,  0.50000 ,  0.00000   &
                    ,-0.50000 ,  0.50000 ,  0.00000   &
                    , 0.00000 ,  0.00000 ,  1.00000/) &
                    ,(/3,3/) )
      real(dp),    parameter  :: Babc(3,3) = RESHAPE( &
!   SG#5, #8, #9, #12, #15.
                  (/  0.50000 ,  0.50000 ,  0.00000   &
                    ,-0.50000 ,  0.50000 ,  0.00000   &
                    , 0.00000 ,  0.00000 ,  1.00000/) &
                    ,(/3,3/) )
!     real(dp),    parameter  :: Aabc(3,3) = RESHAPE( &
!   SG#38-41.
!                 (/  0.00000 , -0.50000 , -0.50000   &
!                   , 0.00000 ,  0.50000 , -0.50000   &
!                   , 1.00000 ,  0.00000 ,  0.00000/) &
!                   ,(/3,3/) )
      real(dp),    parameter  :: Aabc(3,3) = RESHAPE( &
!   SG#38-41.
                  (/  1.00000 ,  0.00000 ,  0.00000   &
                    , 0.00000 ,  0.50000 ,  0.50000   &
                    , 0.00000 , -0.50000 ,  0.50000/) &
                    ,(/3,3/) )
     !real(dp),    parameter  :: Rabc(3,3) = RESHAPE( &
     !            (/    ds2   , -0.500_dp,   inv3     &
     !              ,   ds2   ,  0.500_dp,   inv3     &
     !              ,  -ds1   ,  0.000_dp,   inv3  /) &
     !              ,(/3,3/) )
      real(dp),    parameter  :: Rabc(3,3) = RESHAPE( &
                  (/  dinv3   ,  inv3    ,   inv3     &
                    , -inv3   ,  inv3    ,   inv3     &
                    , -inv3   ,-dinv3    ,   inv3  /) &
                    ,(/3,3/) )
      real(dp),    parameter  :: Fabc(3,3) = RESHAPE( &
                  (/  0.00000 ,  0.50000 ,  0.50000   &
                    , 0.50000 ,  0.00000 ,  0.50000   &
                    , 0.50000 ,  0.50000 ,  0.00000/) &
                    ,(/3,3/) )
      real(dp),    parameter  :: Iabc(3,3) = RESHAPE( &
                  (/ -0.50000 ,  0.50000 ,  0.50000   &
                    , 0.50000 , -0.50000 ,  0.50000   &
                    , 0.50000 ,  0.50000 , -0.50000/) &
                    ,(/3,3/) )
!
      public :: mpntgrp
      public :: dumptableofIrs
      public :: getkid
      public :: wrtirb
!
CONTAINS

      SUBROUTINE MPNTGRP(IORD,IZ,IIZ,TAU,DZ2,SU2,sgn)
      integer  , intent(in)    ::    IORD
      integer  , intent(inout) ::    IZ(3,3,NSYM), IIZ(3,3,NSYM)
      real(dp) , intent(inout) ::    TAU(3,NSYM)
      real(dp) , intent(inout) ::    DZ2(3,3,NSYM)
      complex(dp),intent(inout)::    SU2(2,2,NSYM)        
      integer  , intent(in)    ::    sgn

      integer       ::    IZt(3,3,IORD), IIZt(3,3,IORD)
      real(dp)      ::    TAUt(3,IORD)
      real(dp)      ::    DZ2t(3,3,IORD)
      complex(dp)   ::    SU2t(2,2,IORD)

      integer       ::    series(NSYM)
      real(dp)      ::    dr(3),DTEST,tmp33(3,3)
      integer       ::    Numk,antiss,tnir,wi

      character*1   ::    tag1 
      character*3   ::    tag3, csgn
      character*30  ::    TabIrrepName
      character*40  ::    ListIrrep 
      integer       ::    i,j,k,ierr,itmp,jtmp
      real(dp)      ::    Df(2,4),abcde(5)
      real(dp)      ::    ktmp(3),ttmp(1,3)
      character*5   ::    irtmp
      character*2   ::    nametmp,nametmp2
      character*15  ::    ckpoint
      integer       ::    iir,nele,ikt

      
      CHARACTER(len=180) :: spgpath, spgfile
      !https://fortran-lang.discourse.group/t/conditional-compilation-based-on-environment-variables/2983/9?u=hongyi
      INTEGER            :: spgstat
      
     !**********************************************************************
      IF    (sgn<10 ) THEN; write(csgn,'(I1)') sgn
      ELSEIF(sgn<100) THEN; write(csgn,'(I2)') sgn
      ELSE                ; write(csgn,'(I3)') sgn
      ENDIF
      

!https://www.intel.com/content/www/us/en/develop/documentation/fortran-compiler-oneapi-dev-guide-and-reference/top/language-reference/a-to-z-reference/g-1/get-environment-variable.html
!GET_ENVIRONMENT_VARIABLE
!Intrinsic Subroutine:
!Gets the value of an environment variable.
!CALL GET_ENVIRONMENT_VARIABLE
!(name[,value,length,status,trim_name,errmsg])
!status
!    (Output; optional) Must be a scalar of type integer. If specified, it is assigned a value of 0 if the environment variable exists and either has no value or its value is successfully assigned to value
!    .
!    It is assigned a value of -1 if the value
!    argument is present and has a length less than the significant length of the environment variable value. It is assigned a value of 1 if the environment variable does not exist. For other error conditions, it is assigned a processor-dependent value greater than 2.
!      
!https://fortran-lang.discourse.group/t/conditional-compilation-based-on-environment-variables/2983/9?u=hongyi
call get_environment_variable('IRVSPDATA',spgpath,status=spgstat)
if (spgstat == 1) then
   write(6,*) "Environment variable 'IRVSPDATA' must be provided."
   stop
end if


      spgfile = TRIM(spgpath)//'/kLittleGroups/kLG_'//TRIM(csgn)//'.data'
      WRITE(*,*) "SPGFILE :", spgfile
      open(11,file=spgfile,status='old',form='unformatted')
      open (unit=9, file='Littlegroup_'//trim(csgn)//'.cht',status='unknown')
     !open (unit=99,file='KlistSgroup_'//trim(csgn)//'.cht',status='unknown')
      read(11) DoubNum, spacegroupsymbol
      write(9,*)  spacegroupsymbol

      !https://www.cmt.york.ac.uk/compmag/resources/FortranGuide.pdf
      ! 3.3 Precision
      !       So far we have been using what is known as single precision variables; these
      !       have well defined limits dependent on machine architecture and OS. But what
      !       happens when you want to represent (as we will) much bigger numbers, with
      !       much more precision? The answer is to use double precision variables. My
      !       preffered way of using double precisions is as follows: Declare a double precision
      !       parameter (dp - one point zero ’d’ zero) and then use this in all other declarations
      !       by putting the variable name in brackets after variable type:
      !       integer, parameter :: dp=kind(1.0d0) ! Do this before ANY others
      !       real(dp) :: a, b ! a and b are now double precision
      !       NB: There is a rounding error pitfall to watch out for here. If variable ’a’
      !       is double precision, one must be careful when using it to make use of it being
      !       double. That is:
      !       a = 0.0
      !       May give a to be 0.0000000248192746. The extra precision contains whatever
      !       was previously in memory. Therefore we can make sure it is double precision
      !       zero (0.0000000000000000) by appending a dp tag to the end:
      ! a = 0.0_dp
      ! Note the underscore and then the dp. For clarity, you may call your DP variable
      ! anything, for example:
      ! integer, parameter :: example = kind(1.0d0)
      ! real(example) :: aVariable
      ! aVariable = 0.0_example
      ! Is equivalent to before.
      ! See: SELECTED_REAL_KIND and SELECTED_INT_KIND for even more precisional control.
      
      Kc2p=0.0_dp
      IF(spacegroupsymbol(1:1)=='P')  THEN
        Kc2p=Pabc
      ELSEIF(spacegroupsymbol(1:1)=='C')   THEN
        Kc2p=Cabc
        if (sgn==68) Kc2p=Cabc68
      ELSEIF(spacegroupsymbol(1:1)=='B')   THEN
        Kc2p=Babc
      ELSEIF(spacegroupsymbol(1:1)=='A')   THEN
        Kc2p=Aabc
      ELSEIF(spacegroupsymbol(1:1)=='R')   THEN
        Kc2p=Rabc
      ELSEIF(spacegroupsymbol(1:1)=='F')   THEN
        Kc2p=Fabc
      ELSEIF(spacegroupsymbol(1:1)=='I')   THEN
        Kc2p=Iabc
      ELSE
        STOP "Error in space-group-symbol" 
      ENDIF

      WRITE(9,"(' From conv. to prim. reciprocal space  (DB1)')")
      WRITE(9,'(3(3F16.8,/))') Kc2p
      call invreal33(Kc2p,p2cR)
      WRITE(9,"(' From prim. to conv. reciprocal space  (DR1)')")
      WRITE(9,'(3(3F16.8,/))') p2cR

         !!-----------
         !WRITE(99,"(' A primitive cell in the convention of VASP ')")
         !WRITE(99,'(3(3F16.8,/))') Kc2p
         !WRITE(99,"('The Inverse of the above matrix             ')")
         !WRITE(99,'(3(3F16.8,/))') p2cR
         !WRITE(99,"(35X,A22)") '(conv.) -->    (prim.)'
         !!-----------
      
      SymElemR(:,:,:)=0;SymElemS(2,2,NSYM)=czero
      Df=0.0_dp 
      write(9,600)  ' Ri     taui     spin transf.' 
      do i=1,DoubNum
         read(11) SymElemR(:,:,i),SymElemt(:,i),Df(:,:)
         SymElemS(1,1,i)=cmplx(Df(1,1)*dcos(PI*Df(2,1)),Df(1,1)*dsin(PI*Df(2,1)),dp)
         SymElemS(1,2,i)=cmplx(Df(1,2)*dcos(PI*Df(2,2)),Df(1,2)*dsin(PI*Df(2,2)),dp)
         SymElemS(2,1,i)=cmplx(Df(1,3)*dcos(PI*Df(2,3)),Df(1,3)*dsin(PI*Df(2,3)),dp)
         SymElemS(2,2,i)=cmplx(Df(1,4)*dcos(PI*Df(2,4)),Df(1,4)*dsin(PI*Df(2,4)),dp)
      enddo

      tmp33(:,:)=0.0_dp
      do i=1,DoubNum/2
         write(9,601) i, SymElemR(:,1,i), SymElemt(1,i)
         write(9,602)    SymElemR(:,2,i), SymElemt(2,i), SymElemS(1,1,i), SymElemS(1,2,i)
         write(9,602)    SymElemR(:,3,i), SymElemt(3,i), SymElemS(2,1,i), SymElemS(2,2,i)
        !-- conventional cell -> primitive cell  by zjwang 11.7.2016
          tmp33  (:,:  )=  DBLE(Transpose(SymElemR(:,:,i)))
         SymElemR(:,:,i)=  nint(matmul(matmul(p2cR(:,:), tmp33(:,:) ),Kc2p(:,:)))
         SymElemt(:,i  )=matmul(p2cR,SymElemt(:,i))
        !-- conventional cell -> primitive cell  by zjwang 11.7.2016
      enddo

      !do i = 1, IORD
      !  write(*,*) SymElemR(:,:,i)
      !enddo 
      !write(*,*)
      !do i = 1, IORD
      !  write(*,*) IZ(:,:,i)
      !enddo 

      !------reorder the elements
       IF(IORD/=DoubNum/2) STOP "Elements Error"
       series=0; series(1)=1
       Do i=2,IORD
       Do j=2,IORD
          tmp33(:,:)=IZ(:,:,j)
          IF(tmp33(1,1)==SymElemR(1,1,i).and.tmp33(2,1)==SymElemR(2,1,i).and.tmp33(3,1)==SymElemR(3,1,i) &
       .and. tmp33(1,2)==SymElemR(1,2,i).and.tmp33(2,2)==SymElemR(2,2,i).and.tmp33(3,2)==SymElemR(3,2,i) &
       .and. tmp33(1,3)==SymElemR(1,3,i).and.tmp33(2,3)==SymElemR(2,3,i).and.tmp33(3,3)==SymElemR(3,3,i) ) THEN
          series(i)=j; EXIT
          ENDIF
       ENDDO
          DR(:)=SymElemt(:,i)-TAU(:,j)
          DTEST=DABS(NINT(DR(1))-DR(1))   &
               +DABS(NINT(DR(2))-DR(2))   &
               +DABS(NINT(DR(3))-DR(3))
          IF(DTEST.LE.epsil) THEN
             TAU(:,j)=SymElemt(:,i)
          ELSE
             WRITE(*,'(I3,3F10.5)') I,SymElemt(:,I)
             WRITE(*,'(I3,3F10.5)') J,TAU(:,J)
             WRITE(*,*) "TAU Error !!!"
             WRITE(*,'(A)') " 1. Please check the space group (SG) number by the program 'phonopy':"
             WRITE(*,'(A)') "    ##$ ./phonopy --symmetry --tolerance 0.01 -c POSCAR"
             WRITE(*,'(A)') " 2. Please get 'POSCAR' in the standard (default) setting of the SG: "
             WRITE(*,'(A)') "    a. Download the code 'pos2aBR' below: "
             WRITE(*,'(A)') "       https://github.com/zjwang11/irvsp/blob/master/src_pos2aBR.tar.gz"
             WRITE(*,'(A)') "    b. or you can sent the POSCAR to 'jcgao@iphy.ac.cn' for help."
             WRITE(*,'(A)') "    c. Rerun your DFT calculations with the standard 'POSCAR_std'."
             STOP
          ENDIF
       ENDDO
      !write(*,'(10I5)') series(1:IORD)!;stop
       Do i=1,IORD
          j=series(i) 
           IZt(:,:,i) =  IZ(:,:,j)
          IIZt(:,:,i) = IIZ(:,:,j)
           TAUt( :,i) =  TAU( :,j)
          DZ2t(:,:,i) = DZ2(:,:,j)
          SU2t(:,:,i) = SU2(:,:,j)
       ENDDO
         IZ(:,:,1:IORD) =  IZt(:,:,1:IORD)
        IIZ(:,:,1:IORD) = IIZt(:,:,1:IORD)
         TAU( :,1:IORD) =  TAUt( :,1:IORD)
        DZ2(:,:,1:IORD) = DZ2t(:,:,1:IORD)
        SU2(:,:,1:IORD) = SU2t(:,:,1:IORD)   
       !SU2(:,:,1:IORD)=SymElemS(:,:,1:IORD) ! by zjwang on 12.6.2017 
        CALL getsign(Doubnum,SymElemR,SymElemS,IZ,SU2) !by zjwang on 10.25.2019
      !!---debug
      ! write(5556,'(I5)') Doubnum
      ! write(5556,'(9I5)') SymElemR(:,:,1:Doubnum)
      ! write(5556,'(4F18.12)') SymElemS(:,:,1:Doubnum)
      ! write(5556,'(9I5)') IZ(:,:,1:Doubnum/2)
      ! write(5556,'(4F18.12)') SU2(:,:,1:Doubnum/2)
      ! stop
      !!---debug

     !**********************************************************************
     !--------character tables-----
      nirreps(:)=0
      sirreps(:)=0
      antisym(:)=0
      Herringrule(:,:)=-2

      labels=0
      iir=0;ikt=0
      tableTraces(:,:,:)=czero
      chartTraces(:,:,:)=czero
      coeff_uvw(:,:,:,:)=czero
      factsTraces(:,:,:)='            '
      chkpoint(:)='***************'
      nametmp='  '
      read(11)  Numk,tnir

      DO wi=1,tnir
        read(11) ListIrrep
        nametmp2=nametmp
        read(ListIrrep,*) ktmp(:),antiss,irtmp, itmp,itmp,nametmp,jtmp
        IF(nametmp/=nametmp2) THEN 
           IF(ikt>0) THEN
              nirreps(ikt)=iir
             !call dumptableofIrs(ikt,9)
              call Kreal2string(samplek(:,ikt),ckpoint) 
              ttmp(1,1:3)=samplek(1:3,ikt)
             !WRITE(99,"(1X,I2,A21,3F6.3,A5,3F8.3,I4,A4)") ikt, samplekname(ikt)//' ('//ckpoint//')' & 
             !  , samplek(:,ikt),' --> ',matmul(ttmp(:,:),Kc2p(:,:)),ikt,samplekname(ikt)
           ENDIF
          !------out
           ikt=ikt+1; iir=0
           antisym(ikt)=antiss
        ENDIF
!
        samplek(:,ikt)=ktmp(:)
        samplekname(ikt)=nametmp
        iir=iir+1
        nele=0
        irk:DO j=1,Doubnum
            read(11) itmp,itmp
            IF(itmp==1) THEN
               nele=nele+1
               labels(1,j,iir,ikt)=itmp
               read(11) itmp;labels(2,j,iir,ikt)=itmp
               abcde(:)=0._dp
               IF(itmp==1) THEN
                  read(11) abcde(1:2)
               ELSEIF(itmp==2) THEN
                  read(11) abcde(1:5)
                  coeff_uvw(:,j,iir,ikt)=abcde(3:5)
                  write(factsTraces(j,iir,ikt),"(3F4.1)") abcde(3:5)
               ELSE
                  STOP "Error!" 
               ENDIF
                  chartTraces(j,iir,ikt)=cmplx(abcde(1)*dcos(PI*abcde(2)) &
                                              ,abcde(1)*dsin(PI*abcde(2)),dp)
            ELSEIF(itmp==0) THEN
            ELSE
             STOP "Error!!" 
            ENDIF
        END DO irk
!
        nelelittle (    ikt) = nele
        Irrepsname (iir,ikt) = irtmp
        Herringrule(iir,ikt) = jtmp
      END DO
     !print*,"tnir=",tnir
!
              nirreps(ikt)=iir
             !call dumptableofIrs(ikt,9)
              call Kreal2string(samplek(:,ikt),ckpoint) 
              ttmp(1,1:3)=samplek(1:3,ikt)
             !WRITE(99,"(1X,I2,A21,3F6.3,A5,3F8.3,I4,A4)") ikt, samplekname(ikt)//' ('//ckpoint//')' & 
             !  , samplek(:,ikt),' --> ',matmul(ttmp(:,:),Kc2p(:,:)),ikt,samplekname(ikt)
      IF(ikt/=Numk) STOP"ERROR in little groups of k-points"
      nktp=ikt
!
      close(11)
      close(9)
      close(99)
!
!.....output
!     WRITE(6,*) sgn
!
 600  FORMAT(/,10X,A30)
 601  FORMAT(/,3X,'i=',I2,3X,3I2,F8.3)
 602  FORMAT(10X,3I2,F8.3,2X,'(',2F6.3,')(',2F6.3,')')
      END SUBROUTINE MPNTGRP

      SUBROUTINE getkid(FL3,tkk,ik,GRPNAM,IIZ,NMAT,KV,iR)
      LOGICAL ,INTENT(IN) :: FL3
      real(dp)   , intent(inout) :: tkk(3)
      integer    , intent(out)   :: ik
      CHARACTER*3, intent(out)   :: GRPNAM
      integer    , intent(in)    :: IIZ(3,3,NSYM)
      INTEGER    , INTENT(IN)    :: NMAT
      REAL(DP)   , INTENT(INOUT) :: KV(3,NMAT)
      integer    , intent(OUT)   :: iR

      integer  :: ikt
      integer  :: nkt(NSYM+1),nvar(NSYM+1)
      integer  :: IIZ_(3,3,DoubNum/2+1)

      integer  :: j,ivar,iw
      real(dp) :: kkk(3),tkkc(3)

      integer     :: iir
      real(dp) :: AGW
      complex(dp) :: PHW

       ik=0
       kkk(:)=0._dp
       nkt(:)=0;nvar(:)=4
      !call   getkid2(tkk,ik,ivar)
       IIZ_(:,:,:)=0; IIZ_(:,:,1:DoubNum/2)=IIZ(:,:,1:DoubNum/2)
       IIZ_(1,1,1+DoubNum/2)=-1;IIZ_(2,2,1+DoubNum/2)=-1;IIZ_(3,3,1+DoubNum/2)=-1
       Do iR=1,DoubNum/2+1
          if(iR==DoubNum/2+1 .and. FL3) Cycle
         !kkk(1)=dot_product(tkk(:),IIZ(:,1,iR))
         !kkk(2)=dot_product(tkk(:),IIZ(:,2,iR))
         !kkk(3)=dot_product(tkk(:),IIZ(:,3,iR))
          kkk(:)=matmul(tkk(:),IIZ_(:,:,iR))
          if (kkk(1) < 0.0 .and. abs(kkk(1)) > 1e-6) kkk(1) = kkk(1) + 1.d0
          if (kkk(2) < 0.0 .and. abs(kkk(2)) > 1e-6) kkk(2) = kkk(2) + 1.d0
          if (kkk(3) < 0.0 .and. abs(kkk(3)) > 1e-6) kkk(3) = kkk(3) + 1.d0
          call   getkid2(kkk,ikt,ivar)
         !write(*,*) iR,ikt,ivar
          nkt(iR)=ikt;nvar(iR)=ivar
          IF(ivar==0) EXIT
       ENDDO

       IF(iR==DoubNum/2+2) THEN
       zj:Do iw=1,3 
          Do iR=1,DoubNum/2+1
             ivar=nvar(iR)
             IF(ivar==iw) EXIT zj
          ENDDO
          ENDDO zj
       ENDIF
       ikt=nkt(iR);ivar=nvar(iR)
       !write(*,*) iR,ikt,ivar;stop
       ik=ikt
      IF(ik>nktp) STOP " Nonsymmorphic kpoint is NOT found."
      GRPNAM=samplekname(ik)//' '

!---add by zjwang on 10.18.2019
     !kkk(1)=dot_product(tkk(:),IIZ(:,1,iR))
     !kkk(2)=dot_product(tkk(:),IIZ(:,2,iR))
     !kkk(3)=dot_product(tkk(:),IIZ(:,3,iR))
      kkk(:)=matmul(tkk(:),IIZ_(:,:,iR))
      tkk(:)=kkk(:)
      DO j=1,NMAT
         KKK(:)=KV(:,j)
         KV(:,j)=matmul(KKK(:),IIZ_(:,:,iR))
      ENDDO
     !AGW=PI*dot_product(tkk(:),TAU(:,iR))
     !PHW=CMPLX(DCOS(AGW),DSIN(AGW),DP)
     !DO j=1,NUME
     !   A(:,j)=A(:,j)*PHW
     !   B(:,j)=B(:,j)*PHW
     !ENDDO
!---add by zjwang on 10.18.2019


      tkkc(:)= matmul(tkk(:), p2cR(:,:))
      do iir=1,nirreps(ik)
         do j=1,Doubnum

            IF(labels(1,j,iir,ik)==1) THEN
             IF(labels(2,j,iir,ik)==1) THEN
                tableTraces(j,iir,ik)=chartTraces(j,iir,ik)
             ELSEIF(labels(2,j,iir,ik)==2) THEN
                AGW=0._DP
               !AGW=PI*DOT_PRODUCT(coeff_uvw(1:3,j,iir,ik),uvw(:))
                AGW=PI*DOT_PRODUCT(coeff_uvw(1:3,j,iir,ik),tkkc(:))
                PHW=CMPLX(DCOS(AGW),DSIN(AGW),DP)
                tableTraces(j,iir,ik)=chartTraces(j,iir,ik)*PHW
             ELSE
                STOP "Error" 
             ENDIF
            ENDIF

         enddo
      enddo
      IF(iR==DoubNum/2+1) iR=0
      return
      END SUBROUTINE getkid

      subroutine getkid2(tkk,ik,ivar)
      real(dp)   , intent(in ) :: tkk(3)
      integer    , intent(out) :: ik
      integer    , intent(out) :: ivar

      real(dp)    :: uvw(3)
      integer     :: ikt
      character*15 , allocatable  :: ckpoint(:)
      logical :: is_variable(3)
      integer :: j
      integer :: num_var_tmp, ref_varnum
      integer, allocatable :: num_var(:)
      integer, allocatable :: ind_u(:), ind_v(:), ind_w(:)
      
      real(dp) :: uvw_tmp(3)
      real(dp) :: tkkc(3), refkpoint(3), refkpointp(3)
      real(dp) :: diff
      real(dp), parameter :: tol=1e-5

      character*15 :: chr_tmp
      ivar=-1

      allocate(ckpoint(nktp))
      do ikt=1,nktp
         call Kreal2string(samplek(:,ikt),ckpoint(ikt)) 
      enddo

      allocate(num_var(nktp))
      num_var = 0
      allocate(ind_u(nktp))
      allocate(ind_v(nktp))
      allocate(ind_w(nktp))
      ind_u=0; ind_v=0; ind_w=0

     !print*,"================================================"
     !   write(*,'(20X,A5,10X,A6,20X,A5)') 'conv', 'symbol', 'prim'
     !do ikt=1,nktp
     !   write(*,'(I5,3F9.3,2X,A15,3F12.6)') ikt,samplek(:,ikt),ckpoint(ikt),matmul(samplek(:,ikt),Kc2p(:,:))
     !enddo
     !   write(*,*) 'Caution!!! conv  : 0.333 -> 0.3333333333'
     !   write(*,*) 'Caution!!! symbol: 0.33  -> 0.3333333333'
     ! !--- add something here

      !--- add something below-------
      ! get the number of variables for the reference kpoints
      do ikt = 1, nktp
         num_var_tmp = 0
         ind_u(ikt)   = index(ckpoint(ikt), "  u  ") 
         ind_v(ikt)   = index(ckpoint(ikt), "  v  ") 
         ind_w(ikt)   = index(ckpoint(ikt), "  w  ") 
         if (ind_u(ikt).ne.0) num_var_tmp = num_var_tmp+1
         if (ind_v(ikt).ne.0) num_var_tmp = num_var_tmp+1
         if (ind_w(ikt).ne.0) num_var_tmp = num_var_tmp+1
         num_var(ikt) = num_var_tmp 
      enddo
     !write(*,*) "Number of variables for the reference kpoints"
     !write(*,"(15I3)") num_var

      ! convert the input kpoint to conventional basis
      tkkc = matmul(tkk(:), p2cR(:,:))
     !write(*,"(A35,5X,3F12.6)") "Input kpoint under primitive :", tkk(:)
     !write(*,"(A35,5X,3F12.6)") "Input kpoint under convention:", tkkc(:)

      ! compare the input kpoint with the reference kpoints
      ! the reference kpoints with less variables have a higher priority
      uvw = 0d0 
      ref_varnum = 9999
      do ikt = 1, nktp
         is_variable = .false.
         refkpoint = 9999d0
         uvw_tmp = 0d0 
         if(ckpoint(ikt)(1:5)=='  u  ') then 
            is_variable(1) = .true.
            uvw_tmp(1) = tkkc(1)
            refkpoint(1) = tkkc(1)
            if(ckpoint(ikt)(6:10)=='  u  ') refkpoint(2) = refkpoint(1)
            if(ckpoint(ikt)(6:10)==' 1-u ') refkpoint(2) = 1d0-refkpoint(1)
            if(ckpoint(ikt)(6:10)=='-2u  ') refkpoint(2) = -2d0*refkpoint(1)
            if(ckpoint(ikt)(6:10)==' -u  ') refkpoint(2) = -refkpoint(1)
            if(ckpoint(ikt)(6:10)==' 1+u ') refkpoint(2) = 1d0+refkpoint(1)
            if(ckpoint(ikt)(11:15)=='  u  ') refkpoint(3) = refkpoint(1)
            if(ckpoint(ikt)(11:15)==' 1-u ') refkpoint(3) = 1d0-refkpoint(1)
            if(ckpoint(ikt)(11:15)=='-2u  ') refkpoint(3) = -2d0*refkpoint(1)
            if(ckpoint(ikt)(11:15)==' -u  ') refkpoint(3) = -refkpoint(1)
            if(ckpoint(ikt)(11:15)==' 1+u ') refkpoint(3) = 1d0+refkpoint(1)
         endif
         if(ckpoint(ikt)(6:10)=='  v  ') then 
            is_variable(2) = .true.
            uvw_tmp(2) = tkkc(2)
            refkpoint(2) = tkkc(2)
            if (ckpoint(ikt)(1:5)==' 1-v ') refkpoint(1) = 1d0-refkpoint(2)
            if (ckpoint(ikt)(11:15)==' 1-v ') refkpoint(3) = 1d0-refkpoint(2)
            if (ckpoint(ikt)(11:15)=='  v  ') refkpoint(3) = refkpoint(2)
         endif
         if(ckpoint(ikt)(11:15)=='  w  ') then 
             is_variable(3) = .true.
             uvw_tmp(3) = tkkc(3)
             refkpoint(3) = tkkc(3)
         endif

         ! some exceptions
         ! 0.5 u 0.0  : 195 198 200 201 205
         ! 1+u 1-u 0.0: 197 199 204 206 211 214 217 220 229 230
         if(ckpoint(ikt)(1:5).ne.'  u  '.and.ckpoint(ikt)(6:10).eq.'  u  ') then 
            is_variable(1) = .true.
            uvw_tmp(1) = tkkc(2)
            refkpoint(2) = tkkc(2)
         endif
         if(ckpoint(ikt)(1:5).eq.' 1+u '.and.ckpoint(ikt)(6:10).eq.' 1-u ') then 
            is_variable(1) = .true.
            uvw_tmp(1) = tkkc(1) - 1
            refkpoint(1) = tkkc(1)
            refkpoint(2) = 2-refkpoint(1)
         endif
         ! end exceptions
         
         do j = 1,3
            if (abs(refkpoint(j)-9999d0) < tol) then 
               chr_tmp = adjustr(ckpoint(ikt)((j-1)*5+1:j*5))
               read(chr_tmp, *) refkpoint(j)
               if (abs(refkpoint(j)-0.33) < tol) refkpoint(j) = 1d0/3d0
            endif
         enddo

         refkpointp = matmul(refkpoint,kc2p)

         !do j = 1, 3
         !   if (refkpointp(j) < 0.0) refkpointp(j) = refkpointp(j) + 1.d0
         !enddo 

         diff = abs(mod(refkpointp(1),1d0)-mod(tkk(1),1d0)) + &
                abs(mod(refkpointp(2),1d0)-mod(tkk(2),1d0)) + &
                abs(mod(refkpointp(3),1d0)-mod(tkk(3),1d0))
         diff = abs((refkpointp(1)-floor(refkpointp(1)))-(tkk(1)-floor(tkk(1)))) + &
                abs((refkpointp(2)-floor(refkpointp(2)))-(tkk(2)-floor(tkk(2)))) + &
                abs((refkpointp(3)-floor(refkpointp(3)))-(tkk(3)-floor(tkk(3))))    
         if (diff < tol) then 
             if (num_var(ikt) < ref_varnum) then 
                 ref_varnum = num_var(ikt)
                 ik = ikt 
                 uvw = 0d0
                 do j = 1, 3
                    if(is_variable(j)) uvw(j) = uvw_tmp(j)
                 enddo
             endif
         endif

      enddo
      ivar= num_var(ik)
     !write(*,*) "The id of reference kpoint for the input kpoint is:", ik,num_var(ik)
     !write(*,"(A10,3F12.6)") "uvw(:)=",uvw(:)
     !print*,"================================================"
      deallocate(ckpoint)
      deallocate(num_var)

      deallocate(ind_u)
      deallocate(ind_v)
      deallocate(ind_w)
      end subroutine getkid2


      SUBROUTINE dumptableofIrs(ikt,iout)
      integer, intent(in) :: ikt
      integer, intent(in) :: iout

      character*15  :: ckpoint
      logical       :: ftg
      integer       :: i,j
      
      call Kreal2string(samplek(:,ikt),ckpoint) 
      chkpoint(ikt)=ckpoint
      WRITE(iout,"(/,/,3X,I2,A15,A54,$)") ikt, samplekname(ikt)//' : kname '&
         ,ckpoint//' :  given in the conventional basis'
      WRITE(iout,'(/,3X,I2,A,$)') antisym(ikt) &
          ,' : the existence of antiunitary symmetries. 1-exist; 0-no'
      WRITE(iout,577) 'Reality' 
     !DO I=1,nelelittle(ikt)/2
     !WRITE(iout,579) I 
     !ENDDO
      DO I=1,Doubnum/2
         IF(labels(1,I,1,ikt)==1) WRITE(iout,579) I 
      ENDDO

      do j=1,nirreps(ikt)
         ftg=.false.
         IF(Irrepsname(j,ikt)(1:1)=='-') THEN
            WRITE(iout,576) Herringrule(j,ikt),Irrepsname(j,ikt)(2:5)//' '
         ELSE
            WRITE(iout,576) Herringrule(j,ikt),Irrepsname(j,ikt)
         ENDIF
         DO I=1,Doubnum/2
        !IF(labels(1,I,j,ikt)==1) WRITE(iout,580) tableTraces(I,j,ikt)
         IF(labels(1,I,j,ikt)==1) WRITE(iout,580) chartTraces(I,j,ikt)+cmplx(epsil,epsil,dp)
         IF(labels(1,I,j,ikt)==1.and. labels(2,I,j,ikt)==2 ) ftg=.true.
         ENDDO

         IF(ftg) THEN
                 WRITE(iout,578) '     '
                 DO I=1,Doubnum/2
                    IF(labels(1,I,j,ikt)==1) THEN
                       IF(labels(2,I,j,ikt)==2) THEN
                          WRITE(iout,582) factsTraces(I,j,ikt)
                       ELSE
                          WRITE(iout,582) '          '
                       ENDIF
                    ENDIF
                ENDDO
         ENDIF


         IF(Irrepsname(j+1,ikt)(1:1)=='-'.and.Irrepsname(j,ikt)(1:1)/='-') THEN
            sirreps(ikt)=j 
           !WRITE(iout,'(I3)') j
            WRITE(iout,581) '----'
            do i=1,nelelittle(ikt)/2
            WRITE(iout,582) '------------'
            enddo
         ENDIF
      enddo
 510  FORMAT(/,7X,'WILL BE IMPLEMENTED')

 576  FORMAT(/,3X,I2,2X,A6,$)
 577  FORMAT(/,1X,A7,$)
 578  FORMAT(/,7X,A6,$)
 579  FORMAT(10X,I2,$)
 580  FORMAT(F5.2,SP,F5.2,'i ',$)
 581  FORMAT(/,8X,A4,$)
 582  FORMAT(A12,$)
      RETURN
      END SUBROUTINE dumptableofIrs

      SUBROUTINE Kreal2string(coorkp,kstr) 
      real(dp), intent(in)::  coorkp(3) ! Coordinates 
      character*15,intent(out):: kstr

!real(dp),    parameter  :: su=0.123_dp,s1_u=0.877_dp,s_2u=-0.246_dp
!real(dp),    parameter  :: sv=0.313_dp,s1_v=0.687_dp
!real(dp),    parameter  :: sw=0.427_dp,s_u=-0.123_dp,s1u = 1.123_dp
      ! u     ->  0.123  : su
      ! v     ->  0.313  : sv
      ! w     ->  0.427  : sw
      ! 1 - u ->  0.877  : s1_u
      ! -2 u  -> -0.246  : s_2u
      ! -u    -> -0.123  : s_u
      ! 1 + u ->  1.123  : s1u
      ! 1 - v ->  0.687  : s1_v
      !------------codes for u,v,w

      if     ( abs( coorkp(1)-su   )<epsil) then;kstr(1:5)="  u  "
      elseif ( abs( coorkp(1)-sv   )<epsil) then;kstr(1:5)="  v  "
      elseif ( abs( coorkp(1)-sw   )<epsil) then;kstr(1:5)="  w  "
      elseif ( abs( coorkp(1)-s1_u )<epsil) then;kstr(1:5)=" 1-u "
      elseif ( abs( coorkp(1)-s_2u )<epsil) then;kstr(1:5)="-2u  "
      elseif ( abs( coorkp(1)-s_u  )<epsil) then;kstr(1:5)=" -u  "
      elseif ( abs( coorkp(1)-s1u  )<epsil) then;kstr(1:5)=" 1+u "
      elseif ( abs( coorkp(1)-s1_v )<epsil) then;kstr(1:5)=" 1-v "
      else
       write(kstr(1:5),501) coorkp(1)
      endif

      if     ( abs( coorkp(2)-su   )<epsil) then;kstr(6:10)="  u  "
      elseif ( abs( coorkp(2)-sv   )<epsil) then;kstr(6:10)="  v  "
      elseif ( abs( coorkp(2)-sw   )<epsil) then;kstr(6:10)="  w  "
      elseif ( abs( coorkp(2)-s1_u )<epsil) then;kstr(6:10)=" 1-u "
      elseif ( abs( coorkp(2)-s_2u )<epsil) then;kstr(6:10)="-2u  "
      elseif ( abs( coorkp(2)-s_u  )<epsil) then;kstr(6:10)=" -u  "
      elseif ( abs( coorkp(2)-s1u  )<epsil) then;kstr(6:10)=" 1+u "
      elseif ( abs( coorkp(2)-s1_v )<epsil) then;kstr(6:10)=" 1-v "
      else
       write(kstr(6:10),501) coorkp(2)
      endif

      if     ( abs( coorkp(3)-su   )<epsil) then;kstr(11:15)="  u  "
      elseif ( abs( coorkp(3)-sv   )<epsil) then;kstr(11:15)="  v  "
      elseif ( abs( coorkp(3)-sw   )<epsil) then;kstr(11:15)="  w  "
      elseif ( abs( coorkp(3)-s1_u )<epsil) then;kstr(11:15)=" 1-u "
      elseif ( abs( coorkp(3)-s_2u )<epsil) then;kstr(11:15)="-2u  "
      elseif ( abs( coorkp(3)-s_u  )<epsil) then;kstr(11:15)=" -u  "
      elseif ( abs( coorkp(3)-s1u  )<epsil) then;kstr(11:15)=" 1+u "
      elseif ( abs( coorkp(3)-s1_v )<epsil) then;kstr(11:15)=" 1-v "
      else
       write(kstr(11:15),501) coorkp(3)
      endif

 501  FORMAT(F5.2)
      END SUBROUTINE Kreal2string

      SUBROUTINE wrtirb(NUME,FL,FGT,EE,XM,NE,LKG,IKG,IK,KPH1,nmin,iR)
      LOGICAL ,INTENT(IN) :: FL(FLMAX),FGT
      INTEGER ,INTENT(IN) :: NUME,IKG,IK
      INTEGER ,INTENT(IN) :: NE,LKG(NSYM),nmin,iR
      REAL(dp),INTENT(IN) :: EE(NUME)
      COMPLEX(dp),INTENT(IN)::XM(NSYM,NUME),KPH1(NSYM)

      COMPLEX(dp)::KPH(NSYM)
      real(dp),    parameter  :: TOL  = 0.05_DP
      LOGICAL  :: FG1
      REAL(dp) :: DSUM
      COMPLEX(dp) :: X,XR(MAXIR),ZIR
      INTEGER :: IE,ND,IG,I1,J1,I2,J2,IS
      INTEGER :: NIR,NIR1,NIR2
      INTEGER :: LCNB
      character*5   ::   CHND
     
!**********************************************************************
!
!.....calc. the characters for all energies and tranformations
     !IE=1
      IE=nmin
       WRITE(66,"(/,A4,$)") samplekname(IK)//': '
       IF(.not. IR) WRITE(67,"(/,I3,$)") -IK
       IF( IR) WRITE(67,"(/,I3,$)") IK
       IF(.not. IR) WRITE(68,"(/,I3,$)") -IK
       IF( IR) WRITE(68,"(/,I3,$)") IK
       LCNB=0
      DO WHILE(IE.LE.NE)
      ND=1
      DO WHILE((IE+ND).LE.NE)
      IF((EE(IE+ND)-EE(IE)).LT.TOLDG) THEN
        ND=ND+1
      ELSE
        EXIT 
      ENDIF
      END DO
      WRITE(6,570) IE,ND,EE(IE)

      IF(ND.GT.MAXDG) GOTO 100

      IF(ND<10) THEN; write(CHND,'(A1,I1,A3)') '(',ND,'); '
      ELSE          ; write(CHND,'(A1,I2,A2)') '(',ND,');'
      ENDIF
      LCNB=LCNB+ND
!.....if 'Cornwell condition' is NOT satisfied
!
      IF(.NOT.FGT) THEN
      WRITE(624,'(I3,I3,F12.6,$)') IE,ND,EE(IE)
      DO IG=1,IKG
         X=XM(LKG(IG),IE)
         WRITE(6,580) DREAL(X)+epsil,DIMAG(X)+epsil
         WRITE(624,'(2F12.6,$)') DREAL(X*KPH1(LKG(IG)))+epsil &
                                ,DIMAG(X*KPH1(LKG(IG)))+epsil
         XR(IG)=X
      ENDDO
         WRITE(624,*)

!.....for all classes
      IF(.NOT.FL(2)) THEN
        IS=1
        NIR1=1
        NIR2=sirreps(IK)
      ELSE
        IS=2
        NIR1=sirreps(IK)+1
        NIR2=nirreps(IK)
      ENDIF
      NIR=NIR2-NIR1+1
!
      XR(:)=conjg(XR(:));KPH(:)=conjg(KPH1(:))  ! added by zjwang at Sep29 2018
!--------for one irrep
      DSUM=0._DP
         DO I1=NIR1,NIR2 
      FG1=.TRUE.
      DO IG=1,IKG
         DSUM=ABS(XR(IG)*KPH(LKG(IG))-tableTraces(LKG(IG),I1,IK))
        !PRINT*,XR(IG),I1,IG,DSUM
        !PRINT*,KPH(LKG(IG)),tableTraces(LKG(IG),I1,IK)
         IF(DSUM.GT.(TOL*DBLE(NIR))) THEN
           FG1=.FALSE.;EXIT
         ENDIF
      ENDDO
         IF(FG1) EXIT
         ENDDO
      IF(FG1) THEN
       WRITE(6,"('=',A4,$)") Irrepsname(I1,IK)(IS:IS+3) 
       WRITE(66,"(A4,A5,$)") Irrepsname(I1,IK)(IS:IS+3) ,CHND
       WRITE(67,"(I3,$)") I1
       WRITE(68,"(I3,$)") I1
       GOTO 200
      ENDIF
!
!--------for two irreps
      LP:DO I1=NIR1,NIR2 
         DO J1=I1,NIR2 
      FG1=.TRUE.
      DO IG=1,IKG
         DSUM=ABS(XR(IG)*KPH(LKG(IG))-tableTraces(LKG(IG),I1,IK)-tableTraces(LKG(IG),J1,IK))
         IF(DSUM.GT.(TOL*DBLE(NIR))) THEN
           FG1=.FALSE.;EXIT
         ENDIF
      ENDDO
         IF(FG1) EXIT LP
         ENDDO
         ENDDO LP
      IF(FG1) THEN
       IF(Irrepsname(I1,IK)(IS+3:IS+3)==' ') THEN
       WRITE(6,"('=',A4,'+ ',A4,$)") Irrepsname(I1,IK)(IS:IS+3),Irrepsname(J1,IK)(IS:IS+3) 
       ELSE
       WRITE(6,"('=',A4,' + ',A4,$)") Irrepsname(I1,IK)(IS:IS+3),Irrepsname(J1,IK)(IS:IS+3) 
       ENDIF
       WRITE(66,"(2A4,A5,$)") Irrepsname(I1,IK)(IS:IS+3),Irrepsname(J1,IK)(IS:IS+3),CHND 
       WRITE(67,"(2I3,$)") I1,J1
       WRITE(68,"(2I3,$)") I1,J1
       GOTO 200
      ENDIF

!--------for three irreps
      LP3:DO I1=NIR1,NIR2 
         DO J1=I1,NIR2 
          DO I2=J1,NIR2 
      FG1=.TRUE.
      DO IG=1,IKG
         DSUM=ABS(XR(IG)*KPH(LKG(IG))-tableTraces(LKG(IG),I1,IK)-tableTraces(LKG(IG),J1,IK) &
                                     -tableTraces(LKG(IG),I2,IK) )
         IF(DSUM.GT.(TOL*DBLE(NIR))) THEN
           FG1=.FALSE.;EXIT
         ENDIF
      ENDDO
         IF(FG1) EXIT LP3
         ENDDO
         ENDDO
         ENDDO LP3
      IF(FG1) THEN
       IF(Irrepsname(I1,IK)(IS+3:IS+3)==' ') THEN
       WRITE(6,"('=',A4,'+ ',A4,'+ ',A4,$)") Irrepsname(I1,IK)(IS:IS+3),Irrepsname(J1,IK)(IS:IS+3) &
                    ,Irrepsname(I2,IK)(IS:IS+3)
       ELSE
       WRITE(6,"('=',A4,' + ',A4,' + ',A4,$)") Irrepsname(I1,IK)(IS:IS+3),Irrepsname(J1,IK)(IS:IS+3) &
                    ,Irrepsname(I2,IK)(IS:IS+3)
       ENDIF
       WRITE(66,"(3A4,A5,$)") Irrepsname(I1,IK)(IS:IS+3),Irrepsname(J1,IK)(IS:IS+3) &
                    ,Irrepsname(I2,IK)(IS:IS+3),CHND
       WRITE(67,"(3I3,$)") I1,J1,I2
       WRITE(68,"(3I3,$)") I1,J1,I2
       GOTO 200
      ENDIF

!--------for four irreps
      LP4:DO I1=NIR1,NIR2 
         DO J1=I1,NIR2 
          DO I2=J1,NIR2 
         DO J2=I2,NIR2 
      FG1=.TRUE.
      DO IG=1,IKG
         DSUM=ABS(XR(IG)*KPH(LKG(IG))-tableTraces(LKG(IG),I1,IK)-tableTraces(LKG(IG),J1,IK) &
                                     -tableTraces(LKG(IG),I2,IK)-tableTraces(LKG(IG),J2,IK))
         IF(DSUM.GT.(TOL*DBLE(NIR))) THEN
           FG1=.FALSE.;EXIT
         ENDIF
      ENDDO
         IF(FG1) EXIT LP4
         ENDDO
         ENDDO
         ENDDO
         ENDDO LP4
      IF(FG1) THEN
       IF(Irrepsname(I1,IK)(IS+3:IS+3)==' ') THEN
       WRITE(6,"('=',A4,'+ ',A4,'+ ',A4,'+ ',A4,$)") Irrepsname(I1,IK)(IS:IS+3),Irrepsname(J1,IK)(IS:IS+3) &
                    ,Irrepsname(I2,IK)(IS:IS+3),Irrepsname(J2,IK)(IS:IS+3) 
       ELSE
       WRITE(6,"('=',A4,' + ',A4,' + ',A4,' + ',A4,$)") Irrepsname(I1,IK)(IS:IS+3),Irrepsname(J1,IK)(IS:IS+3) &
                    ,Irrepsname(I2,IK)(IS:IS+3),Irrepsname(J2,IK)(IS:IS+3) 
       ENDIF
       WRITE(66,"(4A4,A5,$)") Irrepsname(I1,IK)(IS:IS+3),Irrepsname(J1,IK)(IS:IS+3) &
                          ,Irrepsname(I2,IK)(IS:IS+3),Irrepsname(J2,IK)(IS:IS+3),CHND
       WRITE(67,"(4I3,$)") I1,J1,I2,J2
       WRITE(68,"(4I3,$)") I1,J1,I2,J2
       GOTO 200
      ENDIF

!--------can NOT find irreps
      WRITE(6,"(A2,$)") '??'
      WRITE(66,"(A3,A5,$)") '?? ',CHND
      WRITE(67,"(A3,$)") ' ??'
      WRITE(68,"(A3,$)") ' ??'
!     STOP "ERROR in irvsp2"
!
!
 200   WRITE(6,*)
      ENDIF

 100  IE=IE+ND
      END DO

      IF(LCNB<10) THEN; write(CHND,'(A1,I1,A3)') '[',LCNB,']  '
      ELSEIF(LCNB<100) THEN; write(CHND,'(A1,I2,A2)') '[',LCNB,'] '
      ELSE            ; write(CHND,'(A1,I3,A1)') '[',LCNB,']'
      ENDIF
       WRITE(66,"(A5,$)") CHND

      RETURN
 570  FORMAT(I3,I3,F10.6,$)
 580  FORMAT(F5.2,SP,F5.2,'i ',$)
      END 
!
end module nonsymm
