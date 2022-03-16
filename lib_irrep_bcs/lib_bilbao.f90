module bilbao 

    use lib_comms
    implicit none 

    real(dp),   public,  save :: Kc2p(3,3), p2cR(3,3)

contains

subroutine bilbao_read(sgn)

    !!! This subroutine read symmetry operations and 
    !!! character tables from Bilbao.
    !!! Compare the input symmetry operations with Bilbao. 

    integer,          intent(in) :: sgn

    character(len=3)   :: csgn
    character(len=180) :: spgpath, spgfile
    character(len=10)  :: symbol_sg

    real(dp)           :: Df(2,4), abcde(5), ktmp(3), ttmp(1,3)
    real(dp)           :: tmp33(3,3) 

    integer            :: Numk, antiss, tnir, wi 

    integer            :: i, j, itmp, jtmp 
    integer            :: iir, nele, ikt 
    character(len=5)   :: irtmp 
    character(len=2)   :: nametmp, nametmp2 
    character(len=15)  :: ckpoint 
    character(len=40)  :: ListIrrep

    ! bilbao table file
    integer,     parameter   :: bb = 11
    ! output file : operations in conventional basis
    integer,     parameter   :: op = 9
    ! output file : special k points
    integer,     parameter   :: sk = 00


    if     (sgn < 10)  then; write(csgn,'(I1)') sgn
    elseif (sgn < 100) then; write(csgn,'(I2)') sgn
    else                   ; write(csgn,'(I3)') sgn
    endif 

    !spgpath = '/storagehome/jcgao/soft/irvsp/src_irvsp4/src_lib'
#ifdef IRVSPDATA 
    call get_environment_variable('IRVSPDATA',spgpath)
#else
    write(6,*) "Environment variable 'IRVSPDATA' must be provided "
    write(6,*) "Please run the following commands to make the library:"
    write(6,*) "./configure.sh"
    write(6,*) "source ~/.bashrc"
    write(6,*) "make lib"
    stop
#endif 
    spgfile = trim(spgpath)//'/kLittleGroups/kLG_'//trim(csgn)//'.data'
    write(*,*) "SPGFILE :", trim(adjustl(spgfile)) 

    open(unit=bb, file=spgfile, status='old', form='unformatted')
    open(unit=sk, file='SGklist_'//trim(csgn)//'.cht', status='unknown')

    read(bb) num_doub_sym, symbol_sg

    if     (symbol_sg(1:1) == 'P') then
        Kc2p = Pabc
    elseif (symbol_sg(1:1) == 'C') then 
        Kc2p = Cabc
        if (sgn==68) Kc2p = Cabc68
    elseif (symbol_sg(1:1) == 'B') then
        Kc2p = Babc
    elseif (symbol_sg(1:1) == 'A') then
        Kc2p = Aabc
    elseif (symbol_sg(1:1) == 'R') then
        Kc2p = Rabc
    elseif (symbol_sg(1:1) == 'F') then
        Kc2p = Fabc
    elseif (symbol_sg(1:1) == 'I') then
        Kc2p = Iabc
    else
        stop "Error in space-group-symbol"
    endif 
    call invreal33(Kc2p, p2cR)

    rot_bilbao = 0
    tau_bilbao = 0.d0
    SU2_bilbao = 0.d0
    do i = 1, num_doub_sym 
        read(bb) rot_bilbao(:,:,i), tau_bilbao(:,i), Df(:,:)
        SU2_bilbao(1,1,i)=cmplx(Df(1,1)*dcos(PI*Df(2,1)),Df(1,1)*dsin(PI*Df(2,1)),dp)
        SU2_bilbao(1,2,i)=cmplx(Df(1,2)*dcos(PI*Df(2,2)),Df(1,2)*dsin(PI*Df(2,2)),dp)
        SU2_bilbao(2,1,i)=cmplx(Df(1,3)*dcos(PI*Df(2,3)),Df(1,3)*dsin(PI*Df(2,3)),dp)
        SU2_bilbao(2,2,i)=cmplx(Df(1,4)*dcos(PI*Df(2,4)),Df(1,4)*dsin(PI*Df(2,4)),dp)
    enddo 

    ! write operations under conventional basis
    open(unit=op, file='SGoperation_'//trim(csgn)//'.cht', status='unknown')
    write(op,*) symbol_sg 
    write(op, "(' From conv. to prim. reciprocal space  (DB1)')")
    write(op, "(3(3F16.8,/))") Kc2p 
    write(op, "(' From prim. to conv. reciprocal space  (DR1)')")
    write(op, "(3(3F16.8,/))") p2cR
    do i = 1, num_doub_sym 
        write(op,601) i, rot_bilbao(:,1,i), tau_bilbao(1,i)
        write(op,602)    rot_bilbao(:,2,i), tau_bilbao(2,i), SU2_bilbao(1,1,i), SU2_bilbao(1,2,i)
        write(op,602)    rot_bilbao(:,3,i), tau_bilbao(3,i), SU2_bilbao(2,1,i), SU2_bilbao(2,2,i)
    enddo 
    close(op)

    ! conventional cell -> primitive cell
    invrot_bilbao = 0
    do i = 1, num_doub_sym 
        tmp33(:,:) = dble(transpose(rot_bilbao(:,:,i)))
        rot_bilbao(:,:,i) = nint(matmul(matmul(p2cR, tmp33), Kc2p(:,:))) 
        tau_bilbao(:,i) = matmul(p2cR, tau_bilbao(:,i))
        call invmati(rot_bilbao(:,:,i), invrot_bilbao(:,:,i))
    enddo 


    ! read character tables
    nirreps(:)=0
    sirreps(:)=0
    antisym(:)=0
    Herringrule(:,:)=-2

    labels=0
    iir=0;ikt=0
    tableTraces(:,:,:)=cmplx_0 
    chartTraces(:,:,:)=cmplx_0 
    coeff_uvw(:,:,:,:)=cmplx_0
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
            WRITE(99,"(1X,I2,A21,3F6.3,A5,3F8.3,I4,A4)") ikt, samplekname(ikt)//' ('//ckpoint//')' & 
              , samplek(:,ikt),' --> ',matmul(ttmp(:,:),Kc2p(:,:)),ikt,samplekname(ikt)
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
      irk:DO j=1,num_doub_sym
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
            WRITE(99,"(1X,I2,A21,3F6.3,A5,3F8.3,I4,A4)") ikt, samplekname(ikt)//' ('//ckpoint//')' & 
              , samplek(:,ikt),' --> ',matmul(ttmp(:,:),Kc2p(:,:)),ikt,samplekname(ikt)
    IF(ikt/=Numk) STOP"ERROR in little groups of k-points"
    num_ktype=ikt


 600  format(/,10X,A30)
 601  format(/,3X,'i=',I2,3X,3I2,F8.3)
 602  format(10X,3I2,F8.3,2X,'(',2F6.3,')(',2F6.3,')')

end subroutine bilbao_read 


subroutine bilbao_getconjg(ind_ope, &
                           num_littg_input, littg_input, &
                           num_littg_bilbao, littg_bilbao, &
                           kopeconjg_bb2input)

    integer, intent(in)  :: ind_ope

    integer, intent(in)  :: num_littg_input
    integer, intent(in)  :: littg_input(MAXSYM)

    integer, intent(in)  :: num_littg_bilbao
    integer, intent(in)  :: littg_bilbao(MAXSYM)

    integer, intent(out) :: kopeconjg_bb2input(MAXSYM)

    integer :: iope, kope, jope, mope  
    integer :: tmprot(3,3)

    kopeconjg_bb2input = 0
    if (num_littg_input/=num_littg_bilbao/2) stop "little group error"

    do iope = 1, num_littg_bilbao/2 
        kope = littg_bilbao(iope)
        tmprot = matmul(invrot_bilbao(:,:,ind_ope), &
                 matmul(rot_bilbao(:,:,kope),rot_bilbao(:,:,ind_ope)))
        do jope = 1, num_littg_input 
            mope = littg_input(jope)
            if (abs(rot_bilbao(1,1,mope)-tmprot(1,1)) + & 
                abs(rot_bilbao(1,2,mope)-tmprot(1,2)) + & 
                abs(rot_bilbao(1,3,mope)-tmprot(1,3)) + & 
                abs(rot_bilbao(2,1,mope)-tmprot(2,1)) + & 
                abs(rot_bilbao(2,2,mope)-tmprot(2,2)) + & 
                abs(rot_bilbao(2,3,mope)-tmprot(2,3)) + & 
                abs(rot_bilbao(3,1,mope)-tmprot(3,1)) + & 
                abs(rot_bilbao(3,2,mope)-tmprot(3,2)) + & 
                abs(rot_bilbao(3,3,mope)-tmprot(3,3)) == 0) then
                
                kopeconjg_bb2input(mope) = kope 
                exit
            endif 
        enddo 
        !if (kopeconjg_bb2input(mope) == 0) stop "Conjg not found"

    enddo 

end subroutine bilbao_getconjg 


subroutine bilbao_reorder(num_sym, &
                          rot_input, tau_input, SU2_input, SO3_input, &
                          rot_reorder, tau_reorder, SU2_reorder, SO3_reorder)

    integer,     intent(in)  :: num_sym 
    integer,     intent(in)  :: rot_input(3,3,num_sym)
    real(dp),    intent(in)  :: tau_input(3,num_sym)
    complex(dp), intent(in)  :: SU2_input(2,2,num_sym)
    real(dp),    intent(in)  :: SO3_input(3,3,num_sym)

    integer,     intent(out) :: rot_reorder(3,3,MAXSYM)
    real(dp),    intent(out) :: tau_reorder(3,MAXSYM)
    complex(dp), intent(out) :: SU2_reorder(2,2,MAXSYM)
    real(dp),    intent(out) :: SO3_reorder(3,3,MAXSYM)

    integer     :: i, j
    real(dp)    :: diff(3), dtest

    ! reorder the input operations according to bilbao
    if (num_sym/=num_doub_sym/2) stop "Element Error"
    reorder = 0 
    do i = 1, num_sym 
        do j = 1, num_sym 
            if (      rot_input(1,1,j)==rot_bilbao(1,1,i) &
                .and. rot_input(1,2,j)==rot_bilbao(1,2,i) &
                .and. rot_input(1,3,j)==rot_bilbao(1,3,i) &
                .and. rot_input(2,1,j)==rot_bilbao(2,1,i) &
                .and. rot_input(2,2,j)==rot_bilbao(2,2,i) &
                .and. rot_input(2,3,j)==rot_bilbao(2,3,i) &
                .and. rot_input(3,1,j)==rot_bilbao(3,1,i) &
                .and. rot_input(3,2,j)==rot_bilbao(3,2,i) &
                .and. rot_input(3,3,j)==rot_bilbao(3,3,i) ) then 
                reorder(i) = j
                exit 
            endif 
        enddo 
        diff(:) = tau_bilbao(:,i) - tau_input(:,j)
        dtest = dabs(nint(diff(1))-diff(1)) &
               +dabs(nint(diff(2))-diff(2)) &
               +dabs(nint(diff(3))-diff(3))
        if (dtest > epsil) then
            write(*,'(I3,3F10.5)') i, tau_bilbao(:,i) 
            write(*,'(I3,3F10.5)') j, tau_input(:,j)
            write(*,*) "TAU Error !!!"
            write(*,'(A)') " 1. Please check the space group number;"
            write(*,'(A)') " 2. Please check the standard (default) settings of the space group: "
            write(*,'(A)') "    a. unique axis b (cell choice 1) for space groups within the monoclinic system."
            write(*,'(A)') "    b. obverse triple hexagonal unit cell for R space groups."
            write(*,'(A)') "    c. the origin choice two - inversion center at (0,0,0) - for the centrosymmetric space groups."
            write(*,'(A)') " ### or you can sent the POSCAR to zjwang11@hotmail.com for help."
            stop
        else
            ! in case tau_input = 0,0,-0.5  tau_bilbao = 0,0,0.5, there may be difference in the generic k point
            !tau_reorder(:,i) = tau_bilbao(:,i)
            difftau_inputbilb(:,i) = diff(:)
        endif 
    enddo   

    rot_reorder = 0
    tau_reorder = 0.d0
    SU2_reorder = 0.d0
    SO3_reorder = 0.d0
    do i = 1, num_sym
        j = reorder(i)
        rot_reorder(:,:,i) = rot_input(:,:,j)
        tau_reorder(:,i)   = tau_input(:,j)
        SU2_reorder(:,:,i) = SU2_input(:,:,j)
        SO3_reorder(:,:,i) = SO3_input(:,:,j)
    enddo 

end subroutine bilbao_reorder 



subroutine bilbao_getkid(k_input, invrot, &
                         ind_ktype, kname, ind_rot, &
                         timerev_k)

    real(dp)        , intent(inout) :: k_input(3)
    integer         , intent(in)    :: invrot(3,3,MAXSYM)
    integer         , intent(out)   :: ind_ktype
    character(len=3), intent(out)   :: kname
    integer         , intent(out)   :: ind_rot 
    logical         , intent(out)   :: timerev_k

    integer     :: ind_ktype_tmp 
    integer     :: nkt(num_doub_sym), nvar(num_doub_sym)

    integer     :: iR,j,ivar,iw
    real(dp)    :: k_trial(3), k_input_conv(3)

    integer     :: iir
    real(dp)    :: AGW
    complex(dp) :: PHW

    timerev_k = .false.

    ind_ktype  = 0
    nkt(:)     = 0
    nvar(:)    = 0
    k_trial(:) = 0._dp
    do iR=1, num_doub_sym/2 + 1
        if (iR < num_doub_sym/2 + 1) then 
            k_trial(:) = matmul(k_input(:),invrot(:,:,iR))
            !if (k_trial(1) < 0.0 .and. abs(k_trial(1)) > 1e-6) k_trial(1) = k_trial(1) + 1.d0
            !if (k_trial(2) < 0.0 .and. abs(k_trial(2)) > 1e-6) k_trial(2) = k_trial(2) + 1.d0
            !if (k_trial(3) < 0.0 .and. abs(k_trial(3)) > 1e-6) k_trial(3) = k_trial(3) + 1.d0
            call getkid2(k_trial,ind_ktype_tmp, ivar)
        else 
            k_trial(:) = -k_input(:)
            call getkid2(k_trial, ind_ktype_tmp, ivar)
        endif 
        nkt(iR)  = ind_ktype_tmp
        nvar(iR) = ivar
        if (ivar == 0) exit 
    enddo 

    if (iR == num_doub_sym/2 + 1) then 
        timerev_k = .true.
    endif 

    if (iR == num_doub_sym/2 + 2) then 
        zj:DO iw = 1,3 
            do iR = 1, num_doub_sym/2
                ivar = nvar(iR)
                if (ivar == iw) exit zj
            enddo
        ENDDO zj
    endif 
    
    ind_rot   = iR
    ind_ktype = nkt(iR)
    ivar      = nvar(iR)

    if (ind_ktype > num_ktype) stop " Nonsymmorphic kpoint is NOT found."
    kname = samplekname(ind_ktype)//' '

    k_trial(1) = dot_product(k_input(:), invrot(:,1,iR))
    k_trial(2) = dot_product(k_input(:), invrot(:,2,iR))
    k_trial(3) = dot_product(k_input(:), invrot(:,3,iR))
    if (timerev_k) k_trial(:) = -k_input(:)
    k_input(:) = k_trial(:)
     
    k_input_conv(:) = matmul(k_input(:), p2cR(:,:))
    do iir = 1,nirreps(ind_ktype)
        do j = 1,num_doub_sym

            if(labels(1,j,iir,ind_ktype) == 1) then
                
                if(labels(2,j,iir,ind_ktype) == 1) then
                    tableTraces(j,iir,ind_ktype)=chartTraces(j,iir,ind_ktype)
                elseif(labels(2,j,iir,ind_ktype) == 2) then
                    AGW=0._DP
                    AGW=PI*DOT_PRODUCT(coeff_uvw(1:3,j,iir,ind_ktype),k_input_conv(:))
                    PHW=CMPLX(DCOS(AGW),DSIN(AGW),DP)
                    tableTraces(j,iir,ind_ktype)=chartTraces(j,iir,ind_ktype)*PHW
                else
                    stop "Error" 
                endif
            endif

        enddo
    enddo

end subroutine bilbao_getkid 


      subroutine getkid2(k_input,ik,ivar)
      real(dp)   , intent(in ) :: k_input(3)
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
      real(dp) :: k_input_conv(3), refkpoint(3), refkpointp(3)
      real(dp) :: diff, diff1, diff2, diff3
      real(dp), parameter :: tol=1e-5

      character*15 :: chr_tmp
      ivar=-1

      allocate(ckpoint(num_ktype))
      do ikt=1, num_ktype
         call Kreal2string(samplek(:,ikt),ckpoint(ikt)) 
      enddo

      

      allocate(num_var(num_ktype))
      num_var = 0
      allocate(ind_u(num_ktype))
      allocate(ind_v(num_ktype))
      allocate(ind_w(num_ktype))
      ind_u=0; ind_v=0; ind_w=0

     !print*,"================================================"
     !   write(*,'(20X,A5,10X,A6,20X,A5)') 'conv', 'symbol', 'prim'
     !do ikt=1,num_ktype
     !   write(*,'(I5,3F9.3,2X,A15,3F12.6)') ikt,samplek(:,ikt),ckpoint(ikt),matmul(samplek(:,ikt),Kc2p(:,:))
     !enddo
     !   write(*,*) 'Caution!!! conv  : 0.333 -> 0.3333333333'
     !   write(*,*) 'Caution!!! symbol: 0.33  -> 0.3333333333'
     ! !--- add something here

      !--- add something below-------
      ! get the number of variables for the reference kpoints
      do ikt = 1, num_ktype
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
      k_input_conv = matmul(k_input(:), p2cR(:,:))
      !write(*,"(A35,5X,3F12.6)") "Input kpoint under primitive :", k_input(:)
      !write(*,"(A35,5X,3F12.6)") "Input kpoint under convention:", k_input_conv(:)
      !PAUSE

      ! compare the input kpoint with the reference kpoints
      ! the reference kpoints with less variables have a higher priority
      uvw = 0d0 
      ref_varnum = 9999
      do ikt = 1, num_ktype
         is_variable = .false.
         refkpoint = 9999d0
         uvw_tmp = 0d0 
         if(ckpoint(ikt)(1:5)=='  u  ') then 
            is_variable(1) = .true.
            uvw_tmp(1) = k_input_conv(1)
            refkpoint(1) = k_input_conv(1)
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
            uvw_tmp(2) = k_input_conv(2)
            refkpoint(2) = k_input_conv(2)
            if (ckpoint(ikt)(1:5)==' 1-v ') refkpoint(1) = 1d0-refkpoint(2)
            if (ckpoint(ikt)(11:15)==' 1-v ') refkpoint(3) = 1d0-refkpoint(2)
            if (ckpoint(ikt)(11:15)=='  v  ') refkpoint(3) = refkpoint(2)
         endif 
         if(ckpoint(ikt)(11:15)=='  w  ') then 
             is_variable(3) = .true.
             uvw_tmp(3) = k_input_conv(3)
             refkpoint(3) = k_input_conv(3)
         endif 

         ! some exceptions
         ! 0.5 u 0.0  : 195 198 200 201 205
         ! 1+u 1-u 0.0: 197 199 204 206 211 214 217 220 229 230
         if(ckpoint(ikt)(1:5).ne.'  u  '.and.ckpoint(ikt)(6:10).eq.'  u  ') then 
            is_variable(1) = .true.
            uvw_tmp(1) = k_input_conv(2)
            refkpoint(2) = k_input_conv(2)
         endif 
         if(ckpoint(ikt)(1:5).eq.' 1+u '.and.ckpoint(ikt)(6:10).eq.' 1-u ') then 
            is_variable(1) = .true.
            uvw_tmp(1) = k_input_conv(1) - 1
            refkpoint(1) = k_input_conv(1)
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

         !diff = abs(mod(refkpointp(1),1d0)-mod(k_input(1),1d0)) + &
         !       abs(mod(refkpointp(2),1d0)-mod(k_input(2),1d0)) + &
         !       abs(mod(refkpointp(3),1d0)-mod(k_input(3),1d0))
         diff = abs((refkpointp(1)-floor(refkpointp(1)))-(k_input(1)-floor(k_input(1)))) + &
                abs((refkpointp(2)-floor(refkpointp(2)))-(k_input(2)-floor(k_input(2)))) + &
                abs((refkpointp(3)-floor(refkpointp(3)))-(k_input(3)-floor(k_input(3))))    
         !write(*,*) "now string", ckpoint(ikt)
         !write(*,"(A35,5X,3F12.6)") "ref   kpoint under primitive :", refkpointp(:)
         !write(*,"(A35,5X,3F12.6)") "Input kpoint under primitive :", k_input(:)
         !write(*,"(A35,5X,3F12.6)") "ref   kpoint under prim mod 1:", mod(refkpointp(1),1d0),mod(refkpointp(2),1d0),mod(refkpointp(3),1d0) 
         !write(*,"(A35,5X,3F12.6)") "Input kpoint under prim mod 1:", mod(k_input(1),1d0),mod(k_input(1),1d0),mod(k_input(1),1d0)
         !write(*,"(A35,5X,3F12.6)") "diff                         :", diff
         !PAUSE
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


subroutine Dmatrix(ih, ik, il, rnx, rny, rnz, rtheta, rphi, degree, twoja1, Dmat)
!!ref    : https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
!!        s, px, py, pz, xy, yz, zx, x2-y2, 3z2-r2
!!        0,  1, -1,  0; -2, -1,   1,   2,  0
!!        fxyz, f5x3-xr2, f5y3-yr2, f5z3-zr2, fx(y2-z2), fy(z2-x2), fz(x2-y2)
!!TABLE I: http://journals.aps.org/prb/pdf/10.1103/PhysRevB.79.045107
  implicit none
  integer, intent(in), optional :: ih, ik, il
  real(dp), intent(in), optional :: rnx, rny, rnz
  real(dp), intent(in), optional :: rtheta, rphi
  real(dp), intent(in) ::  degree
  integer, intent(in) :: twoja1  ! equal to 2*J + 1
  ! Actually, this rotational Dmat should be a real matrix.
  complex(dp),dimension(twoja1,twoja1),intent(out) :: Dmat
  complex(dp),dimension(twoja1,twoja1)             :: Dtmp
  ! the key parameter for Dmatrix
  real(dp) :: nx, ny, nz, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4, nx5,ny5,nz5,nx6,ny6,nz6
  real(dp) ::  omega
  ! complex(dp),dimension(aint(J)*2+1,aint(J)*2+1) :: Jx,Jy,Jz
  ! integer :: m
  complex(dp),parameter:: ci = (0.d0,1.d0)
  real*8     ,parameter:: PI = 3.141592653589793238462643383279d0
  real(dp) :: tmp
  real(dp) :: sqrt3, sqrt15 
  real(dp) :: cosOmega_2, sinOmega_2
  real(dp) :: cos1omega, sin1omega, cos2omega, sin2omega, cos3omega, sin3omega, cos4omega, sin4omega
  real(dp) :: cos5omega, sin5omega, cos6omega, sin6omega
  real(dp) :: cos3omega_2, sin3omega_2, cos5omega_2, sin5omega_2 
  complex(dp) :: exp2omega, exp_2omega, exp4omega, exp_4omega, expomega_2, exp_omega_2
  complex(dp) :: exp3omega, exp_3omega, exp6omega, exp_6omega, expomega, exp_omega, exp5omega, exp_5omega  
  sqrt3 = dsqrt(3.0d0)
  sqrt15 = dsqrt(15.0d0)
  if(present(ih).and.present(ik).and.present(il))then
    tmp = dsqrt(dble(ih)*dble(ih) + dble(ik)*dble(ik) + dble(il)*dble(il))
    nx = dble(ih)/tmp
    ny = dble(ik)/tmp
    nz = dble(il)/tmp
  else if(present(rtheta).and.present(rphi))then
    nx = dsin(rtheta)*dcos(rphi)
    ny = dsin(rtheta)*dsin(rphi)
    nz = dcos(rtheta)
  else if(present(rnx).and.present(rny).and.present(rnz))then
    tmp = dsqrt(rnx*rnx + rny*rny + rnz*rnz )
    if(abs(tmp-1.d0).gt.1.d-6) then
      write(0,'(A,F16.10)')"WARNING!!! The nx, ny, nz input is not normailized: ",tmp
      write(0,*)"The program will normalize them now."
    end if
      nx = rnx/tmp
      ny = rny/tmp
      nz = rnz/tmp
  else
    write(*,*)" Please choose one of the three type of input: "
    write(*,*)" 1. h k l:"
    write(*,*)" 2. nx ny nz:"
    write(*,*)" 3. theta phi:"
  end if
  omega = dble(nint(degree))*PI/180.d0
  IF(twoja1 == 5) omega=omega/2.d0
  cosOmega_2 = dcos(omega/2.d0)
  sinOmega_2 = dsin(omega/2.d0)
  cos3omega_2 = dcos(3.d0*omega/2.d0)
  sin3omega_2 = dsin(3.d0*omega/2.d0)
  cos5omega_2 = dcos(5.d0*omega/2.d0)
  sin5omega_2 = dsin(5.d0*omega/2.d0)
  cos1omega = dcos(omega)
  sin1omega = dsin(omega)
  cos2omega = dcos(2.d0*omega)
  sin2omega = dsin(2.d0*omega)
  cos3omega = dcos(3.d0*omega)
  sin3omega = dsin(3.d0*omega)
  cos4omega = dcos(4.d0*omega)
  sin4omega = dsin(4.d0*omega)
  cos5omega = dcos(5.d0*omega)
  sin5omega = dsin(5.d0*omega)
  cos6omega = dcos(6.d0*omega)
  sin6omega = dsin(6.d0*omega)
  expomega = cmplx(cos1omega, sin1omega)
  exp_omega = cmplx(cos1omega, -sin1omega)
  exp2omega = cmplx(cos2omega,sin2omega)
  exp_2omega = cmplx(cos2omega,-sin2omega)
  exp4omega = cmplx(cos4omega,sin4omega)
  exp_4omega = cmplx(cos4omega,-sin4omega)
  exp3omega = cmplx(cos3omega,sin3omega)
  exp_3omega = cmplx(cos3omega,-sin3omega)
  exp5omega = cmplx(cos5omega, sin5omega)
  exp_5omega = cmplx(cos5omega, -sin5omega)
  exp6omega = cmplx(cos6omega,sin6omega)
  exp_6omega = cmplx(cos6omega,-sin6omega)
  expomega_2 = cmplx(cosOmega_2, sinOmega_2)
  exp_omega_2 = cmplx(cosOmega_2, -sinOmega_2)
  nx2 = nx**2
  ny2 = ny**2
  nz2 = nz**2
  nx3 = nx**3
  ny3 = ny**3
  nz3 = nz**3
  nx4 = nx**4
  ny4 = ny**4
  nz4 = nz**4
  nx5 = nx**5
  ny5 = ny**5
  nz5 = nz**5
  nx6 = nx**6
  ny6 = ny**6
  nz6 = nz**6
  if(twoja1 == 2) then
    Dmat(1,1) = cmplx(cosOmega_2, -nz*sinOmega_2)
    Dmat(2,2) = cmplx(cosOmega_2, nz*sinOmega_2)
    Dmat(1,2) = cmplx(-ny, -nx)*sinOmega_2
    Dmat(2,1) = cmplx(ny, -nx)*sinOmega_2
  else if(twoja1 == 3) then
    Dmat(1,1) = nx**2+(1-nx**2)*cos1omega
    Dmat(2,2) = ny**2+(1-ny**2)*cos1omega
    Dmat(3,3) = nz**2+(1-nz**2)*cos1omega
    Dmat(1,2) = nx*ny*(1-cos1omega) - nz*sin1omega
    Dmat(2,1) = nx*ny*(1-cos1omega) + nz*sin1omega
    Dmat(1,3) = nx*nz*(1-cos1omega) + ny*sin1omega
    Dmat(3,1) = nx*nz*(1-cos1omega) - ny*sin1omega
    Dmat(2,3) = ny*nz*(1-cos1omega) - nx*sin1omega
    Dmat(3,2) = ny*nz*(1-cos1omega) + nx*sin1omega
  else if(twoja1 == 5) then
    Dmat(1,1) = 3.d0*ny**2*nx**2 + &
                0.5d0*(exp_4omega + exp4omega)*(1-ny**2)*(1-nx**2) + &
                0.5d0*(exp_2omega + exp2omega)*(nz**2*(1-nz**2)+(nx**2-ny**2)**2)
    Dmat(2,2) = 3.d0*ny**2*nz**2 + &
                0.5d0*(exp_4omega + exp4omega)*(1-ny**2)*(1-nz**2) + &
                0.5d0*(exp_2omega + exp2omega)*(nx**2*(1-nx**2)+(ny**2-nz**2)**2)
    Dmat(3,3) = 3.d0*nz**2*nx**2 + &
                0.5d0*(exp_4omega + exp4omega)*(1-nz**2)*(1-nx**2) + &
                0.5d0*(exp_2omega + exp2omega)*(ny**2*(1-ny**2)+(nx**2-nz**2)**2)
    Dmat(4,4) = 0.75d0*(nx**2-ny**2)**2 + &
                0.125d0*(exp_4omega + exp4omega)*(4*nz**2+(nx**2-ny**2)**2) + &
                0.5d0*(exp_2omega + exp2omega)*(nz**2*(1-nz**2) + 4*ny**2*nx**2)
    Dmat(5,5) = 0.25d0*(1-3*nz**2)**2 + &
                1.5d0*(exp_2omega + exp2omega)*nz**2*(1-nz**2) + &
                0.375d0*(exp_4omega + exp4omega)*(1-nz**2)**2

    Dmat(1,2) = -2.d0*ci*sin1omega*(ci*ny**3*cos1omega + ci*ny*(1-ny**2)*cos3omega + &
                                      nx*nz*( 3.d0*ci*ny**2*sin1omega+ci*(1-ny**2)*sin3omega))
    Dmat(2,1) =  2.d0*ci*sin1omega*(ci*ny**3*cos1omega + ci*ny*(1-ny**2)*cos3omega + &
                                      nx*nz*(-3.d0*ci*ny**2*sin1omega-ci*(1-ny**2)*sin3omega))
    Dmat(2,3) = -2.d0*ci*sin1omega*(ci*nz**3*cos1omega + ci*nz*(1-nz**2)*cos3omega + &
                                      nx*ny*( 3.d0*ci*nz**2*sin1omega+ci*(1-nz**2)*sin3omega))
    Dmat(3,2) =  2.d0*ci*sin1omega*(ci*nz**3*cos1omega + ci*nz*(1-nz**2)*cos3omega + &
                                      nx*ny*(-3.d0*ci*nz**2*sin1omega-ci*(1-nz**2)*sin3omega))
    Dmat(1,3) =  2.d0*ci*sin1omega*(ci*nx**3*cos1omega + ci*nx*(1-nx**2)*cos3omega + &
                                      ny*nz*(-3.d0*ci*nx**2*sin1omega-ci*(1-nx**2)*sin3omega))
    Dmat(3,1) = -2.d0*ci*sin1omega*(ci*nx**3*cos1omega + ci*nx*(1-nx**2)*cos3omega + &
                                      ny*nz*( 3.d0*ci*nx**2*sin1omega+ci*(1-nx**2)*sin3omega))

    Dmat(1,4) = -ci*sin1omega*(ci*nz*(3.d0-nz**2)*cos1omega + ci*nz*(1+nz**2)*cos3omega - &
                                 4.d0*ci*nx*ny*(ny**2-nx**2)*sin1omega**3)
    Dmat(4,1) =  ci*sin1omega*(ci*nz*(3.d0-nz**2)*cos1omega + ci*nz*(1+nz**2)*cos3omega + &
                                 4.d0*ci*nx*ny*(ny**2-nx**2)*sin1omega**3)
    Dmat(1,5) = -2.d0*sqrt3*sin1omega**2*(nx*ny*(1.d0-nz**2) + nx*ny*(1+nz**2)*cos2omega - &
                                               nz*(ny**2-nx**2)*sin2omega)
    Dmat(5,1) = -2.d0*sqrt3*sin1omega**2*(nx*ny*(1.d0-nz**2) + nx*ny*(1+nz**2)*cos2omega + &
                                               nz*(ny**2-nx**2)*sin2omega)

    Dmat(2,4) = 1.5d0*(nx**2 - ny**2)*ny*nz + ny*nz*(2.d0*ny**2 - 2.d0*nx**2 - 1.d0)*cos2omega + &
                0.5d0*ny*nz*(2.d0*nx**2 + nz**2 + 1.d0)*cos4omega - &
                nx*( 2.d0*ny**2 - nz**2 + (1.d0 - 2.d0*ny**2 + nz**2)*cos2omega)*sin2omega
    Dmat(3,4) = 1.5d0*(nx**2 - ny**2)*nx*nz + nx*nz*(4.d0*ny**2 + 2.d0*nz**2 - 1.d0)*cos2omega - &
                0.5d0*nx*nz*(2.d0*ny**2 + nz**2 + 1.d0)*cos4omega - &
                ny*( 2.d0*nx**2 - nz**2 + (-1.d0 + 2.d0*ny**2 + 3.d0*nz**2)*cos2omega)*sin2omega
    Dmat(4,2) = 1.5d0*(nx**2 - ny**2)*ny*nz + ny*nz*(2.d0*ny**2 - 2.d0*nx**2 - 1.d0)*cos2omega + &
                0.5d0*ny*nz*(2.d0*nx**2 + nz**2 + 1.d0)*cos4omega + &
                nx*( 2.d0*ny**2 - nz**2 + (1.d0 - 2.d0*ny**2 + nz**2)*cos2omega)*sin2omega
    Dmat(4,3) = 1.5d0*(nx**2 - ny**2)*nx*nz + nx*nz*(4.d0*ny**2 + 2.d0*nz**2 - 1.d0)*cos2omega - &
                0.5d0*nx*nz*(2.d0*ny**2 + nz**2 + 1.d0)*cos4omega + &
                ny*( 2.d0*nx**2 - nz**2 + (-1.d0 + 2.d0*ny**2 + 3.d0*nz**2)*cos2omega)*sin2omega

    Dmat(2,5) = 2.d0*ci*sqrt3*(nz**2 + (1 - nz**2)*cos2omega)*sin1omega*ci*(nx*cos1omega - ny*nz*sin1omega)
    Dmat(3,5) =-2.d0*ci*sqrt3*(nz**2 + (1 - nz**2)*cos2omega)*sin1omega*ci*(ny*cos1omega + nx*nz*sin1omega)
    Dmat(5,3) = 2.d0*ci*sqrt3*(nz**2 + (1 - nz**2)*cos2omega)*sin1omega*ci*(ny*cos1omega - nx*nz*sin1omega)
    Dmat(5,2) =-2.d0*ci*sqrt3*(nz**2 + (1 - nz**2)*cos2omega)*sin1omega*ci*(nx*cos1omega + ny*nz*sin1omega)

    Dmat(4,5) = -sqrt3*sin1omega**2*(nx**4 - ny**4 + (nx**2 - ny**2)*(1.d0 + nz**2)*cos2omega - &
                                          4.d0*nx*ny*nz*sin2omega)
    Dmat(5,4) = -sqrt3*sin1omega**2*(nx**4 - ny**4 + (nx**2 - ny**2)*(1.d0 + nz**2)*cos2omega + &
                                          4.d0*nx*ny*nz*sin2omega)
  else if (twoja1 == 7) then 
    Dmat(1,1) = 0.25d0*exp_3omega*(60d0*exp3omega*nx2*ny2*nz2 + 3d0*(1d0-nx2)*(1d0-ny2)*(1d0-nz2) + &
                3d0*exp6omega*(1d0-nx2)*(1d0-ny2)*(1d0-nz2) + &
                5d0*exp2omega*(nx4*(1d0-nx2)+(1d0-nx2)*ny2*nz2 + nx2*(ny4-6d0*ny2*nz2+nz4)) + &
                5d0*exp4omega*(nx4*(1d0-nx2)+(1d0-nx2)*ny2*nz2 + nx2*(ny4-6d0*ny2*nz2+nz4)) + &
                2d0*expomega *(nx6-nx4*(1d0-nx2)+(1d0-nx2)*(ny2-nz2)**2 - nx2*(ny4-3d0*ny2*nz2+nz4)) + &
                2d0*exp5omega*(nx6-nx4*(1d0-nx2)+(1d0-nx2)*(ny2-nz2)**2 - nx2*(ny4-3d0*ny2*nz2+nz4))   )
    Dmat(2,1) = -2d0*sqrt15*(nx2+(1d0-nx2)*cos1omega)*sinOmega_2**2*((1d0-nx2)*ny*nz+(1d0+nx2)*ny*nz*cos1omega+nx*(nz2-ny2)*sin1omega)
    Dmat(3,1) = -2d0*sqrt15*(ny2+(1d0-ny2)*cos1omega)*sinOmega_2**2*((1d0-ny2)*nx*nz+(1d0+ny2)*nx*nz*cos1omega+ny*(nx2-nz2)*sin1omega)
    Dmat(4,1) = -2d0*sqrt15*(nz2+(1d0-nz2)*cos1omega)*sinOmega_2**2*((1d0-nz2)*nx*ny+(1d0+nz2)*nx*ny*cos1omega+nz*(ny2-nx2)*sin1omega)
    Dmat(5,1) = 0.25d0*ci*exp_omega_2*(expomega-1d0)*(2d0*nx3*(2d0*nx2+5d0*(1d0-nx2))*cosOmega_2 + &
                                                      nx*(4d0*nx4+5d0*(1d0-nx2)**2)*cos3omega_2 + &
                                                      6d0*nx3*ny2*cos5omega_2 + 3d0*nx*ny4*cos5omega_2 + 6d0*nx3*nz2*cos5omega_2 + &
                                                      6d0*nx*ny2*nz2*cos5omega_2 + 3d0*nx*nz4*cos5omega_2 + &
                                                      30d0*nx2*ny*nz*(nz2-ny2)*sinOmega_2 + &
                                                      5d0*(1d0-3d0*nx2)*ny*nz*(nz2-ny2)*sin3omega_2 + &
                                                      3d0*ny*nz*(ny4-nz4)*sin5omega_2 )
    Dmat(6,1) = -0.5d0*sinOmega_2*(2d0*ny3*(5d0*nx2+2d0*ny2+5d0*nz2)*cosOmega_2 + &
                                            ny*(4d0*ny4+5d0*(1d0-ny2)**2)*cos3omega_2 + &
                                            6d0*ny3*nx2*cos5omega_2 + 3d0*ny*nx4*cos5omega_2 + 6d0*nx2*ny*nz2*cos5omega_2 + &
                                            6d0*ny3*nz2*cos5omega_2 + 3d0*ny*nz4*cos5omega_2 + &
                                            30d0*ny2*nx*nz*(nx2-nz2)*sinOmega_2 + &
                                            5d0*(1d0-3d0*ny2)*nx*nz*(nx2-nz2)*sin3omega_2 + &
                                            3d0*nx*nz*(nz4-nx4)*sin5omega_2 )
    Dmat(7,1) = -0.5d0*sinOmega_2*(2d0*nz3*(2d0*nz2+5d0*(1d0-nz2))*cosOmega_2 + &
                                   nz*(4d0*nz4+5d0*(1d0-nz2)**2)*cos3omega_2 + &
                                   6d0*nz3*nx2*cos5omega_2 + 3d0*nz*nx4*cos5omega_2 + 6d0*nx2*nz*ny2*cos5omega_2 + &
                                   6d0*nz3*ny2*cos5omega_2 + 3d0*nz*ny4*cos5omega_2 + &
                                   30d0*nz2*nx*ny*(nx2-ny2)*sinOmega_2 + &
                                   5d0*(1d0-3d0*nz2)*nx*ny*(ny2-nx2)*sin3omega_2 + &
                                   3d0*nx*ny*(nx4-ny4)*sin5omega_2) 
    Dmat(1,2) = -2d0*sqrt15*(nx2+(1d0-nx2)*cos1omega)*sinOmega_2**2*((1d0-nx2)*ny*nz+(1d0+nx2)*ny*nz*cos1omega+nx*(ny2-nz2)*sin1omega)
    Dmat(2,2) = 1d0/16d0*exp_3omega*( 3d0*exp2omega*(1d0-5d0*nx2)**2*(1d0-nx2) + 3d0*exp4omega*(1d0-5d0*nx2)**2*(1d0-nx**2) + &
                                       30d0*expomega*nx2*(1d0-nx2)**2 + 30d0*exp5omega*nx2*(1d0-nx2)**2 + &
                                       5d0*(1d0-nx2)**3 + 5d0*exp6omega*(1d0-nx2)**3 + 4d0*exp3omega*nx2*(2d0*nx2-3d0*(1d0-nx2))**2 )
    Dmat(3,2) = -0.25d0*sinOmega_2*(2d0*nz*(6d0*nx4-3d0*nx2*ny2+6d0*ny4+nz4+7d0*nz2*(1d0-nz2))*cosOmega_2 + &
                                    5d0*nz*(9d0*nx2*ny2+nz4+nz2*(1d0-nz2))*cos3omega_2 - 15d0*nx2*ny2*nz*cos5omega_2 + &
                                    5d0*nx2*nz3*cos5omega_2 + 5d0*ny2*nz3*cos5omega_2 + 5d0*nz5*cos5omega_2 + &
                                    2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*ny*(3d0*nx2-2d0*ny2+3d0*nz2)*sinOmega_2 + &
                                    5d0*nx*ny*(5d0*nx2*ny2-3d0*nz4-3d0*nz2*(1d0-nz2))*sin3omega_2 + &
                                    5d0*nx*ny*(-nx2*ny2+3d0*nz4+3d0*nz2*(1-nz2))*sin5omega_2 )
    Dmat(4,2) = 0.25d0*sinOmega_2*(2d0*ny*(6d0*nx4+nx2*(7d0*ny2-3d0*nz2)+(1d0-nx2)*(ny2+6d0*nz2))*cosOmega_2 + &
                                   5d0*ny*((1-nx2)*ny2+nx2*(ny2+9d0*nz2))*cos3omega_2 - 15d0*nx2*ny*nz2*cos5omega_2 + &
                                   5d0*nx2*ny3*cos5omega_2 + 5d0*ny5*cos5omega_2 + 5d0*ny3*nz2*cos5omega_2 + &
                                   (-2d0)*nx*(2d0*nx2-3d0*(1d0-nx2))*nz*(-2d0*nz2+3d0*(1d0-nz2))*sinOmega_2 + &
                                   5d0*nx*nz*(3d0*(1d0-nx2)*ny2+nx2*(3d0*ny2-5d0*nz2))*sin3omega_2 + &
                                   5d0*nx*nz*(-3d0*(1d0-nx2)*ny2 + nx2*(-3d0*ny2+nz2))*sin5omega_2 )
    Dmat(5,2) = -sqrt15*(nx2+(1d0-nx2)*cos1omega)*sinOmega_2**2*(ny4-nz4+(1d0+nx2)*(ny2-nz2)*cos1omega-4d0*nx*ny*nz*sin1omega)
    Dmat(6,2) = 0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*(2d0*nz*(2d0*nx4+3d0*nx2*ny2+nz4+nz2*(1d0-nz2))*cosOmega_2 + &
                                                              nz*(-nx2*(ny2-7d0*nz2)+(1d0-nx2)*(2d0*ny2+nz2))*cos3omega_2 + &
                                                              3d0*nx2*ny2*nz*cos5omega_2 + 2d0*ny4*nz*cos5omega_2 + &
                                                            (-1d0)*nx2*nz3*cos5omega_2 + 3d0*ny2*nz3*cos5omega_2 + nz5*cos5omega_2 + &
                                                             2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*ny*(nx2-nz2)*sinOmega_2 + &
                                                             nx*ny*(-2d0*ny4-ny2*nz2+nz4+nx2*(3d0*ny2+11d0*nz2))*sin3omega_2 + &
                                                             nx*ny*(2d0*ny4+ny2*nz2-nz4+nx2*(ny2-3d0*nz2))*sin5omega_2 )
    Dmat(7,2) = 0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*(2d0*ny*(2d0*nx4+nx2*ny2+ny4+(3d0*nx2+ny2)*nz2)*cosOmega_2 + &
                                                              ny*(-nx2*(nz2-7d0*ny2)+(1d0-nx2)*(2d0*nz2+ny2))*cos3omega_2 + &
                                                              3d0*nx2*ny*nz2*cos5omega_2 + 2d0*ny*nz4*cos5omega_2 + &
                                                            (-1d0)*nx2*ny3*cos5omega_2 + 3d0*ny3*nz2*cos5omega_2 + ny5*cos5omega_2 + &
                                                             2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*nz*(ny2-nx2)*sinOmega_2 + &
                                                            (-1d0)*nx*nz*((1d0-nx2)*(ny2-2d0*nz2)+nx2*(11d0*ny2+3d0*nz2))*sin3omega_2 + &
                                                             nx*nz*(3d0*nx2*ny2+ny4-2d0*nz4-nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(1,3) = -2d0*sqrt15*(ny2+(1d0-ny2)*cos1omega)*sinOmega_2**2*(nx*(1d0-ny2)*nz+nx*(1d0+ny2)*nz*cos1omega+ny*(nz2-nx2)*sin1omega)
    Dmat(2,3) = 0.25d0*sinOmega_2*(2d0*nz*(6d0*nx4-3d0*nx2*ny2+6d0*ny4+nz4+7d0*nz2*(1d0-nz2))*cosOmega_2 + &
                                   5d0*nz*(9d0*nx2*ny2+nz4+nz2*(1d0-nz2))*cos3omega_2 - 15d0*nx2*ny2*nz*cos5omega_2 + &
                                   5d0*nz3*cos5omega_2 - & 
                                   2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*ny*(3d0*nx2-2d0*ny2+3d0*nz2)*sinOmega_2 + &
                                   5d0*nx*ny*(-5d0*nx2*ny2+3d0*nz4+3d0*nz2*(1d0-nz2))*sin3omega_2 + &
                                   5d0*nx*ny*(nx2*ny2-3d0*nz4-3d0*nz2*(1-nz2))*sin5omega_2 )
    Dmat(3,3) = 1d0/16d0*exp_3omega*(3d0*exp2omega*(1d0-5d0*ny2)**2*(1d0-ny2) + 3d0*exp4omega*(1d0-5d0*ny2)**2*(1d0-ny2) + &
                                     30d0*expomega*ny2*(1d0-ny2)**2 + 30d0*exp5omega*ny2*(1d0-ny2)**2 + 5d0*(1d0-ny2)**3 + &
                                     5d0*exp6omega*(1d0-ny2)**3 + 4d0*exp3omega*ny2*(3d0*nx2-2d0*ny2+3d0*nz2)**2 )
    Dmat(4,3) = 0.125d0*(4d0*ny*nz*(3d0*(3d0*nx4+nx2*(1d0-nx2)-2d0*ny4 + ny2*nz2 - 2d0*nz4) -20d0*ny2*nz2*cos1omega - &
                                    5d0*(3d0*nx4+3d0*nx2*(1d0-nx2)-ny2*nz2)*cos2omega)*sinOmega_2**2 - &
                                    2d0*nx*(nx4+7d0*nx2*(1d0-nx2)+6d0*ny4-33d0*ny2*nz2+6d0*nz4+60d0*ny2*nz2*cos1omega + &
                                    5d0*(nx4+nx2*(1d0-nx2)-3d0*ny2*nz2)*cos2omega) * sin1omega )
    Dmat(5,3) = 1d0/16d0*sqrt15*exp_3omega*(-2d0*(expomega-1d0)**2*(1d0+exp4omega)*nx5*ny + &
                                            2d0*ci*nx4*nz*(-1d0+exp4omega*(2d0*cos2omega-1d0)) + &
                                        8d0*exp3omega*nx3*ny*(-ny2+3d0*nz2+4d0*ny2*cos1omega+(1d0-nx2)*cos2omega)*sinOmega_2**2 - &
                 8d0*exp3omega*nx*ny*(-2d0*ny4+ny2*nz2-3d0*nz4+nz2*(-8d0*ny2*cos1omega + (3d0*ny2+nz2)*cos2omega))*sinOmega_2**2- &
                 4d0*exp3omega*nx2*nz*(5d0*ny2+nz2-4d0*ny2*cos1omega+3d0*(1d0-nx2)*cos2omega)*sin1omega - &
                 4d0*exp3omega*nz*(8d0*ny2*nz2*cos1omega + (ny2-nz2)*(2d0*ny2-nz2-nz2*cos2omega))*sin1omega )
    Dmat(6,3) = sqrt15*(ny2+(1d0-ny2)*cos1omega)*sinOmega_2**2 * (nx4-nz4+(1+ny2)*(nx2-nz2)*cos1omega+4d0*nx*ny*nz*sin1omega)
    Dmat(7,3) = 0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)* &
                ( 2d0*nx*cosOmega_2 * ( nx4 - 3d0*nx2*ny2 + 2d0*ny4 + nz2*(1d0+4d0*ny2-nz2) + 4d0*ny2*(2d0*nx2-nz2)*cos1omega + &
                                      (nx2*(nx2-ny2)+2d0*nz4+3d0*nz2*(1d0-nz2))*cos2omega ) + &
                  ny*nz * (2d0*(nx2-ny2)*(3d0*nx2-2d0*ny2+3d0*nz2)*sinOmega_2 + &
                          (nx4 + 3d0*ny2*nz2 - 2d0*nz4 + nx2*(11d0*ny2-nz2))*sin3omega_2 + &
                          (2d0*nz4 + nz2*(1d0-nz2) - nx2*(1d0+2d0*ny2-nz2))*sin5omega_2 ))
    Dmat(1,4) = -2d0*sqrt15*(nz2+(1d0-nz2)*cos1omega)*sinOmega_2**2*(nx*ny*(1d0-nz2)+nx*ny*(1d0+nz2)*cos1omega+(nx2-ny2)*nz*sin1omega)
    Dmat(2,4) = -0.25d0*sinOmega_2*(2d0*ny*(6d0*nx4+nx2*(7d0*ny2-3d0*nz2) + (1d0-nx2)*(ny2+6d0*nz2))*cosOmega_2 + &
                                    5d0*ny*((1d0-nx2)*ny2+nx2*(ny2+9d0*nz2))*cos3omega_2 + 5d0*nx2*ny3*cos5omega_2 + &
                                    5d0*ny5*cos5omega_2 -15d0*nx2*ny*nz2*cos5omega_2 + 5d0*ny3*nz2*cos5omega_2 + &
                                    2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*nz*(-2d0*nz2+3d0*(1d0-nz2))*sinOmega_2 - &
                                    5d0*nx*nz*(3d0*(1d0-nx2)*ny2+nx2*(3d0*ny2-5d0*nz2))*sin3omega_2 + &
                                    5d0*nx*nz*(3d0*(1d0-nx2)*ny2+nx2*(3d0*ny2-nz2))*sin5omega_2 )
    Dmat(3,4) = 0.125d0*(4d0*ny*nz*(3d0*(3d0*nx4+nx2*(1d0-nx2)-2d0*ny4+ny2*nz2-2d0*nz4) - 20d0*ny2*nz2*cos1omega - &
                                     5d0*(3d0*nx4+3d0*nx2*(1d0-nx2)-ny2*nz2)*cos2omega)*sinOmega_2**2 + &
                          2d0*nx*(nx4+7d0*nx2*(1d0-nx2)+6d0*ny4-33d0*ny2*nz2+6d0*nz4+60*ny2*nz2*cos1omega + &
                          5d0*(nx4+nx2*(1d0-nx2)-3d0*ny2*nz2)*cos2omega)*sin1omega )
    Dmat(4,4) = 1d0/16d0*exp_3omega*(3d0*exp2omega*(1d0-5d0*nz2)**2*(1d0-nz2)+3d0*exp4omega*(1d0-5d0*nz2)**2*(1d0-nz2) + &
                                     30d0*expomega*nz2*(1d0-nz2)**2+30d0*exp5omega*nz2*(1d0-nz2)**2 + &
                                     5d0*(1d0+exp6omega)*(1d0-nz2)**3 + 4d0*exp3omega*(-2d0*nz3+3d0*nz*(1d0-nz2))**2) 
    Dmat(5,4) = 0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*(2d0*ny*cosOmega_2*(ny4-3d0*ny2*nz2+2d0*nz4+nx2*(ny2+5d0*nz2) - &
                                        4d0*nz2*(1d0-3d0*ny2-nz2)*cos1omega + (2d0*nx4+3d0*nx2*(1d0-nx2)+ny2*(ny2-nz2))*cos2omega)+&
                                        nx*nz*(2d0*(ny2-nz2)*(3d0-5d0*nz2)*sinOmega_2 + &
                                        (-2d0*nx4+ny4+11d0*ny2*nz2-nx2*(ny2-3d0*nz2))*sin3omega_2 + &
                                        (2d0*nx4+nx2*(1d0-nx2)-ny2*(ny2+3d0*nz2))*sin5omega_2))
    Dmat(6,4) = 0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*(2d0*nx*cosOmega_2*(nx4-3d0*nx2*nz2+2d0*nz4+ny2*(nx2+5d0*nz2) + &
                                        4d0*nz2*(2d0*nx2-ny2)*cos1omega + (nx4+2d0*ny4+3d0*ny2*nz2+nx2*(3d0*ny2-nz2))*cos2omega)+&
                                        ny*nz*(-2d0*(nx2-nz2)*(3d0-5d0*nz2)*sinOmega_2 - &
                                        (nx4-nx2*ny2-2d0*ny4+(11d0*nx2+3d0*ny2)*nz2)*sin3omega_2 + &
                                        (nx4-nx2*ny2-2d0*ny4+nz2*(3d0*nx2-ny2))*sin5omega_2))
    Dmat(7,4) = -sqrt15*(nz2+(1d0-nz2)*cos1omega)*sinOmega_2**2*(nx4-ny4+(nx2-ny2)*(1d0+nz2)*cos1omega-4d0*nx*ny*nz*sin1omega)
    Dmat(1,5) = -0.25d0*ci*exp_omega_2*(expomega-1d0)*(2d0*nx3*(2d0*nx2+5d0*(1d0-nx2))*cosOmega_2 + &
                                                      nx*(4d0*nx4+5d0*(1d0-nx2)**2)*cos3omega_2 + &
                                                      6d0*nx3*ny2*cos5omega_2 + 3d0*nx*ny4*cos5omega_2 + 6d0*nx3*nz2*cos5omega_2 + &
                                                      6d0*nx*ny2*nz2*cos5omega_2 + 3d0*nx*nz4*cos5omega_2 + &
                                                      40d0*nx2*ny*nz*(ny2-nz2)*sinOmega_2**3 + &
                                                      8d0*nz*(ny5-ny*nz4)*(2d0+3d0*cos1omega)*sinOmega_2**3)
    Dmat(2,5) = -sqrt15*(nx2+(1d0-nx2)*cos1omega)*sinOmega_2**2*(ny4-nz4+(1d0+nx2)*(ny2-nz2)*cos1omega+4d0*nx*ny*nz*sin1omega)
    Dmat(3,5) = 0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*&
        (-2d0*nz*cosOmega_2*(2d0*ny4-3d0*ny2*nz2+nz4+nx2*(5d0*ny2+nz2)-4d0*ny2*(1d0-ny2-3d0*nz2)*cos1omega + &
         (2d0*nx4+3d0*nx2*(1d0-nx2)-ny2*nz2+nz4)*cos2omega) + &
         nx*ny*(2d0*(ny2-nz2)*(3d0*nx2-2d0*ny2+3d0*nz2)*sinOmega_2 + (2d0*nx4+nx2*(-3d0*ny2+nz2)-nz2*(11d0*ny2+nz2))*sin3omega_2 + &
         (-2d0*nx4-nx2*(1d0-nx2)+3d0*ny2*nz2+nz4)*sin5omega_2 ))
    Dmat(4,5) = 1d0/16d0*sqrt15*exp_3omega*( 2d0*(expomega-1d0)**2*(1d0+exp4omega)*nx5*nz - &
                                              2d0*ci*nx4*ny*(exp4omega*(2d0*cos2omega-1d0)-1d0) - &
                                    8d0*exp3omega*nx3*nz*(3d0*ny2-nz2+4d0*nz2*cos1omega + (1d0-nx2)*cos2omega)*sinOmega_2**2 + &
                8d0*exp3omega*nx*nz*(-3d0*ny4+ny2*nz2-2d0*nz4+ny2*(-8d0*nz2*cos1omega+(ny2+3d0*nz2)*cos2omega))*sinOmega_2**2 + &
                4d0*exp3omega*nx2*ny*(ny2+5d0*nz2-4d0*nz2*cos1omega+3d0*(1d0-nx2)*cos2omega)*sin1omega + &
                4d0*exp3omega*ny*(8d0*ny2*nz2*cos1omega+(ny2-nz2)*(ny2-2d0*nz2+ny2*cos2omega))*sin1omega )
    Dmat(5,5) = 1d0/16d0*exp_3omega*(60d0*exp3omega*nx2*(ny2-nz2)**2+3d0*(1d0-nx2)*(4d0*nx4+4d0*nx2*(1d0-nx2)+(ny2-nz2)**2) + &
                                    3d0*exp6omega*(1d0-nx2)*(4d0*nx4+4d0*nx2*(1d0-nx2)+(ny2-nz2)**2) + &
                            2d0*expomega*((-2d0*nx3+nx*ny2)**2-2d0*(2d0*nx4+9d0*nx2*ny2-8d0*ny4)*nz2 + nz4*(1d0+15d0*ny2-nz2)) + &
                            2d0*exp5omega*((-2d0*nx3+nx*ny2)**2-2d0*(2d0*nx4+9d0*nx2*ny2-8d0*ny4)*nz2+ nz4*(1d0+15d0*ny2-nz2)) + &
                            5d0*exp2omega*(4d0*nx4*(1d0-nx2)+(1d0-nx2)*(ny2-nz2)**2 - 4d0*nx2*(ny4-6d0*ny2*nz2+nz4)) + &
                            5d0*exp4omega*(4d0*nx4*(1d0-nx2)+(1d0-nx2)*(ny2-nz2)**2 - 4d0*nx2*(ny4-6d0*ny2*nz2+nz4)))
    Dmat(6,5) = -0.25d0*sinOmega_2*( 2*nz*(-15d0*nx2*ny2+nz4-5d0*nz2*(1d0-nz2))*cosOmega_2 + &
                                    nz*(-5d0*(2d0*nx4-5d0*nx2*ny2+2d0*ny4)-3d0*nz4+5d0*nz2*(1d0-nz2))*cos3omega_2 + &
                                    6d0*nx4*nz*cos5omega_2 -3d0*nx2*ny2*nz*cos5omega_2 + 6d0*ny4*nz*cos5omega_2 - &
                                    3d0*nx2*nz3*cos5omega_2 - &
                                    3d0*ny2*nz3*cos5omega_2 -3d0*nz5*cos5omega_2 + 30d0*nx*ny*(nx2-nz2)*(ny2-nz2)*sinOmega_2 + &
                                    5d0*nx*ny*(2d0*nx4+nx2*ny2+2*ny4-3d0*nz4+5d0*nz2*(1d0-nz2))*sin3omega_2 + &
                                    3d0*nx*ny*(2d0*nx4+5d0*nx2*ny2+2d0*ny4+5d0*nz4+5d0*nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(7,5) = -0.25d0*sinOmega_2*( 2d0*ny*(-ny4+5d0*ny2*nz2+5d0*nx2*(ny2+3d0*nz2))*cosOmega_2 + &
                                    ny*(10d0*nx4-5d0*nx2*ny2+3d0*ny4-5d0*(5d0*nx2+ny2)*nz2+10d0*nz4)*cos3omega_2 - &
                                    6d0*nx4*ny*cos5omega_2 + 3d0*nx2*ny*nz2*cos5omega_2 - 6d0*ny*nz4*cos5omega_2 + &
                                    3d0*nx2*ny3*cos5omega_2 + &
                                    3d0*ny3*nz2*cos5omega_2 + 3d0*ny5*cos5omega_2 - 30d0*nx*nz*(nx2-ny2)*(ny2-nz2)*sinOmega_2 + &
                                    5d0*nx*nz*(2d0*nx4+5d0*nx2*ny2+2*nz4-3d0*ny4+nz2*(1d0+4d0*ny2-nz2))*sin3omega_2 + &
                                    3d0*nx*nz*(2d0*nx4+5d0*nx2*ny2+2d0*nz4+5d0*ny4+5d0*nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(1,6) = 0.5d0*sinOmega_2*(2d0*ny3*(5d0*nx2+2d0*ny2+5d0*nz2)*cosOmega_2 + &
                                            ny*(4d0*ny4+5d0*(1d0-ny2)**2)*cos3omega_2 + &
                                            6d0*ny3*nx2*cos5omega_2 + 3d0*ny*nx4*cos5omega_2 + 6d0*nx2*ny*nz2*cos5omega_2 + &
                                            6d0*ny3*nz2*cos5omega_2 + 3d0*ny*nz4*cos5omega_2 + &
                                            30d0*ny2*nx*nz*(nz2-nx2)*sinOmega_2 + &
                                            5d0*(1d0-3d0*ny2)*nx*nz*(nz2-nx2)*sin3omega_2 + &
                                            3d0*nx*nz*(nx4-nz4)*sin5omega_2 )
    Dmat(2,6) = 0.25d0*sqrt15*sinOmega_2*( 2d0*nz*(2d0*nx4+3d0*nx2*ny2+nz4+nz2*(1d0-nz2))*cosOmega_2 + &
                                           nz*(-nx2*(ny2-7d0*nz2)+(1d0-nx2)*(2d0*ny2+nz2))*cos3omega_2 + &
                                           3d0*nx2*ny2*nz*cos5omega_2 + 2d0*ny4*nz*cos5omega_2 - nx2*nz3*cos5omega_2 + &
                                           3d0*ny2*nz3*cos5omega_2 + nz5*cos5omega_2 - &
                                           2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*ny*(nx2-nz2)*sinOmega_2 - &
                                           nx*ny*(-2d0*ny4-ny2*nz2+nz4+nx2*(3d0*ny2+11d0*nz2))*sin3omega_2 + &
                                           nx*ny*(-2d0*ny4-ny2*nz2+nz4-nx2*(ny2-3d0*nz2))*sin5omega_2 )
    Dmat(3,6) = sqrt15*(ny2+(1d0-ny2)*cos1omega)*sinOmega_2**2*(nx4-nz4+(1d0+ny2)*(nx2-nz2)*cos1omega-4d0*nx*ny*nz*sin1omega)
    Dmat(4,6) =-0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*(2d0*nx*cosOmega_2*(nx4-3d0*nx2*nz2+2d0*nz4+ny2*(nx2+5d0*nz2) + &
                                        4d0*nz2*(2d0*nx2-ny2)*cos1omega + (nx4+2d0*ny4+3d0*ny2*nz2+nx2*(3d0*ny2-nz2))*cos2omega)+&
                                        ny*nz*(2d0*(nx2-nz2)*(3d0-5d0*nz2)*sinOmega_2 + &
                                        (nx4-nx2*ny2-2d0*ny4+(11d0*nx2+3d0*ny2)*nz2)*sin3omega_2 - &
                                        (nx4-nx2*ny2-2d0*ny4+nz2*(3d0*nx2-ny2))*sin5omega_2))
    Dmat(5,6) = -0.25d0*sinOmega_2*(-2*nz*(-15d0*nx2*ny2+nz4-5d0*nz2*(1d0-nz2))*cosOmega_2 + &
                                    nz*(5d0*(2d0*nx4-5d0*nx2*ny2+2d0*ny4)+3d0*nz4-5d0*nz2*(1d0-nz2))*cos3omega_2 - &
                                    6d0*nx4*nz*cos5omega_2 +3d0*nx2*ny2*nz*cos5omega_2 - 6d0*ny4*nz*cos5omega_2 + &
                                    3d0*nx2*nz3*cos5omega_2 + &
                                    3d0*ny2*nz3*cos5omega_2 +3d0*nz5*cos5omega_2 + 30d0*nx*ny*(nx2-nz2)*(ny2-nz2)*sinOmega_2 + &
                                    5d0*nx*ny*(2d0*nx4+nx2*ny2+2*ny4-3d0*nz4+5d0*nz4*(1d0-nz2))*sin3omega_2 + &
                                    3d0*nx*ny*(2d0*nx4+5d0*nx2*ny2+2d0*ny4+5d0*nz4+5d0*nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(6,6) = 1d0/16d0*exp_3omega*(60d0*exp3omega*ny2*(1d0-ny2-2d0*nz2)**2 + &
                                     3d0*(1d0-ny2)*(nx4+nx2*(4d0*ny2-2d0*nz2)+(2d0*ny2+nz2)**2) + &
                                     3d0*exp6omega*(1d0-ny2)*(nx4+nx2*(4d0*ny2-2d0*nz2)+(2d0*ny2+nz2)**2) + &
                                2d0*expomega*(nx4*(ny2+16d0*nz2)+(-2d0*ny3+ny*nz2)**2-2d0*nx2*(2d0*ny4+9d0*ny2*nz2-8d0*nz4)) + &
                                2d0*exp5omega*(nx4*(ny2+16d0*nz2)+(-2d0*ny3+ny*nz2)**2-2d0*nx2*(2d0*ny4+9d0*ny2*nz2-8d0*nz4)) + &
                    5d0*(exp2omega+exp4omega)*(nx6-nx4*(4d0*ny2+nz2)+(-2d0*ny2*nz+nz3)**2+nx2*(4d0*ny4+24d0*ny2*nz2-nz4)) )
    Dmat(7,6) = -0.25d0*sinOmega_2*(2d0*nx*(nx4-5d0*nx2*(1d0-nx2)-15d0*ny2*nz2)*cosOmega_2 + &
                                    nx*(-3d0*nx4+5d0*nx2*(1d0-nx2)-5d0*(2d0*ny4-5d0*ny2*nz2+2d0*nz4))*cos3omega_2 + &
                                    cos5omega_2*(-3d0*nx5-3d0*nx3*ny2+6d0*nx*ny4-3d0*nx3*nz2-3d0*nx*ny2*nz2+6d0*nx*nz4) + &
                                    30d0*ny*nz*(nx2-ny2)*(nx2-nz2)*sinOmega_2 + &
                                    5d0*ny*nz*(-3d0*nx4+5d0*nx2*(1d0-nx2)+2d0*ny4+ny2*nz2+2d0*nz4)*sin3omega_2 + &
                                    3d0*ny*nz*(5d0*nx4+5d0*nx2*ny2+2d0*ny4+2d0*nz4+5d0*nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(1,7) =  0.5d0*sinOmega_2*(2d0*nz3*(2d0*nz2+5d0*(1d0-nz2))*cosOmega_2 + &
                                   nz*(4d0*nz4+5d0*(1d0-nz2)**2)*cos3omega_2 + &
                                   6d0*nz3*nx2*cos5omega_2 + 3d0*nz*nx4*cos5omega_2 + 6d0*nx2*nz*ny2*cos5omega_2 + &
                                   6d0*nz3*ny2*cos5omega_2 + 3d0*nz*ny4*cos5omega_2 + &
                                   30d0*nz2*nx*ny*(nx2-ny2)*sinOmega_2 + &
                                   5d0*(1d0-3d0*nz2)*nx*ny*(nx2-ny2)*sin3omega_2 + &
                                   3d0*nx*ny*(ny4-nx4)*sin5omega_2) 
    Dmat(2,7) =-0.125d0*ci*sqrt15*exp_omega_2*(expomega-1d0)*(2d0*ny*(2d0*nx4+nx2*ny2+ny4+(3d0*nx2+ny2)*nz2)*cosOmega_2 + &
                                                              ny*(-nx2*(nz2-7d0*ny2)+(1d0-nx2)*(2d0*nz2+ny2))*cos3omega_2 + &
                                                              3d0*nx2*ny*nz2*cos5omega_2 + 2d0*ny*nz4*cos5omega_2 + &
                                                            (-1d0)*nx2*ny3*cos5omega_2 + 3d0*ny3*nz2*cos5omega_2 + ny5*cos5omega_2 + &
                                                             2d0*nx*(2d0*nx2-3d0*(1d0-nx2))*nz*(nx2-ny2)*sinOmega_2 + &
                                                            nx*nz*((1d0-nx2)*(ny2-2d0*nz2)+nx2*(11d0*ny2+3d0*nz2))*sin3omega_2 + &
                                                             nx*nz*(-3d0*nx2*ny2-ny4+2d0*nz4+nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(3,7) = 1d0/16d0*sqrt15*exp_3omega*(-ci*(exp2omega-1d0)*(exp2omega+1d0)**2*nx5 - &
                                            (expomega-1d0)**2*(1d0-6d0*exp2omega+exp4omega)*nx4*ny*nz + &
                                8d0*exp3omega*nx2*ny*nz*(ny2-3d0*nz2-8d0*ny2*cos1omega+(3d0*ny2-nz2)*cos2omega)*sinOmega_2**2 - &
                            8d0*exp3omega*ny*nz*(ny2*(2d0*ny2-nz2+4d0*nz2*cos1omega)+nz2*(ny2+2d0*nz2)*cos2omega)*sinOmega_2**2 + &
                            4d0*exp3omega*nx3*(-3d0*ny2+nz2+8d0*ny2*cos1omega-(ny2-3d0*nz2)*cos2omega)*sin1omega + &
                            4d0*exp3omega*nx*(ny2*(2d0*ny2+5d0*nz2-4d0*nz2*cos1omega)+nz2*(3d0*ny2+2d0*nz2)*cos2omega)*sin1omega )
    Dmat(4,7) = -sqrt15*(nz2+(1d0-nz2)*cos1omega)*sinOmega_2**2*(nx4-ny4+(nx2-ny2)*(1d0+nz2)*cos1omega+4d0*nx*ny*nz*sin1omega)
    Dmat(5,7) = -0.25d0*sinOmega_2*(-2d0*ny*(-ny4+5d0*ny2*nz2+5d0*nx2*(ny2+3d0*nz2))*cosOmega_2 - &
                                    ny*(10d0*nx4-5d0*nx2*ny2+3d0*ny4-5d0*(5d0*nx2+ny2)*nz2+10d0*nz4)*cos3omega_2 + &
                                    6d0*nx4*ny*cos5omega_2 - 3d0*nx2*ny*nz2*cos5omega_2 + 6d0*ny*nz4*cos5omega_2 - &
                                    3d0*nx2*ny3*cos5omega_2 - &
                                    3d0*ny3*nz2*cos5omega_2 - 3d0*ny5*cos5omega_2 - 30d0*nx*nz*(nx2-ny2)*(ny2-nz2)*sinOmega_2 + &
                                    5d0*nx*nz*(2d0*nx4+5d0*nx2*ny2+2*nz4-3d0*ny4+nz2*(1d0+4d0*ny2-nz2))*sin3omega_2 + &
                                    3d0*nx*nz*(2d0*nx4+5d0*nx2*ny2+2d0*nz4+5d0*ny4+5d0*nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(6,7) = -0.25d0*sinOmega_2*(-2d0*nx*(nx4-5d0*nx2*(1d0-nx2)-15d0*ny2*nz2)*cosOmega_2 + &
                                    nx*(3d0*nx4-5d0*nx2*(1d0-nx2)+5d0*(2d0*ny4-5d0*ny2*nz2+2d0*nz4))*cos3omega_2 + &
                                    cos5omega_2*(3d0*nx5+3d0*nx3*ny2-6d0*nx*ny4+3d0*nx3*nz2+3d0*nx*ny2*nz2-6d0*nx*nz4) + &
                                    30d0*ny*nz*(nx2-ny2)*(nx2-nz2)*sinOmega_2 + &
                                    5d0*ny*nz*(-3d0*nx4+5d0*nx2*(1d0-nx2)+2d0*ny4+ny2*nz2+2d0*nz4)*sin3omega_2 + &
                                    3d0*ny*nz*(5d0*nx4+5d0*nx2*ny2+2d0*ny4+2d0*nz4+5d0*nz2*(1d0-nz2))*sin5omega_2 )
    Dmat(7,7) = 1d0/16d0*exp_3omega*(60d0*exp3omega*nz2*(1d0-2d0*ny2-nz2)**2 + &
                                     2d0*expomega*((nx4-18d0*nx2*ny2+ny4)*nz2+4d0*nz6+16d0*nx2*ny2*(1d0-nz2)-4d0*nz4*(1d0-nz2)) + &
                                     2d0*exp5omega*((nx4-18d0*nx2*ny2+ny4)*nz2+4d0*nz6+16d0*nx2*ny2*(1d0-nz2)-4d0*nz4*(1d0-nz2)) + &
                                     3d0*(1d0+exp6omega)*(1d0-nz2)*(4d0*nz4+4d0*nz2*(1d0-nz2)+(1d0-2d0*ny2-nz2)**2) + &
                        5d0*(exp2omega+exp4omega)*(-4d0*(nx4-6d0*nx2*ny2+ny4)*nz2+4d0*nz4*(1d0-nz2)+(1d0-nz2)*(1d0-2d0*ny2-nz2)**2) )

  else
    write(*,*) " This subroutine can only deal with J = 1/2, 1, 2, 3 !!! "
  end if

    ! using the same order of orbitals as Wannier90 
    ! l=0 : s
    ! l=1 : pz, px, py
    ! l=2 : dz2, dxz, dyz, dx2-y2, dxy
  !if (twoja1 == 3) then
  !  Dtmp(1,:) = Dmat(1,:)
  !  Dmat(1,:) = Dmat(3,:)
  !  Dmat(3,:) = Dtmp(1,:)
  !  Dtmp(:,1) = Dmat(:,1)
  !  Dmat(:,1) = Dmat(:,3)
  !  Dmat(:,3) = Dtmp(:,1)
  !elseif (twoja1 == 5) then 
  !  Dtmp(1,:) = Dmat(1,:)
  !  Dmat(1,:) = Dmat(5,:)
  !  Dmat(5,:) = Dtmp(1,:)
  !  Dtmp(:,1) = Dmat(:,1)
  !  Dmat(:,1) = Dmat(:,5)
  !  Dmat(:,5) = Dtmp(:,1)
  !  Dtmp(2,:) = Dmat(2,:)
  !  Dmat(2,:) = Dmat(3,:)
  !  Dmat(3,:) = Dtmp(2,:)
  !  Dtmp(:,2) = Dmat(:,2)
  !  Dmat(:,2) = Dmat(:,3)
  !  Dmat(:,3) = Dtmp(:,2)
  !endif 

end subroutine Dmatrix

end module bilbao
