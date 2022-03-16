! lib for irvsp and ir2tb
! jcgao95@gmail.com

subroutine irrep_bcs(sgn, num_sym, &
                      rot_input, tau_input, SO3_input, SU2_input, &
                      KKK, WK, kphase, &
                      num_bands, m, n, ene_bands, &
                      spinor, dim_basis, num_basis, &
                      coeffa, coeffb, &
                      G_phase_pw, rot_vec_pw, rot_mat_tb)

    use lib_comms
    use bilbao
    use chrct
    use dump
    implicit none 

    ! space group number
    ! should not be changed for different k
    integer,     intent(in) :: sgn

    ! number of space-group operations (module the integer lattice translations) 
    ! should not be changed for different k
    integer,     intent(in) :: num_sym 

    ! the rotation part of space-group operations with respect to primitive lattice vectors
    ! should not be changed for different k
    integer,     intent(in) :: rot_input(3,3,num_sym)

    ! the translation part of space-group operations with respect to primitive lattice vectors
    ! should not be changed for different k
    real(dp),    intent(in) :: tau_input(3,num_sym)

    ! the rotation part of space-group operations given in Cartesian coordinates
    ! should not be changed for different k
    real(dp),    intent(in) :: SO3_input(3,3,num_sym)

    ! the rotation part of space-group operations given in spin-1/2 space
    ! should not be changed for different k
    complex(dp), intent(in) :: SU2_input(2,2,num_sym)

    ! the sequential number of the given k-point
    integer,     intent(in) :: KKK

    ! the coordinate of the k-point with respect to primitive reciprocal lattice vectors
    real(dp),    intent(in) :: WK(3)

    ! the k-dependent phase factors due to the translation part of space-group operations
    complex(dp), intent(in) :: kphase(num_sym)

    ! the total number of bands
    integer,     intent(in) :: num_bands

    ! the IRs of the set of bands [m,n] are computed
    integer,     intent(in) :: m, n 

    ! the energy of bands at the k-point
    real(dp),    intent(in) :: ene_bands(num_bands) 

    ! set to .true. if the underlying electronic structure calculation has been performed with spinor wavefunctions
    logical,     intent(in) :: spinor

    ! the reserved number of the PW/TB basis
    ! if rot_vec_pw is given, dim_basis should be larger than the PW number of any k-points
    ! if rot_mat_tb is given, one should set dim_basis = num_basis
    integer,     intent(in) :: dim_basis

    ! the number of PW or orthogonal TB basis for the k-point 
    ! (note: the PW basis numbers for different k-points are usually different)
    integer,     intent(in) :: num_basis  

    ! the coefficient of spinor up part of wave functions at the given k-point 
    ! (note: coeffup_basis(1:num_basis,1:num_bands) is nonzero)
    complex(dp), intent(in) :: coeffa(dim_basis, num_bands)

    ! the coefficent of spinor down part of wave functions at the given k-point if spinor is .true.
    ! (note: coeffdn_basis(1:num_basis,1:num_bands) is nonzero)
    complex(dp), intent(in) :: coeffb(dim_basis, num_bands)
 
    ! the phase factor dependent by the PW vectors
    complex(dp), intent(in), optional :: G_phase_pw(dim_basis, num_sym)

    ! the transformation vectors of space-group operations, which send the jth PW to the j'th PW
    integer,     intent(in), optional :: rot_vec_pw(dim_basis, num_sym)

    ! the transformation matrices of space-group operations in orthogonal TB basis
    complex(dp), intent(in), optional :: rot_mat_tb(dim_basis, dim_basis, num_sym)

    ! the index of the type of kpoint in Bilbao
    integer                  :: ind_ktype
    character(len=3)         :: kname 
    ! the index of rotation which translate input k to Bilabo k
    integer                  :: ind_rot 

    ! reordered symmetry operations according to Bilbao
    integer                  :: rot_reorder(3,3,MAXSYM)
    integer                  :: invrot_reorder(3,3,MAXSYM)
    real(dp)                 :: convrot_reorder(3,3,MAXSYM)
    real(dp)                 :: tau_reorder(3,MAXSYM)
    complex(dp)              :: SU2_reorder(2,2,MAXSYM)
    real(dp)                 :: SO3_reorder(3,3,MAXSYM)
    ! SU2 matrix after get sign
    complex(dp)              :: SU2_reorder_getsignk(2,2,MAXSYM)

    ! little group of input k 
    integer                  :: littg_input(MAXSYM)
    integer                  :: num_littg_input 
    ! little group of bilbao k
    integer                  :: littg_bilbao(MAXSYM)
    integer                  :: num_littg_bilbao

    ! bilbao k
    real(dp)                 :: WK_getkid(3)

    ! Do we need time reversal to find the correct Bilbao k ?
    ! Useful for some space groups without inversion
    ! for example 121, K->KA
    logical                  :: timerev_k 

    ! kopeconjg_bb2inp(2)=9 means the conjg operator of ope2 of WK
    ! is ope9 in Bilbao 
    integer                  :: kopeconjg_bb2inp(MAXSYM)

    ! calculated character tables, without outer k phase
    complex(dp)              :: chrct_set(num_bands, MAXSYM)

    ! degenerate informations of input bands
    integer                  :: deg_set(num_bands)
    integer                  :: numdeg_set(num_bands)

    ! representation names
    character(len=20)        :: reps_set(num_bands)
    integer                  :: numreps_set(num_bands)

    integer  :: i, j

    ! useful when the first time the lib is called
    ! read informations from bilbao table
    ! reorder the input symmetry operations
    ! initialize savedata part 
    if (.not.allocated(reorder)) then 
        allocate(reorder(MAXSYM))
        call bilbao_read(sgn)
        call bilbao_reorder(num_sym, &
                            rot_input, tau_input, SU2_input, SO3_input, &
                            rot_reorder, tau_reorder, SU2_reorder, SO3_reorder)
        !call getsign(num_doub_sym, rot_bilbao, SU2_bilbao, &
        !                           rot_reorder, SU2_reorder)
        
        
        
        invrot_reorder  = 0
        convrot_reorder = 0.d0 
        do i = 1, num_sym 
            call invmati(rot_reorder(:,:,i), invrot_reorder(:,:,i))
            convrot_reorder(:,:,i) = matmul(Kc2p, &
                                     matmul(dble(rot_reorder(:,:,i)),p2cR))
        enddo
        call dump_opes(MAXSYM, num_sym, &
                       rot_reorder, invrot_reorder, tau_reorder, & 
                       SO3_reorder, SU2_reorder) 

        ! initialize
        save_kcount = 0
        allocate(save_chrct(num_bands, MAXKSAVE))
        allocate(save_numdeg(num_bands, MAXKSAVE))
        allocate(save_numrep(num_bands, MAXKSAVE))
        allocate(save_ktype(MAXKSAVE))
        save_chrct = -7
        save_numdeg = 0
        save_numrep = 0
        save_ktype = 0
    endif 

    save_kcount = save_kcount + 1

    ! get little group of input k
    call kgroup(MAXSYM, num_sym, invrot_reorder, &
                WK, littg_input, num_littg_input)
    
    call difftauphase(WK, num_littg_input, littg_input)

    ! get bilbao k and little group of bilbao k
    WK_getkid = WK 
    call bilbao_getkid(WK_getkid, invrot_reorder, &
                       ind_ktype, kname, ind_rot, &
                       timerev_k)
    if (save_kcount <= MAXKSAVE) then 
        save_ktype(save_kcount) = ind_ktype 
        if (timerev_k) save_ktype(save_kcount) = -ind_ktype 
    endif 
    call kgroup(MAXSYM, num_doub_sym, invrot_bilbao, &
                WK_getkid, littg_bilbao, num_littg_bilbao)

    ! get the conjugate relation of little group of input k and bilbao k
    call bilbao_getconjg(ind_rot, &
                         num_littg_input, littg_input, &
                         num_littg_bilbao, littg_bilbao, &
                         kopeconjg_bb2inp)

    ! fix the sign problem of SU2 matrix
    SU2_reorder_getsignk = SU2_reorder 
    call getsign_littg(num_doub_sym, rot_bilbao, SU2_bilbao, &
                                num_littg_bilbao, littg_bilbao, &
                       num_sym, rot_reorder, SU2_reorder_getsignk, &
                                num_littg_input, littg_input, &
                                kopeconjg_bb2inp)

    ! dump to file
    call dump_kinformation(KKK, WK, kname, num_littg_input, sgn, timerev_k)
    call dump_bilbaochrct(ind_ktype, timerev_k,6)
    call dump_littg(num_littg_input, littg_bilbao, &
                    rot_reorder, tau_reorder, WK)

    ! calculate characters
    if (present(rot_mat_tb)) then 
        call irvsp_chrct(spinor, num_sym, SU2_reorder_getsignk, &
                         num_littg_input, littg_input, &
                         num_bands, m, n, ene_bands, &
                         dim_basis, num_basis, &
                         coeffa, coeffb, &
                         chrct_set, deg_set, numdeg_set, &
                         map_mat=rot_mat_tb) 
    else if (present(rot_vec_pw)) then 
        call irvsp_chrct(spinor, num_sym, SU2_reorder_getsignk, &
                         num_littg_input, littg_input, &
                         num_bands, m, n, ene_bands, &
                         dim_basis, num_basis, &
                         coeffa, coeffb, &
                         chrct_set, deg_set, numdeg_set, &
                         map_vec=rot_vec_pw, G_phase_pw=G_phase_pw) 
    else
        stop "The basis transformation under symmetry operations are not given"
    endif 

    ! compare to bilbao characters to get representation names
    call irvsp_reps(spinor, num_sym, num_littg_input, littg_input, &
                    ind_ktype, timerev_k, &
                    num_bands, m, n, &
                    chrct_set, deg_set, numdeg_set, kphase, &
                    kopeconjg_bb2inp, &
                    reps_set, numreps_set)

    call dump_reps(num_bands, m, n, &
                   numdeg_set, ene_bands, &
                   num_littg_input, littg_input, littg_bilbao, &
                   kopeconjg_bb2inp, &
                   chrct_set, reps_set, numreps_set)

    call dump_save(sgn, m, n)

end subroutine irrep_bcs


! setup for plane wave part of irrep_bcs
! prepare kphase, Gphase_pw, rot_vec_pw
subroutine pw_setup(WK, lattice, &
                    num_sym, det, angle, axis, tau, &
                    dim_basis, num_basis, Gvec, &
                    rot, SO3, SU2, &
                    kphase, Gphase_pw, rot_vec_pw)

    use lib_comms
    use bilbao
    implicit none 

    ! the coordinates of the k-point with respect to primitive reciprocal lattice vectors
    real(dp),    intent(in)  :: WK(3)

    ! the primitive lattice vectors in Cartesian coordinates
    ! |  t1x, t2x, t3x  |
    ! |  t1y, t2y, t3y  |
    ! |  t1z, t2z, t3z  | 
    real(dp),    intent(in)  :: lattice(3,3)

    ! the number of space-group operations 
    integer,     intent(in)  :: num_sym

    ! the determination of the rotation part of space-group operations
    real(dp),    intent(in)  :: det(num_sym)
    
    ! the rotation angles of space-group operations
    real(dp),    intent(in)  :: angle(num_sym)

    ! the rotation axis of space-group operations in Cartesian coordinates
    real(dp),    intent(in)  :: axis(3, num_sym)

    ! the translation part of space-group operations with respect to primitive lattice vectors
    real(dp),    intent(in)  :: tau(3,num_sym)

    ! the reserved number of the PW_basis (dim_basis >= num_basis)
    integer,     intent(in)  :: dim_basis 

    ! the number of PW for the k-point (note: num_basis for different k-points are usually different)
    integer,     intent(in)  :: num_basis 

    ! the plane-wave G-vector with respected to reciprocal lattice vectors
    integer,     intent(in)  :: Gvec(3, dim_basis)

    ! the rotation part of space-group operations with respect to primitive lattice vectors
    integer,     intent(out) :: rot(3,3,num_sym)

    ! the rotation part of space-group operations in Cartesian coordinates
    real(dp),    intent(out) :: SO3(3,3,num_sym)

    ! the rotation part of space-group operations in spin-1/2 space
    complex(dp), intent(out) :: SU2(2,2,num_sym)

    ! the k-dependent phase factors due to the translation part of space-group operations
    complex(dp), intent(out) :: kphase(num_sym)

    ! the phase factor dependent by the PW vectors
    complex(dp), intent(out) :: Gphase_pw(dim_basis, num_sym)

    ! the transformation vectors of Rs, which send the jth PW to the j'th PW
    integer,     intent(out) :: rot_vec_pw(dim_basis, num_sym)


    ! parameter used inside the library 
    integer :: i, j, irot

    real(dp) :: br2(3,3), br4(3,3)
    
    complex(dp) :: so3tmp(3,3), su2tmp(2,2)

    integer  :: num_litt_group 
    integer  :: litt_group(num_sym)

    integer  :: ind, ikv 
    real(dp) :: RWK(3), RKV(3)
    real(dp) :: diff, ang
    real(dp) :: KV(3,dim_basis)

    integer  :: invrot(3,3,num_sym)


    ! useful when the first time the library is called 
    if (.not.allocated(rot_save)) then 
        allocate(SO3_save(3,3,num_sym))
        allocate(SU2_save(2,2,num_sym))
        allocate(rot_save(3,3,num_sym))
        allocate(invrot_save(3,3,num_sym))

        br2 = lattice 
        call invreal33(br2, br4) 
        
        do irot = 1, num_sym 
            call Dmatrix(rnx=axis(1,irot),rny=axis(2,irot),rnz=axis(3,irot),degree=angle(irot),twoja1=3,Dmat=so3tmp)
            call Dmatrix(rnx=axis(1,irot),rny=axis(2,irot),rnz=axis(3,irot),degree=angle(irot),twoja1=2,Dmat=su2tmp)
            SU2_save(:,:,irot) = su2tmp
            if (det(irot) <1.d-2) then 
                SO3_save(:,:,irot) = real(-so3tmp) 
            else 
                SO3_save(:,:,irot) = real(so3tmp)
            endif 
            rot_save(:,:,irot) = nint(matmul(matmul(br4, SO3_save(:,:,irot)), br2))
            call invmati(rot_save(:,:,irot), invrot_save(:,:,irot))
        enddo 
    endif

    
    SO3 = SO3_save
    SU2 = SU2_save
    rot = rot_save 
    invrot = invrot_save
    
    ! get little group 
    num_litt_group = 0
    litt_group = 0
    do irot = 1, num_sym 
        do j = 1, 3
            RWK(j) = dot_product(WK(:), invrot(:,j,irot)) - WK(j)
        enddo 
        diff = dabs(nint(RWK(1))-RWK(1)) &
              +dabs(nint(RWK(2))-RWK(2)) &
              +dabs(nint(RWK(3))-RWK(3))
        if (diff < 0.1e-4) then 
            num_litt_group = num_litt_group + 1
            litt_group(num_litt_group) = irot 
        endif 
    enddo 

    ! get kphase
    kphase = 0.d0 
    do irot = 1, num_litt_group 
        ang = -2.d0*PI*( WK(1)*tau(1,litt_group(irot)) &
                          +WK(2)*tau(2,litt_group(irot)) &
                          +WK(3)*tau(3,litt_group(irot)) )
        kphase(litt_group(irot)) = cmplx(dcos(ang), dsin(ang))
    enddo 

    ! get Gphase and rot_vec_pw
    KV = 0.d0 
    do ikv = 1, num_basis 
        KV(:,ikv) = dble(Gvec(:,ikv)) + WK(:)
    enddo 


    Gphase_pw = 0.d0
    rot_vec_pw = 0
    do irot = 1, num_litt_group 
        do ikv = 1, num_basis 
            
            do j = 1, 3
                RKV(j) = KV(1,ikv)*dble(invrot(1,j,litt_group(irot))) &
                        +KV(2,ikv)*dble(invrot(2,j,litt_group(irot))) &
                        +KV(3,ikv)*dble(invrot(3,j,litt_group(irot)))
            enddo 

            ind = 1
            do while ( (abs(RKV(1)-KV(1,ind)) &
                       +abs(RKV(2)-KV(2,ind)) &
                       +abs(RKV(3)-KV(3,ind))) > 1.d-3 .and. (ind <= num_basis) )
                ind = ind + 1
            enddo 

            if (ind == num_basis + 1) then 
                write(6,*) "cannot find (k+G)inv(Ri)"
                stop
            endif 

            ang = 0.d0 
            do j = 1, 3
                ang = ang - 2.d0*PI*tau(j, litt_group(irot))*(KV(j,ind)-WK(j))
            enddo 

            Gphase_pw(ikv, litt_group(irot))  = cmplx(dcos(ang), dsin(ang))
            rot_vec_pw(ikv, litt_group(irot)) = ind 

        enddo 
    enddo     


end subroutine pw_setup


! setup for Tight-binding part of irrep_bcs
! prepare kphase, rot_mat_tb
subroutine tb_setup(WK, lattice, &
                    num_sym, det, angle, axis, tau, &
                    num_atom, atom_position, &
                    dim_basis, num_basis, angularmom, orbt, &
                    rot, SO3, SU2, &
                    kphase, rot_mat_tb)

    use lib_comms
    use bilbao
    implicit none 

    ! the coordinates of the k-point with respected to primitive reciprocal lattice vectors
    real(dp),    intent(in)  :: WK(3)

    ! the primitive lattice vectors in Cartesian coordinates
    ! |  a1, a2, a3  |
    ! |  b1, b2, b3  |
    ! |  c1, c2, c3  | 
    real(dp),    intent(in)  :: lattice(3,3)

    ! the number of space-group operations 
    integer,     intent(in)  :: num_sym

    ! the determination of rotation part of space-group operations
    real(dp),    intent(in)  :: det(num_sym)
    
    ! the rotation angles of space-group operations
    real(dp),    intent(in)  :: angle(num_sym)

    ! the rotation axis of space-group operations in Cartesian coordinates
    real(dp),    intent(in)  :: axis(3, num_sym)

    ! the translation part of space-group operations with respect to primitive lattice vectors
    real(dp),    intent(in)  :: tau(3,num_sym)

    ! the number of atoms in the TB Hamiltonian
    integer,     intent(in)  :: num_atom

    ! the coordinates of atoms with respect to primitive lattice vectors
    real(dp),    intent(in)  :: atom_position(3, num_atom)

    ! the reserved number of the TB basis (dim_basis = num_basis)
    integer,     intent(in)  :: dim_basis 

    ! the number of orthogonal local orbitals for the k-point
    integer,     intent(in)  :: num_basis 

    ! the local orbital information on each atom. Detailed explainations can be found in Table S3
    ! could be 1,3,4,5,6,7,8,9 for s,p,s+p,d,s+d,f,p+d,s+p+d 
    integer,     intent(in)  :: angularmom(num_atom)

    ! the convention of the local orbitals on each atom
    ! if orbt=1, local orbitals are in the order of Table S3
    ! if orbt=2, local orbitals are in the order as implemented in Wannier90
    integer,     intent(in)  :: orbt

    ! the rotation part of space-group operations with respect to primitive lattice vectors
    integer,     intent(out) :: rot(3,3,num_sym)

    ! the rotation part of space-group operations in Cartesian coordinates
    real(dp),    intent(out) :: SO3(3,3,num_sym)

    ! the rotation part of space-group operations in spin-1/2 space
    complex(dp), intent(out) :: SU2(2,2,num_sym)

    ! the k-dependent phase factors due to the translation part
    complex(dp), intent(out) :: kphase(num_sym)

    ! the transformation matrices of Rs in the orthogonal TB basis
    complex(dp), intent(out) :: rot_mat_tb(dim_basis, dim_basis, num_sym)


    ! parameter used inside the library 
    integer :: i, j, irot, iatom, jatom 

    real(dp) :: br2(3,3), br4(3,3)
    
    complex(dp) :: rotmt(3,3), crotmt(3,3), srotmt(2,2), protmt(3,3), drotmt(5,5), frotmt(7,7)
    complex(dp) :: cmat3(3,3), cmat5(5,5), cmat7(7,7)

    integer  :: num_litt_group 
    integer  :: litt_group(num_sym)

    integer  :: ind, ikv 
    real(dp) :: RWK(3), RKV(3), RPOS(3), diffr(3)
    real(dp) :: diff, ang
    real(dp) :: KV(3,dim_basis)

    integer  :: invrot(3,3,num_sym)

    integer :: norb_i, norb_j, startorb_i, startorb_j
    complex(dp) :: phk 


    ! useful when the first time the library is called 
    if (.not.allocated(rot_save)) then 
        allocate(SO3_save(3,3,num_sym))
        allocate(SU2_save(2,2,num_sym))
        allocate(rot_save(3,3,num_sym))
        allocate(invrot_save(3,3,num_sym))
        allocate(rot_orb(mix2l, mix2l, num_atom, num_sym))
        allocate(startorb(num_atom))
        allocate(rot_atom(num_atom, num_sym))
        allocate(rot_phivec(3,num_atom,num_sym))

        br2 = lattice 
        call invreal33(br2, br4)  
        
        do irot = 1, num_sym 
            call Dmatrix(rnx=axis(1,irot),rny=axis(2,irot),rnz=axis(3,irot),degree=angle(irot),twoja1=3,Dmat=crotmt)
            if (det(irot) <1.d-2) then 
                SO3_save(:,:,irot) = real(-crotmt) 
            else 
                SO3_save(:,:,irot) = real(crotmt)
            endif 
            rot_save(:,:,irot) = nint(matmul(matmul(br4, SO3_save(:,:,irot)), br2))
            call invmati(rot_save(:,:,irot), invrot_save(:,:,irot))
            protmt = crotmt
            if (orbt == 2) then
                protmt = matmul(protmt, pwann)
                cmat3  = transpose(pwann)
                protmt = matmul(cmat3, protmt) 
            endif 

            call Dmatrix(rnx=axis(1,irot),rny=axis(2,irot),rnz=axis(3,irot),degree=angle(irot),twoja1=5,Dmat=drotmt)
            if (orbt == 2) then 
                drotmt = matmul(drotmt, dwann)
                cmat5 = transpose(dwann)
                drotmt = matmul(cmat5, drotmt)
            endif 

            call Dmatrix(rnx=axis(1,irot),rny=axis(2,irot),rnz=axis(3,irot),degree=angle(irot),twoja1=2,Dmat=frotmt)
            if (orbt == 2) then 
                frotmt = matmul(frotmt, fwann)
                cmat7 = transpose(fwann)
                frotmt = matmul(cmat7, frotmt)
            endif

            call Dmatrix(rnx=axis(1,irot),rny=axis(2,irot),rnz=axis(3,irot),degree=angle(irot),twoja1=2,Dmat=srotmt)
            SU2_save(:,:,irot) = srotmt
            
            do iatom = 1, num_atom
                if     (angularmom(iatom) == 1) then 
                    rot_orb(1,1,iatom,irot) = 1.d0 
                elseif (angularmom(iatom) == 3) then 
                    rot_orb(1:3,1:3,iatom,irot) = real(protmt, dp)
                    if (det(irot) < 1.d-2) rot_orb(1:3,1:3,iatom,irot) = -real(protmt, dp)
                elseif (angularmom(iatom) == 4) then 
                    rot_orb(1,1,iatom,irot) = 1.d0 
                    rot_orb(2:4,2:4,iatom,irot) = real(protmt, dp)
                    if (det(irot) < 1.d-2) rot_orb(2:4,2:4,iatom,irot) = -real(protmt, dp)
                elseif (angularmom(iatom) == 5) then 
                    rot_orb(1:5,1:5,iatom,irot) = real(drotmt, dp)
                elseif (angularmom(iatom) == 6) then 
                    rot_orb(1,1,iatom,irot) = 1.d0
                    rot_orb(2:6,2:6,iatom,irot) = real(drotmt, dp)
                elseif (angularmom(iatom) == 7) then 
                    rot_orb(1:7,1:7,iatom,irot) = real(frotmt, dp)
                    if (det(irot) < 1.d-2) rot_orb(1:7,1:7,iatom,irot) = -real(frotmt, dp)
                elseif (angularmom(iatom) == 8) then 
                    rot_orb(1:3,1:3,iatom,irot) = real(protmt, dp)
                    rot_orb(4:8,4:8,iatom,irot) = real(drotmt, dp)
                    if (det(irot) < 1.d-2) rot_orb(1:3,1:3,iatom,irot) = -real(protmt, dp)
                elseif (angularmom(iatom) == 9) then 
                    rot_orb(1,1,iatom,irot) = 1.d0 
                    rot_orb(2:4,2:4,iatom,irot) = real(protmt, dp)
                    rot_orb(5:9,5:9,iatom,irot) = real(drotmt, dp)
                    if (det(irot) < 1.d-2) rot_orb(2:4,2:4,iatom,irot) = -real(protmt, dp)
                elseif (angularmom(iatom) ==12) then 
                    rot_orb(1:5,1:5,iatom,irot) = real(drotmt, dp)
                    rot_orb(6:12,6:12,iatom,irot) = real(frotmt, dp)
                    if (det(irot) < 1.d-2) rot_orb(6:12,6:12,iatom,irot) = -real(frotmt, dp)
                else
                    write(*,*) "Please set the correct angularmom for each atom"
                    stop
                endif 
            enddo 
        enddo

        startorb = 0
        do iatom = 1, num_atom
            if (iatom /= 1) startorb(iatom) = startorb(iatom-1) + angularmom(iatom-1)
        enddo 

        rot_atom = 0
        do irot = 1, num_sym
            do iatom = 1, num_atom 
                RPOS(:) = matmul(rot_save(:,:,irot), atom_position(:,iatom)) + tau(:,irot)
                do jatom = 1, num_atom
                    diffr = RPOS(:) - atom_position(:,jatom)
                    if (abs(diffr(1)-nint(diffr(1))) .lt. 0.1d-2 .and. &
                        abs(diffr(2)-nint(diffr(2))) .lt. 0.1d-2 .and. &
                        abs(diffr(3)-nint(diffr(3))) .lt. 0.1d-2 ) then 
                        exit
                    endif 
                enddo
                if (jatom == num_atom + 1) then 
                    write(*,*) "rotation on atoms error", iatom
                    stop 
                endif 
                rot_atom(iatom, irot) = jatom
                rot_phivec(:,iatom,irot) = -matmul(invrot_save(:,:,irot),tau(:,irot)) &
                                           +matmul(invrot_save(:,:,irot),nint(diffr(:))) 
            enddo 
        enddo 

    endif

    SO3 = SO3_save
    SU2 = SU2_save
    rot = rot_save 
    invrot = invrot_save
    
    ! get little group 
    num_litt_group = 0
    litt_group = 0
    do irot = 1, num_sym 
        do j = 1, 3
            RWK(j) = dot_product(WK(:), invrot(:,j,irot)) - WK(j)
        enddo 
        diff = dabs(nint(RWK(1))-RWK(1)) &
              +dabs(nint(RWK(2))-RWK(2)) &
              +dabs(nint(RWK(3))-RWK(3))
        if (diff < 0.1e-4) then 
            num_litt_group = num_litt_group + 1
            litt_group(num_litt_group) = irot 
        endif 
    enddo 

    ! get kphase
    kphase = 0.d0 
    do irot = 1, num_litt_group
        RWK = matmul(WK(:), invrot(:,:,litt_group(irot))) 
        ang = -2.d0*PI*( RWK(1)*tau(1,litt_group(irot)) &
                        +RWK(2)*tau(2,litt_group(irot)) &
                        +RWK(3)*tau(3,litt_group(irot)) )
        kphase(litt_group(irot)) = cmplx(dcos(ang), dsin(ang))
    enddo 

    ! get rot_mat_tb
    rot_mat_tb = 0.d0
    do irot = 1, num_litt_group
        do iatom = 1, num_atom 
            jatom = rot_atom(iatom, litt_group(irot))
            norb_i = angularmom(iatom)
            norb_j = angularmom(jatom)
            startorb_i = startorb(iatom)
            startorb_j = startorb(jatom)
            if (norb_i /= norb_j) then 
                write(*,*) "error symmetry"
                stop
            endif 
            phk = -dot_product(WK(:),rot_phivec(:,iatom,litt_group(irot)))
            rot_mat_tb((/1:norb_j/)+startorb_j, (/1:norb_i/)+startorb_i,litt_group(irot)) = &
                rot_orb(1:norb_j, 1:norb_i, iatom, litt_group(irot))*exp(2.d0*PI*cmplx_i*phk)
        enddo 
    enddo 

end subroutine tb_setup
