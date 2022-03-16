module chrct 

    use lib_comms
    implicit none 

contains

subroutine irvsp_chrct(spinor, num_sym, SU2, &
                       num_litt_group, litt_group, &
                       num_bands, bot_band, top_band, ene_bands, &
                       dim_basis, num_basis, & 
                       coeffup_basis, coeffdn_basis, & 
                       chrct_sum, deg_set, numdeg_set, &
                       map_mat, map_vec, G_phase_pw)

    !!! This subroutine calculates the characters
    !!! from bot_band to top_band on one kpoints for all operations.

    logical,     intent(in)  :: spinor 
    integer,     intent(in)  :: num_sym
    complex(dp), intent(in)  :: SU2(2,2,MAXSYM)

    integer,     intent(in)  :: num_litt_group 
    integer,     intent(in)  :: litt_group(MAXSYM)

    integer,     intent(in)  :: bot_band, top_band 
    integer,     intent(in)  :: num_bands 
    real(dp),    intent(in)  :: ene_bands(num_bands)

    integer,     intent(in)  :: dim_basis, num_basis 
    complex(dp), intent(in)  :: coeffup_basis(dim_basis, num_bands)
    complex(dp), intent(in)  :: coeffdn_basis(dim_basis, num_bands)
    
    complex(dp), intent(out) :: chrct_sum(num_bands, MAXSYM)
    integer,     intent(out) :: deg_set(num_bands)
    integer,     intent(out) :: numdeg_set(num_bands)

    complex(dp), intent(in), optional :: map_mat(dim_basis,dim_basis,num_sym)
    integer,     intent(in), optional :: map_vec(dim_basis,num_sym)
    complex(dp), intent(in), optional :: G_phase_pw(dim_basis,num_sym)
    
    integer     :: i, j, j2, iband, iope, kope 
    integer     :: num_deg, cnt_deg_set 
    complex(dp) :: csum
    integer  :: m
    complex(dp) :: gphase

    complex(dp) :: coeffup_norm(dim_basis, num_bands)
    complex(dp) :: coeffdn_norm(dim_basis, num_bands)

    complex(dp) :: coeffup_tmp(dim_basis, MAXDG)
    complex(dp) :: coeffdn_tmp(dim_basis, MAXDG)


    chrct_sum = cmplx(7.777777_dp, 7.777777_dp, dp)

    ! normalization
    coeffup_norm = 0.d0
    coeffdn_norm = 0.d0 
    do j = 1, num_bands
        csum = cmplx_0
        csum = dot_product(coeffup_basis(:,j), coeffup_basis(:,j))
        if (spinor) csum = csum + dot_product(coeffdn_basis(:,j),coeffdn_basis(:,j))
        csum = zsqrt(csum)
        coeffup_norm(:,j) = coeffup_basis(:,j)/csum
        if (spinor) coeffdn_norm(:,j) = coeffdn_basis(:,j)/csum 
    enddo 
    !coeffup_norm = coeffup_basis
    !coeffdn_norm = coeffdn_basis 

    ! calculate characters
    deg_set     = 0
    cnt_deg_set = 0
    iband = bot_band 
    DO WHILE (iband <= top_band)

        ! find the degenerate sets
        num_deg = 1
        cnt_deg_set = cnt_deg_set + 1
        if (iband < top_band) then
            do while (ene_bands(iband+num_deg)-ene_bands(iband) <= TOLDG)
                num_deg = num_deg + 1
                if ((iband+num_deg) >  top_band) exit 
            enddo 
        endif 

        do i = iband, iband+num_deg-1 
            deg_set(i)    = cnt_deg_set
            numdeg_set(i) = num_deg
            if (save_kcount <= MAXKSAVE) save_numdeg(i,save_kcount) = num_deg 
        enddo 

        if (num_deg > MAXDG) then 
            write(6, '(A,I3)') "WARNING: num_deg is larger than MAXDG for the band:", iband
            iband = iband + num_deg 
            cycle 
        endif 

        ! for each operations in little group
        do iope = 1, num_litt_group 
            kope = litt_group(iope)
           
            ! rotation on spinor part
            if (spinor) then 
                j2 = 0
                coeffup_tmp = cmplx_0
                coeffdn_tmp = cmplx_0 
                do j = iband, iband + num_deg - 1
                    j2 = j2 + 1
                    do i = 1, num_basis 
                        coeffup_tmp(i,j2) = (SU2(1,1,kope)*coeffup_norm(i,j)) + &
                                            (SU2(1,2,kope)*coeffdn_norm(i,j))
                        coeffdn_tmp(i,j2) = (SU2(2,1,kope)*coeffup_norm(i,j)) + &
                                            (SU2(2,2,kope)*coeffdn_norm(i,j))
                    enddo 
                enddo 
                                 
            else
                j2 = 0
                coeffup_tmp = cmplx_0
                coeffdn_tmp = cmplx_0
                do j = iband, iband + num_deg - 1
                    j2 = j2 + 1
                    do i = 1, num_basis 
                        coeffup_tmp(i,j2) = coeffup_norm(i,j)
                    enddo 
                enddo 

            endif 

    
            ! rotation effect on the other part encoded in map_basis
            ! which should be provided in the main code
            chrct_sum(iband, kope) = cmplx_0
            do j = 0, num_deg - 1
                csum = cmplx_0
                if (present(map_mat)) then 
                    csum = dot_product(coeffup_norm(:,iband+j), &
                            matmul(map_mat(:,:,reorder(kope)), coeffup_tmp(:,1+j))) + &
                           dot_product(coeffdn_norm(:,iband+j), &
                            matmul(map_mat(:,:,reorder(kope)), coeffdn_tmp(:,1+j))) 
                else if (present(map_vec)) then 
                    do i = 1, num_basis  
                        m = map_vec(i,reorder(kope))
                        gphase = G_phase_pw(i,reorder(kope))
                        csum = csum + dconjg(coeffup_norm(m,iband+j))*coeffup_tmp(i,1+j)*gphase + &
                                      dconjg(coeffdn_norm(m,iband+j))*coeffdn_tmp(i,1+j)*gphase 
                    enddo 
                else
                    stop "The basis transformation under symmetry operations are not given"
                endif 
            
                chrct_sum(iband, kope) = chrct_sum(iband, kope) + csum 
            enddo 

            ! set all degenerate states the same character
            do j = iband+1, iband+num_deg-1
                chrct_sum(j, kope) = chrct_sum(iband, kope)
            enddo 

        enddo ! iope

        iband = iband + num_deg 

    ENDDO ! while  

end subroutine irvsp_chrct


subroutine irvsp_reps(spinor, num_sym, num_litt_group, litt_group, &
                      ind_ktype, timerev_k, &
                      num_bands, bot_band, top_band, &
                      chrct_set, deg_set, numdeg_set, kphase, &
                      kopeconjg_bb2input, &
                      reps_set, numreps_set)
    
    logical, intent(in)            :: spinor

    integer, intent(in)            :: num_sym
    integer, intent(in)            :: num_litt_group
    integer, intent(in)            :: litt_group(MAXSYM)

    integer, intent(in)            :: ind_ktype
    logical, intent(in)            :: timerev_k 
    
    integer, intent(in)            :: num_bands
    integer, intent(in)            :: bot_band, top_band
    
    complex(dp), intent(in)        :: chrct_set(num_bands, MAXSYM)
    
    ! number of bands of one degenerate set 
    integer,     intent(in)        :: deg_set(num_bands)
    integer,     intent(in)        :: numdeg_set(num_bands)
    complex(dp), intent(in)        :: kphase(num_sym)

    integer, intent(in)            :: kopeconjg_bb2input(num_sym)
    
    character(len=20), intent(out) :: reps_set(num_bands)

    ! number of representations of one degenerate set 
    integer,           intent(out) :: numreps_set(num_bands)
    
    logical     :: match 

    integer     :: iope, kope, cope, iband, ir1, ir2, nir, is  
    integer     :: i, j, i1, j1, i2, j2
    integer     :: maxdeg_table
    complex(dp) :: chrct_phase(num_bands, num_sym)
    real(dp)    :: diff  

    reps_set = '***'
    numreps_set = 0

    if (spinor) then 
        is  = 2
        ir1 = sirreps(ind_ktype) + 1
        ir2 = nirreps(ind_ktype)
    else 
        is  = 1
        ir1 = 1
        ir2 = sirreps(ind_ktype)
    endif 
    nir = ir2 - ir1 + 1

    chrct_phase = conjg(chrct_set(:,:))
    do iope = 1, num_litt_group
        kope = litt_group(iope)
        chrct_phase(:,kope) = chrct_phase(:,kope)*conjg(kphase(reorder(kope)))*conjg(kphase_difftau(kope))
    enddo 

    maxdeg_table = 0
    do i = ir1, ir2 
        if (maxdeg_table<nint(real(tableTraces(1,i,ind_ktype)))) then 
            maxdeg_table=nint(real(tableTraces(1,i,ind_ktype)))
        endif 
    enddo 

    iband = bot_band 
    do while (iband <= top_band) 

        if (nint(real(chrct_set(iband,1))) > MAXIRSUM*nint(real(maxdeg_table))) then 
            reps_set(iband) = '??'
            do i = iband + 1, iband + numdeg_set(iband) - 1
                reps_set(i) = '??'
            enddo 
            iband = iband + numdeg_set(iband) 
            cycle 
        endif 

        reps_set(iband) = '###'
    
        diff = 0.d0

        ! for one irrep
        do i1 = ir1, ir2 
            match = .true. 
            do iope = 1, num_litt_group 
                kope = litt_group(iope)
                cope = kopeconjg_bb2input(kope)
                if (.not.timerev_k) then 
                    diff = abs(chrct_phase(iband,kope)-tableTraces(cope,i1,ind_ktype))
                else
                    diff = abs(chrct_phase(iband,kope)-dconjg(tableTraces(cope,i1,ind_ktype)))
                endif ! timerev_k
                if (diff > (dble(nir)*TOLREP)) then 
                    match = .false.
                    exit
                endif 
            enddo 
            if (match) exit 
        enddo 
        if (match) then 
            if (.not.timerev_k) then 
                reps_set(iband) = Irrepsname(i1,ind_ktype)(is:is+3)
            else 
                reps_set(iband) = Irrepsname(i1,ind_ktype)(is:is)//'A'//Irrepsname(i1,ind_ktype)(is+1:is+3)
            endif ! timerev_k 
            numreps_set(iband) = 1
            if (save_kcount <= MAXKSAVE) then 
                save_chrct(iband, save_kcount) = i1 
                save_numrep(iband, save_kcount) = 1
            endif 
        endif 

        ! for two irreps
        if (.not.match) then 
    two:do i1 = ir1, ir2
        do j1 = i1, ir2 
            match = .true.
            do iope = 1, num_litt_group 
                kope = litt_group(iope)
                cope = kopeconjg_bb2input(kope)
                if (.not.timerev_k) then 
                    diff = abs(chrct_phase(iband,kope)-tableTraces(cope,i1,ind_ktype) &
                                                      -tableTraces(cope,j1,ind_ktype))
                else 
                    diff = abs(chrct_phase(iband,kope)-dconjg(tableTraces(cope,i1,ind_ktype)) &
                                                      -dconjg(tableTraces(cope,j1,ind_ktype)))

                endif ! timerev_k
                if (diff > (dble(nir)*TOLREP)) then 
                    match = .false.
                    exit
                endif 
            enddo 
            if (match) exit two 
        enddo 
        enddo two 
        if (match) then 
            if (.not.timerev_k) then 
                reps_set(iband)(1:5)  = Irrepsname(i1,ind_ktype)(is:is+3)
                reps_set(iband)(6:10) = Irrepsname(j1,ind_ktype)(is:is+3)
            else 
                reps_set(iband)(1:5)  = Irrepsname(i1,ind_ktype)(is:is)//'A'//Irrepsname(i1,ind_ktype)(is+1:is+3)
                reps_set(iband)(6:10) = Irrepsname(j1,ind_ktype)(is:is)//'A'//Irrepsname(j1,ind_ktype)(is+1:is+3)

            endif ! timerevk
            numreps_set(iband) = 2
            if (save_kcount <= MAXKSAVE) then 
                save_chrct(iband,  save_kcount)  = i1
                save_chrct(iband+1,save_kcount)  = j1
                save_numrep(iband, save_kcount)  = 2
                save_numrep(iband+1,save_kcount) = 2
            endif 
        endif 
        endif 

        ! for three irreps
        if (.not.match) then 
    thr:do i1 = ir1, ir2
        do j1 = i1, ir2 
        do i2 = j1, ir2  
            match = .true.
            do iope = 1, num_litt_group 
                kope = litt_group(iope)
                cope = kopeconjg_bb2input(kope)
                if (.not.timerev_k) then 
                    diff = abs(chrct_phase(iband,kope)-tableTraces(cope,i1,ind_ktype) &
                                                      -tableTraces(cope,j1,ind_ktype) &
                                                      -tableTraces(cope,i2,ind_ktype))
                else
                    diff = abs(chrct_phase(iband,kope)-dconjg(tableTraces(cope,i1,ind_ktype)) &
                                                      -dconjg(tableTraces(cope,j1,ind_ktype)) &
                                                      -dconjg(tableTraces(cope,i2,ind_ktype)))
                endif ! timerev_k
                if (diff > (dble(nir)*TOLREP)) then 
                    match = .false.
                    exit
                endif 
            enddo 
            if (match) exit thr
        enddo 
        enddo 
        enddo thr  
        if (match) then 
            if (.not.timerev_k) then 
                reps_set(iband)(1:5)  = Irrepsname(i1,ind_ktype)(is:is+3)
                reps_set(iband)(6:10) = Irrepsname(j1,ind_ktype)(is:is+3)
                reps_set(iband)(11:15)= Irrepsname(i2,ind_ktype)(is:is+3)
            else 
                reps_set(iband)(1:5)  = Irrepsname(i1,ind_ktype)(is:is)//'A'//Irrepsname(i1,ind_ktype)(is+1:is+3)
                reps_set(iband)(6:10) = Irrepsname(j1,ind_ktype)(is:is)//'A'//Irrepsname(j1,ind_ktype)(is+1:is+3)
                reps_set(iband)(11:15)= Irrepsname(i2,ind_ktype)(is:is)//'A'//Irrepsname(i2,ind_ktype)(is+1:is+3)

            endif ! timerev_k
            numreps_set(iband) = 3
            if (save_kcount <= MAXKSAVE) then 
                save_chrct(iband, save_kcount) = i1
                save_chrct(iband+1,save_kcount) = j1
                save_chrct(iband+2,save_kcount) = i2
                save_numrep(iband, save_kcount) = 3
                save_numrep(iband+1, save_kcount) = 3
                save_numrep(iband+2, save_kcount) = 3
            endif 
        endif 
        endif 

        ! for four irreps
        if (.not.match) then 
    fou:do i1 = ir1, ir2
        do j1 = i1, ir2 
        do i2 = j1, ir2
        do j2 = i2, ir2 
            match = .true.
            do iope = 1, num_litt_group 
                kope = litt_group(iope)
                cope = kopeconjg_bb2input(kope)
                if (.not.timerev_k) then 
                    diff = abs(chrct_phase(iband,kope)-tableTraces(cope,i1,ind_ktype) &
                                                      -tableTraces(cope,j1,ind_ktype) &
                                                      -tableTraces(cope,i2,ind_ktype) &
                                                      -tableTraces(cope,j2,ind_ktype))
                else 
                    diff = abs(chrct_phase(iband,kope)-dconjg(tableTraces(cope,i1,ind_ktype)) &
                                                      -dconjg(tableTraces(cope,j1,ind_ktype)) &
                                                      -dconjg(tableTraces(cope,i2,ind_ktype)) &
                                                      -dconjg(tableTraces(cope,j2,ind_ktype)))
                endif ! timerev_k
                if (diff > (dble(nir)*TOLREP)) then 
                    match = .false.
                    exit
                endif 
            enddo 
            if (match) exit fou 
        enddo 
        enddo 
        enddo 
        enddo fou
        if (match) then 
            if (.not.timerev_k) then 
                reps_set(iband)(1:5)  = Irrepsname(i1,ind_ktype)(is:is+3)
                reps_set(iband)(6:10) = Irrepsname(j1,ind_ktype)(is:is+3)
                reps_set(iband)(11:15)= Irrepsname(i2,ind_ktype)(is:is+3)
                reps_set(iband)(16:20)= Irrepsname(j2,ind_ktype)(is:is+3)
            else
                reps_set(iband)(1:5)  = Irrepsname(i1,ind_ktype)(is:is)//'A'//Irrepsname(i1,ind_ktype)(is+1:is+3)
                reps_set(iband)(6:10) = Irrepsname(j1,ind_ktype)(is:is)//'A'//Irrepsname(j1,ind_ktype)(is+1:is+3)
                reps_set(iband)(11:15)= Irrepsname(i2,ind_ktype)(is:is)//'A'//Irrepsname(i2,ind_ktype)(is+1:is+3)
                reps_set(iband)(16:20)= Irrepsname(j2,ind_ktype)(is:is)//'A'//Irrepsname(j2,ind_ktype)(is+1:is+3)

            endif 
            numreps_set(iband) = 4
            if (save_kcount <= MAXKSAVE) then 
                save_chrct(iband, save_kcount)   = i1
                save_chrct(iband+1, save_kcount) = j1
                save_chrct(iband+2, save_kcount) = i2 
                save_chrct(iband+3, save_kcount) = j2 
                save_numrep(iband, save_kcount) = 4
                save_numrep(iband+1, save_kcount) = 4
                save_numrep(iband+2, save_kcount) = 4
                save_numrep(iband+3, save_kcount) = 4
            endif 
        endif 
        endif 
        iband = iband + numdeg_set(iband)

    enddo ! iband 

end subroutine irvsp_reps


subroutine difftauphase(WK, num_littg, littg)

    real(dp), intent(in) :: WK(3)
    integer,  intent(in) :: num_littg
    integer,  intent(in) :: littg(MAXSYM)

    real(dp) :: ang
    integer  :: irot 

    kphase_difftau = 0.d0 
    do irot = 1, num_littg  
        ang = -2.d0*PI*( WK(1)*difftau_inputbilb(1,littg(irot)) &
                        +WK(2)*difftau_inputbilb(2,littg(irot)) &
                        +WK(3)*difftau_inputbilb(3,littg(irot)) )
        kphase_difftau(littg(irot)) = cmplx(dcos(ang), dsin(ang))
    enddo 
end subroutine difftauphase

end module chrct
