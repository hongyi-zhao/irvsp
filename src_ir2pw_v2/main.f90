program irvsp 

    use comms 
    use init
    use wave_data
    implicit none 

    integer            :: kkk

    ! command argument 
    integer            :: narg, iarg, lens, stat 
    character(len=100) :: arg, cmd 

    integer,     allocatable :: rot(:,:,:)
    real(dp),    allocatable :: SO3(:,:,:)
    complex(dp), allocatable :: SU2(:,:,:)
    complex(dp), allocatable :: kphase(:)
    complex(dp), allocatable :: Gphase(:,:)
    integer,     allocatable :: tilte_vec(:,:)

interface
subroutine irrep_bcs(sgn, num_sym, &
                      rot_input, tau_input, SO3_input, SU2_input, &
                      KKK, WK, kphase, &
                      num_bands, m, n, ene_bands, &
                      spinor, dim_basis, num_basis, &
                      coeffa, coeffb, &
                      G_phase_pw, rot_vec_pw, rot_mat_tb)

    integer, parameter :: dp = 8
    integer,     intent(in) :: sgn
    integer,     intent(in) :: num_sym 
    integer,     intent(in) :: rot_input(3,3,num_sym)
    real(dp),    intent(in) :: tau_input(3,num_sym)
    real(dp),    intent(in) :: SO3_input(3,3,num_sym)
    complex(dp), intent(in) :: SU2_input(2,2,num_sym)

    integer,     intent(in) :: KKK 
    real(dp),    intent(in) :: WK(3)
    complex(dp), intent(in) :: kphase(num_sym)

    integer,     intent(in) :: num_bands, m, n 
    real(dp),    intent(in) :: ene_bands(num_bands) 

    logical,     intent(in) :: spinor
    integer,     intent(in) :: dim_basis, num_basis  
    complex(dp), intent(in) :: coeffa(dim_basis, num_bands)
    complex(dp), intent(in) :: coeffb(dim_basis, num_bands)
    complex(dp), intent(in), optional :: G_phase_pw(dim_basis, num_sym)
    integer,     intent(in), optional :: rot_vec_pw(dim_basis, num_sym)
    complex(dp), intent(in), optional :: rot_mat_tb(dim_basis, dim_basis, num_sym)
end subroutine irrep_bcs
end interface 

    sgn = 0
    bot_band = 0
    top_band = 0
    call get_command(cmd)
    write(*,*) 'Current command : ', trim(cmd) 
    narg = command_argument_count()
    write(*,*) 'Argument count : ', narg 

    if (narg == 0) then
        write(*,"('Please input the correct space group number (sgn) &
                   by the command below:')")
        write(*,"('###$: irvsp -sg $sgn')")
        stop
    else
        iarg = 1
        do while (.true.)
            call get_command_argument(iarg, arg, lens, stat)
            if (len_trim(arg) == 0) exit 
            if (trim(arg) == '-sg') then 
                iarg = iarg + 1
                call get_command_argument(iarg, arg, lens, stat)
                read(arg, *) sgn 
            elseif (trim(arg) == '-nb') then 
                iarg = iarg + 1
                call get_command_argument(iarg, arg, lens, stat)
                read(arg, *) bot_band
                iarg = iarg + 1
                call get_command_argument(iarg, arg, lens, stat)
                read(arg, *) top_band
            else
                write(*,"('Please input the correct space group number (sgn) &
                           by the command below : ')")
                write(*, "('###$: irvsp -sg $sgn')")
                stop
            endif 

            iarg = iarg + 1
        enddo 
    endif 

    if (sgn == 0 .or. sgn > 230) then
        write(*,"('Please input the correct space group number (sgn) &
                   by the command below:')")
        write(*,"('###$: irvsp -sg $sgn')")
        stop
    endif 

    !
    call read_outcar() 

    call setarray()
    do kkk = 1, nspin*num_k

        allocate (rot(3,3,num_sym))                      ; rot=0
        allocate (SO3(3,3,num_sym))                      ; SO3=0.d0
        allocate (SU2(2,2,num_sym))                      ; SU2=0.d0 
        allocate (kphase(num_sym))                       ; kphase=0.d0
        allocate (tilte_vec(max_plane,num_sym))          ; tilte_vec=0
        allocate (Gphase(max_plane,num_sym))             ; gphase=0.d0 
    
        call read_wavecar(kkk)

        call pw_setup(WK, lattice, &
                      num_sym, det_input, angle_input, axis_input, tau_input, &
                      max_plane, ncnt, igall, &
                      rot, SO3, SU2, &
                      kphase, Gphase, tilte_vec)

        call irrep_bcs(sgn, num_sym, &
                        rot, tau_input, SO3, SU2,&
                        kkk, WK, kphase, &
                        num_bands, bot_band, top_band, EE, &
                        isSpinor, max_plane, ncnt, &
                        coeffa, coeffb, &
                        G_phase_pw=Gphase, rot_vec_pw=tilte_vec)
        
        deallocate(rot)
        deallocate(SO3)
        deallocate(SU2)
        deallocate(kphase)
        deallocate(tilte_vec)
        deallocate(Gphase)

    enddo 
    call downarray()

    WRITE(6,*) 
    WRITE(6,*) "*****************************************************"
    WRITE(6,*) 
    WRITE(6,*) "TOTAL END"


end program irvsp 
