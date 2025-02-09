module init 

    use comms
    implicit none 

    public :: setarray
    public :: downarray
    public :: read_outcar

contains


subroutine setarray()

    integer :: irot 

    allocate (EE(num_bands))                         ; EE=0.d0
    allocate (igall(3,max_plane))                    ; igall=0
    allocate (KV(3,max_plane))                       ; KV=0.d0
    allocate (coeffa(max_plane,num_bands))           ; coeffa=0.d0
    allocate (coeffb(max_plane,num_bands))           ; coeffb=0.d0

    allocate(det_input(num_sym))
    allocate(angle_input(num_sym))
    allocate(axis_input(3,num_sym))
    allocate(tau_input(3,num_sym))
    
    do irot = 1, num_sym 
      det_input(irot) = det_read(irot)
      angle_input(irot) = angle_read(irot)
      axis_input(:,irot) = axis_read(:,irot)
      tau_input(:,irot) = tau_read(:,irot)
    enddo 

end subroutine

subroutine downarray()
    if(allocated(EE    ))    deallocate (EE    )
    if(allocated(igall ))    deallocate (igall )
    if(allocated(KV    ))    deallocate (KV    )
    if(allocated(coeffa))    deallocate (coeffa)
    if(allocated(coeffb))    deallocate (coeffb)

end subroutine


subroutine read_outcar  

    character(len=120) :: chaps 
    character(len=5)   :: chtp5
    character(len=15)  :: chtp15
    character(len=35)  :: chtp35
    character(len=10)  :: title 

    integer            :: itmp
    real(dp)           :: rtmp(8)

    integer            :: i, j, i1, i2, irot 

    integer            :: ierr 

    integer, parameter :: outcar = 1001


    open(unit=outcar, file='OUTCAR', form='formatted', status='old')

    ! read symmetry opeartions
    isInv = .false.
    isSymmorphic = .true.
    do 
        read(outcar, "(A90)") chaps
        chtp5 = chaps(1:5)
        if (chtp5 == 'Space') exit 
    enddo
    read(outcar, *)
    do irot = 1, MAXSYM
        read(outcar, *, iostat=ierr) itmp, rtmp
        if (ierr /= 0) exit 
        num_sym = irot 
        det_read(irot) = rtmp(1)
        angle_read(irot) = -rtmp(2)
        axis_read(1,irot) = rtmp(3)
        axis_read(2,irot) = rtmp(4)
        axis_read(3,irot) = rtmp(5)
        tau_read(:,irot) = rtmp(6:8)
        if (abs(rtmp(1)+1.d0+rtmp(2)) < 1.d-5) isInv = .true.
        if (abs(rtmp(6))+abs(rtmp(7))+abs(rtmp(8)) > 1.d-5) isSymmorphic = .false.
    enddo 

    ! read number of kpoints, number of energy bands
    do 
        read(outcar, "(A120)") chaps
        chtp15 = chaps(4:11)
        if (chtp15 == 'k-points       ') exit 
    enddo 
    chtp15 = chaps(32:38)  ; read(chtp15, *) num_k
    chtp15 = chaps(105:115); read(chtp15, *) num_bands
    if (bot_band == 0) bot_band = 1
    if (top_band == 0) top_band = num_bands

    ! read title
    do
        read(outcar, "(A90)") chaps
        chtp15 = chaps(2:7)
        if (chtp15 == 'SYSTEM         ') exit 
    enddo 
    read(chaps, *) title, title, title 

    ! read nspin and soc
    do 
        read(outcar, "(A90)") chaps 
        chtp5 = chaps(4:8)
        if (chtp5 == 'ISPIN') exit
    enddo 
    read(chaps, "(A15I5)") chtp15, nspin 
    read(outcar, "(A90)") chaps
    read(outcar, "(A15L5)") chtp15, isSpinor
    isSpinPola = .false.
    if (nspin == 2) isSpinPola = .true. 
    isComplexWF = .true. 
    if (.not.isSpinor .and. isInv) isComplexWF = .false. 

    ! read lattice information, convert SO3 to rot
    do 
        read(outcar,"(A90)") chaps
        chtp5 = chaps(7:11)
        if (chtp5 == 'direc') exit 
    enddo 
    do i = 1, 3
        !read(outcar, *) br2(:,i), (br4(i,j), j=1,3)
        read(outcar,'(3X,F13.9,F13.9,F13.9,3X,F13.9,F13.9,F13.9)') br2(:,i),(br4(i,j),j=1,3)
    enddo 
    lattice = br2

    ! read max number of planewave
    do 
        read(outcar, "(A90)") chaps
        chtp15 = chaps(1:15)
        if (chtp15 == ' maximum number') exit
    enddo 
    read(chaps, "(A35I10)") chtp35, max_plane

    close(outcar)

    write(6, 529) title 

    if (     isSymmorphic) write(6,'(A19)',advance='NO') ' Symmorphic crystal'
    if (.not.isSymmorphic) write(6,'(A23)',advance='NO') ' Non-symmorphic crystal'

    if (     isInv) write(6,'(A24)') ' with inversion symmetry'
    if (.not.isInv) write(6,'(A27)') ' without inversion symmetry'

    if (     isComplexWF) write(6, '(A23)') ' Complex eigenfunctions'
    if (.not.isComplexWF) write(6, '(A20)') ' Real eigenfunctions'

    if (     isSpinor) write(6, '(A45)') ' Spin-orbit eigenfunctions (->time inversion)'
    if (.not.isSpinor) write(6, '(A29)') ' No spin-orbit eigenfunctions'

    if (     isSpinPola) write(6, '(A18)') 'Spin-polarization'
    if (.not.isSpinPola) write(6, '(A21)') ' No spin-polarization'

      WRITE(6,590)
      WRITE(6,592)
      WRITE(6,595) ((br2(I1,I2),I2=1,3),I1=1,3)
      WRITE(6,591)
      WRITE(6,596)  (br4( 1,I2),I2=1,3),' : g1/2pi'
      WRITE(6,596)  (br4( 2,I2),I2=1,3),' : g2/2pi'
      WRITE(6,596)  (br4( 3,I2),I2=1,3),' : g3/2pi'
   
      RETURN
!529  FORMAT(1X,A10,/,1X,A4,' lattice')
 529  FORMAT(1X,A10,/)
 590  FORMAT(//,' Transformations:',/, &
             ' Direct lattice vectors in Cartesian coord. system (BR2)')
 591  FORMAT(' Reciprocal lattice vectors in Cartesian coord. system (BR4)')
 592  FORMAT('        t1 :            t2 :            t3 : ')
 595  FORMAT(3(3F16.8,/))
 596  FORMAT((3F16.8),A9)

end subroutine read_outcar 

end module init
