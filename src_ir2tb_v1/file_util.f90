    !=========+=========+=========+=========+=========+=========+=========+
    subroutine get_key_para_int(keyword,nfile,intn)
    character(15),intent(inout):: keyword
    integer      , intent(in)  :: nfile 
    integer      , intent(out) :: intn
    character(15)  :: keyword2

    character(50)  :: chartmp

    integer       :: ierr
    integer       :: i,j
       rewind(nfile)
       do while (.true.)
           read(nfile, "(A)", iostat=ierr) chartmp
           if (ierr /= 0)  exit
           do i=1,50
             if(chartmp(i:i)==":".or.chartmp(i:i)=="=") exit
           enddo
           if (i==51) cycle
           read(chartmp(2:i-1), *)  keyword2
           if(keyword2==keyword) then
          !write(6,*) chartmp
           read(chartmp(i+1:50), *)  intn
           exit
           endif
       enddo ! over while (.true.) loop

    end subroutine get_key_para_int
    !=========+=========+=========+=========+=========+=========+=========+

    !=========+=========+=========+=========+=========+=========+=========+
    subroutine get_key_para_cht(keyword,nfile,chat)
    character(15),intent(inout):: keyword
    integer      , intent(in)  :: nfile 
    character(15), intent(inout) :: chat
    character(15)  :: keyword2

    character(50)  :: chartmp

    integer       :: ierr
    integer       :: i,j
       rewind(nfile)
       do while (.true.)
           read(nfile, "(A)", iostat=ierr) chartmp
           if (ierr /= 0)  exit
           do i=1,50
             if(chartmp(i:i)==":".or.chartmp(i:i)=="=") exit
           enddo
           if (i==51) cycle
           read(chartmp(2:i-1), *)  keyword2
           if(keyword2==keyword) then
          !write(6,*) chartmp
           read(chartmp(i+1:50), *)  chat
           exit
           endif
       enddo ! over while (.true.) loop

    end subroutine get_key_para_cht
    !=========+=========+=========+=========+=========+=========+=========+

    !=========+=========+=========+=========+=========+=========+=========+
    subroutine get_key_para_rel(keyword,nfile,relt)
    character(15),intent(inout):: keyword
    integer      , intent(in)  :: nfile 
    real(8), intent(out) :: relt
    character(15)  :: keyword2

    character(50)  :: chartmp

    integer       :: ierr
    integer       :: i,j
       rewind(nfile)
       do while (.true.)
           read(nfile, "(A)", iostat=ierr) chartmp
           if (ierr /= 0)  exit
           do i=1,50
             if(chartmp(i:i)==":".or.chartmp(i:i)=="=") exit
           enddo
           if (i==51) cycle
           read(chartmp(2:i-1), *)  keyword2
           if(keyword2==keyword) then
          !write(6,*) chartmp
           read(chartmp(i+1:50), *)  relt
           exit
           endif
       enddo ! over while (.true.) loop

    end subroutine get_key_para_rel
    !=========+=========+=========+=========+=========+=========+=========+

    !=========+=========+=========+=========+=========+=========+=========+
    subroutine get_key_para_vec(keyword,nfile,nv,vect)
    character(15),intent(inout):: keyword
    integer      , intent(in)  :: nfile 
    integer      , intent(in)  :: nv
    real(8), intent(out) :: vect(nv)
    character(15)  :: keyword2

    character(50)  :: chartmp

    integer       :: ierr
    integer       :: i,j
       rewind(nfile)
       do while (.true.)
           read(nfile, "(A)", iostat=ierr) chartmp
           if (ierr /= 0)  exit
           do i=1,50
             if(chartmp(i:i)==":".or.chartmp(i:i)=="=") exit
           enddo
           if (i==51) cycle
           read(chartmp(2:i-1), *)  keyword2
           if(keyword2==keyword) then
          !write(6,*) chartmp
           read(chartmp(i+1:50), *)  vect(:)
           exit
           endif
       enddo ! over while (.true.) loop

    end subroutine get_key_para_vec
    !=========+=========+=========+=========+=========+=========+=========+



    !=========+=========+=========+=========+=========+=========+=========+
    subroutine get_key_para_loc(keyword,nfile)
    character(15),intent(inout) :: keyword
    integer      , intent(in)   :: nfile 
   !character(15), intent(inout):: chat
    character(15)  :: keyword2

    character(50)  :: chartmp

    integer       :: ierr
    integer       :: i,j
       rewind(nfile)
       do while (.true.)
           read(nfile, "(A)", iostat=ierr) chartmp
          !if (ierr /= 0)  exit
           do i=1,50
             if(chartmp(i:i)==":".or.chartmp(i:i)=="=") exit
           enddo
           if (i==51) cycle
           read(chartmp(2:i-1), *)  keyword2
           if(keyword2==keyword) then
          !write(6,*) chartmp
   !       read(chartmp(i+1:50), *)  chat
           exit
           endif
       enddo ! over while (.true.) loop

    end subroutine get_key_para_loc
    !=========+=========+=========+=========+=========+=========+=========+

    !=========+=========+=========+=========+=========+=========+=========+
    subroutine get_key_para_intct(keyword,nfile,intn)
    character(15),intent(inout):: keyword
    integer      , intent(in)  :: nfile 
    integer      , intent(out) :: intn
    character(15)  :: keyword2

    character(50)  :: chartmp

    integer       :: ierr
    integer       :: i,j
     ! rewind(nfile)
       do while (.true.)
           read(nfile, "(A)", iostat=ierr) chartmp
           if (ierr /= 0)  exit
           do i=1,50
             if(chartmp(i:i)==":".or.chartmp(i:i)=="=") exit
           enddo
           if (i==51) cycle
           read(chartmp(2:i-1), *)  keyword2
           if(keyword2==keyword) then
          !write(6,*) chartmp
           read(chartmp(i+1:50), *)  intn
           exit
           endif
       enddo ! over while (.true.) loop

    end subroutine get_key_para_intct
    !=========+=========+=========+=========+=========+=========+=========+


    !=========+=========+=========+=========+=========+=========+=========+
    subroutine get_key_para_nvec(nfile,nk,kpoints)
    integer      , intent(in)  :: nfile 
    integer      , intent(in)  :: nk
    real(8), intent(out) :: kpoints(3,0:nk)

    integer       :: ierr
    integer       :: i,j
       do j=0,nk
           read(nfile, *, iostat=ierr) kpoints(:,j)
           if (ierr /= 0)  exit
       enddo ! over while (.true.) loop

    end subroutine get_key_para_nvec
    !=========+=========+=========+=========+=========+=========+=========+

