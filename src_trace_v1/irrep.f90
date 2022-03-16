      PROGRAM IRREDUCIBLE_REPRESENTATIONS
      USE SYMM,only:ssym
      USE WAVE_data
      USE STRUCT_data
      implicit none
!-----
!
      integer  :: KKK,nele

      logical  :: exists,FL2
!
      character*90 :: chaps
      character*10 :: chap1,chap2
      character*15 :: chtp15
      real(dp) :: rele

      ! command argument
      integer :: narg,iarg,lens,stat
      character(len=60) :: arg,cmd
      
      nele=0
      
      call get_command(cmd)
      write(*,*) 'Current command : ',cmd
      narg = command_argument_count()
      write(*,*) 'Argument count : ',narg
      IF(narg>1) THEN
       write(*,*)  'ERROR!!! Please input a(one) number  of bands by the following command:'
       STOP '##$ vasp2trace $nele'
      ENDIF
      iarg = 1
      do while(.true.) 
         call get_command_argument(iarg, arg,lens,stat)
         IF(len_trim(arg) == 0) THEN
           exit
         ELSE
           read(arg,*) nele
           write(*,"('The given number of occupied bands : ',I5)") nele
         ENDIF
         iarg = iarg +1
      end do

!-------check files: OUTCAR and WAVECAR!!!
     exists=.false.
     inquire(file="OUTCAR",exist=exists)
     if(exists .eqv. .false.)  then
        write(6,*) " ERROR: NO OUTCAR exists !!! "
        stop
     endif

     OPEN(unit=547,file='OUTCAR',form='formatted',status='old')
     FL2=.FALSE.
     do 
      read(547,"(A90)") chaps
      if(chaps(4:9)=='LNONCO' ) exit
     enddo
      read(547,"(A15L5)") chtp15,FL2
     do 
      read(547,"(A90)") chaps
      if(chaps(4:9)=='NELECT' ) exit
     enddo
     !print*, chaps
      read(chaps,*) chap1,chap2,rele
     CLOSE(unit=547)

     exists=.false.
     inquire(file="WAVECAR",exist=exists)
     if(exists .eqv. .false.)  then
        write(6,*) " ERROR: NO WAVECAR exists !!! "
        stop
     endif
!-------check files: OUTCAR and WAVECAR!!!
!
      IF(nele==0) THEN
       nele=nint(rele)
       IF(.NOT.FL2) THEN
        if (mod(nele,2)==1) write(6,*) " WARNING: An odd number of electrons !!! "
        nele=nele/2
       ENDIF
       WRITE(0,*) 'The total number of occupied bands(nele) was read in OUTCAR !!!'
       WRITE(0,'(A10,I5)') 'nele:',nele
       WRITE(0,*) '### to give a desirable number of occupied bands by the following command:'
       WRITE(0,*) '###$ vasp2trace $nele'
      ENDIF


     OPEN(unit=624,file='trace.txt',form='formatted',status='unknown')
      WRITE(624,'(I3)') nele
      CALL SSYM()
!
      call setarray()
!
!.....loop over k-points
      DO KKK=1,NSPIN*NKPTS
      CALL KPTIN(KKK,nele)
      ENDDO
!
      call downarray()
!
     CLOSE(unit=624)
      STOP 'TOTAL END-----Congratulations!!!'
      END
      
