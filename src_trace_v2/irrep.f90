      PROGRAM IRREDUCIBLE_REPRESENTATIONS
      USE SYMM,only:ssym
      USE WAVE_data
      USE STRUCT_data
      implicit none
!-----
!
      integer  :: KKK,nele,nele_up,nele_dn

      logical  :: exists,FL2
!
      character*90 :: chaps
      character*10 :: chap1,chap2
      character*15 :: chtp15
      real(dp) :: rele

      ! command argument
      integer :: narg,iarg,lens,stat
      character(len=90) :: arg,cmd
      
      nele=0
      
      call get_command(cmd)
      write(*,*) 'Current command : ',trim(cmd)
      narg = command_argument_count()
      write(*,*) 'Argument count : ',narg
      IF(narg==1) THEN
       write(*,*)  'WARNING: This program is ONLY for spin-polarized claculations without SOC!!!'
       write(*,*)  'ERROR!!! Please input the both numbers of bands by the following command:'
       STOP '##$ vasp2trace $nele_up $nele_dn'
      ENDIF
      iarg = 1
      do while(.true.) 
         call get_command_argument(iarg, arg,lens,stat)
         IF(len_trim(arg) == 0) THEN
           exit
         ELSE
           read(arg,*) nele_up; nele=nele_up
           iarg = iarg + 1
           call get_command_argument(iarg, arg,lens,stat)
           read(arg,*) nele_dn
           write(*,"('The given numbers of occupied bands (up- and down-spin, respectively) : ',2I5)") nele_up,nele_dn
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
       nele_up=nele/2
       nele_dn=nele/2
       IF(.NOT.FL2) THEN
        if (mod(nele,2)==1) write(6,*) " WARNING: An odd number of total electrons !!! "
       ENDIF
       WRITE(0,*) 'The total number of occupied bands(nele) was read in OUTCAR !!!'
       WRITE(0,'(2(A10,I5))') 'nele_up:',nele_up,'nele_dn:',nele_dn
       WRITE(0,*) '### to give desirable numbers of occupied bands by the following command:'
       WRITE(0,*) '###$ vasp2trace $nele_up $nele_dn'
      ENDIF

      IF(FL2) THEN
       WRITE(0,*) 'ERROR!!! *** calculation with spin-orbital coupling (SOC) ***'
       WRITE(0,*) 'ERROR!!! Please download the common version of vasp2trace at'
       WRITE(0,*) 'ERROR!!! https://github.com/zjwang11/irvsp/blob/master/src_trace_v1.tar.gz!!!'
       STOP
      ENDIF

     OPEN(unit=624,file='trace_up.txt',form='formatted',status='unknown')
     OPEN(unit=625,file='trace_dn.txt',form='formatted',status='unknown')
      WRITE(624,'(I3)') nele_up
      WRITE(625,'(I3)') nele_dn
      CALL SSYM()
!
      call setarray()
!
!.....loop over k-points
     !DO KKK=1,NSPIN*NKPTS
      DO KKK=1,NKPTS
      CALL KPTIN(KKK,nele_up)
      ENDDO
     CLOSE(unit=624)
     CLOSE(unit=625)
     OPEN(unit=624,file='trace_dn.txt',form='formatted',status='old',access='append')
      DO KKK=NKPTS+1,NSPIN*NKPTS
      CALL KPTIN(KKK,nele_dn)
      ENDDO
!
      call downarray()
!
     CLOSE(unit=624)
      STOP 'TOTAL END-----Congratulations!!!'
      END
      
