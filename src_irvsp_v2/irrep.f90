      PROGRAM IRREDUCIBLE_REPRESENTATIONS
      USE SYMM
      USE WAVE_data
      USE STRUCT_data
      implicit none
!
      integer  :: KKK,nmax
      integer  :: sgn   = 0  !! space group number
      integer  :: ver_n = 0  !! version number
!

     ! command argument
      integer :: narg,iarg,lens,stat
      character(len=50) :: arg,cmd
      
      sgn=0;ver_n=4 ! default version number
      nmax=0
      
      call get_command(cmd)
      write(*,*) 'Current command : ',cmd
      narg = command_argument_count()
      write(*,*) 'Argument count : ',narg
      IF(narg==0) then
          write(*,"('Please input the correct space group number (sgn) and a version number (vn) by the command below:')")
          write(*,"('###$: irvsp -sg $sgn -v $vn')  ")
          STOP
      ELSE
      iarg = 1
      do while(.true.) 
         call get_command_argument(iarg, arg,lens,stat)
         IF(len_trim(arg) == 0)  EXIT
         
         IF(trim(arg)=='-sg') THEN
           iarg = iarg +1; call get_command_argument(iarg, arg,lens,stat) !tmp
           read(arg,*) sgn!; print*, sgn
         ELSEIF(trim(arg)=='-v') THEN
           iarg = iarg +1; call get_command_argument(iarg, arg,lens,stat) !tmp
           read(arg,*) ver_n!; print*, ver_n
         ELSEIF(trim(arg)=='-nb') THEN
           iarg = iarg +1; call get_command_argument(iarg, arg,lens,stat) !tmp
           read(arg,*) nmin!; print*, nmin
           iarg = iarg +1; call get_command_argument(iarg, arg,lens,stat) !tmp
           read(arg,*) nmax!; print*, nmax
         ELSE
          write(*,"('Please input the correct space group number (sgn) and a version number (vn) by the command below:')")
          write(*,"('###$: irvsp -sg $sgn -v $vn -nb $nmin $nmax')  ")
          STOP
         ENDIF
!        
         iarg = iarg +1
      end do
      ENDIF

      IF(sgn==0 .or. sgn>230) THEN
          write(*,"('Please input the CORRECT space group number (sgn) by the following option:')")
          write(*,"('###$: irvsp -sg $sgn ')  ")
          STOP
      ENDIF
      IF(nmin/=1 .or. nmax/=0) THEN
      IF(nmax<0 .or. nmin<1 .or. nmax < nmin) THEN
          write(*,"('Please input the proper range of bands from nmin to nmax:')")
          write(*,"('###$: irvsp -nb $nmin $nmax')  ")
          STOP
      ENDIF
      ENDIF

      IF    (ver_n==4) THEN
          write(*,"('You can choose another version by inputing a version number (nv):')")
          write(*,"('###$: irvsp -v $vn ')")
      ELSEIF(ver_n==1) THEN                            
       write(*,"('This calcalulation is for SG #',I3,' in Version I')  ") sgn   !  - version I  (v1)
      ELSEIF(ver_n==2) THEN                            
       write(*,"('This calcalulation is for SG #',I3,' in Version II') ") sgn   !  - version II (v2)
      ELSEIF(ver_n==3) THEN                            
       write(*,"('This calcalulation is for SG #',I3,' in Version III')") sgn   !  - version III(v3)
      ELSE                                             
          write(*,"('Please by input a VALID version number (nv=1,2,3,4):')")
          write(*,"('###$: irvsp -v $vn ')")
          STOP
      ENDIF
!
     !IF(sgn==0 .or. ver_n==0) THEN
     !   STOP
     !ENDIF
!
      CALL SSYM(sgn,ver_n,nmax)
!
      call setarray()
      
!  
    
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
       OPEN(UNIT=66,file='tqc.txt',form='formatted',status='unknown')
       OPEN(UNIT=68,file='tqc.data',form='formatted',status='unknown')
       WRITE(66,"(A,I3,A,I3,$)") "Computed bands:" ,nmin, ' -',nmax
       WRITE(67,"(3I5,$)") sgn,NKPTS,nmax-nmin+1
       WRITE(68,"(3I5,$)") sgn,NKPTS,nmax-nmin+1
!      ENDIF 

!
!.....loop over k-points
      DO KKK=1,NSPIN*NKPTS
      CALL KPTIN(KKK,nmax)
      ENDDO
!
      call downarray()
      
!      IF(NKPTS<10) THEN
      CLOSE(66)
      CLOSE(68)
      CLOSE(unit=624)
!      ENDIF
      
      
!
      WRITE(6,*) 
      WRITE(6,*) "*****************************************************"
      WRITE(6,*) 
      WRITE(6,*) "TOTAL END"
      STOP 'IRREP END-----Congratulations you finished the calculation.'
      END
      


