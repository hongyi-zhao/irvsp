      PROGRAM IRREDUCIBLE_REPRESENTATIONS
      USE SYMM
      USE WAVE_data
      USE STRUCT_data
      implicit none
!
      integer  :: KKK
      logical  :: exists
!
!---------check files: OUTCAR and WAVECAR!!!
      exists=.false.
      inquire(file="OUTCAR",exist=exists)
      if(exists .eqv. .false.)  then
         write(6,*) " ERROR: NO OUTCAR exists !!! "
         write(0,*) " ERROR: NO OUTCAR exists !!! "
         stop
      endif
!---------check the file: WAVECAR!!!
      exists=.false.
      inquire(file="WAVECAR",exist=exists)
      if(exists .eqv. .false.)  then
         write(6,*) " ERROR: NO WAVECAR exists !!! "
         write(0,*) " ERROR: NO WAVECAR exists !!! "
         stop
      endif
!---------
!
      CALL SSYM()
!
      call setarray()
!
!.....loop over k-points
      DO KKK=1,NSPIN*NKPTS
      CALL KPTIN(KKK)
      ENDDO
!
      call downarray()
!
      WRITE(6,*) 
      WRITE(6,*) "*****************************************************"
      WRITE(6,*) 
      WRITE(6,*) "TOTAL END"
      STOP 'IRREP END-----Congratulations you finished the calculation.'
      END

