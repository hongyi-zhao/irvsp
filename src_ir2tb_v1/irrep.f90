      PROGRAM IRREDUCIBLE_REPRESENTATIONS
      USE SYMM
      USE WAVE_data
      USE STRUCT_data
      implicit none
 
!*******************************************************************
!
      integer  :: KKK
!
      CALL SSYM()
!
      call setarray()
!
!.....loop over k-points
      DO KKK=0,NKPTS
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
      
