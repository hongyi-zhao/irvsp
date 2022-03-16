      SUBROUTINE WRTDATE(IUNIT)
!
      CHARACTER*10     CDATE,CTIME,CZONE
      INTEGER          IUNIT,ICVALUES(8)
!-------------------------------------------------------
!
!.....writes out date and time
      CALL DATE_AND_TIME(CDATE,CTIME,CZONE,ICVALUES)
      WRITE(IUNIT,100) CDATE(1:4),CDATE(5:6),CDATE(7:8), &
                       CTIME(1:2),CTIME(3:4),CTIME(5:6)

      RETURN
 100  FORMAT('### ',A4,'-',A2,'-',A2,2X,A2,':',A2,':',A2,/)
      END
