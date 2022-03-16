module dump

    use lib_comms
    use bilbao 

    real(dp), save :: RAN(4,MAXSYM)
    integer,  save :: ILC(MAXSYM)

    public :: dump_opes
    public :: dump_bilbaochrct

contains



SUBROUTINE dump_opes(NSYM, IORD, IZ,IIZ,TAU,DZ2, SU2_getsign)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      LOGICAL          FG1,FG2
      CHARACTER*6      RNAM(NSYM)
      DIMENSION        IZ(3,3,NSYM),IIZ(3,3,NSYM),TAU(3,NSYM)
      DIMENSION        ANG(3,NSYM)
      COMPLEX*16       SU2(2,2,NSYM),SU2_getsign(2,2,NSYM),SDET
      COMPLEX*16       SX,SY,SZ(2,2),SZT(2,2)
      DIMENSION        IZ2(3,3),IZ3(3,3)
!     DIMENSION        BR1(3,3),DR1(3,3)
      DIMENSION        JIJ(NSYM,NSYM),IUNI(3,3)
      DIMENSION        DTMP(3,3),DZ2(3,3,NSYM)
      DATA             TOLM/1E-3/
!****************************************************************************
!
      DO 2 I1=1,3
        IUNI(I1,I1)=1
        DO 2 I2=1,3
 2        IF(I1.NE.I2) IUNI(I1,I2)=0
!     DO 4 I1=1,NSYM
!     DO 3 J1=1,4
!3       RAN(J1,I1)=0
!     DO 4 I2=1,NSYM
!4       JIJ(I1,I2)=0
      RAN = 0
      JIJ = 0
!
!.....calculate JIJ(i,j)=n, where Rn=Rj*Ri*inv(Rj); 
!.....for all Ri and Rj:
      DO 300 I=1,IORD
      DO 200 J=1,IORD
        JIJ(I,J)=0
!
!.......calculate -inv(Rj)*tj
!        DO 250 I1=1,3
! 250       TAU1(I1)=-IZ(1,I1,J)*TAU(1,J)
!     &              -IZ(2,I1,J)*TAU(2,J)
!     &              -IZ(3,I1,J)*TAU(3,J)
!
!.........calculate Pi*inv(Pj)
          DO 202 I1=1,3
          DO 202 I2=1,3
 202         IZ2(I1,I2)=IZ(I1,1,I)*IIZ(1,I2,J) &
                       +IZ(I1,2,I)*IIZ(2,I2,J) &
                       +IZ(I1,3,I)*IIZ(3,I2,J)
!          DO 252 I1=1,3
! 252         TAU2(I1)=IZ(I1,1,I)*TAU1(1)
!     &               +IZ(I1,2,I)*TAU1(2)
!     &               +IZ(I1,3,I)*TAU1(3)
!     &               +TAU(I1,I)
!
!.........calculate Pj*Pi*inv(Pj)
          DO 204 I1=1,3
          DO 204 I2=1,3
 204         IZ3(I1,I2)=IZ(I1,1,J)*IZ2(1,I2) &
                       +IZ(I1,2,J)*IZ2(2,I2) &
                       +IZ(I1,3,J)*IZ2(3,I2)
!
!          DO 254 I1=1,3
! 254         TAU3(I1)=IZ(I1,1,J)*TAU2(1)
!     &               +IZ(I1,2,J)*TAU2(2)
!     &               +IZ(I1,3,J)*TAU2(3)
!     &               +TAU(I1,J)
!
!.........find the smallest non-primitiv translation of tau3
!          DO 260 I1=1,3
!          DO 260 I2=4,1,-1
!            IF(TAU3(I1).GE. I2) TAU3(I1)=TAU3(I1)-I2
!            IF(TAU3(I1).LE.-I2) TAU3(I1)=TAU3(I1)+I2
! 260      CONTINUE
!
!.........find Rn
          FG1=.TRUE.
          K=0
          DO WHILE((FG1).AND.(K.LT.IORD))
            K=K+1
            IDFF=0
            DO 206 I1=1,3
            DO 206 I2=1,3
 206           IDFF=IDFF+ABS(IZ(I1,I2,K)-IZ3(I1,I2))
            IF(IDFF.EQ.0) THEN
              FG1=.FALSE.
              JIJ(I,J)=K
            ENDIF
          END DO
          IF(FG1) STOP 'rmprop: cannot find Rj*Ri*inv(Rj)'
 200  CONTINUE
 300  CONTINUE
!
!
      WRITE(6,*)
      DO 100 I=1,IORD
!
!.....type of rotations, E,C2,C3,C4,C6
      IDET=IZ(1,1,I)*IZ(2,2,I)*IZ(3,3,I) &
          +IZ(1,2,I)*IZ(2,3,I)*IZ(3,1,I)   &
          +IZ(1,3,I)*IZ(2,1,I)*IZ(3,2,I) &
          -IZ(1,1,I)*IZ(2,3,I)*IZ(3,2,I)        &
          -IZ(1,2,I)*IZ(2,1,I)*IZ(3,3,I) &
          -IZ(1,3,I)*IZ(2,2,I)*IZ(3,1,I)
      IF(IDET.NE.1.AND.IDET.NE.-1) STOP 'classes: |det| not 1'
      IDFF=0
      DO 23 I1=1,3
      DO 23 I2=1,3
         IZ3(I1,I2)=IDET*IZ(I1,I2,I)
 23      IDFF=IDFF+ABS(IZ3(I1,I2)-IUNI(I1,I2))
      IP=1
      DO WHILE(IP.LE.6.AND.IDFF.NE.0) 
        IP=IP+1
        IDFF=0
        DO 24 I1=1,3
        DO 24 I2=1,3
           IZ2(I1,I2)=IDET*IZ(I1,1,I)*IZ3(1,I2) &
                     +IDET*IZ(I1,2,I)*IZ3(2,I2) &
                     +IDET*IZ(I1,3,I)*IZ3(3,I2)
 24        IDFF=IDFF+ABS(IZ2(I1,I2)-IUNI(I1,I2))
        DO 25 I1=1,3
        DO 25 I2=1,3
 25       IZ3(I1,I2)=IZ2(I1,I2)
      END DO 
      IF(IDFF.NE.0) STOP 'classes: cannot find type of rot (1)' 
      IF(IP.EQ.5)   STOP 'classes: cannot find type of rot (2)' 
!
      IP=IP*IDET
      IF(IP.EQ.1) WRITE(RNAM(I)(1:6),'(A6)')       '   E  '
      IF(IP.GT.1) WRITE(RNAM(I)(1:6),'(A3,I1,A2)') '  C',IP,'   '
      IF(IP.EQ.-1)WRITE(RNAM(I)(1:6),'(A6)')       '   I  '
      IF(IP.LT.-1)WRITE(RNAM(I)(1:6),'(A3,I1,A2)') ' IC',-IP,'  '
      ILC(I)=IP
!
!.....rotation in the Cartesian coordinate system
!.....IZ_c = (inv(BR1)~)*IZ*(BR1~)
!     DO 320 I1=1,3
!     DO 320 I2=1,3
!320     DTMP(I1,I2)=DFLOAT(IZ(I1,1,I))*BR1(I2,1) &
!                   +DFLOAT(IZ(I1,2,I))*BR1(I2,2) &
!                   +DFLOAT(IZ(I1,3,I))*BR1(I2,3)
!     DO 321 I1=1,3
!     DO 321 I2=1,3
!321     DZ2(I1,I2,I)=DR1(1,I1)*DTMP(1,I2) &
!                    +DR1(2,I1)*DTMP(2,I2) &
!                    +DR1(3,I1)*DTMP(3,I2)
!
      DET= DZ2(1,1,I)*DZ2(2,2,I)*DZ2(3,3,I) &
          +DZ2(1,2,I)*DZ2(2,3,I)*DZ2(3,1,I) &
          +DZ2(1,3,I)*DZ2(2,1,I)*DZ2(3,2,I) &
          -DZ2(1,1,I)*DZ2(2,3,I)*DZ2(3,2,I) &
          -DZ2(1,2,I)*DZ2(2,1,I)*DZ2(3,3,I) &
          -DZ2(1,3,I)*DZ2(2,2,I)*DZ2(3,1,I)
!
      IF(DABS(DET-DBLE(IDET)).GT.TOLM) &
                  STOP 'rmprop: wrong DETERM'
!
!.....rotation axis 
      CALL RAXIS(IP,DET,DZ2(1,1,I),RAN(1,I))
!
!.....determine the operations acting on the spin components
      CALL SU2OP(JIJ,IORD,DZ2(1,1,I),DET,SU2(1,1,I),ANG(1,I))
!
 100  CONTINUE
!
!.....use trace(u(Pi)) >= 0
      DO 109 I=1,IORD
      SX=SU2(1,1,I)+SU2(2,2,I)
      IF(DREAL(SX).LT.0.D0) THEN
        DO 90 I1=1,2
        DO 90 J1=1,2
 90        SU2(I1,J1,I)=-SU2(I1,J1,I)
        SX=-SX
      ENDIF
      IF(DREAL(SX).LT.TOLM.AND.DREAL(SU2(1,1,I)).LT.TOLM )THEN
        YY1=ABS(DREAL(SU2(1,1,I)))+ABS(DIMAG(SU2(1,1,I)))
        YY2=ABS(DREAL(SU2(1,2,I)))+YY1
        IF(DREAL(SU2(1,1,I)).LT.-TOLM ) THEN
           FG2=.TRUE.
        ELSEIF(DIMAG(SU2(1,1,I)).LT.-TOLM ) THEN
           FG2=.TRUE.
        ELSEIF(YY1.LT.TOLM.AND.DREAL(SU2(1,2,I)).LT.-TOLM ) THEN
           FG2=.TRUE.
        ELSEIF(YY2.LT.TOLM.AND.DIMAG(SU2(1,2,I)).LT.-TOLM ) THEN
           FG2=.TRUE.   
        ELSE
           FG2=.FALSE.   
        ENDIF 
        IF(FG2) THEN
          DO 91 I1=1,2
          DO 91 J1=1,2
 91          SU2(I1,J1,I)=-SU2(I1,J1,I)
          SX=-SX
        ENDIF
      ENDIF
109  CONTINUE
!
!.....use only unbarred operations [u(Pi)=-u(Pi_bar)]
      DO 110 I=1,IORD
      SX=SU2(1,1,I)+SU2(2,2,I)
      SY=SU2(1,2,I)+SU2(2,1,I)
      SDET=SU2(1,1,I)*SU2(2,2,I)-SU2(2,1,I)*SU2(1,2,I)
      IF(DABS(DIMAG(SX))  .GT.TOLM) STOP 'imagina trace(u)'
      IF(DABS(DIMAG(SDET)).GT.TOLM) STOP 'imagina det(u)  '
      DO 120 J=1,IORD         
         K=JIJ(I,J)
         IF(K.LT.I) THEN
!..........inv(Uj)
           SDET =SU2(1,1,J)*SU2(2,2,J) &
                -SU2(2,1,J)*SU2(1,2,J) 
           SZ(1,1)= SU2(2,2,J)/SDET
           SZ(1,2)=-SU2(1,2,J)/SDET
           SZ(2,1)=-SU2(2,1,J)/SDET
           SZ(2,2)= SU2(1,1,J)/SDET
!las           do 888 i1=1,2
!           do 888 i2=1,2
! 888          ss(i1,i2)=sz(i1,i2)
!..........Ui*inv(Uj)
           SZT(1,1)=SU2(1,1,I)*SZ(1,1)+SU2(1,2,I)*SZ(2,1)  
           SZT(1,2)=SU2(1,1,I)*SZ(1,2)+SU2(1,2,I)*SZ(2,2)  
           SZT(2,1)=SU2(2,1,I)*SZ(1,1)+SU2(2,2,I)*SZ(2,1)  
           SZT(2,2)=SU2(2,1,I)*SZ(1,2)+SU2(2,2,I)*SZ(2,2)  
!..........Uj*Ui*inv(Uj)
           SZ(1,1)=SU2(1,1,J)*SZT(1,1)+SU2(1,2,J)*SZT(2,1)  
           SZ(1,2)=SU2(1,1,J)*SZT(1,2)+SU2(1,2,J)*SZT(2,2)  
           SZ(2,1)=SU2(2,1,J)*SZT(1,1)+SU2(2,2,J)*SZT(2,1)  
           SZ(2,2)=SU2(2,1,J)*SZT(1,2)+SU2(2,2,J)*SZT(2,2)  
!
           DFP=0D0
           DFM=0D0
           DO 121 I1=1,2
           DO 121 I2=1,2
              DFP=DFP+ABS(SZ(I1,I2)+SU2(I1,I2,K))
              DFM=DFM+ABS(SZ(I1,I2)-SU2(I1,I2,K))
 121       CONTINUE
           IF(ABS(DFP).GT.TOLM.AND.ABS(DFM).GT.TOLM)  &
               STOP 'rmprop: not correct un=uj*ui*inv(uj)'
           IF(ABS(DFP).LT.TOLM) THEN
              DO 122 I1=1,2
              DO 122 I2=1,2
 122             SU2(I1,I2,I)=-SU2(I1,I2,I)
           ENDIF
!..........check determinates
           SDETK=SU2(1,1,K)*SU2(2,2,K) &
                -SU2(2,1,K)*SU2(1,2,K)
           SDET =SU2(1,1,I)*SU2(2,2,I) &
                -SU2(2,1,I)*SU2(1,2,I)
           IF(ABS(SDETK-SDET).GT.TOLM)  &
                STOP 'rmprop: |u(i)| not |u(k)| despite same class'
         ENDIF 
 120  CONTINUE
 110  CONTINUE
!
!.....output
       WRITE(6,500)
       DO 150 I=1,IORD
       WRITE(6,516)(IZ( 1,I2,I),I2=1,3),TAU(1,I), &
                   (IIZ(1,I2,I),I2=1,3), &
                   (DZ2(1,I2,I),I2=1,3),ANG(1,I)
       DO 410 I1=2,3                                                 
  410     WRITE(6,515)(IZ(I1,I2,I),I2=1,3),TAU(I1,I), &
            (IIZ(I1,I2,I),I2=1,3), &
            (DZ2(I1,I2,I),I2=1,3),ANG(I1,I),(SU2_getsign(I1-1,I2,I),I2=1,2)
       WRITE(6,517) I
       IF(ILC(I).GE.1) THEN
         IF(ILC(I).EQ. 1) WRITE(6,518) '(E)','unity op.   '
         IF(ILC(I).EQ. 2) WRITE(6,519)  ILC(I),180
         IF(ILC(I).EQ. 3) WRITE(6,519)  ILC(I),120
         IF(ILC(I).EQ. 4) WRITE(6,521)  ILC(I), 90
         IF(ILC(I).EQ. 6) WRITE(6,521)  ILC(I), 60
       ELSE
         IF(ILC(I).EQ.-1) WRITE(6,518) '(I)','inverse op. '
         IF(ILC(I).EQ.-2) WRITE(6,525) -ILC(I),180
         IF(ILC(I).EQ.-3) WRITE(6,525) -ILC(I),120
         IF(ILC(I).EQ.-4) WRITE(6,527) -ILC(I), 90
         IF(ILC(I).EQ.-6) WRITE(6,527) -ILC(I), 60
       ENDIF
       IF(RAN(4,I).LT.-0.5) WRITE(6,530) 'clockwise '
       IF(RAN(4,I).GT.0.5)  WRITE(6,531) 'counterclockwise '
       WRITE(6,532) (RAN(I1,I),I1=1,3)
       WRITE(6,*)
  150  CONTINUE

       WRITE(6,*)
       WRITE(6,501)
       WRITE(6,*)
!
!      WRITE(6,570) &
!        'Note:                                                                         ',&
!        '(1) Labelling of IRs depends on atom positions.                               ',&
!        '(2) Labelling of IRs in the character table of Koster et al. as well as of    ',&
!        'Altmann et al. (see below) is defined by the chosen symmetry axes and for     ',&
!        'each point group. Hence, in order to compare the labelling of the present     ',&
!        'IRs with the literature one has to determine the relation, for each k-point,  ',&
!        'between their symmetry axes and the present ones.                             ',&
!        'That is, please check the group elements of the different classes.            '
!
!
 500  FORMAT(/,'SYMMETRY OPERATIONS Pi={Ri|taui+tm}',/, &
            '  Ri     taui ', &
            '  inv(Ri) Ri(Cartesian coord)', &
            '  Eulers angles  spin transf.'/) 
 501  FORMAT(80('*'))
 515  FORMAT(3I2,F8.3,2X,3I2,2X,3F6.3,2X,F6.3,2X,2('(',2F6.3,')'))
 516  FORMAT(3I2,F8.3,2X,3I2,2X,3F6.3,2X,F6.3)
 517  FORMAT(I4,$)
 518  FORMAT(1X,A3,/,A12) 
 519  FORMAT(1X,'(C',I1,')', /,I3,'-degree rotation') 
 521  FORMAT(1X,'(C',I1,')', /,I2,'-degree rotation')
 525  FORMAT(1X,'(IC',I1,')',/,I3,'-degree rotation times inversion')
 527  FORMAT(1X,'(IC',I1,')',/,I2,'-degree rotation times inversion')

 530  FORMAT(A10,$)
 531  FORMAT(A17,$)
 532  FORMAT('rotation through (',F6.3,2(',',F6.3),')')

 555  FORMAT(/)
! 570  FORMAT(/,8(A78,/))
!
end subroutine dump_opes


subroutine dump_kinformation(KKK, WK, kname, num_littg_input, sgn, timerev_k)

    implicit none 

    integer,          intent(in)  :: KKK 
    real(dp),         intent(in)  :: WK(3)
    character(len=3), intent(in)  :: kname 
    integer,          intent(in)  :: num_littg_input 
    integer,          intent(in)  :: sgn 
    logical,          intent(in)  :: timerev_k
    
    write(6,*)
    if (.not.timerev_k) then 
        write(6,"('knum =',I3,4X,'kname= ',A2)") KKK, trim(adjustl(kname)) 
    else 
        write(6,"('knum =',I3,4X,'kname= ',A2)") KKK, trim(adjustl(kname))//'A'
    endif 
    write(6,"('k =',3F9.6)") WK(:)
    write(6,*)
    if (.not.timerev_k) then 
        write(6,"(7X,'The k-point name is ',A3)") kname 
    else 
        write(6,"(7X,'The k-point name is ',A3)") trim(kname)//'A'
    endif 

    if (num_littg_input < 10) then 
        write(6,"(7X,I1,' symmetry operations (module lattice translations) in space group ', I3)") num_littg_input, sgn 
    else 
        write(6,"(7X,I2,' symmetry operations (module lattice translations) in space group ', I3)") num_littg_input, sgn 
    endif 
    write(6,"(7X,'We do NOT classify the elements into calsses.')")
    write(6,"(7X,'Tables can be found on website: http://www.cryst.ehu.es/.')")

end subroutine dump_kinformation 


subroutine dump_bilbaochrct(ikt,timerev_k, iout)

      implicit none 
      integer, intent(in) :: ikt
      logical, intent(in) :: timerev_k 
      integer, intent(in) :: iout

      character*15  :: ckpoint
      logical       :: ftg
      integer       :: i,j
      real(dp)      :: tmp_samplek(3)
      
      if (.not.timerev_k) then 
          call Kreal2string(samplek(:,ikt),ckpoint) 
          WRITE(iout,"(/,3X,I2,A15,A54,$)") ikt, samplekname(ikt)//' : kname '&
         ,ckpoint//' :  given in the conventional basis'
      else 
          tmp_samplek(:) = -samplek(:,ikt)
          call Kreal2string(tmp_samplek(:),ckpoint) 
          WRITE(iout,"(/,/,3X,I2,A15,A54,$)") ikt, trim(samplekname(ikt))//'A'//' : kname '&
         ,ckpoint//' :  given in the conventional basis'
      endif 

      WRITE(iout,'(/,3X,I2,A,$)') antisym(ikt) &
          ,' : the existence of antiunitary symmetries. 1-exist; 0-no'
      WRITE(iout,577) 'Reality' 
     !DO I=1,nelelittle(ikt)/2
     !WRITE(iout,579) I 
     !ENDDO
      DO I=1,num_doub_sym/2
         IF(labels(1,I,1,ikt)==1) WRITE(iout,579) I 
      ENDDO

      do j=1,nirreps(ikt)
         ftg=.false.

         IF (.not.timerev_k) then 
         IF(Irrepsname(j,ikt)(1:1)=='-') THEN
            WRITE(iout,576) Herringrule(j,ikt),Irrepsname(j,ikt)(2:5)//' '
         ELSE
            WRITE(iout,576) Herringrule(j,ikt),Irrepsname(j,ikt)
         ENDIF
         DO I=1,num_doub_sym/2
        !IF(labels(1,I,j,ikt)==1) WRITE(iout,580) tableTraces(I,j,ikt)
         IF(labels(1,I,j,ikt)==1) WRITE(iout,580) chartTraces(I,j,ikt)+cmplx(epsil,epsil,dp)
         IF(labels(1,I,j,ikt)==1.and. labels(2,I,j,ikt)==2 ) ftg=.true.
         ENDDO

         IF(ftg) THEN
                 WRITE(iout,578) '     '
                 DO I=1,num_doub_sym/2
                    IF(labels(1,I,j,ikt)==1) THEN
                       IF(labels(2,I,j,ikt)==2) THEN
                          WRITE(iout,582) factsTraces(I,j,ikt)
                       ELSE
                          WRITE(iout,582) '          '
                       ENDIF
                    ENDIF
                ENDDO
         ENDIF

         ELSE ! timerev_k
         IF(Irrepsname(j,ikt)(1:1)=='-') THEN
            WRITE(iout,576) Herringrule(j,ikt),Irrepsname(j,ikt)(2:2)//'A'//Irrepsname(j,ikt)(3:5)
         ELSE
            WRITE(iout,576) Herringrule(j,ikt),Irrepsname(j,ikt)(1:1)//'A'//Irrepsname(j,ikt)(2:4)
         ENDIF
         DO I=1,num_doub_sym/2
        !IF(labels(1,I,j,ikt)==1) WRITE(iout,580) tableTraces(I,j,ikt)
         IF(labels(1,I,j,ikt)==1) WRITE(iout,580) dconjg(chartTraces(I,j,ikt))+cmplx(epsil,epsil,dp)
         IF(labels(1,I,j,ikt)==1.and. labels(2,I,j,ikt)==2 ) ftg=.true.
         ENDDO

         IF(ftg) THEN
                 WRITE(iout,578) '     '
                 DO I=1,num_doub_sym/2
                    IF(labels(1,I,j,ikt)==1) THEN
                       IF(labels(2,I,j,ikt)==2) THEN
                          WRITE(iout,582) factsTraces(I,j,ikt)
                       ELSE
                          WRITE(iout,582) '          '
                       ENDIF
                    ENDIF
                ENDDO
         ENDIF

         ENDIF ! timerev_k


         IF(Irrepsname(j+1,ikt)(1:1)=='-'.and.Irrepsname(j,ikt)(1:1)/='-') THEN
            sirreps(ikt)=j 
           !WRITE(iout,'(I3)') j
            WRITE(iout,581) '----'
            do i=1,nelelittle(ikt)/2
            WRITE(iout,582) '------------'
            enddo
         ENDIF
      enddo

 510  FORMAT(/,7X,'WILL BE IMPLEMENTED')

 576  FORMAT(/,3X,I2,2X,A6,$)
 577  FORMAT(/,1X,A7,$)
 578  FORMAT(/,7X,A6,$)
 579  FORMAT(10X,I2,$)
 580  FORMAT(F5.2,SP,F5.2,'i ',$)
 581  FORMAT(/,8X,A4,$)
 582  FORMAT(A12,$)
end subroutine dump_bilbaochrct


subroutine dump_reps(num_bands, bot_band, top_band, &
                     numdeg_set, ene_bands, &
                     num_littg_input, littg_input, littg_bilbao, &
                     kopeconjg_bb2inp, &
                     chrct_set, reps_set, numreps_set)

    implicit none 

    integer,           intent(in)  :: num_bands, bot_band, top_band
    integer,           intent(in)  :: numdeg_set(num_bands) 
    real(dp),          intent(in)  :: ene_bands(num_bands)

    integer,           intent(in)  :: num_littg_input
    integer,           intent(in)  :: littg_input(MAXSYM)
    integer,           intent(in)  :: littg_bilbao(MAXSYM)

    integer,           intent(in)  :: kopeconjg_bb2inp(MAXSYM)

    complex(dp),       intent(in)  :: chrct_set(num_bands, MAXSYM)
    character(len=20), intent(in)  :: reps_set(num_bands)
    integer,           intent(in)  :: numreps_set(num_bands)

    integer :: i, j, k, m, n

    j = bot_band
    do while (j <= top_band)
        write(6, "(2I3, F10.6, $)") j, numdeg_set(j), ene_bands(j)
        do i = 1, num_littg_input
            m = littg_bilbao(i)
            do n = 1, MAXSYM
                if (kopeconjg_bb2inp(n) == m) exit 
            enddo 
            write(6, "(F5.2,SP,F5.2,'i ',$)") dreal(chrct_set(j,n))+epsil, &
                                               dimag(chrct_set(j,n))+epsil 
        enddo 
        if (numreps_set(j) == 1) then
            write(6, "('=',A4)") reps_set(j)(1:4)
        else if (numreps_set(j) == 2) then
            write(6, "('=',A4,' + ',A4)") reps_set(j)(1:4), &
                                          reps_set(j)(6:9)
        else if (numreps_set(j) == 3) then
            write(6, "('=',A4,' + ',A4,' + ',A4)") reps_set(j)(1:4), &
                                                   reps_set(j)(6:9), &
                                                   reps_set(j)(11:14)
        else if (numreps_set(j) == 4) then 
            write(6, "('=',A4,' + ',A4,' + ',A4, ' + ',A4)") reps_set(j)(1:4), &
                                                             reps_set(j)(6:9), &
                                                             reps_set(j)(11:14), &
                                                             reps_set(j)(16:19)
        else if (numreps_set(j) == 0) then 
            write(6, "('=',A2)") "??"
        endif 

        j = j + numdeg_set(j)
    enddo 

    write(6,*)
    write(6,"(80('*'))")
    write(6,*)

end subroutine dump_reps 


subroutine dump_save(sgn, bot_band, top_band)

    implicit none 

    integer, intent(in) :: sgn 
    integer, intent(in) :: bot_band, top_band 

    integer :: ik, ib, ir   
    integer :: irtmp, iktmp

    IF (save_kcount <= MAXKSAVE) then  

    open(unit=67, file='fort.67')
    write(67,'(3I5)') sgn, save_kcount, top_band-bot_band+1
    do ik = 1, save_kcount 
        write(67,'(I3,$)') save_ktype(ik)
        ib = bot_band
        do while (ib <= top_band)
            if (save_numrep(ib, ik) ==0) then 
                write(67,'(A3,$)') ' ??'
            else 
                do ir = 1, save_numrep(ib, ik)
                    write(67, '(I3,$)') save_chrct(ib+ir-1,ik)
                enddo
            endif 
            ib = ib + save_numdeg(ib, ik)
        enddo 
        write(67,*)
    enddo 
    close(67)
    open(unit=68, file='tqc.data')
    write(68,'(3I5)') sgn, save_kcount, top_band-bot_band+1
    do ik = 1, save_kcount 
        write(68,'(I3,$)') save_ktype(ik)
        ib = bot_band
        do while (ib <= top_band)
            if (save_numrep(ib, ik) ==0) then 
                write(68,'(A3,$)') ' ??'
            else 
                do ir = 1, save_numrep(ib, ik)
                    write(68, '(I3,$)') save_chrct(ib+ir-1,ik)
                enddo
            endif 
            ib = ib + save_numdeg(ib, ik)
        enddo 
        write(68,*)
    enddo 
    close(68)

    open(unit=66, file='tqc.txt')
    write(66,'(A,I3,A,I3)') 'Computed bands:', bot_band, ' -', top_band
    do ik = 1, save_kcount 
        if (save_ktype(ik) > 0) then 
            write(66, '(A4,$)') samplekname(save_ktype(ik))//': '
        else
            write(66, '(A4,$)') samplekname(-save_ktype(ik))(1:1)//'A'//': '
        endif 
        ib = bot_band
        do while (ib <= top_band)
            if (save_numrep(ib, ik) ==0) then 
                write(66,'(A3,$)') ' ??'
                if (save_numdeg(ib,ik) < 10) write(66,'(A1,I1,A3,$)') '(',save_numdeg(ib,ik),'); '
                if (save_numdeg(ib,ik) >=10) write(66,'(A1,I2,A3,$)') '(',save_numdeg(ib,ik),'); '
            else 
                do ir = 1, save_numrep(ib, ik)
                    irtmp = save_chrct(ib+ir-1,ik)
                    if (save_ktype(ik) > 0) then 
                        iktmp = save_ktype(ik)
                        if (Irrepsname(irtmp,iktmp)(1:1)=='-') then 
                            write(66,'(A4,$)') Irrepsname(irtmp,iktmp)(2:5)
                        else
                            write(66,'(A4,$)') Irrepsname(irtmp,iktmp)(1:4)
                        endif 
                    else 
                        iktmp =-save_ktype(ik)
                        if (Irrepsname(irtmp,iktmp)(1:1)=='-') then 
                            write(66,'(A4,$)') Irrepsname(irtmp,iktmp)(2:2)//'A'//Irrepsname(irtmp,iktmp)(3:4)
                        else
                            write(66,'(A4,$)') Irrepsname(irtmp,iktmp)(1:1)//'A'//Irrepsname(irtmp,iktmp)(2:3)
                        endif 
                    endif 
                enddo
                if (save_numdeg(ib,ik) < 10) write(66,'(A1,I1,A3,$)') '(',save_numdeg(ib,ik),'); '
                if (save_numdeg(ib,ik) >=10) write(66,'(A1,I2,A3,$)') '(',save_numdeg(ib,ik),'); '
            endif 
            ib = ib + save_numdeg(ib, ik)
        enddo 
        write(66,'(A1,I3,A1)') '[',top_band-bot_band+1,']'
    enddo 
    close(66)

    ENDIF 

end subroutine dump_save 


subroutine dump_littg(IKG, LKG, IZ, TAU, SK)
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
      INTEGER     :: IKG
      INTEGER     :: LKG(MAXSYM)

      INTEGER     :: IZ(3,3,MAXSYM)
      COMPLEX(DP) :: PH(MAXSYM)

      DIMENSION        SK(3),TAU(3,MAXSYM)
      DATA             TOL/5.0E-6/

      WRITE(6,516)

      DO IG=1,IKG 
        !WRITE(6,520) IG,LKG(IG)
         I=LKG(IG)
       IF(ILC(I).GE.1) THEN
         IF(ILC(I).EQ. 1) WRITE(6,545) 'E'    ,LKG(IG)
         IF(ILC(I).EQ. 2) WRITE(6,546)  ILC(I),LKG(IG)
         IF(ILC(I).EQ. 3) WRITE(6,546)  ILC(I),LKG(IG)
         IF(ILC(I).EQ. 4) WRITE(6,546)  ILC(I),LKG(IG)
         IF(ILC(I).EQ. 6) WRITE(6,546)  ILC(I),LKG(IG)
       ELSE
         IF(ILC(I).EQ.-1) WRITE(6,545) 'I'    ,LKG(IG)
         IF(ILC(I).EQ.-2) WRITE(6,547) -ILC(I),LKG(IG)
         IF(ILC(I).EQ.-3) WRITE(6,547) -ILC(I),LKG(IG)
         IF(ILC(I).EQ.-4) WRITE(6,547) -ILC(I),LKG(IG)
         IF(ILC(I).EQ.-6) WRITE(6,547) -ILC(I),LKG(IG)
       ENDIF

      DO 443 I1=8,19
 443     WRITE(6,'(A,$)') ' '
         ARG=-2*PI* ( SK(1)*TAU(1,LKG(IG)) & 
                     +SK(2)*TAU(2,LKG(IG)) &
                     +SK(3)*TAU(3,LKG(IG)) )
         PH(LKG(IG))= DCMPLX(DCOS(ARG),DSIN(ARG))
         WRITE(6,525) DREAL(PH(LKG(IG))), &
                      DIMAG(PH(LKG(IG)))
         WRITE(6,531) (RAN(I2,LKG(IG)), I2=1,3)
      ENDDO
!
      WRITE(6,538)

      do IG=2,IKG
      WRITE(6,579) LKG(IG)
      enddo
      WRITE(6,*)
!
!
      RETURN
 500  FORMAT(/,7X,'Due to time-reversal symmetry, the following pairs', &
             /,7X,'of irreducible representations are degenerate:', &
             /,7X,'',$)
 502  FORMAT(I6,3F10.6,       ' # k number, k-vector')
 503  FORMAT(3X,A3,I3,I10,17X,  ' # pnt-grp, no symm.ops, no classes')
 507  FORMAT(/,1X,A6,1X,I2,$)
 508  FORMAT(/,1X,A6,/,I3,$)
 510  FORMAT('(',A,' and ',A,'); ',$)
 515  FORMAT(//,'elemt ,symmetry ops, exp(-i*k*taui)'$)!by zjwang 11.9.2016
 516  FORMAT(//,'elemt ,symmetry ops, exp(-i*k*taui),  main axes',$)
!520  FORMAT(/,A6,1X,I2,$)
 520  FORMAT(/,I3,4X,I2,$)
 521  FORMAT(I3,$)
 525  FORMAT(' (',SP,F5.2,F5.2,'i)',$)
 531  FORMAT('  (',F6.3,2(',',F6.3),')',$)
 538  FORMAT(//,'bnd ndg  eigval     E  ',$)
 539  FORMAT(6X,A6,$)
 543  FORMAT(A1,$)
 544  FORMAT(' ',$)
 579  FORMAT(10X,I2,$)
 581  FORMAT(' WARNING: PH=(',2F8.5,') PHold=(',2F8.5,')')

 545  FORMAT(/,A3     ,4X,I2,$)
 546  FORMAT(/,' C',I1,4X,I2,$)
 547  FORMAT(/,'IC',I1,4X,I2,$)

end subroutine dump_littg


end module dump
