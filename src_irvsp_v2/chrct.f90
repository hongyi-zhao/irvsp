subroutine chrct(NMAT,NUME,A,B,SU2,EE,L,PH,XM,FL,NE,NV,LKG,IKG,IORD,KKK,IR)
    use symm, only:FLMAX,MAXIR,NSYM,MAXDG,TOLDG,nmin
    implicit none 

    integer ,  parameter :: dp=kind(1.0d0)

    integer , intent(in) :: NMAT, NUME, NE, NV, IKG, IORD, KKK ,IR
    logical , intent(in) :: FL(FLMAX)
    integer , intent(in) :: LKG(NSYM), L(NSYM,NMAT)
    real(dp), intent(in) :: EE(NUME)
    complex(dp), intent(inout) :: A(NMAT,NUME), B(NMAT,NUME)
    complex(dp), intent(in)    :: SU2(2,2,NSYM),PH(NSYM,NMAT)
    complex(dp), intent(out)   :: XM(NSYM,NUME)


    complex(dp) :: A2(NMAT,MAXDG), B2(NMAT,MAXDG)
    complex(dp) :: CSUM, GAM(MAXDG,MAXDG)

    integer :: i, j, i1,  j2, no,np
    integer :: IE
    integer :: ND, IG 

    complex(dp),parameter :: CZERO=cmplx(0._DP,0._DP,DP)

!**********************************************************************
!
    A2=CZERO 
    B2=CZERO 
    XM = cmplx(1.2345678_dp,1.2345678_dp,dp)
   
!
!.....normalization
    do j = 1, NE 
        csum = czero
        csum = dot_product(A(:,j),A(:,j))
        if (FL(2)) csum = csum+dot_product(B(:,j),B(:,j))
        csum = zsqrt(csum)
        A(:, j) = A(:, j)/csum
       !if (FL(2))  B(:, j) = B(:, j)/csum
        if (FL(2))  then
           B(:, j) = B(:, j)/csum
           if (IR .eq. 0 ) then
               A2(:,1)= conjg(A(:,j))
               B2(:,1)= conjg(B(:,j))
               A( :,j)=-B2(:,1)
               B( :,j)= A2(:,1)
           elseif(IR .ne. 1) then
               A2(:,1)= (SU2(1,1,IR)*A(:,j)+SU2(1,2,IR)*B(:,j))
               B2(:,1)= (SU2(2,1,IR)*A(:,j)+SU2(2,2,IR)*B(:,j))
               A(:, j)= A2(:,1)
               B(:, j)= B2(:,1)
           endif
        else
           if(IR.eq.0)  A(:,j)= conjg(A(:,j))
        endif
       !if ( IR .ne. 0 .and. FL(2) ) then
       !   A2(:,1)= (SU2(1,1,IR)*A(:,j)+SU2(1,2,IR)*B(:,j))
       !   B2(:,1)= (SU2(2,1,IR)*A(:,j)+SU2(2,2,IR)*B(:,j))
       !   A(:, j)= A2(:,1)
       !   B(:, j)= B2(:,1)
       !endif
    enddo 

!
!.....characters for all energies and symm. ops.
   !IE = 1
    IE =nmin
    DO WHILE(IE .LE. NE)
        ND = 1
        IF(IE<NE) THEN
        do while( (EE(IE+ND)-EE(IE)).LT.TOLDG )
            ND = ND + 1
            if((IE+ND).GT.NE) exit
        enddo 
        ENDIF

        if (ND.GT.MAXDG) then 
            WRITE(*,'(A,I3)') "WARNING: ND is larger than MAXDG for the band:",IE
            IE = IE + ND 
            cycle
        endif 
!
!.....for all classes
        do IG = 1, IKG 
!
!.....rotation of spin components; include phase from real space rotation 
        
            if (FL(2)) then 
                j2 = 0
                A2=CZERO 
                B2=CZERO 
                do j = IE, IE+ND-1
                j2 = j2 + 1
                do i = 1, NV 
                    A2(i,j2)=(SU2(1,1,LKG(IG))*A(i,j)+ &
                        SU2(1,2,LKG(IG))*B(i,j))*PH(IG,i)
                    B2(i,j2)=(SU2(2,1,LKG(IG))*A(i,j)+ &
                        SU2(2,2,LKG(IG))*B(i,j))*PH(IG,i)
                enddo 
                enddo 
            else 
                j2 = 0
                A2=CZERO 
                B2=CZERO 
                do j = IE, IE+ND-1
                j2 = j2 + 1
                do i = 1, NV
                    A2(i,j2)=A(i,j)*PH(IG,i)
                enddo 
                enddo 
            endif 
!
!.....character
            XM(LKG(IG), IE) = CZERO 
            GAM(:,:)=CZERO
            do no = 0, ND-1 
            do np = 0, ND-1
               if(no.ne.np) cycle
                CSUM = CZERO 
                do i = 1, NV 
                    CSUM=CSUM+DCONJG(A(L(IG,i),IE+np))*A2(i,1+no) &
                             +DCONJG(B(L(IG,i),IE+np))*B2(i,1+no)
                enddo 

                GAM(np+1,no+1)=CSUM
                IF(no.EQ.np) XM(LKG(IG),IE)=XM(LKG(IG),IE)+CSUM
            enddo 
            enddo 
!
!.....for all degenerate states
            do i1=IE+1,IE+ND-1
                XM(LKG(IG),i1)=XM(LKG(IG),IE)
            enddo 
!
        enddo ! IG
        IE=IE+ND

    ENDDO ! while

!
!.....output (characters for all bands and all symm ops.)
    do IE=1,NE
        write(5,510) KKK,IE,EE(IE)
        write(5,530) (DREAL(XM(I,IE)),DIMAG(XM(I,IE)),I=1,IORD)
    enddo 
!
    return 
 510  FORMAT(2I5,1X,F11.7,$)
 530  FORMAT(16(/,3(2F11.7,2X)))

end subroutine
