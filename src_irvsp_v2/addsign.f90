    subroutine getsign(DoubNum,SymElemR,SymElemS,IZ,SU2)
    implicit none
    integer, parameter   :: dp=8!kind(1.0d0)
    integer     , intent(in)   :: DoubNum
    integer     , intent(in)   :: SymElemR(3,3,DoubNum)
    complex(dp) , intent(in)   :: SymElemS(2,2,DoubNum)
    integer     , intent(in)   :: IZ (3,3,DoubNum)
    complex(dp) , intent(inout):: SU2(2,2,DoubNum)        
    integer :: i, j, ir, jr 
    real(dp) :: s1,s2,s3,s4

    integer :: InvLoc
    !integer, allocatable :: SymElemR(:,:,:), IZ(:,:,:)
    !complex(dp), allocatable :: SymElemS(:,:,:), SU2(:,:,:), DSU2(:,:,:)
    complex(dp), allocatable :: DSU2(:,:,:)
    integer :: tmpdet 
    complex(dp) :: tmpmat(2,2)

    integer, allocatable :: MultbRef(:,:), Multb(:,:)
    integer, allocatable :: DetRef(:), Det(:)

    integer, allocatable :: compare(:,:)

    !> variables used for finding generators
    logical, allocatable :: IsGene(:), ByGene(:)
    integer :: irtmp, ir2 

    !> variables used for assign assumptions
    integer :: AssuGene, AssuNum, ia, ia_tmp, bina_l, gene_cnt 
    integer, allocatable :: SignList(:,:), SignBina(:)
    
    integer :: irCorres, ir2Corres, irPos, SignListAssu_tmp  
    integer, allocatable :: MultbAssu(:,:), SignListAssu(:)
    logical, allocatable :: AssuFail(:)
    logical :: AllDone 

    integer :: ProdRef, ProdCorres, Prod


!> get input from file fort.
   !open(unit=5556, file='fort.5556')
   !read(5556,*) DoubNum 

   !allocate(SymElemR(3,3,DoubNum))
   !allocate(IZ(3,3,DoubNum/2))
   !allocate(SymElemS(2,2,DoubNum))
   !allocate(SU2(2,2,DoubNum/2))

   !read(5556, *) SymElemR(:,:,1:DoubNum)
   !do ir = 1, DoubNum 
   !   do i = 1, 2
   !      read(5556, *) s1,s2,s3,s4
   !      SymElemS(1,i,ir) = dcmplx(s1,s2)
   !      SymElemS(2,i,ir) = dcmplx(s3,s4)
   !   enddo 
   !enddo 
   !read(5556, *) IZ(:,:,1:DoubNum/2)
   !do ir = 1, DoubNum/2  
   !   do i = 1, 2
   !      read(5556, *) s1,s2,s3,s4
   !      SU2(1,i,ir) = dcmplx(s1,s2)
   !      SU2(2,i,ir) = dcmplx(s3,s4)
   !   enddo 
   !enddo 
   !close(5556)
!> end get input


!> enlarge the su2 matrix 
    allocate(DSU2(2,2,DoubNum))
    do ir = 1, DoubNum/2 
       DSU2(:,:,ir) = SU2(:,:,ir)
       DSU2(:,:,DoubNum/2+ir) = -SU2(:,:,ir)
    enddo 
!>


!> calculate the determination of SO(3) matrices
!> as another dimension when calculating multiply table
    allocate(DetRef(DoubNum))
    allocate(Det(DoubNum))

    do ir = 1, DoubNum 
       call detso3(SymElemR(:,:,ir), DetRef(ir))
    enddo 
    do ir = 1, DoubNum/2 
       call detso3(IZ(:,:,ir), Det(ir))
       Det(DoubNum/2+ir) = Det(ir)
    enddo 
!> end the calculation of determination 

    
!> calculate the multiply table
    allocate(MultbRef(DoubNum,DoubNum))
    allocate(Multb(DoubNum,DoubNum))

    do ir = 1, DoubNum 
       do jr = 1, DoubNum 
          tmpmat = matmul(SymElemS(:,:,ir), SymElemS(:,:,jr))
          tmpdet = DetRef(ir)*DetRef(jr)
          call findindex(DoubNum, SymElemS, DetRef, tmpmat, tmpdet, MultbRef(ir,jr))
       enddo 
    enddo

    do ir = 1, DoubNum 
       do jr = 1, DoubNum 
          if (MultbRef(ir,jr) > DoubNum/2) MultbRef(ir,jr) = DoubNum/2 - MultbRef(ir,jr)
       enddo 
    enddo 

   !write(6,*) "Multb for the Reference input is"
   !do ir = 1, DoubNum/2 
   !   write(6,"(50I3)") MultbRef(ir,1:DoubNum/2)
   !enddo

    do ir = 1, DoubNum 
       do jr = 1, DoubNum 
          tmpmat = matmul(DSU2(:,:,ir), DSU2(:,:,jr))
          tmpdet = Det(ir)*Det(jr)
          call findindex(DoubNum, DSU2, Det, tmpmat, tmpdet, Multb(ir,jr))
       enddo 
    enddo 
    
    do ir = 1, DoubNum 
       do jr = 1, DoubNum
          if (Multb(ir,jr) > DoubNum/2) Multb(ir,jr) = DoubNum/2 - Multb(ir,jr)
       enddo 
    enddo 

   !write(6,*) "Multb for the input is"
   !do ir = 1, DoubNum/2 
   !   write(6,"(50I3)") Multb(ir,1:DoubNum/2)
   !enddo 
!> end the calculation of multiply table


!> compare the two multb
    allocate(compare(DoubNum,DoubNum))
    compare = 4
    do ir = 1, DoubNum 
       do jr = 1, DoubNum 
          if (Multb(ir,jr) == MultbRef(ir,jr)) compare(ir,jr) = 0
          if (Multb(ir,jr) ==-MultbRef(ir,jr)) compare(ir,jr) =-1
       enddo 
    enddo 
   !write(6,*) "Compare matrix"
   !do ir = 1, DoubNum/2
   !   write(6,"(50I3)") compare(ir,1:DoubNum/2)
   !enddo 
!>


!> find the location of inversion
    InvLoc = 0
    do ir = 1, DoubNum/2 
       if (SymElemR(1,1,ir)+SymElemR(2,2,ir)+SymElemR(3,3,ir)==-3) InvLoc = ir 
    enddo 
   !write(6,*) "The location of inversion is", InvLoc 
!>


!> find the generators of the single group from the multb
    allocate(IsGene(DoubNum/2))
    allocate(ByGene(DoubNum/2))
    IsGene = .true.
    ByGene = .false.
    !> The choice for generator is abritary, we choose the second 
    !> symmetry operation to be one of the generators,
    !> note that this choise may lead to a larger generator set
    do ir = 2, DoubNum/2 
       if(IsGene(ir)==.true.) then 
          ByGene(ir) = .true.
          do ir2 = 2, DoubNum/2  
             if (ByGene(ir2)) then 
                irtmp = ir 
                do 
                   irtmp = abs(MultbRef(ir2,irtmp))
                   if (ByGene(irtmp) == .true.) exit
                   IsGene(irtmp) = .false.
                   ByGene(irtmp) = .true.
                  !write(6,*) "ir, ir2, irtmp", ir, ir2, irtmp
                enddo 
             endif 
          enddo 
       endif 
    enddo 

   !do ir = 2, DoubNum/2 
   !   if (IsGene(ir)) write(6,*) "Generator:", ir
   !enddo
!>


!> make assumptions on the generators except inversion
    AssuGene = 0 
    do ir = 2, DoubNum/2 
       if (IsGene(ir).and.ir/=InvLoc) AssuGene = AssuGene + 1
    enddo 
    allocate(SignBina(AssuGene)) 
    
    AssuNum = 2**AssuGene 
    allocate(SignList(DoubNum, AssuNum)) 
    SignList = 0

    do ia = 1, AssuNum 
       !> convert the decimal number ia to binary number
       !> which is the sign of the generators 
       !> the binary number SignBina is sorted inversely 
       SignBina = -1 
       bina_l = 1 
       ia_tmp = ia-1
       do 
          if (ia_tmp>0) then 
              SignBina(bina_l) = SignBina(bina_l) + 2*mod(ia_tmp,2)
          else 
              exit 
          endif 
          ia_tmp = ia_tmp/2 
          bina_l = bina_l + 1 
       enddo 
      !write(6,*) "Sign of generator is", SignBina 

       !> assign this assumption to one of the signlist
       gene_cnt = 0
       do ir = 2, DoubNum/2 
          if (IsGene(ir).and.ir/=InvLoc) then 
              gene_cnt = gene_cnt + 1
              SignList(ir,ia) = SignBina(gene_cnt)
          endif 
          if (ir==InvLoc) SignList(ir,ia) = 1 
       enddo 
       !> the double group part have the inverse sign
       do ir = 2, DoubNum/2 
          if (SignList(ir,ia)/=0) SignList(DoubNum/2+ir, ia) = -SignList(ir,ia)
       enddo 
      !write(6,*) "SignList is"
      !write(6,"(50I2)") SignList(1:DoubNum/2,ia)

    enddo ! ia 
!>


!> filling the other positions in the Signlist 
!> to find contradiction in assumptions
    allocate(SignListAssu(DoubNum))
    allocate(AssuFail(AssuNum))
    AssuFail = .false.
    do ia = 1, AssuNum
       SignListAssu = SignList(:,ia)
       
       DO 
          do ir = 1, DoubNum 
             if (SignListAssu(ir)/=0) then 
                irCorres = SignListAssu(ir)*(mod(ir-1, DoubNum/2)+1)
                if (irCorres < 0)  irCorres = abs(irCorres) + DoubNum/2 
                !write(*,*) "ir is", ir, "irCorres is", irCorres 

                do ir2 = 1, DoubNum 
                   if (SignListAssu(ir2)/=0) then 
                      ir2Corres = SignListAssu(ir2)*(mod(ir2-1, DoubNum/2)+1)
                      if (ir2Corres < 0)  ir2Corres = abs(ir2Corres) + DoubNum/2 
                      !write(*,*) "ir2 is", ir2, "ir2Corres is", ir2Corres 

                      SignListAssu_tmp = sign(1, Multb(irCorres, ir2Corres))
                      irPos = MultbRef(ir, ir2)
                      if (irPos < 0) irPos = abs(irPos) + DoubNum/2 
                      if (SignListAssu(irPos) == 0) then 
                          SignListAssu(irPos) = SignListAssu_tmp 
                      else if (SignListAssu(irPos) /= SignListAssu_tmp) then 
                          !write(*,*) "irPos:", irPos, SignListAssu(irPos), SignListAssu_tmp 
                          AssuFail(ia) = .true.
                      endif  
                   endif 

                enddo ! ir2
             endif 
          enddo ! ir

          if (AssuFail(ia)) then 
             !write(6,*) "Assumption ", ia, " failed"
              exit 
          endif

          AllDone = .true. 
          do ir = 1, DoubNum 
             if (SignListAssu(ir) == 0) then 
                 AllDone = .false.
                 exit
             endif 
          enddo 

          if (AllDone) then 
             SignList(:,ia) = SignListAssu 
             exit
          endif 
                
       ENDDO 

    enddo ! ia

   !do ia = 1, AssuNum 
   !   write(6,*) "Assumption is "
   !   write(6,"(50I2)") SignList(:,ia)
   !enddo
!>


!> Check the multb for the right assumptions
    allocate(MultbAssu(DoubNum, DoubNum))
    do ia = 1, AssuNum 
       if (.not.AssuFail(ia)) then 
          
           !> construct the multb under this assumption
           do ir = 1, DoubNum 
              irCorres = SignList(ir,ia)*(mod(ir-1, DoubNum/2)+1)
              if (irCorres < 0)  irCorres = abs(irCorres) + DoubNum/2 
              do ir2 = 1, DoubNum 
                 ir2Corres = SignList(ir2,ia)*(mod(ir2-1, DoubNum/2)+1)
                 if (ir2Corres < 0)  ir2Corres = abs(ir2Corres) + DoubNum/2 
                 
                 ProdRef = MultbRef(ir, ir2)
                 if (ProdRef < 0)  ProdRef = abs(ProdRef) + DoubNum/2 
                 ProdCorres = SignList(ProdRef, ia)*(mod(ProdRef-1, DoubNum/2)+1)
                 Prod = Multb(irCorres, ir2Corres) 
                 if (ProdCorres /= Prod)  AssuFail(ia) = .true. 
              enddo ! ir2
           enddo ! ir
       endif 
    enddo ! ia

   !write(*,*) " "
   !write(*,*) "Further verification"
   !do ia = 1, AssuNum 
   !   if (.not.AssuFail(ia)) then 
   !       write(6,*) "Assumption", ia, " passed"
   !       write(6,"(50I2)") SignList(1:DoubNum/2,ia)
   !   endif 
   !enddo 
!>
    do ia = 1, AssuNum 
       if (.not.AssuFail(ia)) then 
           do ir = 1, DoubNum/2
              jr=SignList(ir,ia)
              if(jr==-1) then
                 SU2(:,:,ir)=-SU2(:,:,ir)
                 SU2(:,:,ir+DoubNum/2)=-SU2(:,:,ir+DoubNum/2)
              endif
           enddo
           return
       endif 
    enddo 
    if(ia==AssuNum+1)STOP "ERROR: The symmetries are not consisent with the space group"
    return
    END subroutine getsign



subroutine detso3(A,det)
    implicit none 
    integer, intent(in) :: A(3,3)
    integer, intent(out) :: det 

    det = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) + &
          A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3)) + &
          A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

end subroutine 


subroutine findindex(DoubNum, SymElemS, DetRef, tmpmat, tmpdet, ind)
    implicit none 
    integer, parameter :: dp = kind(1.0d0)

    integer, intent(in) :: DoubNum 
    complex(dp), intent(in) :: SymElemS(2,2,DoubNum)
    integer, intent(in) :: DetRef(DoubNum)
    complex(dp), intent(in) :: tmpmat(2,2)
    integer, intent(in) :: tmpdet 
    integer, intent(out) :: ind 

    integer :: ir, i, j
    logical :: EqualDet(DoubNum)
    real(dp) :: DiffBar 

    EqualDet = .false. 
    do ir = 1, DoubNum 
       if (DetRef(ir)==tmpdet) EqualDet(ir) = .true.
    enddo 

    do ir = 1, DoubNum 
       DiffBar = 0d0 
       if (EqualDet(ir)) then 
           do i = 1, 2
              do j = 1, 2
                 DiffBar = DiffBar + abs(SymElemS(i,j,ir)-tmpmat(i,j))
              enddo 
           enddo 
           if (DiffBar < 1e-4)  ind = ir 
       endif 
    enddo 

    end subroutine 
