subroutine kgroup(dim_sym, num_sym, invrot, veck, litt_group, num_litt_group)

    !!! This subroutine find the little group G(k) of wave vector k
    !!! \hat{R}\vec{k} = \vec{k} + \vec{Km}
    !!! Note that \hat{R}\vec{k} = (k1,k2,k3)Z^{-1}(g1,g2,g3)^T
    !!! Thus the input rotation matrix is inv(Ri)

    implicit none 

    integer, intent(in)  :: dim_sym, num_sym 
    integer, intent(in)  :: invrot(3,3,dim_sym)
    real(8), intent(in)  :: veck(3)

    integer, intent(out) :: litt_group(dim_sym)
    integer, intent(out) :: num_litt_group 

    integer :: i, j
    real(8) :: Rveck(3), diff 
    real(8), parameter   :: tolk = 0.1E-4 

    litt_group = 0
    num_litt_group = 0
    do i = 1, num_sym 
        do j = 1, 3
            Rveck(j) = dot_product(veck(:),invrot(:,j,i))-veck(j)
        enddo 
        diff = dabs(nint(Rveck(1))-Rveck(1)) &
              +dabs(nint(Rveck(2))-Rveck(2)) &
              +dabs(nint(Rveck(3))-Rveck(3))  
        if (diff <= tolk) then 
            num_litt_group = num_litt_group + 1
            litt_group(num_litt_group) = i
        endif 
    enddo 
    return 

end subroutine kgroup 
