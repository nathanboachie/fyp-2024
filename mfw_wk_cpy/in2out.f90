subroutine in2out(px,py,pz,n)
    real :: xx(26), yy(26), zz(26)
    real :: px(n), py(n), pz(n)
    integer :: i
    integer :: inout = 0
    real :: nearest_dist
    
    ! Open the CSV file
    open(unit=10, file='/home/nkb20/FYP/mfw_code/heli-fw/mfw_wk_cpy/fuselage_array.csv', status='old', action='read')
    
    ! Read data from CSV file
    do i = 1, 26
        ! Read data from each line of the CSV file
        read(10, *) xx(i), yy(i), zz(i)
    end do

    close(10)

    do i = 1, n
        call pnpoly3d(px(i), py(i), pz(i), xx, yy, zz, 26, inout)
        if (inout == -1) then
        else if (inout == 0) then
            ! If on the edge, increase coordinates by 20%
            px(i) = px(i) * 1.2
            py(i) = py(i) * 1.2
            pz(i) = pz(i) * 1.2
        else if (inout == 1) then
            ! Map to nearest coordinate
            nearest_dist = HUGE(1.0) ! Initialize with a large value
            do j = 1, 26
                dist = sqrt((px(i)-xx(j))**2 + (py(i)-yy(j))**2 + (pz(i)-zz(j))**2)
                if (dist < nearest_dist) then
                    nearest_dist = dist
                    nearest_index = j
                end if
            end do
            ! Increase coordinates by 20% (radially)
            px(i) = xx(nearest_index) * 1.2
            py(i) = yy(nearest_index) * 1.2
            pz(i) = zz(nearest_index) * 1.2
        end if
    end do
    
end subroutine in2out
