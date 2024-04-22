!Reverse coordinate transformation

subroutine ml_inv(vx, vy, vz,n)
    implicit none
    integer :: n
    real(8) :: vx(n), vy(n), vz(n)
    real(8), parameter :: pi = 3.141592653589793238462643383279502884d0
    real(8), dimension(3,3) :: inv_TS_matrix
    real(8) :: vx_tmp, vy_tmp, vz_tmp
    integer :: i

    ! Inverse of shaft transformation coordinates
    inv_TS_matrix = reshape((/ &
    1.0d0, 0.0d0, 0.0d0, &
    0.0d0, cos(3*pi/180.0d0), sin(3*pi/180.0d0), &
    0.0d0, -sin(3*pi/180.0d0), cos(3*pi/180.0d0) & 
    /), shape(inv_TS_matrix))


    do i = 1, n
        vx_tmp = inv_TS_matrix(1,1)*vx(i) + inv_TS_matrix(1,2)*vy(i) + inv_TS_matrix(1,3)*vz(i)
        vy_tmp = inv_TS_matrix(2,1)*vx(i) + inv_TS_matrix(2,2)*vy(i) + inv_TS_matrix(2,3)*vz(i)
        vz_tmp = inv_TS_matrix(3,1)*vx(i) + inv_TS_matrix(3,2)*vy(i) + inv_TS_matrix(3,3)*vz(i)

        ! Putting back into arrays
        vx(i) = vx_tmp
        vy(i) = vy_tmp
        vz(i) = vz_tmp
    end do 
end subroutine ml_inv
