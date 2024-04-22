!Transform coordinates from bc system into mach line coordinate system

subroutine ml_coord(Px, Py, Pz, n, px_ml, py_ml, pz_ml)
    integer, intent(in) :: n
    real(8), intent(in) :: Px(n), Py(n), Pz(n)
    real(8), parameter :: pi = 3.141592653589793238462643383279502884d0
    real(8), dimension(3,3) :: TS_matrix 
    real(8), dimension(n), intent(out) :: px_ml, py_ml, pz_ml
    real(8) :: x_tmp, y_tmp, z_tmp
    integer :: i

    ! Shaft transformation coordinates
    TS_matrix = reshape([1.0d0, 0.0d0, 0.0d0, &
                        0.0d0, cos(3*pi/180.0d0), -sin(3*pi/180.0d0), &
                        0.0d0, sin(3*pi/180.0d0), cos(3*pi/180.0d0)], [3, 3])

    do i = 1, n
        x_tmp = TS_matrix(1,1)*Px(i) + TS_matrix(1,2)*Py(i) + TS_matrix(1,3)*Pz(i)
        y_tmp = TS_matrix(2,1)*Px(i) + TS_matrix(2,2)*Py(i) + TS_matrix(2,3)*Pz(i)
        z_tmp = TS_matrix(3,1)*Px(i) + TS_matrix(3,2)*Py(i) + TS_matrix(3,3)*Pz(i)

        ! Store transformed coordinates into new arrays
        px_ml(i) = x_tmp
        py_ml(i) = y_tmp
        pz_ml(i) = z_tmp
    end do 
end subroutine ml_coord
