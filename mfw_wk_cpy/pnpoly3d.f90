SUBROUTINE pnpoly3d(px, py, pz, xx, yy, zz, n, inout)
    REAL :: px, py, pz
    REAL :: xx(n), yy(n), zz(n)
    INTEGER :: n, inout
    REAL :: min_range, max_range

    ! Determine which coordinate to delete for projection
    ! Find the range of each coordinate
    min_range = MINVAL(xx) - MAXVAL(xx)
    max_range = MINVAL(yy) - MAXVAL(yy)
    z_range = MINVAL(zz) - MAXVAL(zz)
    ! Determine which coordinate has the smallest range
    if (min_range < max_range .and. min_range < z_range) then
        ! Delete the x-coordinate for projection
        call pnpoly(py, pz, yy, zz, n, inout)
    else if (max_range < z_range) then
        ! Delete the y-coordinate for projection
        call pnpoly(px, pz, xx, zz, n, inout)
    else
        ! Delete the z-coordinate for projection
        call pnpoly(px, py, xx, yy, n, inout)
    end if

    ! Perform the 2D point-in-polygon test on the projected polygon
    
END SUBROUTINE pnpoly3d
