
subroutine wake_pertub(pcx_ml, pcy_ml, pcz_ml, n, vcx, vcy, vcz)
    implicit none

    integer, parameter :: max_points = 1000
    real(8), intent(in) :: pcx_ml(n), pcy_ml(n), pcz_ml(n)
    integer, intent(in) :: n
    character(len=255) :: machline_command
    character(len=200) :: file_output_velocity
    integer, parameter :: num_columns_ml = 20
    integer :: i
    real(8), dimension(n, num_columns_ml) :: wake_data
    real(8), intent(out) :: vcx(n), vcy(n), vcz(n)
    character(len=200) :: csv_file
    character(len=100) :: header

    ! Part 1: Read the test points into the program
    csv_file = '/home/nkb20/FYP/mfw_code/heli-fw/mfw_wk_cpy/MachLine/' // &
               'studies/helicopter-body-fyp/tests/test-points/' // &
               'collocation_points.csv'
    header = 'x,y,z'
    open(unit=10, file=csv_file, status='replace', action='write')
    write(10, '(A)') header
    do i = 1, n
        write(10, '(3(F12.6, A1))') pcx_ml(i), ",", pcy_ml(i), ",", pcz_ml(i)
    end do
    close(unit=10)

    ! Part 2: Call Machline program 
    machline_command = "/home/nkb20/FYP/mfw_code/heli-fw/mfw_wk_cpy/MachLine/machline.exe " // &
                   & "/home/nkb20/FYP/mfw_code/heli-fw/mfw_wk_cpy/MachLine/" // &
                   & "studies/helicopter-body-fyp/tests/fyp.json"

    CALL SYSTEM(machline_command)

    ! Part 3: Read the output from Machline
    file_output_velocity = '/home/nkb20/FYP/mfw_code/heli-fw/mfw_wk_cpy/MachLine/' // &
                            'studies/helicopter-body-fyp/tests/test-points/' // &
                            'ind_velocities.csv'
    open(unit=20, file=file_output_velocity, status='old', action='read')
    ! reading entire csv as a matrix, disregarding labels 
    read(20, *)
    do i = 1, n
        read(20, *) wake_data(i, :)
    end do
    close(20)

    vcx = wake_data(:, 18)
    vcy = wake_data(:, 19)
    vcz = wake_data(:, 20)
end subroutine wake_pertub
