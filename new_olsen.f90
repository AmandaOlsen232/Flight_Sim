module olsen_m
implicit none
! integer, parameter :: dp = selected_real_kind(16,307)
real, parameter :: PI=3.1415926535897932384626433832795
real, parameter :: TOLERANCE = 1.0e-12
real, parameter :: g_ssl = 9.80665 ! [m/s^2]
real, parameter :: REz = 6356766 ! [m]
real, parameter :: ft_to_m = 0.3048 ![m]
real, parameter :: R_to_K = 1.8 ![K]
real, parameter :: Pen_to_PSI = 47.880258 ! [N/m^2]
real, parameter :: densen_to_densSI = 515.379 ! [kg/m^3]
real, parameter :: Po = 101325 ! [N/m^2]
real, parameter :: R = 287.0528 ! [J/kg-K]
real, parameter :: gamma = 1.4
real, dimension(9), parameter :: Zi = [0., 11000., 20000., 32000., 47000., 52000., 61000., 79000., 90000.] ! [m]
real, dimension(8), parameter :: Zi_1 = [11000., 20000., 32000., 47000., 52000., 61000., 79000., 90000.] ! [m]
real, dimension(8), parameter :: Ti = [288.150, 216.650, 216.650, 228.650, 270.650, 270.650, 252.650, 180.650] ! [K]
real, dimension(8), parameter :: Ti_p = [-6.5, 0.0, 1.0, 2.8, 0.0, -2.0, -4.0,0.0]/1000 ! [K/km]
real, dimension(8), parameter :: P_list = [101325.000000000, &
                                           22632.0318222212, &
                                           5474.87352827083, &
                                           868.014769086723, &
                                           110.905588989225, &
                                           59.0005242789244, &
                                           18.2099249050177, &
                                           1.03770045489203]



contains
function quat_mult(A, B) result(answer)
    ! 1.8.3 Write a function named quat_mult that accepts two arrays of length
    ! 4 representing 2 quaternions and returns the quaternion product
    implicit none
    real, dimension(4), intent(in) :: A, B
    real :: A0, Ax, Ay, Az, B0, Bx, By, Bz
    real, dimension(4) :: answer
    A0 = A(1)
    Ax = A(2)
    Ay = A(3)
    Az = A(4)

    B0 = B(1)
    Bx = B(2)
    By = B(3)
    Bz = B(4)

    ! answer(1) = A0*B0 - Ax*Bx - Ay*By - Az*Bz
    ! answer(2) = A0*Bx + Ax*B0 + Ay*Bz - Az*By
    ! answer(3) = A0*By - Ax*Bz + Ay*B0 + Az*Bx
    ! answer(4) = A0*Bz + Ax*By - Ay*Bx + Az*B0
    answer = [A0*B0 - Ax*Bx - Ay*By - Az*Bz, &
              A0*Bx + Ax*B0 + Ay*Bz - Az*By, &
              A0*By - Ax*Bz + Ay*B0 + Az*Bx, &
              A0*Bz + Ax*By - Ay*Bx + Az*B0]
    ! write(*,*) answer
end function quat_mult


function quat_base_to_dependent(v1, e) result(v2)
    !1.8.4 Write a function named quat_base_to_dependent to compute a vector
    ! transformation using the Euler-Rodrigues quaternion. The function
    ! should accept an array of length three containing the three
    ! components of a vector in the base coordinate system as well as
    ! an array of length four containing the four components of the
    ! quaternion that represent the orientation of the dependent
    ! coordinate system. This function should return an array of length
    ! three containing the vector components in the dependent
    ! coordinate system.
    implicit none
    real, dimension(3), intent(in) :: v1 !array of length three containing the three components of a vector in the base coordinate system
    real, dimension(4), intent(in) :: e !array of length four containing the four components of the quaternion that represent the orientation of the dependent coordinate system
    real, dimension(4) :: e_neg, T, v1_mult, v2_four
    real, dimension(3) :: v2

    e_neg = [e(1), &
            -e(2), &
            -e(3), &
            -e(4)]

    v1_mult = [0., &
               v1(1), &
               v1(2), &
               v1(3)]

    T = quat_mult(v1_mult, e) !eq. 1.5.6
    v2_four = quat_mult(e_neg, T) !eq. 1.5.7

    v2 = [v2_four(2), &
          v2_four(3), &
          v2_four(4)]

end function quat_base_to_dependent

function quat_dependent_to_base(v2, e) result(v1)
    !1.8.5 Write a function named quat_dependent_to_base to compute an inverse
    ! transformation using the Euler-Rodrigues quaternion. The function should
    ! accept an array of length three containing the three components of a
    ! vector in the dependent coordinate system as well as an array of length
    ! four containing the four components of the quaternion that represent the
    ! orientation of the dependent coordinate system. This function should return
    ! an array of length three containing the vector components in the base
    ! coordinate system.
    implicit none
    real, dimension(3), intent(in) :: v2
    real, dimension(4), intent(in) :: e
    real, dimension(4) :: e_neg, T, v2_mult, v1_four
    real, dimension(3) :: v1

    e_neg = [e(1), &
            -e(2), &
            -e(3), &
            -e(4)]

    v2_mult = [0., &
               v2(1), &
               v2(2), &
               v2(3)]

    T = quat_mult(v2_mult, e_neg) !eq. 1.5.5
    v1_four = quat_mult(e, T) !eq. 1.5.5
    v1 = [v1_four(2), &
          v1_four(3), &
          v1_four(4)]
end function quat_dependent_to_base