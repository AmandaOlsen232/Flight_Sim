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
real, dimension(8), parameter :: Ti_p = [-6.5, 0.0, 1.0, 2.8, 0.0, -2.0, -4.0, 0.0]/1000 ! [K/km]
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
    !    1.8.3 Write a function named quat_mult that accepts two arrays of length 
    !   4 representing 2 quaternions and returns the quaternion product
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
    !   transformation using the Euler-Rodrigues quaternion. The function
    !   should accept an array of length three containing the three 
    !   components of a vector in the base coordinate system as well as
    !   an array of length four containing the four components of the
    !   quaternion that represent the orientation of the dependent 
    !   coordinate system. This function should return an array of length
    !   three containing the vector components in the dependent 
    !   coordinate system.
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
    !   transformation using the Euler-Rodrigues quaternion. The function should
    !   accept an array of length three containing the three components of a 
    !   vector in the dependent coordinate system as well as an array of length 
    !   four containing the four components of the quaternion that represent the 
    !   orientation of the dependent coordinate system. This function should return
    !   an array of length three containing the vector components in the base
    !   coordinate system.
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

subroutine quat_norm(Q)
    !1.8.6 Write a function named quat_norm that accepts an array of length four 
    !   containing a quaternion, and returns the normalized quaternion. For this
    !   computation, use the exact solution given in Eq. (1.5.12). You may wish to
    !   overwrite the quaternion that was passed in to the function with the 
    !   normalized solution.
    implicit none
    real, dimension(4), intent(inout) :: Q
    ! real, dimension(4) :: Q_norm
    real :: mag

    mag = sqrt(Q(1)*Q(1) + Q(2)*Q(2) + Q(3)*Q(3) + Q(4)*Q(4)) !eq. 1.5.12
    Q = Q/mag

end subroutine quat_norm

function euler_to_quat(euler) result(e)
    !1.8.7 Write a function called euler_to_quat that accepts an array of length
    !   three containing the three Euler angles and returns an array of length 
    !   four containing values for the equivalent quaternion components.
    implicit none
    real, dimension(3), intent(in) :: euler
    real, dimension(4) :: e 
    real :: phi_rad, theta_rad, psi_rad, cphi, cpsi, ctheta, sphi, spsi, stheta

    phi_rad = Euler(1)
    theta_rad = Euler(2)
    psi_rad = Euler(3)

    cphi = cos(phi_rad*0.5)
    cpsi = cos(psi_rad*0.5)
    ctheta = cos(theta_rad*0.5)
    sphi = sin(phi_rad*0.5)
    spsi = sin(psi_rad*0.5)
    stheta = sin(theta_rad*0.5)

    e = [cphi*ctheta*cpsi + sphi*stheta*spsi, &
         sphi*ctheta*cpsi- cphi*stheta*spsi, &
         cphi*stheta*cpsi + sphi*ctheta*spsi, &
         cphi*ctheta*spsi - sphi*stheta*cpsi]

end function euler_to_quat

function quat_to_euler(e) result(euler)
    !1.8.8 Write a function named quat_to_euler that accepts an array of length 
    !   four containing the four components of a quaternion and returns an array 
    !   of length three containing values for the equivalent Euler angles.
    implicit none
    real, dimension(4) :: e
    real, dimension(3) :: euler
    real :: e0, ex, ey, ez, psi

    e0 = e(1)
    ex = e(2)
    ey = e(3)
    ez = e(4)
    psi = 0.

    if (abs(e0*ey - ex*ez - 0.5) < TOLERANCE) then
        euler = [2*asin(ex/cos(PI/4)) + psi, &
                 PI*0.5, &
                 psi]
    end if

    if ((e0*ey - ex*ez + 0.5) < TOLERANCE) then
        euler = [2*asin(ex/cos(PI/4)) - psi, &
                 -PI*0.5, &
                 psi]
    else
    euler = [atan2(2*(e0*ex + ey*ez), (e0*e0 + ez*ez - ex*ex - ey*ey)), &
             asin(2*(e0*ey - ex*ez)), &
             atan2(2*(e0*ez + ex*ey), (e0*e0 + ex*ex - ey*ey - ez*ez))]
    end if

end function quat_to_euler

!___________________________________________________________________________________
function gravity_SI(H) result(g)
    !3.13.1 Write a function to compute the gravity as a function of altitude
    !from Eq. (3.2.1). The input to the function call should be the geometric
    !altitude in meters, and the function should return the gravity in SI
    !units of m/s^2. Test your code against the solutions given in Table. B.1.1
    implicit none
    real :: H, g

    g = g_ssl*(REz/(REz+H))**2
end function gravity_SI

function gravity_English(H) result(g)
    !3.13.2 Write a function to compute the gravity as a function of geometric 
    !altitude in English units. The input to the function call should be the 
    !geometric altitude in feet, and the function should return the gravity 
    !in English units of ft/s^2. Test your code against the solutions given in 
    !Table. B.1.1
    implicit none
    real :: H, g 
    real :: g_ssl_english, REz_english

    g_ssl_english = g_ssl/ft_to_m
    REz_english = REz/ft_to_m
    g = g_ssl_english*(REz_english/(REz_english+H))**2
end function gravity_English

subroutine std_atm_SI(H, Z, T, Pz, rho, a)
    !3.13.3 Write a function to compute standard-atmospheric properties in SI units as
    !a function of geometric altitude. Given a geometric altitude in meters, your
    !code must return the geopotential altitude (m), temperature (K), pressure 
    !(N/m^2), density (kg/m^3), and speed of sound (m/s).
    implicit none
    real, intent(inout) :: Z, T, Pz, rho, a
    real, intent(in) :: H
    real :: Pi
    real, dimension(5) :: std_atm
    integer :: i

    !geopotential altitude
    Z = REz*H/(REz + H)
 
    !temperature
    do i=1, size(Zi)
        if (abs(Z - Zi(i)) < TOLERANCE) then
            exit
        else if (Z >= Zi(i) .and. Z < Zi(i+1)) then
            exit
        end if
    end do

    T = Ti(i) + (Ti_p(i))*(Z - Zi(i))
    
    !Pressure
    Pi = P_list(i)

    ! do k=1,i-1
    !     if (abs(Ti_p(k)) < TOLERANCE) then
    !         Pi = Pi*exp(-g_ssl*(Zi(k+1) - Zi(k))/(R*Ti(k)))
    !     else 
    !         Pi = Pi*((Ti(k) + (Ti_p(k))*(Zi(k+1) - Zi(k)))/Ti(k))**(-g_ssl/(R*Ti_p(k)))
    !     end if
    ! end do

    if (abs(Ti_p(i)) < TOLERANCE) then
        Pz = Pi*exp(-g_ssl*(Z - Zi(i))/(R*Ti(i)))
    else 
        Pz = Pi*((Ti(i) + (Ti_p(i))*(Z - Zi(i)))/Ti(i))**(-g_ssl/(R*Ti_p(i)))
    end if
    
    !density
    rho = Pz/(R*T)

    !speed of sound 
    a = sqrt(gamma*R*T)

    std_atm(1) = Z
    std_atm(2) = T
    std_atm(3) = Pz
    std_atm(4) = rho
    std_atm(5) = a 

end subroutine std_atm_SI

subroutine vary_alt_SI(Z, T, P, rho, a)
    !3.13.4 Use your code developed for Problem 3.13.3 to compute the atmospheric
    !properites at altitudes varying from 0 to 90,000 m in increments of
    !5000 m. Compare your results to the solutions given in Appendix B
    implicit none 
    real, intent(inout) :: Z, T, P, rho, a
    real, dimension(:), allocatable :: geometric_alt
    integer :: i,k
    k=1
    geometric_alt = [ (i, i=0, 90000, 5000)]

    do k=1, size(geometric_alt)
        call std_atm_SI(geometric_alt(k), Z, T, P, rho, a)
        write(*,'(7ES25.11)') geometric_alt(k),Z,T,P,rho,a
    end do
end subroutine vary_alt_SI 

subroutine std_atm_English(H, Z, T, P, rho, a)
    !3.13.5 Write a function to compute standard atmospheric properites in English
    !units as a funciton of geometric altitude. Given a geometric altitude in feet,
    !your code must return the geopotential altitude (ft), temperature (R), pressure
    !(lbf/ft^2), density (slugs/ft^3), and speed of sound (ft/s). The easiest way to
    !do this is to create a wrapper to the function written for Problem 3.13.3. The
    !wrapper would accept the altitude in feet, convert the altitude to meters, call
    !the function in SI units, and convert the results to English units.
    implicit none
    real, intent(in) :: H
    real :: H_SI
    real, intent(inout) :: Z, T, P, rho, a

    H_SI = H*ft_to_m
    call std_atm_SI(H_SI, Z, T, P, rho, a)
    Z = Z/ft_to_m
    T = T*R_to_K
    P = P/Pen_to_PSI
    rho = rho/densen_to_densSI
    a = a/ft_to_m
end subroutine std_atm_English

subroutine vary_alt_english(Z, T, P, rho, a)
    !3.13.6 Use your code developed for Problem 3.13.5 to compute the atmospheric
    !properites at altitudes varying from 0 to 200,000 ft in increments of
    !10000 ft. Compare your results to the solutions given in Appendix B
    implicit none 
    real, intent(inout) :: Z, T, P, rho, a
    real, dimension(:), allocatable :: geometric_alt
    integer :: i,k
    k=1
    geometric_alt = [ (i, i=0, 200000, 10000)]

    do k=1, size(geometric_alt)
        call std_atm_english(geometric_alt(k), Z, T, P, rho, a)
        write(*,'(7ES25.11)') geometric_alt(k),Z,T,P,rho,a
    end do
end subroutine vary_alt_english 

! function runge_kutta(t0, y0, del_t) result(y1)
!     !to = time at begining of step
!     !y0 = function value at beginning of step
!     !del_t = step size 
!     !y = function value at end of step
!     implicit none 
!     real, intent(in) :: t0, y0, del_t 
!     real :: y1
!     real :: k1, k2, k3, k4 

!     k1 = differential_equations(t0, y0)
!     k2 = differential_equations(t0 + del_t/2, y0 + k1*del_t/2)
!     k3 = differential_equations(t0 + del_t/2, y0 + k2*del_t/2)
!     k4 = differential_equations(t0 + del_t, y0 + k3*del_t)
    
!     y1 = y0 + (del_t/6)*(k1 + 2*k2 + 2*k3 + k4)
!     write(*,*) t0, y0, k1, k2, k3, k4, y1
! end function runge_kutta

! function differential_equations(t, y) result(dy_dt)
!     !t = time 
!     !y = value of y at time t 
!     !dy_dt = change in y with respecdt to t 
!     implicit none 
!     real, intent(in) :: t, y
!     real :: dy_dt 

!     dy_dt = 1 + tan(y)
! end function differential_equations

! subroutine test_main(t0, y0, del_t, y1, it) 
!     implicit none 
!     integer, intent(in) :: it
!     integer :: i 
!     real :: t0, y0, del_t
!     real, intent(inout) :: y1

!     do i=1, it
!         y1 = runge_kutta(t0, y0, del_t)
!         t0 = t0 + del_t 
!         y0 = y1 
!     end do
! end subroutine test_main 

function runge_kutta(t0, y0, del_t) result(y1)
    !to = time at begining of step
    !y0 = function value at beginning of step
    !del_t = step size 
    !y = function value at end of step
    implicit none 
    real, intent(in) :: t0, del_t
    real, dimension(:), intent(in) :: y0
    real, dimension(:), allocatable :: y1, k1, k2, k3, k4

    k1 = differential_equations(t0, y0)
    k2 = differential_equations(t0 + del_t/2, y0 + k1*del_t/2)
    k3 = differential_equations(t0 + del_t/2, y0 + k2*del_t/2)
    k4 = differential_equations(t0 + del_t, y0 + k3*del_t)
    
    y1 = y0 + (del_t/6)*(k1 + 2*k2 + 2*k3 + k4)
    write(*,*) t0, y0
end function runge_kutta

function differential_equations(t, y) result(dy_dt)
    !t = time 
    !y = value of y at time t 
    !dy_dt = change in y with respecdt to t 
    implicit none 
    real, intent(in) :: t
    real, dimension(:), intent(in) :: y
    real, dimension(:), allocatable :: dy_dt
    real :: x, z 

    x = y(1)
    z = y(2) 
    dy_dt = [t + z**2*sin(x), t + x*cos(z)]

end function differential_equations

! subroutine test_main(t0, y0, del_t, y1, it) 
!     implicit none 
!     integer, intent(in) :: it
!     integer :: i 
!     real :: t0, y0, del_t
!     real, intent(inout) :: y1

!     do i=1, it
!         y1 = runge_kutta(t0, y0, del_t)
!         t0 = t0 + del_t 
!         y0 = y1 
!     end do
! end subroutine test_main 

end module olsen_m

program main
use olsen_m

real :: t0, del_t
real, dimension(2) :: y0, dy_dt

t0 = 0.0
y0 = [0, 0]
del_t = 0.025
dy_dt = runge_kutta(t0, y0, del_t)

write(*,*) dy_dt

end program main
