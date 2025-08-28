module my_functions
implicit none
! integer, parameter :: dp = selected_real_kind(16,307)

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

function quat_base_to_dependent(v1, E) result(v2)
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
    real, dimension(4), intent(in) :: E !array of length four containing the four components of the quaternion that represent the orientation of the dependent coordinate system
    real, dimension(4) :: e_rod, e_rod_neg, T, v1_mult, v2_four
    real, dimension(3) :: v2
    real :: pi=3.1415926535897932385
    real :: theta_deg, theta_rad, Ex, Ey, Ez

    theta_deg = E(1)
    theta_rad = theta_deg*pi/180
    Ex = E(2)
    Ey = E(3)
    Ez = E(4)

    
    e_rod = [cos(theta_rad/2), &  !eq. 1.5.1
            Ex*sin(theta_rad/2), &
            Ey*sin(theta_rad/2), &
            Ez*sin(theta_rad/2)]
    
    e_rod_neg = [e_rod(1), &
                 -e_rod(2), &
                 -e_rod(3), &
                 -e_rod(4)]

    v1_mult = [0., &
               v1(1), &
               v1(2), &
               v1(3)]
   
    T = quat_mult(v1_mult, e_rod) !eq. 1.5.6
    v2_four = quat_mult(e_rod_neg, T) !eq. 1.5.7

    v2 = [v2_four(2), &
          v2_four(3), &
          v2_four(4)]

end function quat_base_to_dependent

function quat_dependent_to_base(v2, E) result(v1)
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
    real, dimension(4), intent(in) :: E
    real, dimension(3) :: v1

    real :: theta_deg, theta_rad, Ex, Ey, Ez
    real, dimension(4) :: e_rod, e_rod_neg, v2_mult, T, v1_four
    real :: pi=3.1415926535897932385

    theta_deg = E(1)
    theta_rad = theta_deg*pi/180
    Ex = E(2)
    Ey = E(3)
    Ez = E(4)

    
    e_rod = [cos(theta_rad/2), &  !eq. 1.5.1
            Ex*sin(theta_rad/2), &
            Ey*sin(theta_rad/2), &
            Ez*sin(theta_rad/2)]

    e_rod_neg = [e_rod(1), &
                 -e_rod(2), &
                 -e_rod(3), &
                 -e_rod(4)]

    v2_mult = [0., &
               v2(1), &
               v2(2), &
               v2(3)]

    T = quat_mult(v2_mult, e_rod_neg) !eq. 1.5.5
    v1_four = quat_mult(e_rod, T) !eq. 1.5.5
    v1 = [v1_four(2), &
          v1_four(3), &
          v1_four(4)] 
end function quat_dependent_to_base

function quat_norm(Q) result(Q_norm)
    !1.8.6 Write a function named quat_norm that accepts an array of length four 
    !   containing a quaternion, and returns the normalized quaternion. For this
    !   computation, use the exact solution given in Eq. (1.5.12). You may wish to
    !   overwrite the quaternion that was passed in to the function with the 
    !   normalized solution.
    implicit none
    real, dimension(4), intent(in) :: Q
    real, dimension(4) :: Q_norm
    real :: mag

    mag = sqrt(Q(1)*Q(1) + Q(2)*Q(2) + Q(3)*Q(3) + Q(4)*Q(4)) !eq. 1.5.12
    Q_norm = Q/mag

end function quat_norm

function euler_to_quat(Euler) result(e)
    !1.8.7 Write a function called euler_to_quat that accepts an array of length
    !   three containing the three Euler angles and returns an array of length 
    !   four containing values for the equivalent quaternion components.
    implicit none
    real, dimension(3), intent(in) :: Euler
    real, dimension(4) :: e 
    real :: phi_rad, theta_rad, psi_rad
    real :: pi=3.1415926535897932385

    phi_rad = Euler(1)*pi/180
    theta_rad = Euler(2)*pi/180
    psi_rad = Euler(3)*pi/180

    e = [cos(phi_rad/2)*cos(theta_rad/2)*cos(psi_rad/2) + sin(phi_rad/2)*sin(theta_rad/2)*sin(psi_rad/2), &
         sin(phi_rad/2)*cos(theta_rad/2)*cos(psi_rad/2) - cos(phi_rad/2)*sin(theta_rad/2)*sin(psi_rad/2), &
         cos(phi_rad/2)*sin(theta_rad/2)*cos(psi_rad/2) + sin(phi_rad/2)*cos(theta_rad/2)*sin(psi_rad/2), &
         cos(phi_rad/2)*cos(theta_rad/2)*sin(psi_rad/2) - sin(phi_rad/2)*sin(theta_rad/2)*cos(psi_rad/2)]

end function euler_to_quat

! function quat_to_euler() result()
!     !1.8.8 Write a function named quat_to_euler that accepts an array of length 
!     !   four containing the four components of a quaternion and returns an array 
!     !   of length three containing values for the equivalent Euler angles.
!     implicit none


! end function quat_to_euler

end module my_functions

program Flight_Sim
use my_functions

implicit none
real, dimension(4) :: test

test = euler_to_quat([45., 45., 45.])
write(*,*) test
end program Flight_Sim