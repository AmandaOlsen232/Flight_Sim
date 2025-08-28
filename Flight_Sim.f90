module my_functions
implicit none
! integer, parameter :: dp = selected_real_kind(16,307)

contains
!Write a function named quat_mult that accepts two arrays of length 
!   4 representing 2 quaternions and returns the quaternion product
function quat_mult(A, B) result(answer)
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


!Write a function named quat_base_to_dependent to compute a vector
!   transformation using the Euler-Rodrigues quaternion. The function
!   should accept an array of length three containing the three 
!   components of a vector in the base coordinate system as well as
!   an array of length four containing the four components of the
!   quaternion that represent the orientation of the dependent 
!   coordinate system. This function should return an array of length
!   three containing the vector components in the dependent 
!   coordinate system.
function quat_base_to_dependent(base, E) result(v2)
    
    implicit none
    real, dimension(3), intent(in) :: base !array of length three containing the three components of a vector in the base coordinate system
    real, dimension(4), intent(in) :: E !array of length four containing the four components of the quaternion that represent the orientation of the dependent coordinate system
    real, dimension(4) :: e_rod, T
    real, dimension(3) :: v2
    real :: pi=3.1415926535897932385
    real :: theta_deg, theta_rad, Ex, Ey, Ez, vx1, vy1, vz1, erod0, erodx, erody, erodz, T0, Tx, Ty, Tz

    theta_deg = E(1)
    theta_rad = theta_deg*pi/180
    Ex = E(2)
    Ey = E(3)
    Ez = E(4)
    
    vx1 = base(1)
    vy1 = base(2)
    vz1 = base(3)

    e_rod = [cos(theta_rad/2), &
            Ex*sin(theta_rad/2), &
            Ey*sin(theta_rad/2), &
            Ez*sin(theta_rad/2)]
    
    T = [-vx1*erodx - vy1*erody - vz1*erodz, &
         vx1*erod0 + vy1*erodz - vz1*erody, &
         -vx1*erodz + vy1*erod0 + vz1*erodz, &
          vx1*erody - vy1*erodx + vz1*erod0]
    
    T0 = T(1)
    Tx = T(2)
    Ty = T(3)
    Tz = T(4)

    v2 = [erod0*Tx - erodx*T0 - erody*Tz + erodz*Ty, &
          erod0*Ty + erodx*Tz - erody*T0 - erodz*Tx, &
         erod0*Tz - erodx*Ty - erody*Tx - erodz*T0]    

end function quat_base_to_dependent


end module my_functions

program Flight_Sim
use my_functions

implicit none
real, dimension(3) :: test

test = quat_base_to_dependent([1., 2., 3.], [45., 2., 3., 4.])
write(*,*) test
end program Flight_Sim