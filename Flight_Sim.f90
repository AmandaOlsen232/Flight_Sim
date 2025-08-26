module my_functions
implicit none
integer, parameter :: dp = selected_real_kind(16,307)

contains
!Write a function named quat_mult that accepts two arrays of length 
!   4 representing 2 quaternions and returns the quaternion product
function quat_mult(A, B) result(answer)
    implicit none

    real(dp), dimension(4), intent(in) :: A, B
    real(dp) :: A0, Ax, Ay, Az, B0, Bx, By, Bz
    real(dp), dimension(4) :: answer

    A0 = A(1)
    Ax = A(2)
    Ay = A(3)
    Az = A(4)

    B0 = B(1)
    Bx = B(2)
    By = B(3)
    Bz = B(4)

    answer(1) = A0*B0 - Ax*Bx - Ay*By - Az*Bz
    answer(2) = A0*Bx + Ax*B0 + Ay*Bz - Az*By
    answer(3) = A0*By - Ax*Bz + Ay*B0 + Az*Bx
    answer(4) = A0*Bz + Ax*By - Ay*Bx + Az*B0

    write(*,*) answer
    
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
! function quat_base_to_dependent() result()

! end function quat_base_to_dependent

end module my_functions

program Flight_Sim
use my_functions

implicit none
real(dp), dimension(4) :: test

test = quat_mult([1._dp, 2._dp, 3._dp, 4._dp], [1._dp, 2._dp, 3._dp, 4._dp])
end program Flight_Sim