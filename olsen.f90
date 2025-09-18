module olsen_m
use quat_math
use atmosphere
implicit none
real :: Weight=0.006, & ! [lbf]
        Diameter=1., & ! [ft]
        Thickness=0.00131, & ! [ft]
        alt=0., & ! [ft]
        vel=200., &! [ft/s]
        alpha=15., & ! [deg]
        beta=10., &
        P(3)=[1.,2.,3.] ! [deg]
contains 

function runge_kutta(t0, y0, del_t) result(y1)
    implicit none
    real, intent(in) :: t0, del_t
    real, dimension(:), intent(in) :: y0
    real, dimension(:), allocatable :: y1, k1, k2, k3, k4
    k1 = differential_equations(t0, y0)
    k2 = differential_equations(t0 + del_t/2, y0 + k1*del_t/2)
    k3 = differential_equations(t0 + del_t/2, y0 + k2*del_t/2)
    k4 = differential_equations(t0 + del_t, y0 + k3*del_t)
    y1 = y0 + (del_t/6)*(k1 + 2*k2 + 2*k3 + k4)

end function runge_kutta

subroutine mass_inertia(t, y, M, I)
    implicit none 
    real, intent(in) :: t
    real, dimension(:), intent(in) :: y 
    real, intent(inout) :: M 
    real, dimension(:,:), intent(inout), allocatable :: I
    real :: r1, r2, g 

    allocate(I(3,3))
    g = gravity_English(0.)
    M = Weight/g 
    r2 = Diameter/2. 
    r1 = r2 - Thickness
    I = M*(2./5.)*((r2**5 - r1**5)/(r2**3 - r1**3)) * reshape([1.0, 0.0, 0.0, &
                                                             0.0, 1.0, 0.0, &
                                                             0.0, 0.0, 1.0], shape=[3,3])
    write(*,*) M, M*(2./5.)*((r2**5 - r1**5)/(r2**3 - r1**3))
end subroutine mass_inertia

subroutine pseudo_aerodynamics(t, y, Fb, Mb)
    implicit none 
    real, intent(in) :: t 
    real, dimension(:), intent(in) :: y 
    real, dimension(3), intent(inout) :: Fb, Mb 
    real, dimension(3,3) :: R 
    real :: eo, ex, ey, ez 
    real :: Vb(3), Mc(3)

    Mc = 0.

    eo = y(10)
    ex = y(11)
    ey = y(12)
    ez = y(13)
    R = reshape([ex**2+eo**2-ey**2-ez**2, 2*(ex*ey-ez*eo), 2*(ex*ez+ey*eo), &
                 2*(ex*ey+ez*eo), ey**2+eo**2-ex**2-ez**2, 2*(ey*ez-ex*eo), &
                 2*(ex*ez-ey*eo), 2*(ey*ez+ex*eo), ez**2+eo**2-ex**2-ey**2], shape=[3,3])
    write(*,*) R
    Vb =vel*[cos(alpha*PI/180)*cos(beta*PI/180), sin(beta*PI/180), sin(alpha*PI/180)*cos(beta*PI/180)]
    write(*,*) Vb
end subroutine pseudo_aerodynamics

function inv_3d(A) result(A_inv)
    implicit none 
    real, dimension(3,3), intent(in) :: A
    real, dimension(3,3) :: A_inv, A_tilde 
    real :: det_A  
    det_A = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) &
            - A(1,3)*A(2,2)*A(3,1) - A(1,2)*A(2,1)*A(3,3) - A(1,1)*A(2,3)*A(3,2)

    A_tilde = reshape([A(2,2)*A(3,3)-A(2,3)*A(3,2), A(1,3)*A(3,2)-A(1,2)*A(3,3), A(1,2)*A(2,3)-A(1,3)*A(2,2), &
                       A(2,3)*A(3,1)-A(2,1)*A(3,3), A(1,1)*A(3,3)-A(1,3)*A(3,1), A(1,3)*A(2,1)-A(1,1)*A(2,3), &
                       A(2,1)*A(3,2)-A(2,2)*A(3,1), A(1,2)*A(3,1)-A(1,1)*A(3,2), A(1,1)*A(2,2)-A(1,2)*A(2,1)], shape=[3,3])

    A_inv = A_tilde/det_A
end function inv_3d

function matmul_3d_vec(A,B) result(C)
    implicit none 
    real, dimension(3,3), intent(in) :: A
    real, dimension(3), intent(in) :: B
    real, dimension(3) :: C 

    C = [A(1,1)*B(1) + A(1,2)*B(2) + A(1,3)*B(3), &
         A(2,1)*B(1) + A(2,2)*B(2) + A(2,3)*B(3), &
         A(3,1)*B(1) + A(3,2)*B(2) + A(3,3)*B(3)]
end function matmul_3d_vec

function differential_equations(t, y_array) result(dy_dt)
    implicit none
    real, intent(in) :: t
    real, dimension(:), intent(in) :: y_array
    real, dimension(:) :: dy_dt(13), one(3), two(3), three(4), four(4)
    real :: g, Mass, u, v, w, p, q, r, x, y, z, eo, ex, ey, ez, Fxb, Fyb, Fzb, Ixxb, Ixyb, Ixzb, Iyyb, Iyzb, Izzb, Mxb, Myb, Mzb
    real, dimension(:,:), allocatable :: I
    real, dimension(3) :: F, M

    !y_array = [u, v, w, p, q, r, x, y, z, eo, ex, ey, ez]

    !uncomment for 5.9.4
    call pseudo_aerodynamics(t, y_array, F, M)
    call mass_inertia(t, y_array, Mass, I)
    g = gravity_English(alt)

    Fxb = F(1)
    Fyb = F(2) 
    Fzb = F(3) 
    Mxb = M(1)
    Myb = M(2)
    Mzb = M(3) 

    Ixxb = I(1,1)
    Ixyb = I(2,1)
    Ixzb = I(3,1)
    Iyyb = I(2,2)
    Iyzb = I(3,2)
    Izzb = I(3,3)
    u = y_array(1) 
    v = y_array(2)
    w = y_array(3)
    p = y_array(4)
    q = y_array(5)
    r = y_array(6)
    x = y_array(7)
    y = y_array(8)
    z = y_array(9)
    eo = y_array(10)
    ex = y_array(11)
    ey = y_array(12)
    ez= y_array(13)
    
    one = [(1/Mass)*Fxb + g*2*(ex*ez - ey*eo) + (r*v - q*w), &
             (1/Mass)*Fyb + g*2*(ey*ez + ex*eo) + (p*w - r*u), &
             (1/Mass)*Fzb + g*(ez**2 + eo**2 - ex**2 - ey**2) + (q*u - p*v)]

    two = [Mxb + (Iyyb - Izzb)*q*r + Iyzb*(q**2 - r**2) + Ixzb*p*q - Ixyb*p*r, &
           Myb + (Izzb - Ixxb)*p*r + Ixzb*(r**2 - p**2) + Ixyb*q*r - Iyzb*p*q, &
           Mzb + (Ixxb - Iyyb)*p*q + Ixyb*(p**2 - q**2) + Iyzb*p*r - Ixzb*q*r]
    two = matmul_3d_vec(inv_3d(I), two)

    three = quat_mult([eo, ex, ey, ez], quat_mult([0., u, v, w], [eo, -ex, -ey, -ez])) !add wind vec 

    four = 0.5*[-ex*p -ey*q -ez*r, &
                eo*p -ez*q +ey*r, &
                ez*p +eo*q -ex*r, &
                -ey*p +ex*q +eo*r]

    dy_dt = [one(1), one(2), one(3), two(1), two(2), two(3), three(2), three(3), three(4), four(1), four(2), four(3), four(4)]
end function differential_equations

function simulation_main(t0, tf, del_t, y0) result(y1)
    implicit none 
    real :: t0
    real, intent(in) :: tf, del_t
    real, dimension(:) :: y0
    real, dimension(:), allocatable :: y1

    do while (t0 < tf)
        write(*,*) y0
        y1 = runge_kutta(t0, y0, del_t)
        call quat_norm(y1(10:13))
        t0 = t0 + del_t
        y0 = y1
        
    end do
end function simulation_main

end module olsen_m


program main
use olsen_m


real :: t0, tf, del_t 
real, dimension(13) :: y0, test 
! real :: M 
real, dimension(:,:), allocatable :: F, M
!y_array = [u, v, w, p, q, r, x, y, z, eo, ex, ey, ez]
! t0 = 0.0
! tf = 10.0 
! del_t = 1.
! y0 = [50., 0., 0., 0., 0., 0., 0., 0., -20., 1., 0., 0., 0.]
! call quat_norm(y0(10:13))
! test = simulation_main(t0, tf, del_t, y0)
call pseudo_aerodynamics(0.1, [0.,0.], F, M)

end program main
