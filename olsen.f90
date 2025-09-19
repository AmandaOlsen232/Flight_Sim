module olsen_m
use quat_math
use atmosphere
implicit none
real :: Weight=0.006, & ! [lbf]
        Diameter=0.13084, & ! [ft]
        R2 = 0.06542, & ! [ft]
        R1 = 0.06411, & ![ft]
        Skin_thickness=0.00131, & ! [ft]
        alt=200., & ! [ft]
        V_0=50., & ! [ft/s]
        alpha = 15. ![deg]

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
    I = M*(2./5.)*((R2**5 - R1**5)/(R2**3 - R1**3)) * reshape([1.0, 0.0, 0.0, &
                                                             0.0, 1.0, 0.0, &
                                                             0.0, 0.0, 1.0], shape=[3,3])

end subroutine mass_inertia

subroutine pseudo_aerodynamics(t, y, Fc, Mc)
    implicit none 
    real, intent(in) :: t 
    real, dimension(:), intent(in) :: y 
    real, intent(inout) :: Fc(3), Mc(3) 

    real :: Vb(3), V_mag, uc(3), volume, rho, Mass, Re, Z, Temp, P, a, mu, CD
    real, allocatable :: I(:,:) 
    
    call std_atm_English(alt, Z, Temp, P, rho, a, mu)
    
    V_mag = norm2(y(1:3))
    uc = y(1:3)/V_mag
    Re = 2*rho*V_mag*R2/mu
    
    if (0<Re .and. Re <=450000) then
        CD = (24/Re) + 6/(1+sqrt(Re)) + 0.4
    else if (450000<Re .and. Re<=560000) then 
        CD = 1.0e+29*Re**(-5.211)
    else if (560000<Re .and. Re<=14000000) then 
        CD = -2.0e-23*Re**3 - 1.0e-16*Re**2 + 9.0e-09*Re + 0.069
    else if (Re > 14000000) then 
        CD = 0.12
    end if 

    Fc = -0.5*rho*(V_mag**2)*PI*(R2**2)*CD*uc
    Mc = 0
    
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
real, dimension(3) :: F, M
!y_array = [u, v, w, p, q, r, x, y, z, eo, ex, ey, ez]
t0 = 0.0
tf = 10.0 
del_t = 0.01
y0 = [50., 0., 0., 0., 0., 0., 0., 0., -200., 1., 0., 0., 0.]
call quat_norm(y0(10:13))
test = simulation_main(t0, tf, del_t, y0)


end program main
