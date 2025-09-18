module atmosphere 
use quat_math
real, parameter :: g_ssl = 9.80665 ! [m/s^2]
real, parameter :: REz = 6356766 ! [m]
real, parameter :: gamma = 1.4
real, parameter :: Ts = 273.15 ! [K]
real, parameter :: mu_s = 1.716e-05 ! [kg/m-s]
real, parameter :: Ks = 110.4  
real, dimension(9), parameter :: Zi = [0., 11000., 20000., 32000., 47000., 52000., 61000., 79000., 90000.] ! [m]
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

    g_ssl_english = g_ssl/0.3048
    REz_english = REz/0.3048
    g = g_ssl_english*(REz_english/(REz_english+H))**2
end function gravity_English

subroutine std_atm_SI(H, Z, T, Pz, rho, a, mu)
    !3.13.3 Write a function to compute standard-atmospheric properties in SI units as
    !a function of geometric altitude. Given a geometric altitude in meters, your
    !code must return the geopotential altitude (m), temperature (K), pressure 
    !(N/m^2), density (kg/m^3), and speed of sound (m/s).
    implicit none
    real, intent(inout) :: Z, T, Pz, rho, a, mu
    real, intent(in) :: H
    real :: Pi
    real, dimension(6) :: std_atm
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

    if (abs(Ti_p(i)) < TOLERANCE) then
        Pz = Pi*exp(-g_ssl*(Z - Zi(i))/(287.0528*Ti(i)))
    else 
        Pz = Pi*((Ti(i) + (Ti_p(i))*(Z - Zi(i)))/Ti(i))**(-g_ssl/(287.0528*Ti_p(i)))
    end if
    
    !density
    rho = Pz/(287.0528*T)

    !speed of sound 
    a = sqrt(gamma*287.0528*T)

    !viscosity 
    mu = mu_s*((Ts + Ks)/(T + Ks))*(T/Ts)**(3./2.)

    std_atm(1) = Z
    std_atm(2) = T
    std_atm(3) = Pz
    std_atm(4) = rho
    std_atm(5) = a 
    std_atm(6) = mu 

end subroutine std_atm_SI

subroutine vary_alt_SI(Z, T, P, rho, a, mu)
    !3.13.4 Use your code developed for Problem 3.13.3 to compute the atmospheric
    !properites at altitudes varying from 0 to 90,000 m in increments of
    !5000 m. Compare your results to the solutions given in Appendix B
    implicit none 
    real, intent(inout) :: Z, T, P, rho, a, mu
    real, dimension(:), allocatable :: geometric_alt
    integer :: i,k
    k=1
    geometric_alt = [ (i, i=0, 90000, 5000)]

    do k=1, size(geometric_alt)
        call std_atm_SI(geometric_alt(k), Z, T, P, rho, a, mu)
        write(*,'(7ES25.11)') geometric_alt(k),Z,T,P,rho, mu
    end do
end subroutine vary_alt_SI 

subroutine std_atm_English(H, Z, T, P, rho, a, mu)
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
    real, intent(inout) :: Z, T, P, rho, a, mu

    H_SI = H*0.3048
    call std_atm_SI(H_SI, Z, T, P, rho, a, mu)
    Z = Z/0.3048
    T = T*1.8
    P = P/47.880258
    rho = rho/515.379
    a = a/0.3048
    mu = (mu/14.5939029)*0.3048
end subroutine std_atm_English

subroutine vary_alt_english(Z, T, P, rho, a, mu)
    !3.13.6 Use your code developed for Problem 3.13.5 to compute the atmospheric
    !properites at altitudes varying from 0 to 200,000 ft in increments of
    !10000 ft. Compare your results to the solutions given in Appendix B
    implicit none 
    real, intent(inout) :: Z, T, P, rho, a, mu
    real, dimension(:), allocatable :: geometric_alt
    integer :: i,k
    k=1
    geometric_alt = [ (i, i=0, 200000, 10000)]

    do k=1, size(geometric_alt)
        call std_atm_english(geometric_alt(k), Z, T, P, rho, a, mu)
        write(*,'(7ES25.11)') geometric_alt(k),Z,T,P,rho,a,mu
    end do
end subroutine vary_alt_english 

end module atmosphere