program quadratic_equation
 
    implicit none
    real :: a, b, c     ! coefficients
    real :: d           ! discriminant
 
    ! read a, b, c
    read(*, *) a, b, c
 
    ! calculate discriminant
    d = b**2 - 4*a*c
 
    if (d > 0) then
        ! two distinct solutions
        write(*, *) (-b + sqrt(d)) / (2*a), &
                    (-b - sqrt(d)) / (2*a)
    else if (d == 0) then
        ! one solution
        write(*, *) -b / (2*a)
    else
        ! no real solution
        write(*, *) 'no real solution'
    end if
 
end program quadratic_equation