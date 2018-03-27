subroutine forces(xf,grav,fict)
  
! Calculate both real (ie, gravity and whatnot) 
! and fictitious (ie, coriolis and centrifugal) forces.
!-------------------------------------------------------------

! GLOBALS
use global
use sweeps

IMPLICIT NONE

! LOCALS
integer :: n
real :: xf0, sinxf0
real, dimension(maxsweep) :: grav, fict, xf

!--------------------------------------------------------------

xf0 = xf0
 
if (sweep == 'x') then

   do n = nmin-4, nmax+5
      grav(n) = - GMP / xf(n)**2

      fict(n) = (w(n)*w(n)+v(n)*v(n))/xf(n)
   enddo

else if (sweep == 'y') then

   do n = nmin-4, nmax+5
      grav(n) = 0.0d0

      fict(n) = v(n)*v(n)*cos(xf(n))/(radius*sin(xf(n))) - u(n)*w(n)/radius
   enddo

else

   do n = nmin-4, nmax+5
      grav(n) = 0.0d0
      
      fict(n) = - u(n)*w(n)/radius*cos(theta) - u(n)*v(n)/radius*stheta 
   enddo

endif

return
end
