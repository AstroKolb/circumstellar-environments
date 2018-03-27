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
real :: xf0, sinxf0, radS
real, dimension(maxsweep) :: grav, fict, xf

!--------------------------------------------------------------
 
if (sweep == 'x') then

   if (yin) then
      do n = nmin-4, nmax+5
         radS = sqrt(xf(n)**2 + sep**2 - 2.0*xf(n)*sep*stheta*cphi)

         grav(n) = - GMP / xf(n)**2                &
                   - GMS / radS**3 * (xf(n)-sep*stheta*cphi)

         fict(n) = + (w(n)*w(n)+v(n)*v(n))/xf(n)   &
                   + omega**2*xf(n)*stheta**2      &
                   - omega**2*rcm*stheta*cphi      &
                   + omega*2.0*stheta*w(n)         &
                   + 0.0
      enddo
   else
      do n = nmin-4, nmax+5
         radS = sqrt(xf(n)**2 + sep**2 + 2.0*xf(n)*sep*stheta*cphi)

         grav(n) = - GMP / xf(n)**2                               &
                   - GMS / radS**3 * (xf(n)+sep*stheta*cphi)

         fict(n) = + (w(n)*w(n)+v(n)*v(n))/xf(n)                  &
                   + omega**2*xf(n)*(ctheta**2*sphi**2+cphi**2)   &
                   + omega**2*rcm*stheta*cphi                     &
                   - omega*2.0*(sphi*ctheta*w(n)-cphi*v(n))       &
                   + 0.0
      enddo
   endif     

else if (sweep == 'y') then

   if (yin) then
      do n = nmin-4, nmax+5
         stheta = sin(xf(n))
         ctheta = cos(xf(n))
         radS = sqrt(radius**2 + sep**2 - 2.0*radius*sep*stheta*cphi)

         grav(n) = + 0.0d0                                     &
                   + GMS / radS**3 * sep*ctheta*cphi

         fict(n) = + v(n)*v(n)*cos(xf(n))/(radius*sin(xf(n)))  &
                   - u(n)*w(n)/radius                          &
                   + omega**2*radius*stheta*ctheta             &
                   - omega**2*rcm*ctheta*cphi                  &
                   + omega*2.0*ctheta*v(n)                     &
                   + 0.0
      enddo
   else
      do n = nmin-4, nmax+5
         stheta = sin(xf(n))
         ctheta = cos(xf(n))
         radS = sqrt(radius**2 + sep**2 + 2.0*radius*sep*stheta*cphi)

         grav(n) = + 0.0d0                                     &
                   - GMS / radS**3 * sep*ctheta*cphi

         fict(n) = + v(n)*v(n)*cos(xf(n))/(radius*sin(xf(n)))  &
                   - u(n)*w(n)/radius                          &
                   - omega**2*radius*sphi**2*ctheta*stheta     &
                   + omega**2*rcm*ctheta*cphi                  &
                   - omega*2.0*(cphi*w(n)-sphi*stheta*v(n))    &
                   + 0.0
      enddo
   endif

else

   if (yin) then
      do n = nmin-4, nmax+5
         sphi = sin(xf(n))
         cphi = cos(xf(n))
         radS = sqrt(radius**2/stheta**2 + sep**2 - 2.0*radius*sep*cphi)

         grav(n) = + 0.0d0                                  &
                   - GMS / radS**3 * sep*sphi
         
         fict(n) = - u(n)*w(n)/radius*cos(theta)            &
                   - u(n)*v(n)/radius*stheta                &
                   + 0.0                                    &
                   + omega**2*rcm*sphi                      &
                   - omega*2.0*(stheta*v(n)+ctheta*w(n))    &
                   + 0.0
      enddo
   else
      do n = nmin-4, nmax+5
         sphi = sin(xf(n))
         cphi = cos(xf(n))
         radS = sqrt(radius**2/stheta**2 + sep**2 + 2.0*radius*sep*cphi)

         grav(n) = + 0.0d0                                           &
                   + GMS / radS**3 * sep*sphi
         
         fict(n) = - u(n)*w(n)/radius*cos(theta)                     &
                   - u(n)*v(n)/radius*stheta                         &
                   - omega**2*radius*cphi*sphi                       &
                   - omega**2*rcm*sphi                               &
                   - omega*2.0*(stheta*sphi*w(n)-ctheta*sphi*v(n))   &
                   + 0.0
      enddo
   endif

endif

return
end
