subroutine sweepx2(xexpand)

! This subroutine performs 1D hydro sweeps in the X direction, looping over j by k rows.
!-----------------------------------------------------------------------

! GLOBALS
use zone
use global
use sweeps

! LOCALS
integer :: i, j, k, n, m, myj, myk
real :: rad, xexpand, costhetan

!-----------------------------------------------------------------------

sweep  = 'x'
ngeom  = ngeomx
nmin   = 7
nmax   = imax + 6
radius = 1.0d0

do i = 1, imax
  n = i + 6
  xa0(n) = zxa(i) * (1.0 + xexpand)
  dx0(n) = zdx(i) * (1.0 + xexpand)
enddo

! Now Loop over each row...

do k = 1, ks
 myk  = k+krow*ks
 cphi = cos(zzc(myk))
 sphi = sin(zzc(myk))

 do j = 1, js
   myj    = j+jcol*js
   ctheta = cos(zyc(myj))
   stheta = sin(zyc(myj))

   if (yin) then
     costhetan = ctheta
   else
     costhetan = stheta * sphi
   endif

   ! Put state variables into 1D arrays, padding with 6 ghost zones
   do i = 1, isy
    do m = 1, npey
     n = i + isy*(m-1) + 6
     r  (n) = recv2(1,k,i,j,m)
     p  (n) = recv2(2,k,i,j,m)
     u  (n) = recv2(3,k,i,j,m)
     v  (n) = recv2(4,k,i,j,m)
     w  (n) = recv2(5,k,i,j,m)
     f  (n) = recv2(6,k,i,j,m)
     c  (n) = recv2(7,k,i,j,m)
    enddo
   enddo

   do i = 1,imax
     n = i + 6
     xa (n) = zxa(i)
     dx (n) = zdx(i)
     p  (n) = max(smallp,p(n))
     e  (n) = p(n)/(r(n)*gamm)+0.5d0*(u(n)**2+v(n)**2+w(n)**2)
   enddo

   ! Set boundary conditions, compute volume elements, evolve flow, then remap

   do n = 1, 6
    dx (nmin-n)= dx (nmin)
    dx0(nmin-n)= dx0(nmin)
    xa (nmin-n)= xa (nmin-n+1) - dx (nmin-n)
    xa0(nmin-n)= xa0(nmin-n+1) - dx0(nmin-n)

    rad = xa(nmin-n) + 0.5d0*dx(nmin-n)

    u (nmin-n) = uin(n)
    v (nmin-n) = 0.0
    w (nmin-n) = 0.0

    r (nmin-n) = mdot / (4.0*pi*u(nmin-n)*rad**2)
    p (nmin-n) = r(nmin-n)*kB*Tmp/mp

    e (nmin-n) = p(nmin-n)/(gamm*r(nmin-n))+0.5d0*(u(nmin-n)**2+v(nmin-n)**2+w(nmin-n)**2)
    f (nmin-n) = 0.0
    c (nmin-n) = 0.0

    dx (nmax+n)= dx (nmax)
    dx0(nmax+n)= dx0(nmax)
    xa (nmax+n)= xa (nmax+n-1) + dx (nmax+n-1)
    xa0(nmax+n)= xa0(nmax+n-1) + dx0(nmax+n-1)

    rad = xa(nmax+n) + 0.5d0*dx(nmax+n)

    u (nmax+n) = u(nmax)
    v (nmax+n) = 0.0
    w (nmax+n) = 0.0

    r (nmax+n) = r(nmax)
    p (nmax+n) = p(nmax)
    e (nmax+n) = p(nmax+n)/(gamm*r(nmax+n))+0.5d0*(u(nmax+n)**2+v(nmax+n)**2+w(nmax+n)**2)
    f (nmax+n) = 0.0
    c (nmax+n) = 0.0
   enddo

   if (ongrid(myj,myk) == 0) call ppmlr 

   ! Put updated values back into 3D arrays, dropping ghost zones
   do i = 1, imax
     n = i + 6
     zro(i,j,k) = r(n)
     zpr(i,j,k) = p(n)
     zux(i,j,k) = u(n)
     zuy(i,j,k) = v(n)
     zuz(i,j,k) = w(n)
     zfl(i,j,k) = f(n)
     zcl(i,j,k) = c(n)
   enddo

 enddo
enddo
      
do i = 1, imax
  n = i + 6
  zxa(i) = xa0(n)
  zdx(i) = dx0(n)
  zxc(i) = xa0(n) + 0.5*dx0(n)
enddo
 
return
end

