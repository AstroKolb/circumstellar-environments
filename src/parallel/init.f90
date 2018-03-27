subroutine init

! Ejecta-Driven SNR with exponential density profile
! 07oct11 blondin
!=======================================================================

! GLOBALS
use global
use zone

include 'mpif.h'

! LOCALS
integer :: i, j, k, mpierr
integer novery, noverz
real :: xmin, xmax, ymin, ymax, zmin, zmax, dely, delz
real :: ridt, rodt, xvel, yvel, zvel, width, widthz, widthy
real :: zoom, rad, Z1
real :: tmpu, RHS

!--------------------------------------------------------------------------------
! Set up grid geometry: This Yin-Yang code is strictly 3D spherical

ndim   = 3
ngeomx = 2 
ngeomy = 4
ngeomz = 5

novery = 6      ! define zone overlap
noverz = 6

pi   = 2.0*asin(1.0)

dely =     pi/2.0/float(jmax-novery)*float(novery)/2.0
delz = 3.0*pi/2.0/float(kmax-noverz)*float(noverz)/2.0

ymin =      pi/4.0 - dely
ymax =  3.0*pi/4.0 + dely
zmin = -3.0*pi/4.0 - delz
zmax =  3.0*pi/4.0 + delz

! set time and cycle counters
time   = 0.0
timep  = 0.
timem  = 0.
ncycle = 0
ncycp  = 0
ncycd  = 0
ncycm  = 0
nfile  = 1000

!======================================================================
! Set up parameters for exponential ejecta; initial scaling determined by time
time   = 0.0

gam    = 1.001!5.0/3.0
gamm   = gam - 1.0d0


GMP = 3.12 * GM
GMS = 9.0 * GM
sep = 5.e3 * Rsun

Tmp  = 3.0e3
mdot = 1.0e-6  * Msun / 3.15e7

xmin = 1000*Rsun
xmax = 50000*Rsun


! parker wind calculations
cs2   = gam * kB * Tmp / mp
uc    = cs2**0.5
rc    = GMP / (2.0 * cs2)
rho   = (mdot*Msun/3.15e7)/(4.0*pi*rc**2*uc)
capI  = rho*uc*rc**2


! rotational parameters
omega = ((GMP+GMS)/sep**3)**0.5
rcm   = sep*GMS/(GMP+GMS)
opd   = 2.0*pi/omega


Z1 = (xmax/xmin)**(1.0/float(imax)) - 1.0

zoom = Z1

!=======================================================================
! Set up grid coordinates and parabolic coefficients

call grid(imax,xmin,xmax,zxa,zxc,zdx,zoom)
call grid(jmax,ymin,ymax,zya,zyc,zdy,0.00)
call grid(kmax,zmin,zmax,zza,zzc,zdz,0.00)

!=======================================================================
! Log parameters of problem in history file

if (mype == 0) then
  write (8,"('SNR with exponential ejecta density profile in 3 dimensions.')")
  write (8,"('Grid dimensions: ',i3,' x ',i3,' x ',i3)") imax,jmax,kmax
  write (8,*) 
  write (8,"('   starting at time = ',1pe11.3)") time
  write (8,"('Adiabatic index, gamma = ',f7.3)") gam
  write (8,*) 
endif

! initialize grid:

if (.true.) then
   do k = 1, ks
      do j = 1, js
         do i = 1, imax
            zux(i,j,k) = 1.0e+04
            zuy(i,j,k) = 0.0
            zuz(i,j,k) = 0.0

            zro(i,j,k) = mdot / (4.0*pi*zxc(i)**2*zux(i,j,k))
            zpr(i,j,k) = zro(i,j,k)*kB*Tmp/mp
            zfl(i,j,k) = 0.0
            zcl(i,j,k) = 0.0
         enddo
      enddo
   enddo
endif

do i = 1, 6
   uin(i) = 1.0e+04
enddo

write(*,*) mdot*Msun/3.14e7/4.0/pi, capI



if (.false.) then
   do i = 1, imax

      RHS = 4.0 * (log(zxc(i)/rc) + rc/zxc(i)) - 3.0

      tmpu = uc
      do while (tmpu**2/uc**2 - log(tmpu**2/uc**2) < RHS)
         if (zxc(i) < rc) tmpu = 0.9999*tmpu
         if (zxc(i) > rc) tmpu = 1.0001*tmpu
      enddo

      do k = 1, ks
         do j = 1, js
            zux(i,j,k) = tmpu
            zuy(i,j,k) = 0.0
            zuz(i,j,k) = 0.0

            zro(i,j,k) = mdot / (4.0*pi*zux(i,j,k)*zxc(i)**2)
            zpr(i,j,k) = zro(i,j,k)*kB*Tmp/mp
            zfl(i,j,k) = 0.0
            zcl(i,j,k) = 0.0
         enddo
      enddo
   enddo

   do i = 1, 6
      rad = zxc(1) - zdx(1)*i

      RHS = 4.0 * (log(rad/rc) + rc/rad) - 3.0

      tmpu = uc
      do while (tmpu**2/uc**2 - log(tmpu**2/uc**2) < RHS)
         tmpu = 0.9999*tmpu
      enddo

      uin(i) = tmpu
   enddo
endif



if (.false.) then
   do k = 1, ks
      do j = 1, js
         do i = 1, imax
            zux(i,j,k) = 10.0
            zuy(i,j,k) = 0.0
            zuz(i,j,k) = 0.0
            
            zro(i,j,k) = 1.0 / zxc(i)**2
            zpr(i,j,k) = 1.0e-10 * zro(i,j,k)*zux(i,j,k)**2
            zcl(i,j,k) = 0.0

            zfl(i,j,k) = 0.0
         enddo
      enddo
   enddo
endif


!########################################################################
! Compute Courant-limited timestep

ridt = 0.
do k = 1, ks
 do j = 1, js
  do i = 1, imax
     widthy = zdy(j+jcol*js)
     widthz = zdz(k+krow*ks)
     widthy = widthy*zxc(i)
     widthz = widthz*zxc(i)*sin(zyc(j+jcol*js))
     width  = min(zdx(i),widthy,widthz)
     svel = sqrt(gam*zpr(i,j,k)/zro(i,j,k))/width
     xvel = abs(zux(i,j,k)) / zdx(i)
     yvel = abs(zuy(i,j,k)) / widthy
     zvel = abs(zuz(i,j,k)) / widthz
     ridt = max(xvel,yvel,zvel,svel,ridt)
  enddo
 enddo
enddo

call MPI_ALLREDUCE( ridt, rodt, 1, VH1_DATATYPE, MPI_MAX, MPI_COMM_WORLD, mpierr )
dt = courant / rodt

return
end

!#######################################################################


subroutine grid( nzones, xmin, xmax, xa, xc, dx, zoom )

! Create grid to cover physical size from xmin to xmax
! number of physical grid zones is nzones
!
! xa(1) is left boundary location - at xmin
! xa(nzones+1) is right boundary location - at xmax
!----------------------------------------------------------------------

! LOCALS
integer :: nzones, n
real, dimension(nzones) :: xa, dx, xc
real :: dxfac, xmin, xmax, zoom

!=======================================================================

if (zoom == 0.0d0) then
  dxfac = (xmax - xmin) / real(nzones)
  do n = 1, nzones
    xa(n) = xmin + (n-1)*dxfac
    dx(n) = dxfac
    xc(n) = xa(n) + 0.5*dx(n)
  enddo
else
  xa(1) = xmin
  dx(1) = zoom*xa(1)
  xc(1) = xa(1) + 0.5d0*dx(1)
  do n  = 2, nzones
    xa(n) = xa(n-1)+dx(n-1)
    dx(n) = zoom*xa(n)
    xc(n) = xa(n)+0.5d0*dx(n)
  enddo
endif

return
end
