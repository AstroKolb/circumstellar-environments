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
real :: xmin, xmax, ymin, ymax, zmin, zmax
real :: ridt, rodt, xvel, yvel, zvel, width, widthz, widthy
real :: zoom, rad
real :: perturb

!--------------------------------------------------------------------------------
! Set up grid geometry: This Yin-Yang code is strictly 3D spherical

ndim   = 3
ngeomx = 2 
ngeomy = 4
ngeomz = 5

pi   = 2.0 * asin(1.0)
ymin =  0.245 * pi
ymax =  0.755 * pi
zmin = -0.755 * pi
zmax = +0.755 * pi

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
time   = 1.0e-03
xmax = -1.25*time*log(time)
xmin = xmax * 0.5

gam    = 5.0/3.0

zoom = 0.0
rad = 0.0

gamm   = gam - 1.0d0
!=======================================================================
! Set up grid coordinates and parabolic coefficients

call grid(imax,xmin,xmax,zxa,zxc,zdx)
call grid(jmax,ymin,ymax,zya,zyc,zdy)
call grid(kmax,zmin,zmax,zza,zzc,zdz)

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


do k = 1, ks
 do j = 1, js
  do i = 1, imax/2
    zux(i,j,k) = zxc(i) / time
    zro(i,j,k) = 4.0*sqrt(3.0) * exp(-2.0*sqrt(3.0)*zux(i,j,k)) / time**3
    zpr(i,j,k) = 1.0e-10 * zro(i,j,k)*zux(i,j,k)**2
    zuy(i,j,k) = 0.0
    zuz(i,j,k) = 0.0
    zfl(i,j,k) = 0.0
    zcl(i,j,k) = 0.0
  enddo
  do i = imax/2+1, imax
    zux(i,j,k) = 0.0
    zro(i,j,k) = 1.0
    zpr(i,j,k) = 1.0e-10
    zuy(i,j,k) = 0.0
    zuz(i,j,k) = 0.0
    zfl(i,j,k) = 0.0
    zcl(i,j,k) = 0.0
  enddo
 enddo
enddo

! break the spherical symmetry with some small density wiggles
do k = 1, ks
 do j = 1, js
  perturb = 1.0 + 0.10*cos(20.0*zyc(j+jcol*js))*cos(20.0*zzc(k+krow*ks))
  do i = imax/2-20, imax/2
    zro(i,j,k) = zro(i,j,k) * perturb
  enddo
 enddo
enddo

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


subroutine grid( nzones, xmin, xmax, xa, xc, dx  )

! Create grid to cover physical size from xmin to xmax
! number of physical grid zones is nzones
!
! xa(1) is left boundary location - at xmin
! xa(nzones+1) is right boundary location - at xmax
!----------------------------------------------------------------------

! LOCALS
integer :: nzones, n
real, dimension(nzones) :: xa, dx, xc
real :: dxfac, xmin, xmax

!=======================================================================

dxfac = (xmax - xmin) / real(nzones)
do n = 1, nzones
  xa(n) = xmin + (n-1)*dxfac
  dx(n) = dxfac
  xc(n) = xa(n) + 0.5*dx(n)
enddo

return
end
