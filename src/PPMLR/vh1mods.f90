module global
!=======================================================================
! global variables accessible by (most) anything
!-----------------------------------------------------------------------

 logical :: yin

 character(len=50) :: prefix              ! prefix for output filenames

 integer :: ndim
 integer :: ncycle, ncycp, ncycm, ncycd  ! cycle number
 integer :: nfile                        ! output file marker


 real :: time, dt, timem, timep, svel 
 real :: gam, pi, gamm
 real, parameter :: courant = 0.5        ! timestep fraction of courant limit
 real, parameter :: xwig = 0.00          ! fraction of a zone to wiggle grid for dissipation
 real, parameter :: smallp = 1.0e-45     ! Set small values to prevent divide by zero
 real, parameter :: smallr = 1.0e-45
 real, parameter :: small  = 1.0e-45


 ! constants
 real, parameter :: GRV  = 6.67e-08
 real, parameter :: Msun = 1.99e+33
 real, parameter :: Rsun = 6.95e+10
 real, parameter :: Lsun = 3.89e+33
 real, parameter :: AU   = 1.50e+13
 real, parameter :: mp   = 1.67e-24
 real, parameter :: kB   = 1.67e-16
 real, parameter :: GM   = GRV*Msun

 real :: uinflo, dinflo, vinflo, winflo, pinflo, einflo 
 real :: uotflo, dotflo, votflo, wotflo, potflo, eotflo

 ! parker wind parameters
 real :: Tmp, cs2, uc, rc, rho, capI
 
 ! rotational parameters
 real :: GMP, GMS, sep, mdot
 real :: omega, rcm, opd

 real, dimension(6) :: uin
      
end module 


module sweepsize
!=======================================================================
! Dimension of 1D sweeps.  maxsweep must be as long as the longest of the 
! 3D arrays PLUS the ghost zones:  maxsweep = max(imax,jmax,kmax) + 12
!----------------------------------------------------------------------

 integer, parameter :: maxsweep=1036 

end module sweepsize

module sweeps      
!=======================================================================
!  data structures used in 1D sweeps, dimensioned maxsweep (in sweep_size)
!-----------------------------------------------------------------------

use sweepsize

character(len=1) :: sweep                                    ! direction of sweep: x,y,z
integer :: nmin, nmax, ngeom                                 ! number of first and last real zone  
real, dimension(maxsweep) :: r, p, e, q, u, v, w, c          ! fluid variables
real, dimension(maxsweep) :: xa, xa0, dx, dx0, dvol          ! coordinate values
real, dimension(maxsweep) :: f, flat                         ! flattening parameter
real, dimension(maxsweep,5) :: para                          ! parabolic interpolation coefficients
real :: radius, theta, stheta, sphi, ctheta, cphi

end module sweeps

