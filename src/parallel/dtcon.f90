subroutine dtcon

! set the timestep using various constraints.
! NB - additional constraints may be required depending on the problem!
!-----------------------------------------------------------------------

! GLOBALS
use global
use zone

include 'mpif.h'

! LOCALS
integer :: i, j, k, mpierr
real ::  ridt, rodt, dtx, dt3, xvel, yvel, zvel, widthy, widthz

!------------------------------------------------------------------------
!   Hydro constraint on timestep.  Use R*d(theta) if y geometry is angular
ridt = 0.

do k = 1, ks
 do j = 1, js
  do i = 1, imax
    widthy = zdy(j+jcol*js)
    widthz = zdz(k+krow*ks)
    widthy = widthy*zxc(i)
    widthz = widthz*zxc(i)*sin(zyc(j+jcol*js))
    xvel = abs(zux(i,j,k)) / zdx(i)
    yvel = abs(zuy(i,j,k)) / widthy
    zvel = abs(zuz(i,j,k)) / widthz
!    if (xvel>1e10 .or. yvel>1e10 .or. zvel>1e10) then
!      write(*,*) xvel, yvel, zvel, i, j, k
!      stop
!    endif
    ridt = max(xvel,yvel,zvel,ridt)
  enddo
 enddo
enddo
!write(*,*) svel, ridt

ridt = max(svel,ridt)
call MPI_ALLREDUCE( ridt, rodt, 1, VH1_DATATYPE, MPI_MAX, MPI_COMM_WORLD, mpierr )

dtx  = courant / rodt     ! global time constraint for given courant parameter
dt3  = 1.1 * dt           ! limiting constraint on rate of increase of dt              
dt   = min( dt3, dtx )    ! use smallest required timestep                                        
      
if (time/dt .gt. 1.e6) then   ! if timestep becomes too small, stop the program!
  if (mype.eq.0) then
    write(*,*) 'Timestep has become too small: dt = ',dt
    write(*,*) '                             time = ',time
  endif
  call prin('ABORT')
  call MPI_FINALIZE(mpierr)
  stop
endif

return
end
