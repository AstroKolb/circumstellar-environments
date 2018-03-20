subroutine images

! gather data from processors to make a 2D slice, append data to a netCDF file
! use ncview to animate by scanning through the time dimension
!-----------------------------------------------------------------------------

! GLOBALS

use NETCDF
use global
use zone

include 'mpif.h'

! LOCALS

CHARACTER(LEN=50) :: filename
CHARACTER(LEN=10) :: coord, varname
INTEGER :: ncstat, ncid, frame_number, rank, nvar, natt, icheck, kcheck, count, jmid, yystat(MPI_STATUS_SIZE)
INTEGER :: i, j, imn, ipn, kmn, kpn, jmn, jpn, kc, jc, image1_pe, image2_pe, count1, count2, mpitag
INTEGER :: xDimID, yDimID, tDimID, varID, XScale_varID, YScale_varID, TScale_varID
INTEGER, DIMENSION(3) :: start
REAL :: a1, a2, c1, c2, c3, c4, ytop, ybot, dely, radius, theta, phi, entropy

REAL, DIMENSION(imax,kmax/pez) :: buffone
REAL, DIMENSION(imax,jmax/pey) :: bufftwo
REAL, DIMENSION(imax,kmax) :: dens1
REAL, DIMENSION(imax,jmax) :: dens2

INTEGER, PARAMETER :: iplt = 400
INTEGER, PARAMETER :: jplt = 400  
REAL, DIMENSION(iplt,jplt) :: imagea
REAL, DIMENSION(iplt) :: x
REAL, DIMENSION(jplt) :: y

!======================================================================
varname = 'Density'
image2_pe = (npe/4) + (npey/2)  ! this pe on the yin grid will collect data and write out netcdf file
image1_pe = (npe/2) + image2_pe ! this pe on the yang grid collects data and sends to image2_pe
mpitag = 87
count  = imax * kmax

if (.not.yin) then

 if (jcol == npey/2) then   ! pull together an r-phi (at theta=pi/2) slice on the yang grid
  do k = 1, ks
   do i = 1, imax
     buffone(i,k) = zro(i,1,k)
   enddo
  enddo
  count1 = imax*ks
  call MPI_GATHER(buffone,count1,VH1_DATATYPE,dens1,count1,VH1_DATATYPE,npez/2,MPI_COMM_COL, mpierr)

  ! root sends this data over to image2_pe collecting yin grid data
  if (krow==npez/2) then 
   call MPI_SEND(dens1, count, VH1_DATATYPE, image2_pe, mpitag, MPI_COMM_WORLD, mpierr)
  endif
 endif

else

 if (krow == npez/2) then   ! pull together an r-theta (at phi=0) slice on the yin grid
  do j = 1, js
   do i = 1, imax
     bufftwo(i,j) = zro(i,j,ks/2)
   enddo
  enddo
  count2 = imax*js
  call MPI_GATHER(bufftwo,count2,VH1_DATATYPE,dens2,count2,VH1_DATATYPE,npey/2,MPI_COMM_ROW, mpierr)

  if (jcol == npey/2) then
   call MPI_RECV(dens1, count, VH1_DATATYPE, image1_pe, mpitag, MPI_COMM_WORLD, yystat, mpierr)
  endif
 endif

endif

!======================================================================

if (mype == image2_pe) then     ! only this unique processor has both data slices

  ! Create a SQUARE cartesian grid for interpolated movie frames
  ytop  = zxa(imax) 
  ybot  = -ytop
  dely = (ytop - ybot) / jplt
  do j = 1, jplt
    y(j) = ybot + real(j-1) * dely
  enddo
  do i = 1, iplt
    x(i) = ybot + real(i-1) * dely
  enddo

  do j = 1, jplt
   do i = 1, iplt

     radius = sqrt(y(j)**2 + x(i)**2) + small
     theta  = acos(y(j)/radius)
     phi    = acos(x(i)/radius)
     if (y(j)<0.0) phi = -phi

     if (radius < zxc(1)) then
       entropy  = 0.0
     else if (radius > zxc(imax)) then
       entropy  = 0.0
     else
       imn = 1
       ipn = imax
       do while(ipn-imn > 1)
         ic = (imn + ipn)/2
         if(radius > zxc(ic)) then
           imn = ic
         else
           ipn = ic
         endif
       enddo
       c1 = (radius - zxc(ipn)) / (zxc(imn) - zxc(ipn))
       c2 = (radius - zxc(imn)) / (zxc(ipn) - zxc(imn))

       if (abs(phi) < zzc(kmax)) then  ! interpolate from yang grid

         kmn = 1
         kpn = kmax 
         do while(kpn-kmn > 1)
           kc = (kmn + kpn)/2
           if(phi > zzc(kc)) then
              kmn = kc
           else
              kpn = kc
           endif
         enddo
         c3 = (phi    - zzc(kpn)) / (zzc(kmn) - zzc(kpn))
         c4 = (phi    - zzc(kmn)) / (zzc(kpn) - zzc(kmn))
         a1 = dens1(imn,kmn)*c1 + dens1(ipn,kmn)*c2
         a2 = dens1(imn,kpn)*c1 + dens1(ipn,kpn)*c2
         entropy = a1*c3+a2*c4

       else   ! interpolate from yin grid

         jmn = 1
         jpn = jmax 
         do while(jpn-jmn > 1)
           jc = (jmn + jpn)/2
           if(theta > zyc(jc)) then
              jmn = jc
           else
              jpn = jc
           endif
         enddo
         c3 = (theta  - zyc(jpn)) / (zyc(jmn) - zyc(jpn))
         c4 = (theta  - zyc(jmn)) / (zyc(jpn) - zyc(jmn))
         a1 = dens2(imn,jmn)*c1 + dens2(ipn,jmn)*c2
         a2 = dens2(imn,jpn)*c1 + dens2(ipn,jpn)*c2
         entropy = a1*c3+a2*c4

       endif

     endif

     imagea(i,j) = entropy
 
   enddo
  enddo


  ! create the netCDF filename and try to open it
  filename = trim(prefix) // 'XZ.nc'
  ncstat = nf90_open(filename, NF90_WRITE, ncid)

  if (ncstat == 2) then   ! file does not exist, so create it and start definitions

    ! create the file and define coordinates/variables
    ncstat = nf90_create(filename, NF90_Clobber, ncid)
    if (ncstat /= nf90_NoErr ) print *, ncstat, NF90_STRERROR(ncstat)

    ! Initialize Dimensions
    ncstat = nf90_def_dim(ncid, "x", iplt, xDimID)      
    ncstat = nf90_def_dim(ncid, "y", jplt, yDimID)
    ncstat = nf90_def_dim(ncid, "time" , NF90_UNLIMITED, tDimID)
        if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

    ! define the coordinate scales as 1D variables
    ncstat = nf90_def_var(ncid, "x", nf90_float,(/ xDimID /), XScale_varID)
    ncstat = nf90_def_var(ncid, "y", nf90_float,(/ yDimID /), YScale_varID)
    ncstat = nf90_def_var(ncid, "time", nf90_float,(/ tDimID /), TScale_varID)
        if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

    ! define the simulation variables
    ncstat = nf90_def_var(ncid,varname,nf90_float,(/ xDimID, yDimID, tDimID /), varID)
        if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

    ! take dataset out of definition mode
    ncstat = NF90_ENDDEF(ncid)

    ! write out coordinate arrays
    ncstat = nf90_put_var(ncid, XScale_varID, x)
    ncstat = nf90_put_var(ncid, YScale_varID, y)
        if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

    frame_number = 1

  else

    ! inquire about dimensions to make sure it fits current run
    ncstat = nf90_inquire(ncid, rank, nvar, natt, tDimID)
    ncstat = nf90_inquire_dimension(ncid, 1, coord, icheck)
      if(icheck.ne.iplt) print *, 'Dim 1 does not equal imax', icheck, imax
    ncstat = nf90_inquire_dimension(ncid, 2, coord, kcheck)
      if(kcheck.ne.jplt) print *, 'Dim 2 does not equal kmax', kcheck, kmax
    ncstat = nf90_inquire_dimension(ncid, tDimID, coord, frame_number)
    frame_number = frame_number + 1

  endif

  start = 1
  start(3) = frame_number
  ncstat = nf90_put_var(ncid, 4, imagea, start)
     if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

  start(1) = frame_number
  ncstat = nf90_put_var(ncid, 3, time, start)
     if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

  ncstat = nf90_close(ncid)

endif

return
end


