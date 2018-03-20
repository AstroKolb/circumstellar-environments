subroutine diagnostics_init

   use global
   use zone

   if (mype == 0) open(unit=16, file='output/' // trim(prefix) // '_diagnostics.dat', form='formatted')

return
end

subroutine diagnostics_finalize

   use global
   use zone

   if (mype == 0) close(16)

return
end


subroutine diagnostics

   use global
   use zone

   include 'mpif.h'

   integer :: i, j, k, mpierr
   real :: dVol, mlocal, mtotal

   ! loop through data and output diagnostics
   mlocal = 0.0
   mtotal = 0.0
   do k = 1, ks
      do j = 1, js
         do i = 1, imax
            dVol = zxc(i)**2*sin(zyc(mypey*js+j))*zdx(i)*zdy(mypey*js+j)*zdz(mypez*ks+k)
            mlocal = mlocal + zro(i,j,k)*dVol
         enddo
      enddo
   enddo

   ! combine across processors
   call MPI_ALLREDUCE(mlocal, mtotal, 1, VH1_DATATYPE, MPI_SUM, MPI_COMM_WORLD, mpierr)

   if (mype == 0) write(16,*) time, mtotal


return
end