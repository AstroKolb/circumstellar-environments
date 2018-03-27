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

   integer :: i, j, k, myj, myk, mpierr
   real :: dVol, dArea
   real :: lmass1, lmass2, lmass3, lmass4
   real :: gmass1, gmass2, gmass3, gmass4
   real :: lflux1, lflux2, lflux3, lflux4
   real :: gflux1, gflux2, gflux3, gflux4

   ! loop through data and output diagnostics
   lmass1 = 0.0
   lmass2 = 0.0
   lmass3 = 0.0
   lmass4 = 0.0
   lflux1 = 0.0
   lflux2 = 0.0
   lflux3 = 0.0
   lflux4 = 0.0

   gmass1 = 0.0
   gmass2 = 0.0
   gmass3 = 0.0
   gmass4 = 0.0
   gflux1 = 0.0
   gflux2 = 0.0
   gflux3 = 0.0
   gflux4 = 0.0

   do k = 1, ks
      myk = mypez*ks+k
      do j = 1, js
         myj = mypey*js+j

         if (ongrid(myj,myk) == 0) then

            ! zone 1 mass flux
            i = 7
            dArea = zxc(i)**2*sin(zyc(myj))*zdy(myj)*zdz(myk)
            lflux1 = lflux1 + zro(i,j,k)*zux(i,j,k)*dArea

            ! zone 1 mass
            do i = 7, imax/4
               dVol = zxc(i)**2*sin(zyc(myj))*zdx(i)*zdy(myj)*zdz(myk)
               lmass1 = lmass1 + zro(i,j,k)*dVol
            enddo

            ! zone 2 mass flux
            i = imax/4
            dArea = zxc(i)**2*sin(zyc(myj))*zdy(myj)*zdz(myk)
            lflux2 = lflux2 + zro(i,j,k)*zux(i,j,k)*dArea
            
            ! zone 2 mass
            do i = imax/4, imax/2
               dVol = zxc(i)**2*sin(zyc(myj))*zdx(i)*zdy(myj)*zdz(myk)
               lmass2 = lmass2 + zro(i,j,k)*dVol
            enddo

            ! zone 3 mass flux
            i = imax/2
            dArea = zxc(i)**2*sin(zyc(myj))*zdy(myj)*zdz(myk)
            lflux3 = lflux3 + zro(i,j,k)*zux(i,j,k)*dArea
            
            ! zone 3 mass
            do i = imax/2, 3*imax/4
               dVol = zxc(i)**2*sin(zyc(myj))*zdx(i)*zdy(myj)*zdz(myk)
               lmass3 = lmass3 + zro(i,j,k)*dVol
            enddo

            ! zone 4 mass flux
            i = 3*imax/4
            dArea = zxc(i)**2*sin(zyc(myj))*zdy(myj)*zdz(myk)
            lflux4 = lflux4 + zro(i,j,k)*zux(i,j,k)*dArea
            
            ! zone 4 mass
            do i = 3*imax/4, imax
               dVol = zxc(i)**2*sin(zyc(myj))*zdx(i)*zdy(myj)*zdz(myk)
               lmass4 = lmass4 + zro(i,j,k)*dVol
            enddo

         endif
      enddo
   enddo

   ! combine across processors
   call MPI_ALLREDUCE(lmass1, gmass1, 1, VH1_DATATYPE, MPI_SUM, MPI_COMM_WORLD, mpierr)
   call MPI_ALLREDUCE(lmass2, gmass2, 1, VH1_DATATYPE, MPI_SUM, MPI_COMM_WORLD, mpierr)
   call MPI_ALLREDUCE(lmass3, gmass3, 1, VH1_DATATYPE, MPI_SUM, MPI_COMM_WORLD, mpierr)
   call MPI_ALLREDUCE(lmass4, gmass4, 1, VH1_DATATYPE, MPI_SUM, MPI_COMM_WORLD, mpierr)

   call MPI_ALLREDUCE(lflux1, gflux1, 1, VH1_DATATYPE, MPI_SUM, MPI_COMM_WORLD, mpierr)
   call MPI_ALLREDUCE(lflux2, gflux2, 1, VH1_DATATYPE, MPI_SUM, MPI_COMM_WORLD, mpierr)
   call MPI_ALLREDUCE(lflux3, gflux3, 1, VH1_DATATYPE, MPI_SUM, MPI_COMM_WORLD, mpierr)
   call MPI_ALLREDUCE(lflux4, gflux4, 1, VH1_DATATYPE, MPI_SUM, MPI_COMM_WORLD, mpierr)

   if (mype == 0) write(16,*) time, gmass1, gmass2, gmass3, gmass4, gflux1, gflux2, gflux3, gflux4


return
end