# Specify search path for subroutines that perform 1D ppmlr hydrodynamics
VPATH = ../PPMLR
#
# Tunable parameters (currently set for Cray YMP)
#
# LDR           Program to load the objects into an executable
# LDFLAGS       For definition of lmset and csu on the Cray 2
#               Flags to the loader
# LIBS          A list of libraries to use on the Cray 2
#               Note that the C library (-lc) must ALWAYS be the last library
#               in the load step on both the XMP and the Cray 2
# COMMAND       Name of the executable
# -g  :: needed to turn on debugging
# -Wall :: look for unused variables, uninitialized functions
# -Wextra :: look for uninitialized variables
# -Warray-temporaries :: look for temp arrays (don't use this)
# -Wconversion :: checks for conversions (not necessary)
# -fimplicit-none :: checks for undefined variables (ALWAYS INCLUDE!)
# -fbacktrace :: ???
# -ffree-line-length-0 :: ???
# -ffpe-trap=zero,overflow,underflow :: ???
# -finit-real=nan :: initialize reals as nan (SLOW!)

# 0.67 mins vs 0.33

F90 = mpif90 
#FFLAGS = -c -O3 -fdefault-real-8 -I/common/software/netcdf/include
FFLAGS = -c -O3 -fdefault-real-8 -g -finit-real=nan -fimplicit-none -I/common/software/netcdf/include
LDR=  mpif90
LDRFLAGS= 
LIBS=   -L/common/software/netcdf/lib -lnetcdf
#
# List of objects to build multidimensional VH-1 hydro code:

VHOBJS = vh1mods.o zonemods.o dtcon.o dump.o evolve.o flatten.o forces.o \
         init.o parabola.o ppmlr.o volume.o prin.o remap.o riemann.o \
         states.o sweepx1.o sweepx2.o sweepy.o sweepz.o diagnostics.o \
         corner.o boundaryY.o  boundaryZ.o  boundaryYsetup.o  boundaryZsetup.o \
			vhone.o


# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build the executable
#
.SUFFIXES: .f90 .o

.f90.o:
	$(F90) $(FFLAGS) $<

vhone:	$(VHOBJS)
	$(LDR) $(LDRFLAGS) -o vhone $(VHOBJS) $(LIBS); mv vhone ../..

merge:	merge_netcdf.o
	$(LDR) $(LDRFLAGS) -o merge merge_netcdf.o $(LIBS)

#
#
clean:
	rm -f *.o *.l

clobber:	clean
	rm -f ../../vhone postprocess
	rm -f *.mod


# Dependencies for the object files
dtcon.o:    dtcon.f90    global.mod zone.mod
dump.o:     dump.f90     global.mod zone.mod
init.o:     init.f90     global.mod zone.mod
vhone.o:    vhone.f90    global.mod zone.mod sweepsize.mod 
sweepx1.o:  sweepx1.f90  global.mod zone.mod sweeps.mod sweepsize.mod
sweepx2.o:  sweepx1.f90  global.mod zone.mod sweeps.mod sweepsize.mod
sweepy.o:   sweepy.f90   global.mod zone.mod sweeps.mod sweepsize.mod
sweepz.o:   sweepz.f90   global.mod zone.mod sweeps.mod sweepsize.mod
prin.o:     prin.f90     global.mod zone.mod sweeps.mod sweepsize.mod

diagnostics.o: diagnostics.f90 global.mod zone.mod



evolve.o:   evolve.f90   global.mod sweeps.mod sweepsize.mod
flatten.o:  flatten.f90  global.mod sweeps.mod sweepsize.mod
forces.o:   forces.f90   global.mod sweeps.mod sweepsize.mod
ppmlr.o:    ppmlr.f90    global.mod sweeps.mod sweepsize.mod
remap.o:    remap.f90    global.mod sweeps.mod sweepsize.mod
riemann.o:  riemann.f90  sweepsize.mod
states.o:   states.f90   global.mod sweeps.mod sweepsize.mod
parabola.o: parabola.f90 sweepsize.mod





