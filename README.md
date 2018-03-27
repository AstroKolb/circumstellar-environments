# circumstellar-environments


to do list

1.  verify working code
2.  run with compiler flags
3.  test mass conservation for steady flow, isothermal+adiabatic
4.  edit grid (zoom, overlap)
5.  add primary gravity, test mass conservation for parker wind
6.  add wind profile, test mass conservation (do this in 1D too...)
7.  add pseudoforces, test mass conservation
8.  test min radius
9.  add accretion, test..
10. add cooling, test..
11. fix grid overlap (see step 4 notes)


# 5. Add primary gravity, test with Parker wind

Appears to be working. Small glitch when initialized but seems to smooth itself out

![alt text](http://astro.physics.ncsu.edu/~cekolb/research/circumstellar-environments/img/step5_mass.png)

mass, normalized to last zone

![alt text](http://astro.physics.ncsu.edu/~cekolb/research/circumstellar-environments/img/step5_flux.png)

mass flux, inner boundary doesn't quite agree, but I think it averages.. shown is zone 7, zone 1 is negative, zone 3 is more positive..

tested with debug compiler


# 4. Edit Grid Setup

Fixed a bug in the grid overlap scheme

![alt text](http://astro.physics.ncsu.edu/~cekolb/research/circumstellar-environments/img/overlap.png)

everything appears to be working
 - mass is conserved
 - tested with gamma = 1, 5/3
 - tested with debug compiler

![alt text](http://astro.physics.ncsu.edu/~cekolb/research/circumstellar-environments/img/mass_conservation2.png)
(mass conservation for new grid)

 NOTE!!! setting novery/z = 6 looks like 4 overlap on the z-plane - fix this


# 3. test mass conservation for steady flow

added diagnostics.f90 subroutines
 - outputs to 'output/*prefix*_diagnostics.dat'

using u=10, rho=1/r^2, iso:
![alt text](http://astro.physics.ncsu.edu/~cekolb/research/circumstellar-environments/img/mass_conservation1.png)

 - identical plot for gamma=5/3
 - checked with debug flags


# 2. Run with compiler flags

Many issues found, should be working now
still need to run with -Wall at some point

verified the results are unchanged from initial problem
 - used compiler flags during test


# 1. Verifying working code 
(downloaded from DrB email, some reorganization)

Code is working, appears to generate a SSDW
![alt text](http://astro.physics.ncsu.edu/~cekolb/research/circumstellar-environments/img/initial.png)