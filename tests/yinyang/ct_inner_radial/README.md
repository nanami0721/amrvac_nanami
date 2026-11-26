# Inner radial CT boundary on a Parker wind

Use the same 3D user module with spherical.par and yinyang.par. Build and run
in separate external case directories with AMRVAC_DIR pointing to an external
source copy, using the project memory guard, make -j1 and two MPI ranks.
Repeat each parameter file with drive=0.001 in the radial_test namelist.

The shell starts with radial B=0.02/r² and a subsonic Parker wind. The inner
fluid ghosts follow the Parker profile and the magnetic ghosts use the CT
closure. The outer boundary is continuous. With drive=0 each inner radial
face must retain its initial flux. The manufactured second case imposes the
edge integrals of E=drive*(z_hat cross r_hat) at r=1. Its exact face flux is
the initial flux minus 2*drive*t times the integral of z_hat dot r_hat over
the face solid angle; the test evaluates this integral in each patch basis.

RADIAL_CT reports step, time, face and ghost sample counts, maximum absolute
face-flux error and first-layer divergence, maximum relative ghost density and
internal-energy errors, and minimum first-layer thermal pressure. Each patch
checks every physical inner face, including seam endpoints. The ghost check
includes radial/angular corners on Yin–Yang. The driven case tests Faraday's
law with a prescribed EMF, independent of the velocity boundary condition.
