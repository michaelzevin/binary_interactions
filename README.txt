# Binary BBH Scattering Experiments -- KSPA 2017
Fortran base code, python wrappers, data, documentation, plotting scripts, and paper for binary BBH scattering project #16


# INFO
The main code is split into two files:
1) 'TEST_module_Nbody_AffineTides_solver.f90'
This file contains the main code in one big 'module'.
2) 'TEST_main_Nbody_AffineTides_solver.f90'
This file contains two different interphases that can be used to communicate with the main code. The two interphases are respectively a subroutine and a program.



# INPUT FILE
b1_const_arr = np.array([b1_mass, b1_rad, b1_gas_n, b1_gas_gamma, b1_Mqp, b1_evoTides_yesno, b1_RigidSph_yesno, 0,0,0], dtype='d')
[use_12PN, use_25PN, Identify_3Body_endstate, max_sim_nrsteps, FREE, outputinfo_screenfiles, ...]
   example: nbody_params_arr_1_INT  = np.array([0, use_25PN, 1, 1000000000, 0, 1, 0,0,0,0], dtype='i')
[scale_dt, max_sim_time (set later in code!!!), evolvetides_threshold, ENDbinsingle_threshold, max_simtime_sec (set later in code!!!), IMSbinsingle_threshold, tidaldisrup_threshold, insp_threshold, ...]
   example: nbody_params_arr_2_REAL = np.array([0.01, -1.0, evolvetides_threshold, 0.1, -1.0, 0.01, tidaldisrup_threshold, insp_threshold,0,0], dtype='d')



# RUN
./main_Nbody_AffineTides_solver.exe <input/test_0.txt

# make gif using imagemagick
sub=1; convert -delay 0 $(for i in $(seq 0 `expr $(ls | wc -l) - $sub`); do echo ${i}.png; done) -loop 0 animated.gif

# make movie using ffmpeg
ffmpeg -framerate 24 -start_number 0 -i %04d.png -pix_fmt yuv420p movie.mp4



# COMPILE
gfortran module_Nbody_AffineTides_solver.f90 main_Nbody_AffineTides_solver.f90 -o main_Nbody_AffineTides_solver.exe libs_odepackdp.so gslwr.o -lgsl -lgslcblas -L/opt/local/lib -llapack



# MODULE INFORMATION
SUBROUTINE Calc_binary_info(pos_1, vel_1, mass_1, pos_2, vel_2, mass_2, out_binary_info_arr):
    Takes in positions, velocities, and masses of binary and returns:
        out_binary_info_arr(:) = [E_kin, E_pot, E_tot, a_bin, e_bin]



# Remove redspace in vim for fortran code:
:hi link fortranTab NONE

# Remove first N characters of file
sed -i -r '1s/^.{N}//' file.dat

# Change cluster files to proper format:
sed -i.bak s/#//g *
sed -i.bak '1s/^.//' *



---------------------------------------------------------------------------
Libraries:
---------------------------------------------------------------------------
The code needs the libraries: ODEPACK, LAPACK, GSL
Download and unpack each of these programs. 
The programs must be compiled as shared objects.
Normally this can be done by just adding '-fPIC' to the compile line.
This will be shown below for each of the programs.
---------------------------------------------------------------------------
---------------------------------------------------------------------------
ODEPACK - We use LSODE from ODEPACK to solve the diff equations.
---------------------------------------------------------------------------
Make shared library (libs_odepackdp_JS.so):
gfortran -fPIC -c opkda1.f opkda2.f opkdmain.f
gfortran -fPIC -shared -o libs_odepackdp_JS.so opkda1.o opkda2.o opkdmain.o
You now have libs_odepackdp_JS.so (or whatever name you choose).
I don't know if its necessary to copy this lib .so to the folder where the main code is.
May need to add path to .so file: LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/projects/b1011/mzevin/binary_interactions/main_files/
---------------------------------------------------------------------------
---------------------------------------------------------------------------
LAPACK - used for lin alg operations. We use it for calc matrix eigen vecs, vals, inverse, etc.
---------------------------------------------------------------------------
It is important to compile LAPACK as a shared library. For this we need to add -fPIC
a few places in the make.inc (also see make.inc.example for a template example as described in README).
Do the following:

In the top of make.inc you see:
FORTRAN  = gfortran 
OPTS     = -O2 -frecursive
DRVOPTS  = $(OPTS)
NOOPT    = -O0 -frecursive
LOADER   = gfortran
LOADOPTS =

Change to:
FORTRAN  = gfortran -fPIC
OPTS     = -O2 -fPIC
DRVOPTS  = $(OPTS)
NOOPT    = -O0 -fPIC
LOADER   = gfortran -fPIC
LOADOPTS =

Im not sure if we need the -frecursive. But it works fine without.
Now compile by typing:
make blaslib
make
More info see e.g. https://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=3589 
To see a list of the LAPACK routines: http://www.netlib.org/lapack/double/
---------------------------------------------------------------------------
---------------------------------------------------------------------------
GSL - we use gsl to calc e.g. the self grav elliptical integral.
---------------------------------------------------------------------------
gsl is written in c, so we need a wrapper.
For the function we need, the wrapper must look like this:

#include <gsl/gsl_sf_ellint.h>
void c_wrapper_gsl_sf_ellint_rj_(double* R, double* x, double* y, double* z, double* p, double* err){
        *R = gsl_sf_ellint_RJ(*x, *y, *z, *p, *err);
}

save it as: gslwr.c
Make shared object:
type: gcc -c -fPIC gslwr.c
link to gsl when compiling: gcc -Wall -I/usr/local/include -c example.c if its not working..
---------------------------------------------------------------------------

