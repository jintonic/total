total - Totani's intepolator of the Livermore model

These are [Tomononi Totani](http://tac.astron.s.u-tokyo.ac.jp/~totani/)'s 
Fortran codes to intepolate the supernova neutrino database described in 
his paper [Future Detection of Supernova Neutrino Burst and Explosion Mechanism](http://stacks.iop.org/0004-637X/496/i=1/a=216).

The source codes are slightly modified to allow specifying data directory by an
enviroment variable called ```TOTAL_DATA_DIR```. The program cannot open the
data files wilson-*.dat if the enviroment variable is not set. The source codes
hence have to be compiled with a morden Fortran compiler.

A makefile is provided to create a shared library ```libTOTAL.so``` with
[gfortran](http://gcc.gnu.org/wiki/GFortran).
The library can be installed to ```/prefix/lib/``` by ```make && make install```.
The prefix can be specified in the first line of the Makefile.

#####Original README
```
   <<<<<<<<<<<<<<<  Supernova Neutrino Data Set  >>>>>>>>>>>>>>>>>>>

                                      Tomo Totani
                                      2003, Feb. 14 last update.

This package provides numerical data of supernova neutrino flux, spectrum, and
their time evolution presented in Totani et al. (1998), ApJ, 496, 216.  Please
reference this paper when you write a paper using this data.  Questions should
be directed to Tomo Totani, ttotani@princeton.edu

The package includes:

README               this file
wilson-early.dat     neutrino data for t < 0.2 sec
wilson-late.dat      neutrino data for t > 0.2 sec
wilson-nb.dat        neutrino data during neutronization burst

wilson_NL.f          Fortran subroutine to read these files and interpolate
wilson_nb_NL.f       Fortran subroutine, called in wilson_NL

test.f               A test program showing how to use wilson_NL.f

spline.f             These are subroutines used in wilson_NL.f
splint.f
linear_int.f
nat_cub_sp.f

```
