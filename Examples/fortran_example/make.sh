#!/bin/bash

# GNU
gfortran -o testf -O3 main.f90 -ffree-form -cpp -I../../include -lgfortran `INSTALLED_DIR/bin/npt-config --libs`

# Intel
#ifort -o testc -O3 main.f90 -fpp -I../../include `INSTALLED_DIR/bin/npt-config --libs`

# Sparc
#frtpx -o testf -O3 main.f90 -Cpp -I../../include `INSTALLED_DIR/bin/npt-config --libs` --host=sparc64-unknown-linux-gnu
