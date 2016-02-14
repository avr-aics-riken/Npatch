#!/bin/bash

# GNU
gcc -o testc -O3 main.c -I../../include `INSTALLED_DIR/bin/npt-config --libs`

# Intel
#icc -o testc -O3 main.c -I../../include `INSTALLED_DIR/bin/npt-config --libs`

# Sparc
#FCCpx -o testc -O3 main.c -I../../include `INSTALLED_DIR/bin/npt-config --libs` --host=sparc64-unknown-linux-gnu
