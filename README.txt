2D radial distribution function calculator g(r)
===============================================

Fortran, with a Python wrapper.

To generate the dynamic library, install f2py and run:

$ f2py -c gr.f90 --fcompiler=intelem -m gr_lib

This generates gr_lib.so. If intel compiler is not available, then run:

$ f2py -c gr.f90 -m gr_lib

Now, to generate a 2D g(r) of polymers in the middle of the channel and
analyze bundle structure formation do:

$ python get_2D_gr.py

In order for this to work:

1) Make shure to run in the parent dir of 1_run 2_run ... n_run
2) Make shure "mfa_input" and "system_input" files are present in the working dir
