This is the readme for the model associated with the paper:

Traub RD, Hawkins K, Adams NE, Hall SP, Simon A, Whittington MA (2020)
Layer 4 pyramidal neuron dendritic bursting underlies a post-stimulus
visual cortical alpha rhythm Nature Communications Biology

This code was contributed by R Traub.

The source code is standard fortran with special instructions for the
mpi parallel environment.  Included are 2 pdf files with illustrative
output, an example of the main program (alphayY67.f), a number of
fortran subroutines and one in C, and the “makefile” which does
compilation and linking.

I have run the code on an IBM power chip residing in IBM’s Cognitive
Computing Cluster, at a time when the AIX operating system was in use
(that is IBM’s version of Linux).  The machine model is given in the
ms text, but is probably not relevant.  The CCC has since switched to
Linux.  The code should also run under Linux, but I have not tried it;
however, I have converted similar programs from AIX to Linux without
difficulty; the routines simply need to be compiled with the Fortran
compiler in use on the operating system, and an appropriate makefile
used.  The details of the makefile are going to be system-dependent.

The code does not run on data.  Each run is autonomous, with
parameters chosen by the programmer.  I also do not expect that the
code will run on a desktop, given the need for mpi.
