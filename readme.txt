This code is based on the coed "NEXD-2D". Detailed content can be referred to "https://github.com/seismology-RUB/NEXD-2D"

## Requirements

* GNU Make
* Fortran compiler (sufficient standard)
* LAPACK libraries
* MPI libraries (e.g., OpenMPI) and METIS (http://glaros.dtc.umn.edu/gkhome/metis/metis/overview; tested with version 4.0.3) for parallel applications


## Installation

1. Install all software dependencies.

2. Modify the "Makefile" ( the path to your METIS installation).
		la = -L/es01/software/apps/lapack-3.8.0/lib/ -llapack -lblas
		lib = /es01/paratera/sce2567/metis-4.0.3/libmetis.a 

3. Run command
     ```
     make all
     ```