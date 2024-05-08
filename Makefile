#  General definitions
#
#  Paths
#
bindir = ./bin
obsdir = ./obj
moduledir = ./mod
srcdir = ./src
AR = ar

# Choose your compiler:
F95 = mpif90
#F95 = mpiifort

ifeq ($(notdir $(F95)),g95)
	FFLAGS = -O3 -Wunused -fmod=$(moduledir)
else ifeq ($(F95),mpif90)
	# ==================================================================
	# GFORTRAN COMPILER
	# ==================================================================
	# Choose a suitable set of compiler flags for gfortran. The default 
	# flags are explained in the following...
	# required:
	#   -J$(moduledir)
	#   -fimplicit-none
	#   -ffixed-line-length-none
	#   -ffree-line-length-none
	# optional:
	#   -O3                : The compiler tries to optimize the code to 
	#                        make it faster. Level 3 is the most 
	#                        optimization available.
	#   -Og                : Enables optimizations that do not interfere
	#                        with debugging. It should be the 
	#                        optimization level of choice for the 
	#                        standard edit-compile-debug cycle.
	#   -funroll-all-loops : This might increase the speed of the code 
	#                        as well.
	#   -fbounds-check     : Add a check that the array index is within 
	#                        the bounds of the array every time an array
	#                        element is accessed.
	#                        This substantially slows down a program 
	#                        using it, but is a very useful way to find 
	#                        bugs related to arrays.
	#   -fbacktrace        : If the program crashes, a backtrace will be
	#                        produced if possible.
	#   -Wall              : gfortran will generate warnings about many 
	#                        common sources of bugs.
	#   -Wpedantic         : Generate warnings about language features 
	#                        that are supported by gfortran but are not 
	#                        part of the official Fortran 95 standard. 
	#   -Wextra            : This enables some extra warning flags that 
	#                        are not enabled by -Wall. Note, 
	#                        -Wunused-parameter (part of -Wextra) has a 
	#                        problem to correctly work with constants.h.
	# recommended flags for developers:
#	FFLAGS = -Og -J$(moduledir) -fimplicit-none -ffixed-line-length-none -ffree-line-length-none -fbounds-check -fbacktrace -Wall -Wpedantic -Wextra
	# recommended flags for profiling (CP CP):
#	FFLAGS = -O3 -J$(moduledir) -fimplicit-none -ffixed-line-length-none -ffree-line-length-none -funroll-all-loops -pg -g -no-pie
	# recommended flags for users:
	FFLAGS = -O3 -J$(moduledir) -fimplicit-none -ffixed-line-length-none -ffree-line-length-none -funroll-all-loops
else ifeq ($(F95),mpiifort)
	# ==================================================================
	# INTEL COMPILER (ifort)
	# ==================================================================
	# Choose a suitable set of compiler flags for ifort. The default 
	# flags are explained in the following...
	# required:
	#   -module $(moduledir)
	#   -implicitnone
	#   -132
	# optional:
	#   -O3                : The compiler tries to optimize the code to 
	#                        make it faster. Level 3 is the most 
	#                        optimization available. Level 0 (-O0) is 
	#                        recommended for developers.
	#   -funroll-loops     : This might increase the speed of the code 
	#                        as well.
	#   -nowarn            : Suppresses all warnings.
	# developers:
#	FFLAGS = -O0 -module $(moduledir) -implicitnone -132
	# users:
	FFLAGS = -O3 -module $(moduledir) -implicitnone -132 -funroll-loops -nowarn
endif

obj_mesher = program_mesher.o \
	constantsMod.o \
	parameterMod.o \
	gllMod.o \
	warpfactorMod.o \
	nodesMod.o \
	triTrafoMod.o \
	vtkMod.o \
	jacobiMod.o \
	simplexMod.o \
	vandermondeMod.o \
	dmatricesMod.o \
	geometricFactorsMod.o \
	rosettaGammaMod.o \
	meshMod.o \
	plotMod.o \
	liftMod.o \
	derMod.o \
	normalsMod.o \
	matrixMod.o \
	dtMod.o \
	stfMod.o \
	sourceReceiverMod.o \
	triangleNeighborMod.o\
	externalModelMod.o \
	linearSystemMod.o \
	fileunitMod.o \
	fileParameterMod.o \
	errorMessage.o \
	realloc.o \
	materialsMod.o \
	progressbarMod.o \
	mpiMod.o \
	waveMod.o

obj_solver = program_solver.o \
	mpiincludeMod.o \
	constantsMod.o \
	parameterMod.o \
	gllMod.o \
	warpfactorMod.o \
	nodesMod.o \
	triTrafoMod.o \
	vtkMod.o \
	jacobiMod.o \
	simplexMod.o \
	vandermondeMod.o \
	dmatricesMod.o \
	geometricFactorsMod.o \
	rosettaGammaMod.o \
	meshMod.o \
	plotMod.o \
	liftMod.o \
	derMod.o \
	normalsMod.o \
	timeloopMod.o \
	matrixMod.o \
	waveMod.o \
	dtMod.o \
	stfMod.o \
	sourceReceiverMod.o \
	triangleNeighborMod.o \
	externalModelMod.o \
	mpiMod.o \
	linearSystemMod.o \
	fileunitMod.o \
	fileParameterMod.o\
	pmlMod.o \
	errorMessage.o \
	realloc.o \
	materialsMod.o \
	progressbarMod.o \
	calcFluxMod.o \
	timestampMod.o \
	calendar.o \
	convert_time.o \
	adjointMod.o \
	timeSeries.o\
	dateTime.o\
	recursiveFilterCoefficients.o\
	fourierTransform.o\
	realloc.o\
	filterCoefficientsDecimate.o\
	timeUtils.o\
	collectMovieMod.o \
	mathConstants.o

obj_movie = program_movie.o \
	constantsMod.o \
	parameterMod.o \
	collectMovieMod.o \
	vtkMod.o \
	gllMod.o \
	warpfactorMod.o \
	nodesMod.o \
	triTrafoMod.o \
	jacobiMod.o \
	simplexMod.o \
	vandermondeMod.o \
	dmatricesMod.o \
	geometricFactorsMod.o \
	rosettaGammaMod.o \
	meshMod.o \
	plotMod.o \
	normalsMod.o \
	matrixMod.o \
	dtMod.o \
	liftMod.o \
	sourceReceiverMod.o \
	triangleNeighborMod.o \
	externalModelMod.o \
	linearSystemMod.o \
	fileunitMod.o \
	fileParameterMod.o \
	errorMessage.o \
	realloc.o \
	materialsMod.o \
	progressbarMod.o \
	waveMod.o 

#-------------------------------------------------------
#  Direcory search
#
vpath %.o $(obsdir)
vpath %.f90 $(srcdir) $(srcdir)/include
vpath %.f $(srcdir) $(srcdir)/include
vpath %.c $(srcdir) $(srcdir)/include
#--------------------------------------------------------
#  additional directories to be searched for module or include dependencies
#  default is search in ./ only
#
DEPDIRS = $(srcdir) $(srcdir)/include
#-------------------------------------------------------
#  Implicit rule to compile .o files from .f90 files.
#  Because of vpath, targets and dependencies need not be
#  in the current directory.
#
%.o: %.f90
	$(F95) -c $(FFLAGS) $< -o $(obsdir)/$@
%.o: %.f
	$(F95) -c $(FFLAGS) -fimplicit-none -ffixed-line-length-132 $< -o $(obsdir)/$@
%.o: %.c
	gcc -c $(CFLAGS) $< -o $(obsdir)/$@
#--------------------------------------------------------------
#  Object string for linking:
#  Adds object dir as prefix and removes directory part
#  of $^ (all dependencies)
#
obstring = $(addprefix $(obsdir)/,$(notdir $^))
obstringtest = $(addprefix $(obsdir)/,$(notdir $^))
#
#   End of generic part of Makefile
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  Library paths
#  (specify the path to your METIS library here!)
#
la = -L/es01/software/apps/lapack-3.8.0/lib/ -llapack -lblas
lib = /es01/paratera/sce2567/metis-4.0.3/libmetis.a 
#---------------------------------------------------------
.PHONY: all

#
#  create dependencies on modules 
#  make.incdep is a Makefile because it is included. Such files are first updated
#  before anything else and make starts from anew including all updated makefiles
#
make.incdep:
	./tools/scripts/makeDepFromUseInclude.py $(DEPDIRS) > $@
-include make.incdep
#
#       Targets
#
all: required mesher solver movie 

allclean :
	rm -f $(bindir)/mesher $(bindir)/movie $(moduledir)/*.mod $(obsdir)/*.o out/* make.incdep

clean :
	rm -f $(bindir)/solver $(bindir)/mesher $(bindir)/movie $(moduledir)/*.mod $(obsdir)/*.o make.incdep


dg2d: $(obj_prog)
	$(F95) $(FFLAGS) -o $(bindir)/$@ $(obstring) $(pgplot) $(la) $(lib)

mesher: $(obj_mesher)
	$(F95) $(FFLAGS) -o $(bindir)/$@ $(obstring) $(pgplot) $(la)  $(lib)

solver: $(obj_solver)
	$(F95) $(FFLAGS) -o $(bindir)/$@ $(obstring) $(pgplot) $(la)  $(lib)

movie: $(obj_movie)
	$(F95) $(FFLAGS) -o $(bindir)/$@ $(obstring) $(pgplot) $(la)  $(lib)

bin:
	mkdir -p $(bindir)

mod:
	mkdir -p $(moduledir)

obj:
	mkdir -p $(obsdir)

required: mod obj bin
