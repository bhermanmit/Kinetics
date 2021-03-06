program = cmfd 

templates = $(wildcard templates/*.o)
xml_fort = xml-fortran/xmlparse.o \
           xml-fortran/read_xml_primitives.o \
           xml-fortran/write_xml_primitives.o

#===============================================================================
# Object Files
#===============================================================================

include OBJECTS

#===============================================================================
# User Options
#===============================================================================

COMPILER = petsc
DEBUG    = no
PROFILE  = no
OPTIMIZE = yes
USE_MPI  = no
USE_HDF5 = yes

#===============================================================================
# Add git SHA-1 hash
#===============================================================================

GIT_SHA1 = $(shell git log -1 | head -n 1 | awk '{print $$2}')

#===============================================================================
# GNU Fortran compiler options
#===============================================================================

ifeq ($(COMPILER),gnu)
  F90 = gfortran
  F90FLAGS := -cpp -fbacktrace -DNO_F2008
  LDFLAGS =

  # Debugging
  ifeq ($(DEBUG),yes)
    F90FLAGS += -g -Wall -pedantic -std=f2008 -fbounds-check \
                -ffpe-trap=invalid,zero,overflow,underflow
    LDFLAGS  += -g
  endif

  # Profiling
  ifeq ($(PROFILE),yes)
    F90FLAGS += -pg
    LDFLAGS  += -pg
  endif

  # Optimization
  ifeq ($(OPTIMIZE),yes)
    F90FLAGS += -O3
  endif
endif

#===============================================================================
# Intel Fortran compiler options
#===============================================================================

ifeq ($(COMPILER),intel)
  F90 = ifort
  F90FLAGS := -cpp -warn -assume byterecl -traceback
  LDFLAGS =

  # Debugging
  ifeq ($(DEBUG),yes)
    F90FLAGS += -g -ftrapuv -fp-stack-check -check all -fpe0
    LDFLAGS  += -g
  endif

  # Profiling
  ifeq ($(PROFILE),yes)
    F90FLAGS += -pg
    LDFLAGS  += -pg
  endif

  # Optimization
  ifeq ($(OPTIMIZE),yes)
    F90FLAGS += -O3 -ipo
  endif
endif

#===============================================================================
# PGI compiler options
#===============================================================================

ifeq ($(COMPILER),pgi)
  F90 = pgf90
  F90FLAGS := -Mpreprocess -DNO_F2008 -Minform=inform -traceback
  LDFLAGS =

  # Debugging
  ifeq ($(DEBUG),yes)
    F90FLAGS += -g -Mbounds -Mchkptr -Mchkstk
    LDFLAGS  += -g
  endif

  # Profiling
  ifeq ($(PROFILE),yes)
    F90FLAGS += -pg
    LDFLAGS  += -pg
  endif

  # Optimization
  ifeq ($(OPTIMIZE),yes)
    F90FLAGS += -fast -Mipa
  endif
endif

#===============================================================================
# IBM XL compiler options
#===============================================================================

ifeq ($(COMPILER),ibm)
  F90 = xlf2003
  F90FLAGS := -WF,-DNO_F2008 -O2

  # Debugging
  ifeq ($(DEBUG),yes)
    F90FLAGS += -g -C -qflag=i:i -u
    LDFLAGS  += -g
  endif

  # Profiling
  ifeq ($(PROFILE),yes)
    F90FLAGS += -p
    LDFLAGS  += -p
  endif

  # Optimization
  ifeq ($(OPTIMIZE),yes)
    F90FLAGS += -O3
  endif
endif

#===============================================================================
# Cray compiler options
#===============================================================================

ifeq ($(COMPILER),cray)
  F90 = ftn
  F90FLAGS := -e Z -m 0

  # Debugging
  ifeq ($(DEBUG),yes)
    F90FLAGS += -g -R abcnsp -O0
    LDFLAGS  += -g
  endif
endif

#===============================================================================
# PETSc compiler options
#===============================================================================

ifeq ($(COMPILER),petsc)

  # variables
  PETSC_DIR=/opt/petsc/petsc_gfortran
  PETSC_ARCH=arch-linux2-c-debug
  SLEPC_DIR=/opt/petsc/slepc-3.3-p2

  # include petsc variables
  include ${SLEPC_DIR}/conf/slepc_common

  # standard options
  F90 = ${FC_LINKER} 
  F90FLAGS := -cpp -fbacktrace -DNO_F2008  \
              -I${PETSC_DIR}/include -I${SLEPC_DIR}/include \
              -DMPI -DPETSC
  LDFLAGS = ${SLEPC_LIB}

  # Debugging
  ifeq ($(DEBUG),yes)
    F90FLAGS += -g -Wall -pedantic -fbounds-check \
                -ffpe-trap=invalid,zero,overflow,underflow
    LDFLAGS += -g
  endif

  # Profiling
  ifeq ($(PROFILE),yes)
    F90FLAGS += -pg
    LDFLAGS += -pg
  endif

  # Optimization
  ifeq ($(OPTIMIZE),yes)
    F90FLAGS += -O3
  endif

endif

#===============================================================================
# Miscellaneous compiler options
#===============================================================================

# Use MPI for parallelism

ifeq ($(USE_MPI),yes)
  MPI_ROOT = /opt/mpich2/1.5b2-$(COMPILER)
  F90 = $(MPI_ROOT)/bin/mpif90
  F90FLAGS += -DMPI
endif

# Use HDF5 for output

ifeq ($(USE_HDF5),yes)
  HDF5_ROOT = /opt/hdf5/hdf5-1.8.9/hdf5
  F90FLAGS += -DHDF5 -I${HDF5_ROOT}/include
  LDFLAGS += -L${HDF5_ROOT}/lib ${HDF5_ROOT}/lib/libhdf5hl_fortran.a \
             ${HDF5_ROOT}/lib/libhdf5_hl.a ${HDF5_ROOT}/lib/libhdf5_fortran.a \
             ${HDF5_ROOT}/lib/libhdf5.a -lz -lrt -lm -Wl,-rpath -Wl,${HDF5_ROOT}/lib
endif


#===============================================================================
# Options for IBM Blue Gene/P ANL supercomputer
#===============================================================================

ifeq ($(MACHINE),bluegene)
  F90 = /bgsys/drivers/ppcfloor/comm/xl/bin/mpixlf2003
  F90FLAGS = -WF,-DNO_F2008,-DMPI -O3
endif

#===============================================================================
# Options for Cray XK6 ORNL Titan supercomputer
#===============================================================================

ifeq ($(MACHINE),crayxk6)
  F90 = ftn
  F90FLAGS += -DMPI
endif

#===============================================================================
# Targets
#===============================================================================

all: xml-fortran $(program)
xml-fortran:
	cd xml-fortran; make MACHINE=$(MACHINE) F90=$(F90) F90FLAGS="$(F90FLAGS)"
	cd templates; make F90=$(F90) F90FLAGS="$(F90FLAGS)"
$(program): $(objects)
	$(F90) $(objects) $(templates) $(xml_fort) $(gnuplt) -o $@ $(LDFLAGS)
distclean: clean
	cd xml-fortran; make clean
	cd templates; make clean
clean:
	@rm -f *.o *.mod $(program)
neat:
	@rm -f *.o *.mod

#===============================================================================
# Rules
#===============================================================================

.SUFFIXES: .F90 .o
.PHONY: all xml-fortran clean neat distclean 

%.o: %.F90
	$(F90) $(F90FLAGS) -DGIT_SHA1="\"$(GIT_SHA1)\"" -Ixml-fortran -Itemplates -c $<

#===============================================================================
# Dependencies
#===============================================================================

include DEPENDENCIES
