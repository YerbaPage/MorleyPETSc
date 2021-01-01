# PETSc paths and packages
export PETSC_DIR = /home/daniel/petsc-3.14.2
export PETSC_ARCH = arch-linux-c-debug
petsc.pc := $(PETSC_DIR)/$(PETSC_ARCH)/lib/pkgconfig/petsc.pc
# petsc.pc := /home/daniel/petsc-3.14.2/arch-linux-c-debug/lib/pkgconfig/petsc.pc

# Additional libraries that support pkg-config can be added to the list of PACKAGES below.
PACKAGES := $(petsc.pc)

########################################################################
# Compiler and external dependences
########################################################################
# CC = cc
# FC = ftn
# CC = pgcc
# FC = pgf90
# CC = icc
# FC = ifort -nofor-main
# CC = gcc
# FC = gfortran -fcray-pointer
# FC = gcc
# AR = ar ruc

# FC = gcc
# CC = gcc
AR = ar ruc

# Compiler for PETSc
CC := $(shell pkg-config --variable=ccompiler $(PACKAGES))
CXX := $(shell pkg-config --variable=cxxcompiler $(PACKAGES))
FC := $(shell pkg-config --variable=fcompiler $(PACKAGES))
CFLAGS_OTHER := $(shell pkg-config --cflags-only-other $(PACKAGES))
CFLAGS := $(shell pkg-config --variable=cflags_extra $(PACKAGES)) $(CFLAGS_OTHER)
CXXFLAGS := $(shell pkg-config --variable=cxxflags_extra $(PACKAGES)) $(CFLAGS_OTHER)
FFLAGS := $(shell pkg-config --variable=fflags_extra $(PACKAGES))
CPPFLAGS := $(shell pkg-config --cflags-only-I $(PACKAGES))
LDFLAGS := $(shell pkg-config --libs-only-L --libs-only-other $(PACKAGES))
LDFLAGS += $(patsubst -L%, $(shell pkg-config --variable=ldflag_rpath $(PACKAGES))%, $(shell pkg-config --libs-only-L $(PACKAGES)))
LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lm


########################################################################
# Compiling and linking options
########################################################################
#COPTS     = -g -pedantic -Wall 
#COPTS     = -g -cc=icc -config=icc -O3
#COPTS     = -g -config=icc -O3
#COPTS     = -g -O3

COPTS     = $(BOPT)
CINCLUDES = -I./include
CDEFS     = 
CFLAGS    = $(COPTS) $(CINCLUDES) $(CDEFS)

FOPTS     = $(BOPT)
FINCLUDES = $(CINCLUDES)
FDEFS     = $(CDEFS)
FFLAGS    = $(FOPTS) $(FINCLUDES) $(FDEFS)

AMGLIB = lib/libAMG.a

LINKOPTS  = $(COPTS) #-Mnomain
LIBS      = -lm

CLFLAGS    = $(LINKOPTS) $(LIBS) -lstdc++
FLFLAGS    = $(LINKOPTS) $(LIBS)


########################################################################
# Rules for compiling the source files
########################################################################
.SUFFIXES: .c .for .f

CSRCDIR = ./csrc
FSRCDIR = ./fsrc

FSRC := $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.for))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f))
CSRC := $(foreach dir,$(CSRCDIR),$(wildcard $(CSRCDIR)/*.c))

OBJSF := $(patsubst %.for,%.o,$(FSRC))
OBJSF += $(patsubst %.f,%.o,$(FSRC))
OBJSC := $(patsubst %.c,%.o,$(CSRC))

########################################################################
# List of all programs to be compiled
########################################################################

# Everything
#ALLPROG = $(TESTPROG) ($TESTMATPROG)
ALLPROG = $(AMGLIB) test testmat

# Test for solvers
TESTPROG = $(AMGLIB) main/test.o 

# Test for matrix properties
TESTMATPROG = $(AMGLIB) main/testmat.o


########################################################################
# Link
########################################################################

# Q: 1. LINK.F 哪里来的 2. 下面的编译
########################################################################
# From PETSc
% : %.F90
	$(LINK.F) -o $@ $^ $(LDLIBS)
%.o: %.F90
	$(COMPILE.F) $(OUTPUT_OPTION) $<
% : %.cxx
	$(LINK.cc) -o $@ $^ $(LDLIBS)
%.o: %.cxx
	$(COMPILE.cc) $(OUTPUT_OPTION) $<
########################################################################


all: $(ALLPROG)

Default: test # bugs here?

$(AMGLIB): $(OBJSC) $(OBJSF)
	ranlib $(AMGLIB)

test: 
	$(CC) $(CFLAGS) -c main/test.c  -o main/test.o
	$(FC) main/test.o $(AMGLIB) $(FLFLAGS) -o test.ex

testmat: 
	$(CC) $(CFLAGS) -c main/testmat.c  -o main/testmat.o
	$(FC) main/testmat.o $(AMGLIB) $(FLFLAGS) -o testmat.ex


########################################################################
# Clean up
########################################################################

.PHONY : clean allclean # 伪目标, 不生成文件, 避免和其他 文件重名

clean:
	rm -f main/*.o
	rm -f $(CSRCDIR)/*.o
	rm -f $(FSRCDIR)/*.o
	rm -f lib/*.a

allclean:
	make clean
	rm -f *~
	rm -f main/*.o
	rm -f *.ex

.SUFFIXES: .c .for .f .o

.for.o:
	$(FC) $(FFLAGS) -c $<  -o $@
	$(AR) $(AMGLIB) $@

.f.o:
	$(FC) $(FFLAGS) -c $<  -o $@
	$(AR) $(AMGLIB) $@

.c.o:
	$(CC) $(CFLAGS) -c $<  -o $@
	$(AR) $(AMGLIB) $@
