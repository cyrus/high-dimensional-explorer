# 
# Copyright (C) 2004,2005,2006,2007,2008  Cyrus Shaoul and Geoff Hollis 
#
# This file is part of HiDEx.
#
#    HiDEx is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    HiDEx is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with HiDEx in the COPYING.txt file.
#    If not, see <http://www.gnu.org/licenses/>.


##############################################################################
# Makefile
##############################################################################

# This Makefile is simple, and only capable of compiling HiDEx on the following 3 platforms:
#
# Mac OS X (Darwin), Linux and AIX
#
# For instuctions on how to compile and install 
# 
# If you port HiDEx to another platform, please contribute your additions to the Makefile to 
# the HiDEx maintainers.
#
#  For more info, see: http://www.psych.ualberta.ca/~westburylab/HiDEx/
#

# Change this line to where you would like to install the program. Make sure
# That it is a directory that is listed in your $PATH.

BINDIR=~/bin

# For system-wide usage, use this directory.
# BINDIR=/usr/local/bin

INSTALL=cp

#Was your OS supported? No? Then you get the wimpy default compiler flags:
CC = g++
FLAGS = -Wall -g 

# Shell command to detect Operating System.
UNAME= $(shell uname)

#Mac OS X PPC specific comiler flags. DEFAULT: Configured for OpenMP support.
ifeq ($(UNAME), Darwin)
CC = g++
FLAGS = -Wall -g -gfull -O3 -falign-loops=16 -D NDEBUG -fopenmp
#To remove OpenMP support, use this line:
#FLAGS = -Wall -g -gfull -O3 -falign-loops=16 -D NDEBUG 
LIBS =  -lm 
endif

# Linux Specific compiler Flags. Configured for OpenMP support. Remove OpenMP flags if you don't need parallel support.
ifeq ($(UNAME), Linux)
CC = g++
# Debug only 
# FLAGS = -Wall -g -O0 -Wextra -fopenmp
# Debug only 
FLAGS = -Wall -g -O3 -Wextra -fopenmp -ftree-vectorize -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops -funit-at-a-time -D NDEBUG
endif

#AIX specific compiler flags. Configured for OpenMP support.
ifeq ($(UNAME), AIX)
CC = xlC_r
FLAGS = -g -O5 -qmaxmem=-1 -q64 -qsmp=omp -I/home/cyrus/include/ -qnoinline -qnoipa
endif

# Payloads below
TARGET = hidex

OBJS = \
  SDDB.o \
  ConfigFile.o \
  Exception.o \
  utilities.o\
  Filesystem.o\
  SDDBAccessor.o \
  hidex.o 

default: $(TARGET)


$(TARGET): $(OBJS)
	$(CC) $(FLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

.SUFFIXES:
.SUFFIXES: .o .cpp

.PHONY: clean test

.cpp.o:
	$(CC) $(FLAGS) -o $@ -c $<

# make dependancies for all of the files
depend:
	g++ -MM *.cpp test/*.cpp -I/usr/global/ibm/boost/boost/ > depend

# clean up everything
clean:
	rm -f *.o test/*.o *~ hidex depend hidex-test

#Testing system. Do not use!
#test: hidex-test
#	./hidex-test

#Install
install:
	@echo WARNING: You might need to be root to install this program. 
	@echo Change the BINDIR  in the Makefile if you have trouble.
	@echo
	$(INSTALL) -v hidex $(BINDIR)
	@echo Installed HiDEx in the direectory that you specified. Please make sure
	@echo that this directory is in your PATH.


















# Old compiler settings that were used during the development of this Makefile. Not very interesting to 
# anyone but me. 


#Old Mac Stuff.

# flags to use during compilation
#FLAGS = -O3 -Wall -g -Wconversion -DNDEBUG
#FLAGS = -fast -mdynamic-no-pic -Wall -g -arch ppc64 -L/usr/local/lib 
#-lMallocDebug
# -D_GLIBCXX_DEBUG 
#CC = g++
#FLAGS = -Wall -fast -gfull -arch ppc64 -fno-inline
# -D_GLIBCXX_DEBUG
#FLAGS = -Wall -O3 -ftree-vectorize -g -gfull -fopenmp
#FLAGS = -Wall -g -gfull -O3 -arch ppc64 -pg
#LIBS =  -lm -lboost_filesystem


#boost stuff 
#	g++ -MM *.cpp test/*.cpp -I/usr/global/ibm/boost/test/include/boost-1_33_1/ > depend
#	$(CC) -qmakedep=gcc *.cpp test/*.cpp -I/usr/global/ibm/boost/boost/ > depend


#Linux Stuff
#FLAGS = -Wall -g -O3 -Wextra -fopenmp -ftree-vectorize -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops -funit-at-a-time -D NDEBUG
#FLAGS = -Wall -g -Wextra -D_GLIBCXX_DEBUG -fstack-protector-all
#LIBS =  -lm 
#FLAGS = -O3 -mtune=opteron -fopenmp -ftree-vectorize -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops -funit-at-a-time -D NDEBUG
#FLAGS = -openmp -O3 -g -ipo -i-static -axP -xW -DNDEBUG
#FLAGS = -O0 -g -openmp -i-static 
#FLAGS = -openmp -O3 -ipo -axP -xW -i-static 
#-openmp-report2
#FLAGS = -mp -tp k8-64 -fastsse
#FLAGS = -Wall -g -O3 -Wextra -fno-strict-aliasing -m64 -mtune=opteron -ffast-math 
#FLAGS = -Wall -g -Wextra -D_GLIBCXX_DEBUG -fstack-protector-all -D NDEBUG
#TEST_LIBS = -lboost_unit_test_framework -lboost_filesystem
#LIBS =  -lm -lacml_mp -lpgftnrtl -lm
#LIBS =  -lm 
#-lmkl_em64t -lguide -lpthread -Wl,-rpath,$INTEL_LIB_PATH 


#AIX stuff
#FLAGS = -g -qsmp=omp -O3 -q64 -I/home/cyrus/include/ -qnoinline
#FLAGS = -g -qnoipa -qNOOPTimize -qsmp=noopt -qextchk -qcheck=all -qmaxmem=-1 -q64 -I/home/cyrus/include/
#-qsmp=omp FLAGS
#FLAGS = -g -qsmp=omp -qnoipa -qnooptimize -qextchk -qcheck=all -qmaxmem=-1 -q64 -I/home/cyrus/include/
#FLAGS = -g -O4 -qmaxmem=-1 -qlargepage -qsmp=omp -q64 -qnoinline -I/home/cyrus/include/
#LIBS = -bstatic -lmass -L/home/cyrus/lib/ -bstatic -lboost_filesystem-xlc-mt-d -bdynamic
#LIBS =  -lmass -L/usr/global/ibm/boost/test/lib -lboost_filesystem-xlc-d
#LIBS =  -lmass -L/usr/global/ibm/boost/boost/lib64/ -lboost_filesystem-xlc-mt-d
#-qextchk 
#-bdynamic
#-L/home/cyrus/lib -lboost_filesystem-xlc-d
#-qinfo=all
#-I/usr/global/ibm/boost/boost/
#-qsmp 
#-O5 -qlargepage
#-pg
#-brtl
#-I /home/cyrus/boost/include/boost-1_33_1/ 
#-L/usr/local/lib
#-L/home/cyrus/boost/lib
#-q64
#-qsmp=omp:noauto
#-qinfo=all 	Produces additional warnings, similar to those generated by Lint.
#-qcheck=all 	Compiles the program to check for "null" pointers, object bounds, and division by zero at run time.
#-qflttrap 	Compiles the program to detect floating-point exceptions at run time. The following form of this option causes the program to abort on floating-point overflow or division by zero.
#-bmaxdata:256000000000 not in -q64 mode
#      -qflttrap=overflow:zerodivide:enable 
#print:
#	echo $(TEST_FLAGS)
