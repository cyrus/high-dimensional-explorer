# 
# Copyright (C) 2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014  Cyrus Shaoul and Geoff Hollis 
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
# Mac OS X (Darwin) and Linux 
#
# For instuctions on how to compile and install HiDEx, please see the documentation provided
# in the file called HowToUseHiDEx.txt
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
#CC = g++
CC = g++-mp-4.9
FLAGS = -Wall -g -gfull -O3 -D NDEBUG -fopenmp
LIBS = -L/opt/local/lib -lm -lunistring 
endif

# Linux Specific compiler Flags. Configured for OpenMP support. Remove OpenMP flags if you don't need parallel support.
ifeq ($(UNAME), Linux)
CC = g++
# Debug only 
#FLAGS = -Wall -g -O0 -Wextra -fopenmp
# Normal Flags
FLAGS = -g -O3 -fopenmp -ffast-math -mtune=native -D NDEBUG
LIBS =  -lunistring 
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

# clean up everything
clean:
	rm -f *.o test/*.o *~ hidex 

#Install
install:
	@echo WARNING: You might need to be root to install this program. 
	@echo Change the BINDIR  in the Makefile if you have trouble.
	@echo
	$(INSTALL) -v hidex $(BINDIR)
	@echo Installed HiDEx in the direectory that you specified. Please make sure
	@echo that this directory is in your PATH.
