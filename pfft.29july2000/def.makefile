# Please go through this file and choose values 
# that are consitent with your local system.
#
# PFFT_ROOT must be your working directory.
# It might be a bad idea to use ~'s. A full 
# path starting with / is recommended, e.g:
PFFT_ROOT       = /u0/buchmann/pFFTkernel
#PFFT_ROOT       = /u0/buchmann/pFFTkernel

# If you don't have gmake, then the current hierarchy 
# may not work. "make -f makesimple" may then work
# but it will not give the same nice directory structure.
MAKE		= gmake

SHELL		= /bin/sh
RM		= /bin/rm -f
RMD		= /bin/rmdir
SED        	= /bin/sed
MV	        = /bin/mv
CP	        = /bin/cp

FIND            = /usr/bin/find
XARGS           = /usr/bin/xargs
RANLIB          = /usr/bin/ranlib

TAR		= /bin/tar
GZIP		= /bin/gzip

DEVNULL		= /dev/null

GCC_EXISTS	= true # Use "false" if no gcc!
GCC		= gcc

# Your favourite C compiler goes here:
CC              = cc
#CC              = gcc

# For ix36-linux machines (and others with gcc):
#CC = $(GCC)
# For debugging:
#CDEBUG_FLAGS    = -ggdb 
#CWARNING_FLAGS  = -ansi -pedantic                        \
#                  -W -Wall -Wtraditional -Wmissing-prototypes \
#                  -Wmissing-declarations
# For running and profiling (choose one). These options will 
# be used for all compilation and linking.
COPTIM_FLAGS    = -O2
#COPTIM_FLAGS    = -O3
#COPTIM_FLAGS    = -O3 -Wuninitialized
#COPTIM_FLAGS    = -O3 -pg
#COPTIM_FLAGS    = -pg


## For alpha-linux machines using the Compaq C Compiler (ccc)
## You CAN put stuff linke this e.g. in "opts.alpha-linux" or 
## a similar file corresponding to your operating system.
#CC              = ccc
#COPTIM_FLAGS    = -O4 -tune ev6 -arch ev6
#CDEBUG_FLAGS    = -O0 -g
#CWARNING_FLAGS  = 

# Libraries:
CLIBS = -L/usr/lib -lc -lm 


## For alpha-DEC machines: 
# [Apparently this makefile setup does not work for dec-alphas.
#  The util/config.guess seems to fail. 
#  You may try "make -f makesimple" to get started]
#CC = cc
#COPTIM_FLAGS    = -O3 


# These directories should not need change.
PFFT_SRC_DIR    = $(PFFT_ROOT)/src/pfft
CLAPACK_SRC_DIR = $(PFFT_ROOT)/src/clapack
GENERAT_SRC_DIR = $(PFFT_ROOT)/src/generate
PFFT_TEST_DIR   = $(PFFT_ROOT)/test
PFFT_UTIL       = $(PFFT_ROOT)/util
PFFT_INC_DIR    = $(PFFT_ROOT)/include
PFFT_BIN_ROOT   = $(PFFT_ROOT)/bin
PFFT_BIN_DIR    = $(PFFT_BIN_ROOT)/$(ARCH)

# This is to remove any default dependencies:
.SUFFIXES:


# If the following fails, then please add your architecture 
# explicitly, e.g.  ARCH = ix86-linux
ARCH	:= $(shell $(PFFT_UTIL)/config.guess)
#ARCH = unknown

# Include possible architecture dependent stuff. 
# Uncomment this line only if you know what you are 
# doing and need added flexibility.
#include $(PFFT_ROOT)/opts.$(ARCH)
