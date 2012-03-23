
#! /bin/sh

# For debugging:
#GDB_FLAG    = -g
#MORE_GDB_FLAGS    = -fno-for-scope

# For optimizing
#OPTIM_FLAGS_FOR_C++ = -O2
#OPTIM_FLAGS_FOR_C = -O6


# Skip generation of a temp variable in constructing an object
# CONSTRUCTOR = -felide-constructors

# For profiling
# GPROF_FLAG = -pg

# For compilation warning message:
WARNING_FLAGS  = -Wnon-virtual-dtor -Wno-long-long -Wundef -Wpointer-arith -D_XOPEN_SOURCE=500 -D_BSD_SOURCE -Wcast-align -Wconversion -Wchar-subscripts -O2 -fno-exceptions -fno-check-new -fno-common -fexceptions
#MORE_WARNING_FLAGS  = -Wall -pedantic 

LD_ROOT         = /usr/bin
# if you installed the newer version gcc and g++, it will be in here
#GCC_ROOT         = /usr/local/bin
# the default place for gcc and g++
GCC_ROOT         = /usr/bin
BIN_ROOT        = /bin
SHELL		= $(BIN_ROOT)/sh
RM		= $(BIN_ROOT)/rm -f
RMR             = $(BIN_ROOT)/rm -rf
SED        	= $(BIN_ROOT)/sed
MV	        = $(BIN_ROOT)/mv
RANLIB          = /usr/bin/ranlib
AR              = ar
LD              = $(LD_ROOT)/ld
CC              = $(GCC_ROOT)/gcc
C++		 = $(GCC_ROOT)/g++
# DEBUG_ELEMENT	= -DDEBUG_ELEMENT

ROOT = ..
SRC_INC_DIR = $(ROOT)/src
INC_DIR = $(ROOT)/dummy
LIB_DIR = $(ROOT)/lib
BIN_DIR = $(ROOT)/bin
UTIL_DIR = $(ROOT)/util
DEF_MAKEFILE = $(ROOT)/def.makefile

PFFT_ROOT =  $(ROOT)/pfft
PFFT_INC_DIR    = $(PFFT_ROOT)
PFFT_SRC_DIR    = $(PFFT_ROOT)
PFFT_LIB_DIR    = $(PFFT_ROOT)

SUPER_LU_INC_DIR    = $(LIB_DIR)/SuperLU3.0/
ITPP_INC_DIR    = $(LIB_DIR)/itpp-4.0.1/include
FFTW_INC_DIR    = $(LIB_DIR)/fftw-2.1.5/include
