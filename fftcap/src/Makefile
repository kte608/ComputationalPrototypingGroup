# Makefile for
#  fastcap - capacitance calculation program
#  busgen - bus crossing input file generator
#  capgen - || plate capacitor input file generator
#  cubegen - cubic capacitor input file generator
#  pipedgen - parallelepiped input file generator
#  pyragen - pyramid input file generator

# assumed places for unaliased executables - must be same in all Makefile.*
SHELL = /bin/sh
AWK = awk
RM = /bin/rm

MAKETYPE = other
# normally overidden by master Makefile
TYPE = other

INCL = .
IO = .
DIR = .
MUL = .
ALC = .
COM = .
GEN = .
BIN = ../bin
TOP = ..
FFT = .

OBJS = $(MUL)/mulGlobal.o $(MUL)/mulSetup.o $(MUL)/mulMats.o $(MUL)/mulDo.o \
$(IO)/mulDisplay.o $(MUL)/mulMulti.o \
$(MUL)/mulLocal.o $(DIR)/direct.o $(DIR)/calcp.o $(MUL)/fastcap.o \
$(DIR)/blkDirect.o $(ALC)/uglyalloc.o \
$(IO)/patran.o $(IO)/quickif.o $(MUL)/electric.o \
$(IO)/input.o $(IO)/psMatDisplay.o $(MUL)/capsolve.o $(IO)/savemat_mod.o \
$(IO)/zbufInOut.o $(IO)/zbuf2fastcap.o $(IO)/zbufProj.o $(IO)/zbufSort.o \
$(MUL)/mulGrid.o  nconvolve.o dsvdc.o nfft.o

OBJSH = $(MUL)/mulGlobal.o $(MUL)/mulSetup.o $(MUL)/mulMats.o $(MUL)/mulDo.o \
$(IO)/mulDisplay.o $(MUL)/mulMulti.o \
$(MUL)/mulLocal.o $(DIR)/direct.o $(DIR)/calcp_helm.o $(MUL)/fastcap.o \
$(DIR)/blkDirect.o $(ALC)/uglyalloc.o \
$(IO)/patran.o $(IO)/quickif.o $(MUL)/electric.o \
$(IO)/input.o $(IO)/psMatDisplay.o $(MUL)/capsolve.o $(IO)/savemat_mod.o \
$(IO)/zbufInOut.o $(IO)/zbuf2fastcap.o $(IO)/zbufProj.o $(IO)/zbufSort.o \
$(MUL)/mulGrid.o $(FFT)/rlft3.o $(FFT)/fourn.o

SRCS = $(MUL)/mulGlobal.c $(MUL)/mulSetup.c $(MUL)/mulMats.c $(MUL)/mulDo.c \
$(IO)/mulDisplay.c $(MUL)/mulMulti.c \
$(MUL)/mulLocal.c $(DIR)/direct.c $(DIR)/calcp.c $(MUL)/fastcap.c \
$(DIR)/blkDirect.c $(ALC)/uglyalloc.c \
$(IO)/patran.c $(IO)/quickif.c $(MUL)/electric.c \
$(IO)/input.c $(IO)/psMatDisplay.c $(MUL)/capsolve.c $(IO)/savemat_mod.c \
$(IO)/zbufInOut.c $(IO)/zbuf2fastcap.c $(IO)/zbufProj.c $(IO)/zbufSort.c \
$(MUL)/mulGrid.c $(FFT)/rlft3.c $(FFT)/fourn.c nconvolve.c 

OBJSB = $(GEN)/busgen.o $(GEN)/disrect.o $(COM)/epsilon.o
SRCSB = $(GEN)/busgen.c $(GEN)/disrect.c $(COM)/epsilon.c

OBJSB3d = $(GEN)/bus3d.o $(GEN)/disrect.o $(COM)/epsilon.o
SRCSB3d = $(GEN)/bus3d.c $(GEN)/disrect.c $(COM)/epsilon.c

OBJSP = $(GEN)/pltcapgen.o $(GEN)/disrect.o $(COM)/epsilon.o
SRCSP = $(GEN)/pltcapgen.c $(GEN)/disrect.c $(COM)/epsilon.c

OBJSC = $(GEN)/cubegen.o $(GEN)/disrect.o $(COM)/epsilon.o
SRCSC = $(GEN)/cubegen.c $(GEN)/disrect.c $(COM)/epsilon.c

OBJST = $(GEN)/teragen.o $(GEN)/disrect.o $(GEN)/distri.o $(COM)/epsilon.o
SRCST = $(GEN)/teragen.c $(GEN)/disrect.c $(GEN)/distri.c $(COM)/epsilon.c

OBJSI = $(GEN)/pipedgen.o $(GEN)/disrect.o $(COM)/epsilon.o
SRCSI = $(GEN)/pipedgen.c $(GEN)/disrect.c $(COM)/epsilon.c

OBJSW = $(GEN)/weavegen.o 
SRCSW = $(GEN)/weavegen.c 

LIBS = -lm -lc 

# These compile flags are needed to account for differences in the
#  resource usage system call used to get CPU time and memory use data
#  - also used to select the correct CC macro definition for AIX machines
#
# for generic machine use (timers not supported)
#CFLAGS = -O -DOTHER
# for most machines running 4.2/3 use (only necessary if timers are needed)
CFLAGS = -O -DFOUR 
# for DEC alpha 
#CFLAGS = -O -Olimit 5000 -DFOUR 
# for most machines running System V use (only necessary if timers are needed)
#CFLAGS = -O -DOTHER
# for IBM AIX systems use
#CC = xlc

SO = /usr/opt/XMDLOA250/shlib/libdxml.so 

def:    fftcap

fftcap: $(OBJS)
	$(CC) -o fftcap $(CFLAGS) $(OBJS) $(LIBS)

install: fftcap 
	mv fftcap ../bin 

clean:
	$(RM) -f *.o

$(OBJS) : $(INCL)/mulStruct.h $(INCL)/mulGlobal.h $(INCL)/patran.h \
$(INCL)/resusage.h 
$(OBJSH) : $(INCL)/mulStruct.h $(INCL)/mulGlobal.h $(INCL)/patran.h \
$(INCL)/resusage.h 
$(SRCS) : $(INCL)/mulStruct.h $(INCL)/mulGlobal.h $(INCL)/patran.h \
$(INCL)/resusage.h 
$(IO)/quickif.c $(IO)/quickif.o : $(IO)/quickif.h 



