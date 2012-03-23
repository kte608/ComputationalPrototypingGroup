include ./def.makefile

FASTMAX_ROOT =  $(shell pwd)

PFFT_ROOT =  $(FASTMAX_ROOT)/pfft
PFFT_LIB_DIR = $(PFFT_ROOT)


FASTMAX_SRC_DIR    = $(FASTMAX_ROOT)/src
FASTMAX_UTIL       = $(FASTMAX_ROOT)/util
FASTMAX_LIB_DIR    = $(FASTMAX_ROOT)/lib
FASTMAX_BIN_DIR    = $(FASTMAX_ROOT)/bin

#SUPER_LU_SRC_DIR    = $(FASTMAX_ROOT)/lib/SuperLU_3.0
SUPER_LU_LIB_DIR    = $(FASTMAX_ROOT)/lib/SuperLU_3.0

FFTW_LIB_DIR = $(FASTMAX_ROOT)/lib/fftw-2.1.5/lib
ITPP_LIB_DIR = $(FASTMAX_ROOT)/lib/itpp-4.0.1/lib
LAPACK_LIB_DIR = $(FASTMAX_ROOT)/lib/lapack-3.1.1

ARCH	:= $(shell $(FASTMAX_UTIL)/config.guess)

all:  $(FASTMAX_BIN_DIR)/$(ARCH)
	${MAKE} -C $(FASTMAX_SRC_DIR);
	${MAKE} -C $(PFFT_ROOT);
	$(C++) -lg2c $(GPROF_FLAG)\
	-o $(FASTMAX_BIN_DIR)/$(ARCH)/FastMaxwell -lm \
	$(FASTMAX_SRC_DIR)/$(ARCH)/fastMax.o \
	$(PFFT_LIB_DIR)/$(ARCH)/pfft.o \
	$(FFTW_LIB_DIR)/libfftw.a \
	$(FFTW_LIB_DIR)/librfftw.a \
	$(SUPER_LU_LIB_DIR)/libsuperlu_3.0.a \
	$(ITPP_LIB_DIR)/libitpp.a \
	$(LAPACK_LIB_DIR)/lapack_LINUX.a \
	$(LAPACK_LIB_DIR)/blas_LINUX.a

#$(SUPER_LU_LIB_DIR)/libblas.a \
#
#$(CLAPACK_LIB_DIR)/F2CLIBS/libI77.a \
#$(CLAPACK_LIB_DIR)/F2CLIBS/libF77.a \
#$(FFTW_LIB_DIR)/libfftw.a \
#$(FFTW_LIB_DIR)/librfftw.a \
#$(SUPER_LU_LIB_DIR)/libblas.a\
#$(CLAPACK_LIB_DIR)/blas_LINUX.a \
#$(CLAPACK_LIB_DIR)/lapack_LINUX.a \



$(FASTMAX_BIN_DIR)/$(ARCH): 
	if [ ! -d $(FASTMAX_BIN_DIR)/$(ARCH) ]; then\
	(mkdir $(FASTMAX_BIN_DIR)/$(ARCH); sleep 1;)\
	fi

$(PFFT_LIB_DIR)/$(ARCH): 
	if [ ! -d $(PFFT_LIB_DIR)/$(ARCH) ]; then\
	(mkdir $(PFFT_LIB_DIR)/$(ARCH); sleep 1;)\
	fi

depend:
	$(MAKE) clean;
	${MAKE} -C $(FASTMAX_SRC_DIR) depend; 
	${MAKE} -C $(MESH_SRC_DIR) depend; 
	${MAKE} -C $(FILEIO_SRC_DIR) depend; 

minorclean:
	${RM} *~; \
	${MAKE} -C $(FASTMAX_SRC_DIR) minorclean; 
	${MAKE} -C $(PFFT_ROOT) minorclean; 

clean:	
	${RM} *~; \
	${MAKE} -C $(FASTMAX_SRC_DIR) clean; 
	${MAKE} -C $(PFFT_ROOT) clean; 
	${RMR} $(FASTMAX_BIN_DIR)/$(ARCH);

realclean:
	${RM} *~; \
	${RMR} $(FASTMAX_BIN_DIR)/$(ARCH); \
	${RM} $(PFFT_ROOT)/$(ARCH); \
	${MAKE} -C $(FASTMAX_SRC_DIR) clean; 
	${MAKE} -C $(PFFT_ROOT) clean;






