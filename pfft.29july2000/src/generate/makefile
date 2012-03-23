include ../../def.makefile

MODULE = generators

#FLAGS = -ansi -pedantic -ggdb -pg

CFLAGS = $(CWARNING_FLAGS) $(COPTIM_FLAGS) $(CDEBUG_FLAGS) -c

IFLAGS = -I. -I$(PFFT_INC_DIR) 
# -I$(PFFT_INC_DIR)

LFLAGS = $(CLIBS)

COMMONOBJS = 
CUBEOBJS   = cubegen.o 
SPHEREOBJS = spheregen.o
SPHEREOBJS2 = spheregen2.o 
WEAVEOBJS   = weavegen.o

OBJS= $(COMMONOBJS) $(CUBEOBJS) $(SPHEREOBJS) $(SPHEREOBJS2)  \
      $(WEAVEOBJS)

$(ARCH)/%.o: %.c
	$(CC) $(IFLAGS) $(CFLAGS) $< -o $@


all: $(ARCH) $(OBJS:%=$(ARCH)/%) $(PFFT_BIN_DIR)/cubegen $(PFFT_BIN_DIR)/spheregen $(PFFT_BIN_DIR)/spheregen2 $(PFFT_BIN_DIR)/weavegen


$(PFFT_BIN_DIR)/cubegen: $(COMMONOBJS:%=$(ARCH)/%) $(CUBEOBJS:%=$(ARCH)/%)
	$(CC) $(IFLAGS) \
              -o $(PFFT_BIN_DIR)/cubegen \
	      $(COMMONOBJS:%=$(ARCH)/%) $(CUBEOBJS:%=$(ARCH)/%) \
	      $(LFLAGS) 

$(PFFT_BIN_DIR)/spheregen: $(COMMONOBJS:%=$(ARCH)/%) $(SPHEREOBJS:%=$(ARCH)/%)
	$(CC) $(IFLAGS) \
              -o $(PFFT_BIN_DIR)/spheregen \
              $(COMMONOBJS:%=$(ARCH)/%) $(SPHEREOBJS:%=$(ARCH)/%)\
	      $(LFLAGS) 

$(PFFT_BIN_DIR)/spheregen2: $(COMMONOBJS:%=$(ARCH)/%) $(SPHEREOBJS2:%=$(ARCH)/%)
	$(CC) $(IFLAGS) \
              -o $(PFFT_BIN_DIR)/spheregen2 \
              $(COMMONOBJS:%=$(ARCH)/%) $(SPHEREOBJS2:%=$(ARCH)/%)\
	      $(LFLAGS) 

$(PFFT_BIN_DIR)/weavegen: $(COMMONOBJS:%=$(ARCH)/%) $(WEAVEOBJS:%=$(ARCH)/%)
	$(CC) $(IFLAGS) \
              -o $(PFFT_BIN_DIR)/weavegen \
              $(COMMONOBJS:%=$(ARCH)/%) $(WEAVEOBJS:%=$(ARCH)/%)\
	      $(LFLAGS) 


# Commonly used stuff, same across multiple makefiles.
include $(PFFT_ROOT)/srcmake.incl

include depends
