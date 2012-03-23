include ../../def.makefile

MODULE = pfft

CFLAGS = $(CWARNING_FLAGS) $(COPTIM_FLAGS) $(CDEBUG_FLAGS) -c

IFLAGS = -I. -I$(PFFT_INC_DIR) 
# -I$(PFFT_INC_DIR)

OBJS= main.o                                          \
      stdmath.o elements.o grid.o gridData.o direct.o \
      sparse.o dense.o calcp.o gmres.o nfft.o solve.o \
      bconds.o

OBJNAMES = $(OBJS:%.o=%)


all: $(ARCH) $(OBJS:%=$(ARCH)/%)

$(ARCH)/%.o: %.c
	$(CC) $(IFLAGS) $(CFLAGS) $< -o $@


# Commonly used stuff, same across multiple makefiles.
include $(PFFT_ROOT)/srcmake.incl

include depends

