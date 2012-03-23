include def.makefile

all:
	${MAKE} archs
	${MAKE} dirs
	${MAKE} depend
	${MAKE} objects
	${MAKE} link
	${MAKE} generators

clobber:
	${MAKE} clean
	${MAKE} all

dirs: archs
	@if [ ! -d $(PFFT_ROOT)/bin ]; then \
	  echo 'Creating directory $(PFFT_ROOT)/bin';\
	  mkdir bin;\
	fi
	@if [ ! -d $(PFFT_ROOT)/test ]; then \
	  echo 'Creating directory $(PFFT_ROOT)/test';\
	  mkdir test;\
	fi
	@if [ ! -d $(PFFT_BIN_DIR) ]; then \
	  echo 'Creating directory $(PFFT_BIN_DIR)';\
	  mkdir $(PFFT_BIN_DIR);\
	fi

objects: archs dirs
	@ cd $(CLAPACK_SRC_DIR); ${MAKE}; 
	@ cd $(PFFT_SRC_DIR);    ${MAKE}; 

link: objects
	$(CC) $(CLIBS) $(COPTIM_FLAGS) $(CDEBUG_FLAGS) \
              $(PFFT_SRC_DIR)/$(ARCH)/*.o \
              $(CLAPACK_SRC_DIR)/$(ARCH)/*.o \
              -o $(PFFT_BIN_DIR)/pfft ;
	$(CP) $(PFFT_BIN_DIR)/pfft $(PFFT_TEST_DIR)/.;

depend: archs
	@ cd $(CLAPACK_SRC_DIR); ${MAKE} depend;
	@ cd $(PFFT_SRC_DIR);    ${MAKE} depend;
	@ cd $(GENERAT_SRC_DIR); ${MAKE} depend; 

generators: archs dirs
	@ cd $(GENERAT_SRC_DIR); ${MAKE};

minorclean:
	@ cd $(CLAPACK_SRC_DIR); ${MAKE} minorclean;
	@ cd $(PFFT_SRC_DIR);    ${MAKE} minorclean; 
	@ cd $(GENERAT_SRC_DIR); ${MAKE} minorclean; 
	${RM} $(PFFT_TEST_DIR)/core;
	${RM} *~;
	${RM} $(PFFT_INC_DIR)/*~;

clean:
	@ cd $(CLAPACK_SRC_DIR); ${MAKE} clean;
	@ cd $(PFFT_SRC_DIR);    ${MAKE} clean; 
	@ cd $(GENERAT_SRC_DIR); ${MAKE} clean; 
	${RM} $(PFFT_BIN_DIR)/pfft;
	${RM} $(PFFT_TEST_DIR)/pfft;
	${RM} $(PFFT_TEST_DIR)/core;
	${RM} *~;
	${RM} $(PFFT_INC_DIR)/*~;

realclean:
	@# Making real-clean in each sub-directory
	@ cd $(CLAPACK_SRC_DIR); ${MAKE} realclean; 
	@ cd $(PFFT_SRC_DIR);    ${MAKE} realclean; 
	@ cd $(GENERAT_SRC_DIR); ${MAKE} realclean; 
	${RM} $(PFFT_TEST_DIR)/pfft;
	${RM} $(PFFT_TEST_DIR)/core;
	@echo 'Removing from all sub-dirs all files named *~'
	@$(FIND) . -name '*~' | $(XARGS) -n1 $(RM);
	@# Removing files in architecture-dependent directories:
	-@for i in `cat archs` ; do \
	  if [ -d $(PFFT_BIN_ROOT)/$$i ]; then\
	    cd $(PFFT_BIN_ROOT)/$$i; $(RM) $(EXECS);\
	  fi\
	done
	@# Removing architecture-dependent directories:
	-@for i in `cat archs` ; do \
	  if [ -d $(PFFT_BIN_ROOT)/$$i ]; then\
	    $(RMD) $(PFFT_BIN_ROOT)/$$i;\
	  fi\
	done
	@# Remove any binaries directly in bin:
	cd $(PFFT_BIN_ROOT);${RM} $(EXECS);

tar:
	$(MAKE) realclean;
	$(RM) ./bak/pfft.tar.gz;
	$(TAR) cvf pfft.tar .;
	$(GZIP) -v pfft.tar;
	$(MV) pfft.tar.gz ./bak;
	ls -al ./bak;


archs: always
	@ touch $@;
	@ if \
	  grep $(ARCH) $@ > $(DEVNULL);\
	then\
	  echo Working on architecture: $(ARCH);\
	else\
	  echo $(ARCH) >> $@;\
	  echo "New architecture included in list"; \
	fi

always:

EXECS = pfft cubegen spheregen spheregen2 weavegen
