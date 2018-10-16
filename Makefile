WORKDIR = /home/zetenyi/work
ROOTNAME = pi-N
ROOTDIR = $(WORKDIR)/$(ROOTNAME)
INCLUDE_DIR = $(ROOTDIR)/include
SOURCE_DIR = $(ROOTDIR)/src
LIB_DIR = $(ROOTDIR)/lib
EXEC_DIR = $(ROOTDIR)/bin
DOCS_DIR = $(ROOTDIR)/docs
DEP_DIR = $(ROOTDIR)/.deps
PREPROC_DIR = $(ROOTDIR)/preprocessed

vpath %.h   $(INCLUDE_DIR)
vpath %.hpp $(INCLUDE_DIR)
vpath %.c   $(SOURCE_DIR)
vpath %.cpp $(SOURCE_DIR)
vpath %.o   $(LIB_DIR)
vpath %     $(EXEC_DIR)
vpath %.ii  $(PREPROC_DIR)

SHELL = /bin/bash

#CPP = colorgcc
CPP = g++

INCLUDE_OPTIONS = -iquote$(INCLUDE_DIR)
LIB_OPTIONS = -L$(LIB_DIR)
WARN_OPTION = -Wall  # all warnings
#WARN_OPTION = -w   # no warnings
OPTIMIZE_OPTION = -O3
#OPTIMIZE_OPTION = 
OPTIONS = -std=c++11 $(OPTIMIZE_OPTION) $(WARN_OPTION) $(INCLUDE_OPTIONS) $(LIB_OPTIONS)
#OPTIONS = -std=c++11 $(OPTIMIZE_OPTION) $(WARN_OPTION) $(INCLUDE_OPTIONS) $(LIB_OPTIONS) -static

MAKEDEPEND = gcc $(INCLUDE_OPTIONS) -MM $(CPPFLAGS) -o $(DEP_DIR)/$*.d $<

df = $(DEP_DIR)/$(*F)

%.o : %.cpp
	@$(MAKEDEPEND); \
          cp $(df).d $(df).P; \
          sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
              -e '/^$$/ d' -e 's/$$/ :/' < $(df).d >> $(df).P; \
          rm -f $(df).d
	$(CPP) $(OPTIONS) -c -o $(LIB_DIR)/$@ $<

-include $(DEP_DIR)/*.P

%.ii: %.cpp
	$(CPP) $(OPTIONS) -E $< > $(PREPROC_DIR)/$@

# documentation with Doxygen:
DOXY = /usr/bin/doxygen

#docs: Doxyfile **/*.cpp **/*.h
.PHONY: docs
docs:
	rm -rf $(DOCS_DIR)/*
	$(DOXY) Doxyfile

.PHONY: clean
clean:
	rm -f $(LIB_DIR)/*.o 

piNdilep: units.o utils.o halfint.o Array.o Vectors.o Spinors.o wavefunc.o Config.o RandomNumberGenerator.o Histogram.o piNdilep.o Vrancx.o
	$(CPP) $(OPTIONS) -o $(EXEC_DIR)/$@ $^ -lm

piNdilep_static: units.o utils.o halfint.o Array.o Vectors.o Spinors.o wavefunc.o Config.o RandomNumberGenerator.o Histogram.o piNdilep.o Vrancx.o
	$(CPP) $(OPTIONS) -static -o $(EXEC_DIR)/$@ $^ -lm

