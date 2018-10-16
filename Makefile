INCLUDE_DIR = include
SOURCE_DIR = src
OBJ_DIR = obj
EXEC_DIR = bin
DOCS_DIR = docs
DEP_DIR = .deps
PREPROC_DIR = preprocessed

vpath %.h   $(INCLUDE_DIR)
vpath %.hpp $(INCLUDE_DIR)
vpath %.c   $(SOURCE_DIR)
vpath %.cpp $(SOURCE_DIR)
vpath %.o   $(OBJ_DIR)
vpath %     $(EXEC_DIR)
vpath %.ii  $(PREPROC_DIR)

SHELL = /bin/bash

#CPP = colorgcc
CPP = g++

INCLUDE_OPTIONS = -iquote$(INCLUDE_DIR)
LIB_OPTIONS = -L$(OBJ_DIR)
WARN_OPTION = -Wall  # all warnings
#WARN_OPTION = -w   # no warnings
OPTIMIZE_OPTION = -O3
#OPTIMIZE_OPTION = 
OPTIONS = -std=c++11 $(OPTIMIZE_OPTION) $(WARN_OPTION) $(INCLUDE_OPTIONS) $(LIB_OPTIONS)
#OPTIONS = -std=c++11 $(OPTIMIZE_OPTION) $(WARN_OPTION) $(INCLUDE_OPTIONS) $(LIB_OPTIONS) -static

MAKEDEPEND = gcc $(INCLUDE_OPTIONS) -MM $(CPPFLAGS) -o $(DEP_DIR)/$*.d $<

df = $(DEP_DIR)/$(*F)

%.o : %.cpp | $(OBJ_DIR) $(DEP_DIR) $(EXEC_DIR)
	@$(MAKEDEPEND); \
          cp $(df).d $(df).P; \
          sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
              -e '/^$$/ d' -e 's/$$/ :/' < $(df).d >> $(df).P; \
          rm -f $(df).d
	$(CPP) $(OPTIONS) -c -o $(OBJ_DIR)/$@ $<

-include $(DEP_DIR)/*.P

$(OBJ_DIR):
	mkdir $(OBJ_DIR)

$(DEP_DIR):
	mkdir $(DEP_DIR)

$(EXEC_DIR):
	mkdir $(EXEC_DIR)

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
	rm -f $(OBJ_DIR)/*.o 

piNdilep: units.o utils.o halfint.o Array.o Vectors.o Spinors.o wavefunc.o Config.o RandomNumberGenerator.o Histogram.o piNdilep.o Vrancx.o
	$(CPP) $(OPTIONS) -o $(EXEC_DIR)/$@ $^ -lm

piNdilep_static: units.o utils.o halfint.o Array.o Vectors.o Spinors.o wavefunc.o Config.o RandomNumberGenerator.o Histogram.o piNdilep.o Vrancx.o
	$(CPP) $(OPTIONS) -static -o $(EXEC_DIR)/$@ $^ -lm

