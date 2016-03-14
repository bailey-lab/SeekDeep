
ROOT = $(realpath ./)
ifdef CXXFLAGS 
	ENV_CXXFLAGS := $(CXXFLAGS)
endif
ifneq (,$(wildcard compfile.mk))
COMPFILE=compfile.mk
endif
include $(COMPFILE)
#include 
include $(ROOT)/makefile-common.mk
#os name
UNAME_S := $(shell uname -s)
# header files
HEADERS = $(call rwildcard, src/, *.h) \
	$(call rwildcard, src/, *.hpp)
## Compile Object dir and files
OBJ_DIR = $(addprefix build/, $(addsuffix Build, $(CXXOUTNAME)))
OBJ = $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(call rwildcard, src/, *.cpp))) 
OBJNOMAIN = $(filter-out $(addsuffix /src/main.o, $(OBJ_DIR)), $(OBJ))
## Out binary name
BIN = $(addsuffix $(CXXOUTNAME), bin/)
## Library Names
LIB_DIR=$(ROOT)/lib
LIBNAME = $(addsuffix $(CXXOUTNAME), lib)
DYLIB = $(addprefix $(addsuffix $(LIBNAME), $(LIB_DIR)/), .dylib)
SOLIB = $(addprefix $(addsuffix $(LIBNAME), $(LIB_DIR)/), .so)

##Phony Targets
.PHONY: all
.PHONY: docs
.PHONY: do_preReqs 
.PHONY: sharedLibrary
.PHONY: dyLibrary
.PHONY: clean
.PHONY: install
.PHONY: cpHeaders
.PHONY: cpEtc
.PHONY: unitTest

## unit tert dir 
TESTDIR=test

###compiler options
CXXFLAGS += $(ENV_CXXFLAGS)
#CXXFLAGS += -Wno-missing-braces
COMMON = $(CXXFLAGS) $(CXXOPT) $(COMLIBS)
-include do_preReqs

############ default for calling make with no arguments
all: do_preReqs $(OBJ_DIR) $(BIN) 
ifeq ($(UNAME_S), Darwin)
	scripts/setUpScripts/fixDyLinking_mac.sh bin $(EXT_PATH)
endif

######### docs
docs: docs/Doxygen
	doxygen docs/Doxygen
	
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)
	mkdir -p bin
	mkdir -p lib

# using automatic variables $<: the name of the prerequisite of the rule and
#							$@: the name of the target of the rule 
$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(OBJ_DIR)/$(shell dirname $<)
	$(CXX) -DNOT_HEADER_ONLY $(COMMON) -fPIC -c $< -o $@

$(BIN): $(OBJ) 
	$(CXX) $(CXXFLAGS) $(CXXOPT) -o $@ $^ $(LD_FLAGS) 

############ remove the objects that were dependant on any changed headers and check for compfile.mk
do_preReqs: 
ifndef COMPFILE
	$(error compfile is not set, do either make COMPFILE=aCompfile.mk or create a file called compfile.mk)
endif
	scripts/setUpScripts/rmNeedToRecompile.py -obj $(OBJ_DIR) -src src/

	
############ shared library
sharedLibrary: do_preReqs $(OBJ_DIR) $(SOLIB)
ifeq ($(UNAME_S), Darwin)
	scripts/setUpScripts/fixDyLinking_mac.sh lib $(EXT_PATH)
endif

$(SOLIB): $(OBJNOMAIN)
	$(CXX) $(CXXFLAGS) $(CXXOPT) -shared -o $@ $^ $(LD_FLAGS) 
	
############ dyLibrary
dyLibrary: do_preReqs $(OBJ_DIR) $(DYLIB)
ifeq ($(UNAME_S), Darwin)
	scripts/setUpScripts/fixDyLinking_mac.sh lib $(EXT_PATH)
endif
$(DYLIB): $(OBJNOMAIN)
	$(CXX) $(CXXFLAGS) $(CXXOPT) -dynamiclib -o $@ $^ $(LD_FLAGS) 


############ clean
clean:
	@rm -f $(BIN)
	@rm -rf $(OBJ_DIR)
	
############ install
INSTALL_OUTNAME=$(INSTALL_DIR)/bin/$(CXXOUTNAME) 
INSTALL_DYLIBNAME=$(INSTALL_DIR)/lib/$(LIBNAME).dylib
INSTALL_SHAREDLIBNAME=$(OBJ_DIR) $(INSTALL_DIR)/lib/$(LIBNAME).so 
INSTALL_FILES=$(INSTALL_OUTNAME) $(INSTALL_DYLIBNAME) $(INSTALL_SHAREDLIBNAME)

install: $(INSTALL_DIR) cpHeaders cpEtc do_preReqs $(INSTALL_FILES)
ifeq ($(UNAME_S), Darwin)
	scripts/setUpScripts/fixDyLinking_mac.sh $(INSTALL_DIR)/lib/ $(EXT_PATH)
	scripts/setUpScripts/fixDyLinking_mac.sh $(INSTALL_DIR)/bin/ $(EXT_PATH)
endif

#### install directories set up
$(INSTALL_DIR): 
	@mkdir -p $(INSTALL_DIR)
	@mkdir -p $(INSTALL_DIR)/include
	@mkdir -p $(INSTALL_DIR)/bin
	@mkdir -p $(INSTALL_DIR)/lib
	
#### installing shared library
$(INSTALL_DIR)/lib/$(LIBNAME).so: $(OBJNOMAIN)
	$(CXX) $(CXXFLAGS) $(CXXOPT) -shared -o $(realpath $(INSTALL_DIR))/lib/$(LIBNAME).so $^ $(LD_FLAGS) 

#### installing dynamic library
$(INSTALL_DIR)/lib/$(LIBNAME).dylib: $(OBJNOMAIN)
ifeq ($(UNAME_S), Darwin)
	$(CXX) $(CXXFLAGS) $(CXXOPT) -dynamiclib -o $(realpath $(INSTALL_DIR))/lib/$(LIBNAME).dylib $^ $(LD_FLAGS)
endif
	  

$(INSTALL_DIR)/bin/$(CXXOUTNAME): $(OBJ)
	$(CXX) $(CXXFLAGS) $(CXXOPT) -o $(realpath $(INSTALL_DIR))/bin/$(CXXOUTNAME) $^ $(LD_FLAGS) 

### Run the move headers script and remove the previous includes directory to get rid 
### of any headers that were removed
cpHeaders: $(INSTALL_DIR)
	scripts/setUpScripts/installHeaders.py -src src/ -dest $(INSTALL_DIR)/include/ -rmDir
	
### Copy etc folder if it exists 
cpEtc: $(INSTALL_DIR)
ifneq ("$(wildcard $(PATH_TO_FILE))","")
	scripts/setUpScripts/installEtc.py -etcFolder etc/ -dest $(INSTALL_DIR)/ -rmDir
endif

### Run unit tests if available
unitTest: 
	scripts/setUpScripts/runUnitTest.sh
	

