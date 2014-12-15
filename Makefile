ROOT = $(realpath ./)
LIB_DIR=$(ROOT)/lib

ifneq (,$(wildcard compfile.mk))
COMPFILE=compfile.mk
endif

include $(COMPFILE)

include $(ROOT)/makefile-common.mk


ifndef COMPFILE
	$(error $(COMPFILE) is not set)
endif

UNAME_S := $(shell uname -s)
# from http://stackoverflow.com/a/8654800
HEADERS = $(call rwildcard, src/, *.h) \
	$(call rwildcard, src/, *.hpp)

OBJ_DIR = $(addprefix build/, $(addsuffix Build, $(CXXOUTNAME)))
OBJ = $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(call rwildcard, src/, *.cpp))) 
OBJNOMAIN = $(filter-out $(addsuffix /src/main.o, $(OBJ_DIR)), $(OBJ))



BIN = $(addsuffix $(CXXOUTNAME), bin/)
LIBNAME = $(addsuffix $(CXXOUTNAME), lib)
DYLIB = $(addprefix $(addsuffix $(LIBNAME), $(LIB_DIR)/), .dylib)
SOLIB = $(addprefix $(addsuffix $(LIBNAME), $(LIB_DIR)/), .so)
COMMON = $(CXXFLAGS) $(CXXOPT) $(COMLIBS)




-include do_preReqs
############ main
.PHONY: all
all: do_preReqs $(OBJ_DIR) $(BIN) 
	setUpScripts/fixDyLinking_mac.sh bin $(EXT_PATH)
	
.PHONY: dev
dev: do_preReqs $(OBJ_DIR) $(BIN) $(PROTO) $(BIOALG) $(EULER) $(QTTEST) 
	setUpScripts/fixDyLinking_mac.sh bin $(EXT_PATH)	

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)
	mkdir -p bin
	mkdir -p lib

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(OBJ_DIR)/$(shell dirname $<)
	$(CXX) -DNOT_HEADER_ONLY $(COMMON) -fPIC -c $< -o $@

$(BIN): $(OBJ) 
	$(CXX) $(CXXFLAGS) $(CXXOPT) -o $@ $^ $(LD_FLAGS) 
	
	
############ remove the objects that were dependant the changed headers 
.PHONY: do_script
do_script:
	setUpScripts/rmNeedToRecompile.py -obj $(OBJ_DIR) -src src/

prerequisites: do_script

.PHONY: do_preReqs 
do_preReqs: prerequisites

	
############ shared library
.PHONY: sharedLibrary
sharedLibrary: $(OBJ_DIR) $(SOLIB)
ifeq ($(UNAME_S), Darwin)
	setUpScripts/fixDyLinking_mac.sh lib $(EXT_PATH)
endif

$(SOLIB): $(OBJNOMAIN)
	$(CXX) $(CXXFLAGS) $(CXXOPT) -shared -o $@ $^ $(LD_FLAGS) 
	
############ dylibLibrary
.PHONY: dylibLibrary
dylibLibrary: $(OBJ_DIR) $(DYLIB)
	setUpScripts/fixDyLinking_mac.sh lib $(EXT_PATH)
$(DYLIB): $(OBJNOMAIN)
	$(CXX) $(CXXFLAGS) $(CXXOPT) -dynamiclib -o $@ $^ $(LD_FLAGS) 


############ clean
.PHONY: clean
clean:
	@rm -f $(BIN)
	@rm -rf $(OBJ_DIR)
	
	
############ install
.PHONY: install
install: $(INSTALL_DIR) moveHeaders $(OBJ_DIR) $(INSTALL_DIR)/lib/$(LIBNAME).so $(INSTALL_DIR)/lib/$(LIBNAME).dylib $(INSTALL_DIR)/bin/$(CXXOUTNAME) 
ifeq ($(UNAME_S), Darwin)
	setUpScripts/fixDyLinking_mac.sh $(INSTALL_DIR)/lib/ $(EXT_PATH)
	setUpScripts/fixDyLinking_mac.sh $(INSTALL_DIR)/bin/ $(EXT_PATH)
endif

$(INSTALL_DIR):
	@mkdir -p $(INSTALL_DIR)
	@mkdir -p $(INSTALL_DIR)/include
	@mkdir -p $(INSTALL_DIR)/bin
	@mkdir -p $(INSTALL_DIR)/lib
	
$(INSTALL_DIR)/lib/$(LIBNAME).so: $(OBJNOMAIN)
	$(CXX) $(CXXFLAGS) $(CXXOPT) -shared -o $(realpath $(INSTALL_DIR))/lib/$(LIBNAME).so $^ $(LD_FLAGS) 

$(INSTALL_DIR)/lib/$(LIBNAME).dylib: $(OBJNOMAIN)
ifeq ($(UNAME_S), Darwin)
	$(CXX) $(CXXFLAGS) $(CXXOPT) -dynamiclib -o $(realpath $(INSTALL_DIR))/lib/$(LIBNAME).dylib $^ $(LD_FLAGS)
endif
	  

$(INSTALL_DIR)/bin/$(CXXOUTNAME):$(OBJ)
	$(CXX) $(CXXFLAGS) $(CXXOPT) -o $(realpath $(INSTALL_DIR))/bin/$(CXXOUTNAME) $^ $(LD_FLAGS) 


.PHONY moveHeaders:
moveHeaders: $(INSTALL_DIR)
	setUpScripts/installHeaders.py -src src/ -dest $(INSTALL_DIR)/include/
	
.PHONY unitTest:
unitTest: 
	setUpScripts/runUnitTest.sh $(COMPFILE) $(NCORETEST)
	
	

	