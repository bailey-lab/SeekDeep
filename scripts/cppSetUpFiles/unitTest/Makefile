ROOT = $(realpath ../)
ifdef CXXFLAGS 
	ENV_CXXFLAGS := $(CXXFLAGS)
endif
### compfile should be be in the make command e.g. make COMPFILE=../compfile.mk
### or it is assumed to be 
ifneq (,$(wildcard ../compfile.mk))
COMPFILE=../compfile.mk
endif
include $(COMPFILE)
#make sure catch is included
USE_CATCH=1
include $(ROOT)/makefile-common.mk

#### File targets
BIN = bin/tester
OBJ_DIR = build
TESTSRC = $(realpath ./src)
SRCSRC = $(realpath ../src)
OBJTEST = $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(call rwildcard, $(TESTSRC), *.cpp)))
OBJSRC = $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(call rwildcard, $(SRCSRC), *.cpp)))
OBJSRCNOMAIN = $(filter-out $(addsuffix $(SRCSRC)/main.o, $(OBJ_DIR)/), $(OBJSRC))
#OBJSRCNOMAIN = $(filter-out build//Users/nick/hathaway/external/build/TwoBit/src/main.o, $(OBJSRC))
#
OBJALL = $(OBJTEST) $(OBJSRCNOMAIN)

HEADERS = $(call rwildcard, $(SRCSRC), *.h) \
			$(call rwildcard, $(SRCSRC), *.hpp)
#### Phony targets
.PHONY: all
.PHONY: clean
.PHONY: do_preReqs 


#compiler options
CXXFLAGS += $(ENV_CXXFLAGS)
COMMON = -I$(SRCSRC) $(CXXFLAGS) $(CXXOPT) $(COMLIBS)


############ main

all: do_preReqs $(OBJ_DIR) $(BIN)
ifeq ($(UNAME_S), Darwin)
	../scripts/setUpScripts/fixDyLinking_mac.sh bin $(EXT_PATH)
endif


$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)
	mkdir -p bin

# using automatic variables $<: the name of the prerequisite of the rule and
#							$@: the name of the target of the rule 
$(OBJ_DIR)/%.o: %.cpp $(HEADERS)
	@mkdir -p $(OBJ_DIR)/$(shell dirname $<)
	$(CXX) -DNOT_HEADER_ONLY $(COMMON) -fpic -c $< -o $@


############ remove the objects that were dependant on any changed headers and check for compfile.mk
do_preReqs: 
ifndef COMPFILE
	$(error COMPFILE is not set, do either make COMPFILE=aCompfile.mk or create a file called compfile.mk)
endif
#../scripts/setUpScripts/rmNeedToRecompile.py -obj $(OBJ_DIR) -src $(SRCSRC)
#../scripts/setUpScripts/rmNeedToRecompile.py -obj $(OBJ_DIR) -src $(TESTSRC)

	
$(BIN): $(OBJALL)
	$(CXX) $(CXXFLAGS) $(CXXOPT) -o $@ $^ $(LD_FLAGS)


clean:
	@rm -f $(BIN)
	@rm -rf $(OBJ_DIR)
