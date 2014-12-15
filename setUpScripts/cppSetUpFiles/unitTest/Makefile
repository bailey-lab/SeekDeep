ROOT = $(realpath ../)


include $(COMPFILE)

include $(ROOT)/makefile-common.mk

SRC += -I../src/

BIN = bin/tester

OBJ_DIR = build
TESTSRC = $(realpath ./src)
SRCSRC = $(realpath ../src)
OBJ = $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(call rwildcard, $(TESTSRC), *.cpp))) \
	$(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(call rwildcard, $(SRCSRC)/programs, *.cpp))) \
	$(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(call rwildcard, $(SRCSRC)/seqTools, *.cpp)))
HEADERS = $(call rwildcard, src/, *.h) \
	$(call rwildcard, src/, *.hpp)

COMMON =  $(CXXFLAGS) $(CXXOPT) $(COMLIBS)
#COMMON = $(COMMON_DEBUG)

############ main
.PHONY: all
all: $(OBJ_DIR) $(BIN)
ifeq ($(UNAME_S), Darwin)
	../setUpScripts/fixDyLinking_mac.sh bin $(EXT_PATH)
endif

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)
	mkdir -p bin

$(OBJ_DIR)/%.o: %.cpp $(HEADERS)
	@mkdir -p $(OBJ_DIR)/$(shell dirname $<)
	$(CXX) -DNOT_HEADER_ONLY $(COMMON) -fpic -c $< -o $@


$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) $(CXXOPT) -o $@ $^ $(LD_FLAGS)



.PHONY: clean
clean:
	@rm -f $(BIN)
	@rm -rf $(OBJ_DIR)
