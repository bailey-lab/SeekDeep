CC = gcc-4.9
CXX = g++-4.8
CXXOUTNAME = out
CXXFLAGS = -std=c++1y -Wall
CXXOPT += -O2 -funroll-loops -DNDEBUG  
ifneq ($(shell uname -s),Darwin)
	CXXOPT += -march=native -mtune=native
endif

#debug
CXXDEBUG = -g -gstabs+ 
INSTALL_DIR=./out
EXT_PATH=$(realpath external)
SCRIPTS_DIR=$(realpath scripts)

USE_CPPITERTOOLS = 0
USE_CPPPROGUTILS = 0
USE_ZI_LIB = 0
USE_BOOST = 0
USE_R = 0
USE_BAMTOOLS = 0
USE_CPPCMS = 0
USE_MATHGL = 0
USE_ARMADILLO = 0
USE_MLPACK = 0
USE_LIBLINEAR = 0
USE_PEAR = 0
USE_CURL = 0
USE_GTKMM = 0
USE_BIBCPP = 0
USE_CATCH = 0
