UNAME_S := $(shell uname -s)
LOCAL_PATH = $(EXT_PATH)/local
#LD_FLAGS += 
#defaults for most progjects
LOCALTOOLS = -I$(LOCAL_PATH)
EXTTOOLS = -I$(EXT_PATH)
SRC = -I./src/
COMLIBS += $(LOCALTOOLS) $(EXTTOOLS) $(SRC)

#bibcpp
ifeq ($(USE_BIBCPP),1)
	COMLIBS += -isystem$(LOCAL_PATH)/bibcpp/include
	USE_JSONCPP=1
	USE_BOOST=1
	
	#currently no compiled components so no need for library flags
	#uncomment bellow in the future if there parts of the package need to be compiled
	#LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/bibcpp/lib \
			-L$(LOCAL_PATH)/bibcpp/lib  \
			-lbibcpp
endif

#SeekDeep
ifeq ($(USE_SEEKDEEP),1)
	COMLIBS += -isystem$(LOCAL_PATH)/SeekDeep/include
	ifeq (,$(findstring $(COMLIBS),bibcpp))
		COMLIBS += -isystem$(LOCAL_PATH)/bibcpp/include
		USE_JSONCPP=1
		USE_BOOST=1
	endif
	USE_BAMTOOLS=1
	USE_R=1
	USE_CURL=1
	LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/SeekDeep/lib \
			-L$(LOCAL_PATH)/SeekDeep/lib  \
			-lSeekDeep
endif


#bibseq
ifeq ($(USE_BIBSEQ),1)
	COMLIBS += -isystem$(LOCAL_PATH)/bibseq/include
	ifeq (,$(findstring $(COMLIBS),bibcpp))
		COMLIBS += -isystem$(LOCAL_PATH)/bibcpp/include
		USE_JSONCPP=1
		USE_BOOST=1
	endif
	USE_ARMADILLO=1
	USE_BAMTOOLS=1
	USE_R=1
	USE_CURL=1
	LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/bibseq/lib \
			-L$(LOCAL_PATH)/bibseq/lib  \
			-lbibseq
endif



#jsoncpp
ifeq ($(USE_JSONCPP),1)
	COMLIBS += -isystem$(LOCAL_PATH)/jsoncpp/include
	LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/jsoncpp/lib \
			-L$(LOCAL_PATH)/jsoncpp/lib  \
			-ljsoncpp
endif


#boost
ifeq ($(USE_BOOST),1)
	CXXOPT += -DBOOST_UBLAS_NDEBUG
	COMLIBS += -isystem$(LOCAL_PATH)/boost/include
	LD_FLAGS +=  -Wl,-rpath,$(LOCAL_PATH)/boost/lib \
			-L$(LOCAL_PATH)/boost/lib  \
			-lpthread -lboost_program_options -lboost_system -lboost_thread \
			-lboost_filesystem -lboost_iostreams -lboost_regex -lboost_serialization
endif

#cppcms
ifeq ($(USE_CPPCMS),1)
	COMLIBS += -isystem$(LOCAL_PATH)/cppcms/include
	LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/cppcms/lib \
			-L$(LOCAL_PATH)/cppcms/lib  \
			-lcppcms -lbooster
endif

#armadillo
ifeq ($(USE_ARMADILLO),1)
	COMLIBS += -isystem$(LOCAL_PATH)/armadillo/include
	LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/armadillo/lib \
			-L$(LOCAL_PATH)/armadillo/lib  \
			-larmadillo
endif





#shark
ifeq ($(USE_SHARK),1)
	COMLIBS += -isystem$(LOCAL_PATH)/shark/include
	LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/shark/lib \
		-L$(LOCAL_PATH)/shark/lib  \
		-lshark
endif

#CPPPROGUTILS
ifeq ($(USE_CPPPROGUTILS),1)
	COMLIBS += -I$(LOCAL_PATH)/cppprogutils
endif

#CATCH
ifeq ($(USE_CATCH),1)
	COMLIBS += -I$(LOCAL_PATH)/catch/single_include
endif

#ZI_LIB
ifeq ($(USE_ZI_LIB),1)
	COMLIBS += -I$(LOCAL_PATH)/zi_lib
	#CXXFLAGS += -DZI_USE_OPENMP
	ifeq ($(UNAME_S),Darwin)

	else
		CXXFLAGS += -DZI_USE_OPENMP
    	LD_FLAGS += -lrt
	endif
endif


#bamtools
ifeq ($(USE_BAMTOOLS),1)
	COMLIBS += -I$(LOCAL_PATH)/bamtools/include/bamtools
	LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/bamtools/lib/bamtools \
			-L$(LOCAL_PATH)/bamtools/lib/bamtools\
			-lbamtools
endif
#c url library 
ifeq ($(USE_CURL),1)
	LD_FLAGS += -lcurl
endif

#gtkmm library, should have 3.0 install and .pc file should be in PKG_CONFIG_PATH  
ifeq ($(USE_GTKMM),1)
	LD_FLAGS += `pkg-config gtkmm-3.0 --libs`
	COMLIBS += `pkg-config gtkmm-3.0 --cflags`
endif

#ml_pack
ifeq ($(USE_MLPACK),1)
	ifeq ($(UNAME_S),Darwin)
    	LD_FLAGS += -llapack  -lcblas # non-threaded
	else
   		LD_FLAGS += -llapack -lf77blas -lcblas -latlas # non-threaded
	endif
endif

#qt5
ifeq ($(USE_QT5),1)
	ifeq ($(UNAME_S),Darwin)
		LD_FLAGS += -Wl,-rpath,/usr/local/opt/qt5/lib \
	 				-L/usr/local/opt/qt5/lib \
	 				-lQt5UiTools
    	COMLIBS += -I/usr/local/opt/qt5/include
	endif
endif
ifeq ($(USE_R),1)
	include $(ROOT)/r-makefile-common.mk
endif

ifeq ($(UNAME_S),Darwin)
    #for dylib path fixing in macs, this gets rid of the name_size limit
    LD_FLAGS += -headerpad_max_install_names
endif



# from http://stackoverflow.com/a/18258352
rwildcard=$(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2) $(filter $(subst *,%,$2),$d))
