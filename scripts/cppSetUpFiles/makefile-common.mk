UNAME_S := $(shell uname -s)
LOCAL_PATH = $(EXT_PATH)/local
#LD_FLAGS += 
#defaults for most progjects
LOCALTOOLS = -I$(LOCAL_PATH)
EXTTOOLS = -I$(EXT_PATH)
SRC = -I./src/
COMLIBS += $(LOCALTOOLS) $(EXTTOOLS) $(SRC)



#dlib
ifeq ($(USE_DLIB),1)
	COMLIBS += -isystem$(LOCAL_PATH)/dlib/
	#LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/dlib/lib \
	#		-L$(LOCAL_PATH)/dlib/lib  \
	#		-ldlib
endif

#libsvm
ifeq ($(USE_LIBSVM),1)
	COMLIBS += -isystem$(LOCAL_PATH)/libsvm/
	#LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/dlib/lib \
	#		-L$(LOCAL_PATH)/dlib/lib  \
	#		-ldlib
endif

#TwoBit
ifeq ($(USE_TWOBIT),1)
	COMLIBS += -isystem$(LOCAL_PATH)/TwoBit/include
	LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/TwoBit/lib \
			-L$(LOCAL_PATH)/TwoBit/lib  \
			-lTwoBit
endif

#TwoBit
ifeq ($(USE_TWOBIT),1)
	COMLIBS += -isystem$(LOCAL_PATH)/TwoBit/include
	LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/TwoBit/lib \
			-L$(LOCAL_PATH)/TwoBit/lib  \
			-lTwoBit
endif




#SeekDeep
ifeq ($(USE_SEEKDEEP),1)
	COMLIBS += -isystem$(LOCAL_PATH)/SeekDeep/include
	USE_SEQSERVER=1
	USE_NJHRINSIDE=1
	USE_BIBSEQ=1
	LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/SeekDeep/lib \
			-L$(LOCAL_PATH)/SeekDeep/lib  \
			-lSeekDeep
endif

#SeekDeepDev
ifeq ($(USE_SEEKDEEPDEV),1)
	COMLIBS += -isystem$(LOCAL_PATH)/SeekDeepDev/include
	USE_SEQSERVER=1
	USE_NJHRINSIDE=1
	USE_BIBSEQDEV=1
	LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/SeekDeepDev/lib \
			-L$(LOCAL_PATH)/SeekDeepDev/lib  \
			-lSeekDeepDev
endif

#SeqServer
ifeq ($(USE_SEQSERVER),1)
	COMLIBS += -isystem$(LOCAL_PATH)/seqServer/include
	ifeq ($(USE_BIBSEQDEV),1)
		USE_BIBSEQDEV=1
		USE_BIBSEQ=0
	else
		USE_BIBSEQ=1
	endif
	USE_CPPCMS=1
	LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/seqServer/lib \
			-L$(LOCAL_PATH)/seqServer/lib  \
			-lseqServer
endif

#njhRInside
ifeq ($(USE_NJHRINSIDE),1)
	COMLIBS += -isystem$(LOCAL_PATH)/njhRInside/include
	USE_R=1
	USE_CPPITERTOOLS=1
	LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/njhRInside/lib \
			-L$(LOCAL_PATH)/njhRInside/lib  \
			-lnjhRInside
endif

#bibseq
ifeq ($(USE_BIBSEQ),1)
	COMLIBS += -isystem$(LOCAL_PATH)/bibseq/include
	USE_BIBCPP=1
	USE_ARMADILLO=1
	USE_BAMTOOLS=1
	USE_CURL=1
	LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/bibseq/lib \
			-L$(LOCAL_PATH)/bibseq/lib  \
			-lbibseq
endif

#bibseqDev
ifeq ($(USE_BIBSEQDEV),1)
	COMLIBS += -isystem$(LOCAL_PATH)/bibseqDev/include
	USE_BIBCPPDEV=1
	USE_ARMADILLO=1
	USE_BAMTOOLS=1
	USE_BIBSEQ=0
	#USE_R=1
	USE_CURL=1
	LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/bibseqDev/lib \
			-L$(LOCAL_PATH)/bibseqDev/lib  \
			-lbibseqDev
endif


#bibcpp
ifeq ($(USE_BIBCPP),1)
	COMLIBS += -isystem$(LOCAL_PATH)/bibcpp/include
	USE_JSONCPP=1
	USE_BOOST=1
	LD_FLAGS += -lpthread
	#currently no compiled components so no need for library flags
	#uncomment below in the future if there parts of the package need to be compiled
	#LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/bibcpp/lib \
			-L$(LOCAL_PATH)/bibcpp/lib  \
			-lbibcpp
endif

#bibcppDev
ifeq ($(USE_BIBCPPDEV),1)
	COMLIBS += -isystem$(LOCAL_PATH)/bibcppDev/include
	USE_JSONCPP=1
	USE_BOOST=1
	LD_FLAGS += -lpthread
	#currently no compiled components so no need for library flags
	#uncomment below in the future if there parts of the package need to be compiled
	#LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/bibcppDev/lib \
			-L$(LOCAL_PATH)/bibcppDev/lib  \
			-lbibcppDev
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
			-lboost_system -lboost_filesystem -lboost_iostreams
			#-lpthread -lboost_program_options -lboost_system -lboost_thread \
			#-lboost_filesystem -lboost_iostreams -lboost_regex -lboost_serialization
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
	ifeq ($(UNAME_S),Darwin)

	else
		#CXXFLAGS += -DZI_USE_OPENMP
    	LD_FLAGS += -lrt
	endif
endif


#bamtools
ifeq ($(USE_BAMTOOLS),1)
	COMLIBS += -isystem$(LOCAL_PATH)/bamtools/include/bamtools
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


#mongocxx  
ifeq ($(USE_MONGOCXX),1)
	COMLIBS += -isystem$(LOCAL_PATH)/mongocxx/include/mongocxx/v0.3 \
	-isystem$(LOCAL_PATH)/mongocxx/include/bsoncxx/v0.3 
	LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/mongocxx/lib \
			-L$(LOCAL_PATH)/mongocxx/lib  \
			-lmongocxx -lbsoncxx 
	USE_MONGOC=1
endif



#mongoc  
ifeq ($(USE_MONGOC),1)
	COMLIBS += -isystem$(LOCAL_PATH)/mongoc/include/libbson-1.0 \
	-isystem$(LOCAL_PATH)/mongoc/include/libmongoc-1.0
	LD_FLAGS += -Wl,-rpath,$(LOCAL_PATH)/mongoc/lib \
			-L$(LOCAL_PATH)/mongoc/lib \
			-lssl -lcrypto -lmongoc-1.0 -lbson-1.0
	ifeq ($(UNAME_S),Darwin)

	else
   		LD_FLAGS += -lrt
	endif
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
    #for dylib path fixing in macs, this gets rid of the name_size limit, which why the hell is there a name size limit
    LD_FLAGS += -headerpad_max_install_names
endif



# from http://stackoverflow.com/a/18258352
rwildcard=$(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2) $(filter $(subst *,%,$2),$d))
