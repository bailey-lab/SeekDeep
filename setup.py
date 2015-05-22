#!/usr/bin/env python

# by purcaro@gmail.com

import subprocess, sys, os, argparse
from setUpScripts.utils import Utils
from setUpScripts.genFuncs import genHelper 
from collections import namedtuple
from setUpScripts.color_text import ColorText as CT


BuildPaths = namedtuple("BuildPaths", 'url build_dir build_sub_dir local_dir')

def shellquote(s):
    # from http://stackoverflow.com/a/35857
    return "'" + s.replace("'", "'\\''") + "'"

def isMac():
    return sys.platform == "darwin"

class Paths():
    '''class to hold and setup all the necessary paths for 
    downloading, building, and then installing packages/libraries'''
    def __init__(self, externalLoc):
        self.base_dir = externalLoc; #top dir to hold tars,build, local directories
        
        self.ext_tars = os.path.join(self.base_dir, "tarballs") #location to keep tarballs of programs/libraries downloads
        self.ext_build = os.path.join(self.base_dir, "build") #location for the building of programs/libraries
        self.install_dir = os.path.join(self.base_dir, "local") #location for the final install of programs/libraries
        
        Utils.mkdir(self.ext_tars) #tar storage directory
        Utils.mkdir(self.ext_build) #build directory
        self.paths = {} #dictionary to hold path infos
        self.paths["zi_lib"] = self.__zi_lib()
        self.paths["pstreams"] = self.__pstreams()
        self.paths["cppitertools"] = self.__cppitertools()
        self.paths["cppprogutils"] = self.__cppprogutils()
        self.paths["boost"] = self.__boost()
        self.paths["r"] = self.__r()
        self.paths["cppcms"] = self.__cppcms()
        self.paths["bamtools"] = self.__bamtools()
        self.paths["jsoncpp"] = self.__jsoncpp()
        self.paths["pear"] = self.__pear()
        self.paths["mathgl"] = self.__mathgl()
        self.paths["armadillo"] = self.__armadillo()
        self.paths["mlpack"] = self.__mlpack()
        self.paths["liblinear"] = self.__liblinear()
        self.paths["bibseq"] = self.__bibseq()
        self.paths["bibcpp"] = self.__bibcpp()
        self.paths["seekdeep"] = self.__SeekDeep()
        self.paths["bibseqdev"] = self.__bibseqDev()
        self.paths["bibcppdev"] = self.__bibcppDev()
        self.paths["seekdeepdev"] = self.__SeekDeepDev()
        self.paths["seqserver"] = self.__seqserver()
        self.paths["njhrinside"] = self.__njhRInside()
        self.paths["catch"] = self.__catch()

    def path(self, name):
        '''get path info if it exists'''
        if name in self.paths:
            return self.paths[name]
        raise Exception(name + " not found in paths")

    def __zi_lib(self):
        url = 'https://github.com/weng-lab/zi_lib.git'
        local_dir = os.path.join(self.install_dir, "zi_lib")
        return BuildPaths(url, '', '', local_dir)
    
    def __pstreams(self):
        url = 'http://git.code.sf.net/p/pstreams/code'
        local_dir = os.path.join(self.install_dir, "pstreams")
        return BuildPaths(url, '', '', local_dir)

    def __bamtools(self):
        url = 'https://github.com/pezmaster31/bamtools.git'
        name = "bamtools"
        return self.__package_dirs(url, name)
    
    def __jsoncpp(self):
        url = "https://github.com/open-source-parsers/jsoncpp.git"
        name = "jsoncpp"
        return self.__package_dirs(url, name)

    def __cppitertools(self):
        url = 'https://github.com/ryanhaining/cppitertools.git'
        name = "cppitertools"
        return self.__headerOnly_dirs(url,name)

    def __catch(self):
        url = 'https://github.com/philsquared/Catch.git'
        name = "catch"
        return self.__headerOnly_dirs(url,name)

    def __cppprogutils(self):
        url = 'https://github.com/bailey-lab/cppprogutils.git'
        name = "cppprogutils"
        return self.__headerOnly_dirs(url,name)

    def __pear(self):
        url = "http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.4-src.tar.gz"
        return self.__package_dirs(url, "pear")

    def __r(self):
        #url = "ftp://ftp.stat.math.ethz.ch/Software/R/R-devel.tar.gz"
        #url = "http://cran.r-project.org/src/base/R-3/R-3.2.0.tar.gz"
        #url = "http://cran.r-project.org/src/base/R-3/R-3.1.0.tar.gz"
        #url = "http://cran.r-project.org/src/base/R-3/R-3.1.1.tar.gz"
        #url = "http://cran.r-project.org/src/base/R-3/R-3.1.2.tar.gz"
        url = "http://cran.r-project.org/src/base/R-3/R-3.1.3.tar.gz"
        return self.__package_dirs(url, "R")

    def __boost(self):
        url = "http://sourceforge.net/projects/boost/files/boost/1.57.0/boost_1_57_0.tar.gz"
        return self.__package_dirs(url, "boost")

    def __armadillo(self):
        url = "http://freefr.dl.sourceforge.net/project/arma/armadillo-4.000.2.tar.gz"
        return self.__package_dirs(url, "armadillo")

    def __mlpack(self):
        url = "http://www.mlpack.org/files/mlpack-1.0.8.tar.gz"
        return self.__package_dirs(url, "mlpack")

    def __liblinear(self):
        url = "http://www.csie.ntu.edu.tw/~cjlin/liblinear/liblinear-1.94.tar.gz"
        return self.__package_dirs(url, "liblinear")

    def __cppcms(self):
        url = "http://freefr.dl.sourceforge.net/project/cppcms/cppcms/1.0.5/cppcms-1.0.5.tar.bz2"
        return self.__package_dirs(url, "cppcms")

    def __mathgl(self):
        url = "http://freefr.dl.sourceforge.net/project/mathgl/mathgl/mathgl%202.2.1/mathgl-2.2.1.tar.gz"
        return self.__package_dirs(url, "mathgl")

    def __bibseq(self):
        url = "https://github.com/bailey-lab/bibseq.git"
        name = "bibseq"
        return self.__package_dirs(url, name)
    
    def __bibseqDev(self):
        url = "https://github.com/bailey-lab/bibseqPrivate.git"
        name = "bibseqDev"
        return self.__package_dirs(url, name)  
      
    def __SeekDeep(self):
        url = "https://github.com/bailey-lab/SeekDeep.git"
        name = "SeekDeep"
        return self.__package_dirs(url, name)
    
    def __SeekDeepDev(self):
        url = "https://github.com/bailey-lab/SeekDeepPrivate.git"
        name = "SeekDeepDev"
        return self.__package_dirs(url, name)
    
    def __seqserver(self):
        url = "https://github.com/nickjhathaway/seqServer.git"
        name = "seqserver"
        return self.__package_dirs(url, name)
    
    def __njhRInside(self):
        url = "https://github.com/nickjhathaway/njhRInside.git"
        name = "njhRInside"
        return self.__package_dirs(url, name)

    def __bibcppDev(self):
        url = "https://github.com/umass-bib/bibcppDev.git"
        name = "bibcppDev"
        return self.__package_dirs(url, name)
    
    def __bibcpp(self):
        url = "https://github.com/umass-bib/bibcpp.git"
        name = "bibcpp"
        return self.__package_dirs(url, name)

    def __package_dirs(self, url, name):
        '''set up dirs for libraries that need some installing as well, 
        will only work for compresed files like .tar.gz or .tar.bz2 and for git repo 
        currently'''
        build_dir = os.path.join(self.ext_build, name)
        fn = os.path.basename(url)
        fn_noex = fn.replace(".tar.gz", "").replace(".tar.bz2", "").replace(".git", "")
        build_sub_dir = os.path.join(build_dir, fn_noex)
        local_dir = os.path.join(self.install_dir, name)
        return BuildPaths(url, build_dir, build_sub_dir, local_dir)
    
    def __headerOnly_dirs(self, url, name):
        '''set up for header only libraries, these just need
         the header copied no need for build_dir build_sub_dir '''
        local_dir = os.path.join(self.install_dir, name)
        return BuildPaths(url, "", "", local_dir)

class Setup:
    def __init__(self, args):
        self.extDirLoc = ""
        #if no compile file set up and assume external is next to setup.py
        if not args.compfile:
            self.extDirLoc = os.path.abspath(os.path.join(os.path.dirname(__file__), "external"))
        else:
            self.extDirLoc = os.path.abspath(self.parseForExtPath(args.compfile[0]))
        self.paths = Paths(self.extDirLoc)
        self.args = args
        self.setUps = {}
        self.setUpsNeeded = []
        self.installed = []
        self.failedInstall = []
        self.CC = ""
        self.CXX = ""
        self.externalLoc = ""
        self.__initSetUpFuncs()
        self.__processArgs()
        
    def setup(self):
        # if no compfile need to determine compiler 
        if not self.args.compfile:
            self.CC = genHelper.determineCC(self.args)
            self.CXX = genHelper.determineCXX(self.args)

        for set in self.setUpsNeeded:
            if not set in self.setUps.keys():
                print CT.boldBlack( "Unrecognized option ") + CT.boldRed(set)
            else:
                self.__setup(set, self.setUps[set])

        for p in self.installed:
            print p, CT.boldGreen("installed")

        for p in self.failedInstall:
            print  p, CT.boldRed("failed to install")

    def __initSetUpFuncs(self):
        self.setUps = {"zi_lib": self.zi_lib,
                       "boost": self.boost,
                       "cppitertools": self.cppitertools,
                       "catch": self.catch,
                       "cppprogutils": self.cppprogutils,
                       "r": self.r,
                       "bamtools": self.bamtools,
                       "cppcms": self.cppcms,
                       "mathgl": self.mathgl,
                       "armadillo": self.armadillo,
                       "mlpack": self.mlpack,
                       "liblinear": self.liblinear,
                       "pear": self.pear,
                       "bibseq": self.bibseq,
                       "seekdeep": self.SeekDeep,
                       "bibcpp": self.bibcpp,
                       "bibseqdev": self.bibseqDev,
                       "seekdeepdev": self.SeekDeepDev,
                       "bibcppdev": self.bibcppDev,
                       "seqserver": self.seqserver,
                       "njhrinside": self.njhRInside,
                       "jsoncpp": self.jsoncpp,
                       "pstreams": self.pstreams
                       }
    def printAvailableSetUps(self):
        self.__initSetUpFuncs()
        print "Available installs:"
        installs = self.setUps.keys()
        installs.sort()
        for set in installs:
            print set

    def __processArgs(self):
        if self.args.libs:
            inLibs = self.args.libs.split(",")
            for lib in inLibs:
                self.setUpsNeeded.append(lib.lower())
        if self.args.compfile:
            self.parseSetUpNeeded(self.args.compfile[0])
            self.parserForCompilers(self.args.compfile[0])

    def parseForExtPath(self, fn):
        args = self.parseCompFile(fn)
        if "EXT_PATH" in args:
            extPath = args["EXT_PATH"].strip()
            extPath = extPath.replace("$(realpath", "")
            extPath = extPath.replace(")", "")
            extPath = extPath.strip()
        else:
            print "did not find external folder location; assuming ./external"
            extPath = "./external"
        return extPath

    def parseSetUpNeeded(self, fn):
        args = self.parseCompFile(fn)
        for k,v in args.iteritems():
            if k.startswith("USE_"):
                if '1' == v:
                    self.setUpsNeeded.append(k[4:].lower())

    def parseCompFile(self, fn):
        ret = {}
        with open(fn) as f:
            for line in f:
                if '=' in line:
                    toks = line.split('=')
                    ret[toks[0].strip()] = toks[1].strip()
        return ret

    def parserForCompilers(self, fn):
        args = self.parseCompFile(fn)
        if 'CC' in args:
            self.CC = args['CC']
        if 'CXX' in args:
            self.CXX = args['CXX']
            if "clang" in self.CXX:
                self.args.clang = True

    def __path(self, name):
        return self.paths.path(name)

    def __setup(self, name, builder_f):
        if os.path.exists(self.__path(name).local_dir):
            print name, CT.boldGreen("found at ") + CT.boldBlack(self.__path(name).local_dir)
        else:
            print name, CT.boldRed("NOT"), "found; building..."
            try:
                builder_f()
                self.installed.append(name)
            except:
                print "failed to install " + name
                self.failedInstall.append(name)

    def showDefaultExample(self):
        print """
Need to supply compfile to parse for needed libraries and compilers"
by giving -compfile"

example:

python ./setUpScripts/generateCompFile.py -outFilename compfile.mk \
-externalLoc external \
-CC gcc -CXX g++ \
-outname seqTools \
-installName bibseq \
-neededLibs zi_lib,cppitertools,cppprogutils,boost,R,bamtools,pear,curl

python ./setup.py -compfile compfile.mk

make COMPFILE=compfile.mk -j {num_cores}
"""



    def num_cores(self):
        c = Utils.num_cores()
        if c > 8:
            return c/2
        if 1 == c:
            return 1
        return c - 1

    def __build(self, i, cmd):
        print "\t Getting file..."
        fnp = Utils.get_file_if_size_diff(i.url, self.paths.ext_tars)
        Utils.clear_dir(i.build_dir)
        Utils.untar(fnp, i.build_dir)
        try:
            Utils.run_in_dir(cmd, i.build_sub_dir)
        except:
            Utils.rm_rf(i.local_dir)
            sys.exit(1)


    def __buildNjhProject(self,i):
        cmd = """
        python ./configure.py -CC {CC} -CXX {CXX} -externalLibDir {external} -prefix {localTop} 
        && python ./setup.py -compfile compfile.mk 
        && make -j {num_cores} && make install""".format(localTop=shellquote(self.paths.install_dir),
                                                          num_cores=self.num_cores(), CC=self.CC, CXX=self.CXX,
                                                           external=self.extDirLoc)
        cmd = " ".join(cmd.split())
        self.__buildFromGit(i, cmd)
    
    
    def __buildFromGit(self, i, cmd):
        if os.path.exists(i.build_dir):
            print "pulling from {url}".format(url=i.url)
            pCmd = "git pull"
            try:
                Utils.run_in_dir(pCmd, i.build_dir)
            except:
                print "failed to pull from {url}".format(url=i.url)
                sys.exit(1)
        else:
            print "cloning from {url}".format(url=i.url)
            cCmd = "git clone {url} {d}".format(url=i.url, d=i.build_dir)
            try:
                print self.paths.ext_build
                Utils.run_in_dir(cCmd, self.paths.ext_build)
            except:
                print "failed to clone from {url}".format(url=i.url)
                sys.exit(1)
        try:
            Utils.run_in_dir(cmd, i.build_dir)
        except:
            Utils.rm_rf(i.local_dir)
            sys.exit(1)
    
    def __git(self, i):
        cmd = "git clone {url} {d}".format(url=i.url, d=shellquote(i.local_dir))
        Utils.run(cmd)
    
    def installRPackageSource(self, sourceFile):
        i = self.__path("r")
        for pack in sourceFile.split(","):
            if isMac():
                cmd = """echo 'install.packages(\"{SOURCEFILE}\", repos = NULL, type="source")' | $({local_dir}/R.framework/Resources/bin/R RHOME)/bin/R --slave --vanilla
                """.format(local_dir=shellquote(i.local_dir).replace(' ', '\ '),SOURCEFILE = pack )
            else:
                cmd = """echo 'install.packages(\"{SOURCEFILE}\", repos = NULL, type="source")' | $({local_dir}/bin/R RHOME)/bin/R --slave --vanilla
                """.format(local_dir=shellquote(i.local_dir).replace(' ', '\ '),SOURCEFILE = pack )
            print CT.boldBlack(cmd)
            cmd = " ".join(cmd.split())
            Utils.run(cmd)

    def installRPackageName(self, packageName):
        i = self.__path("r")
        for pack in packageName.split(","):
            if isMac():
                cmd = """echo 'install.packages(c(\"{PACKAGENAME}\"), repos=\"http://cran.us.r-project.org\")' | $({local_dir}/R.framework/Resources/bin/R RHOME)/bin/R --slave --vanilla
                """.format(local_dir=shellquote(i.local_dir).replace(' ', '\ '), PACKAGENAME = pack )
            else:
                cmd = """echo 'install.packages(\"{PACKAGENAME}\", repos=\"http://cran.us.r-project.org\")'  | $({local_dir}/bin/R RHOME)/bin/R --slave --vanilla
                """.format(local_dir=shellquote(i.local_dir).replace(' ', '\ '),PACKAGENAME = pack )
            print CT.boldBlack(cmd)
            cmd = " ".join(cmd.split())
            Utils.run(cmd)

    def boost(self):
        i = self.__path("boost")
        if self.args.clang:
             if isMac():
                cmd = """./bootstrap.sh --with-toolset=clang --prefix={local_dir} --with-libraries=date_time,filesystem,iostreams,math,program_options,random,regex,serialization,signals,system,test,thread,log  &&  ./b2  -d 2 toolset=clang cxxflags=\"-stdlib=libc++ -std=c++14\" linkflags=\"-stdlib=libc++\" -j {num_cores} install &&  install_name_tool -change libboost_system.dylib {local_dir}/lib/libboost_system.dylib {local_dir}/lib/libboost_thread.dylib&&  install_name_tool -change libboost_system.dylib {local_dir}/lib/libboost_system.dylib {local_dir}/lib/libboost_filesystem.dylib""".format(local_dir=shellquote(i.local_dir).replace(' ', '\ '), num_cores=self.num_cores())
                #cmd = """wget https://github.com/boostorg/atomic/commit/6bb71fdd.diff && wget https://github.com/boostorg/atomic/commit/e4bde20f.diff&&  wget https://gist.githubusercontent.com/philacs/375303205d5f8918e700/raw/d6ded52c3a927b6558984d22efe0a5cf9e59cd8c/0005-Boost.S11n-include-missing-algorithm.patch&&  patch -p2 -i 6bb71fdd.diff&&  patch -p2 -i e4bde20f.diff&&  patch -p1 -i 0005-Boost.S11n-include-missing-algorithm.patch&&  echo "using clang;  " >> tools/build/v2/user-config.jam&&  ./bootstrap.sh --with-toolset=clang --prefix={local_dir} --with-libraries=date_time,filesystem,iostreams,math,program_options,random,regex,serialization,signals,system,test,thread,log  &&  ./b2  -d 2 toolset=clang cxxflags=\"-stdlib=libc++\" linkflags=\"-stdlib=libc++\" -j {num_cores} install &&  install_name_tool -change libboost_system.dylib {local_dir}/lib/libboost_system.dylib {local_dir}/lib/libboost_thread.dylib&&  install_name_tool -change libboost_system.dylib {local_dir}/lib/libboost_system.dylib {local_dir}/lib/libboost_filesystem.dylib""".format(local_dir=shellquote(i.local_dir).replace(' ', '\ '), num_cores=self.num_cores())
             else:
                cmd = """./bootstrap.sh --with-toolset=clang --prefix={local_dir}  --with-libraries=date_time,filesystem,iostreams,math,program_options,random,regex,serialization,signals,system,test,thread,log &&  ./b2  -d 2 toolset=clang cxxflags=\"-std=c++14\" -j {num_cores} install""".format(local_dir=shellquote(i.local_dir).replace(' ', '\ '), num_cores=self.num_cores())
                #cmd = """wget https://github.com/boostorg/atomic/commit/6bb71fdd.diff && wget https://github.com/boostorg/atomic/commit/e4bde20f.diff&&  wget https://gist.githubusercontent.com/philacs/375303205d5f8918e700/raw/d6ded52c3a927b6558984d22efe0a5cf9e59cd8c/0005-Boost.S11n-include-missing-algorithm.patch&&  patch -p2 -i 6bb71fdd.diff&&  patch -p2 -i e4bde20f.diff&&  patch -p1 -i 0005-Boost.S11n-include-missing-algorithm.patch&&  echo "using clang;  " >> tools/build/v2/user-config.jam&&  ./bootstrap.sh --with-toolset=clang --prefix={local_dir}  --with-libraries=date_time,filesystem,iostreams,math,program_options,random,regex,serialization,signals,system,test,thread,log &&  ./b2  -d 2 toolset=clang -j {num_cores} install""".format(local_dir=shellquote(i.local_dir).replace(' ', '\ '), num_cores=self.num_cores())
        elif self.CXX == "g++-4.8":
            if isMac():
                cmd = """echo "using gcc : 4.8 : g++-4.8 ; " >> tools/build/v2/user-config.jam && ./bootstrap.sh --prefix={local_dir} --with-libraries=date_time,filesystem,iostreams,math,program_options,random,regex,serialization,signals,system,test,thread,log && ./b2 -d 2 toolset=darwin-4.8 -j {num_cores} install && install_name_tool -change libboost_system.dylib {local_dir}/lib/libboost_system.dylib {local_dir}/lib/libboost_thread.dylib&&  install_name_tool -change libboost_system.dylib {local_dir}/lib/libboost_system.dylib {local_dir}/lib/libboost_filesystem.dylib""".format(local_dir=shellquote(i.local_dir).replace(' ', '\ '), num_cores=self.num_cores())
            else:
                cmd = """./bootstrap.sh --with-toolset=gcc-4.8 --prefix={local_dir} --with-libraries=date_time,filesystem,iostreams,math,program_options,random,regex,serialization,signals,system,test,thread,log  && ./b2 -d 2 toolset=gcc-4.8 -j {num_cores} install""".format(local_dir=shellquote(i.local_dir).replace(' ', '\ '), num_cores=self.num_cores())
        elif self.CXX == "g++-4.9":
            if isMac():
                cmd = """echo "using gcc : 4.9 : g++-4.9 ; " >> tools/build/v2/user-config.jam && ./bootstrap.sh --prefix={local_dir} --with-libraries=date_time,filesystem,iostreams,math,program_options,random,regex,serialization,signals,system,test,thread,log && ./b2 -d 2 toolset=darwin-4.9 -j {num_cores} install && install_name_tool -change libboost_system.dylib {local_dir}/lib/libboost_system.dylib {local_dir}/lib/libboost_thread.dylib&&  install_name_tool -change libboost_system.dylib {local_dir}/lib/libboost_system.dylib {local_dir}/lib/libboost_filesystem.dylib""".format(local_dir=shellquote(i.local_dir).replace(' ', '\ '), num_cores=self.num_cores())
            else:
                cmd = """./bootstrap.sh --with-toolset=gcc --prefix={local_dir} --with-libraries=date_time,filesystem,iostreams,math,program_options,random,regex,serialization,signals,system,test,thread,log  && ./b2 -d 2 -j {num_cores} install""".format(local_dir=shellquote(i.local_dir).replace(' ', '\ '), num_cores=self.num_cores())
        elif self.CXX == "g++":
            if isMac():
                cmd = """echo "using gcc : 4.9 : g++ ; " >> tools/build/v2/user-config.jam && ./bootstrap.sh --prefix={local_dir} --with-libraries=date_time,filesystem,iostreams,math,program_options,random,regex,serialization,signals,system,test,thread,log && ./b2 -d 2 toolset=darwin-4.9 -j {num_cores} install && install_name_tool -change libboost_system.dylib {local_dir}/lib/libboost_system.dylib {local_dir}/lib/libboost_thread.dylib&&  install_name_tool -change libboost_system.dylib {local_dir}/lib/libboost_system.dylib {local_dir}/lib/libboost_filesystem.dylib""".format(local_dir=shellquote(i.local_dir).replace(' ', '\ '), num_cores=self.num_cores())
            else:
                cmd = """./bootstrap.sh --with-toolset=gcc --prefix={local_dir} --with-libraries=date_time,filesystem,iostreams,math,program_options,random,regex,serialization,signals,system,test,thread,log  && ./b2 cxxflags=\"-std=clib\" -d 2 toolset=gcc -j {num_cores} install""".format(local_dir=shellquote(i.local_dir).replace(' ', '\ '), num_cores=self.num_cores())
        cmd = " ".join(cmd.split())
        self.__build(i, cmd)

    def pear(self):
        i = self.__path("pear")
        cmd = """
            echo $(pwd) && ./configure --prefix={local_dir}/../../../ 
            && make -j {num_cores} && make install""".format(local_dir=shellquote(i.local_dir).replace(' ', '\ '),
                                                              num_cores=self.num_cores())
        cmd = " ".join(cmd.split())
        self.__build(i, cmd)

    def r(self):
        i = self.__path("r")
        if isMac():
            cmd = """
                ./configure --prefix={local_dir} --enable-R-shlib --with-x=no CC={CC} CXX={CXX} OBJC={CC}
                && make -j {num_cores}
                && make install
                && echo 'install.packages(c(\"gridExtra\", \"ape\", \"ggplot2\", \"seqinr\",\"Rcpp\", \"RInside\", \"devtools\"),
                 repos=\"http://cran.us.r-project.org\")' | $({local_dir}/R.framework/Resources/bin/R RHOME)/bin/R --slave --vanilla
                """.format(local_dir=shellquote(i.local_dir).replace(' ', '\ '), num_cores=self.num_cores(), CC=self.CC, CXX=self.CXX)
        else:
            cmd = """
                ./configure --prefix={local_dir} --enable-R-shlib --with-x=no CC={CC} CXX={CXX} OBJC={CC}
                && make -j {num_cores}
                && make install
                && echo 'install.packages(c(\"gridExtra\", \"ape\", \"ggplot2\", \"seqinr\",\"Rcpp\", \"RInside\",\"devtools\"),
                 repos=\"http://cran.us.r-project.org\")' | $({local_dir}/bin/R RHOME)/bin/R --slave --vanilla
            """.format(local_dir=shellquote(i.local_dir).replace(' ', '\ '), num_cores=self.num_cores(), CC=self.CC, CXX=self.CXX)
        cmd = " ".join(cmd.split())
        self.__build(i, cmd)

    def bamtools(self):
        i = self.__path('bamtools')
        cmd = """mkdir -p build && cd build && CC={CC} CXX={CXX} cmake -DCMAKE_INSTALL_PREFIX:PATH={local_dir} .. 
        && make -j {num_cores} install""".format(
            local_dir=shellquote(i.local_dir), num_cores=self.num_cores(), CC=self.CC, CXX=self.CXX)
        cmd = " ".join(cmd.split())
        self.__buildFromGit(i, cmd)

    def bibcpp(self):
        i = self.__path('bibcpp')
        self.__buildNjhProject(i)
    
    def bibcppDev(self):
        i = self.__path('bibcppdev')
        self.__buildNjhProject(i)

    def bibseq(self):
        i = self.__path('bibseq')
        self.__buildNjhProject(i)
        
    def bibseqDev(self):
        i = self.__path('bibseqdev')
        self.__buildNjhProject(i)
        
    def SeekDeep(self):
        i = self.__path('seekdeep')
        self.__buildNjhProject(i)
    
    def SeekDeepDev(self):
        i = self.__path('seekdeepdev')
        self.__buildNjhProject(i)
        
    def seqserver(self):
        i = self.__path('seqserver')
        self.__buildNjhProject(i)
        
    def njhRInside(self):
        i = self.__path('njhrinside')
        self.__buildNjhProject(i)
    
    def jsoncpp(self):
        i = self.__path('jsoncpp')
        cmd = """mkdir -p build && cd build && CC={CC} CXX={CXX} cmake -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_EXE_LINKER_FLAGS=-fPIC -DCMAKE_INSTALL_PREFIX:PATH={local_dir} .. 
        && make -j {num_cores} install""".format(
            local_dir=shellquote(i.local_dir), num_cores=self.num_cores(),CC=self.CC, CXX=self.CXX)
        cmd = " ".join(cmd.split())
        self.__buildFromGit(i, cmd)
    
    def cppcms(self):
        i = self.__path('cppcms')
        cmd = "mkdir -p build && cd build && CC={CC} CXX={CXX} cmake -DCMAKE_INSTALL_PREFIX:PATH={local_dir} .. && make -j {num_cores} install".format(local_dir=shellquote(i.local_dir), num_cores=self.num_cores(), CC=self.CC, CXX=self.CXX)  
        if(sys.platform == "darwin"):
            cmd += " && install_name_tool -change libbooster.0.dylib {local_dir}/lib/libbooster.0.dylib {local_dir}/lib/libcppcms.1.dylib".format(local_dir=shellquote(i.local_dir), num_cores=self.num_cores())
        cmd = " ".join(cmd.split())
        self.__build(i, cmd)

    def armadillo(self):
        i = self.__path('armadillo')
        cmd = "mkdir -p build && cd build && CC={CC} CXX={CXX} cmake -DCMAKE_INSTALL_PREFIX:PATH={local_dir} .. && make -j {num_cores} install".format(
            local_dir=shellquote(i.local_dir), num_cores=self.num_cores(), CC=self.CC, CXX=self.CXX)
        cmd = " ".join(cmd.split())
        self.__build(i, cmd)

    def liblinear(self):
        i = self.__path('liblinear')
        cmd = "make && mkdir -p {local_dir} && cp predict train {local_dir}".format(
            local_dir=shellquote(i.local_dir))
        cmd = " ".join(cmd.split())
        self.__build(i, cmd)

    def mlpack(self):
        i = self.__path('mlpack')
        armadillo_dir = shellquote(i.local_dir).replace("mlpack", "armadillo")
        boost_dir = shellquote(i.local_dir).replace("mlpack", "boost")
        cmd = """
        mkdir -p build
        && cd build
        && CC={CC} CXX={CXX} cmake -D DEBUG=OFF -D PROFILE=OFF
         -D ARMADILLO_LIBRARY={armadillo_dir}/lib/libarmadillo.so.4.0.2
         -D ARMADILLO_INCLUDE_DIR={armadillo_dir}/include/
         -D CMAKE_INSTALL_PREFIX:PATH={local_dir} ..
         -DBoost_NO_SYSTEM_PATHS=TRUE -DBOOST_INCLUDEDIR={boost}/include/ -DBOOST_LIBRARYDIR={boost}/lib/
        && make -j {num_cores} install
        """.format(local_dir=shellquote(i.local_dir),
           armadillo_dir=armadillo_dir,
           num_cores=self.num_cores(),
           boost=boost_dir, CC=self.CC, CXX=self.CXX)
        cmd = " ".join(cmd.split('\n'))
        self.__build(i, cmd)

    def mathgl(self):
        i = self.__path('mathgl')
        if (self.args.clang):
            cmd = """mkdir -p build && cd build && CC={CC} CXX={CXX} cmake -DCMAKE_INSTALL_PREFIX:PATH={local_dir} -Denable-pthread=ON -Denable-openmp=OFF .. 
            && make -j {num_cores} install""".format(
            local_dir=shellquote(i.local_dir), num_cores=self.num_cores(), CC=self.CC, CXX=self.CXX)
        else:
            cmd = """mkdir -p build && cd build && CC={CC} CXX={CXX} cmake -DCMAKE_INSTALL_PREFIX:PATH={local_dir}  .. 
            && make -j {num_cores} install""".format(
            local_dir=shellquote(i.local_dir), num_cores=self.num_cores(), CC=self.CC, CXX=self.CXX)
        cmd = " ".join(cmd.split())
        self.__build(i, cmd)

    def zi_lib(self):
        self.__git(self.__path('zi_lib'))
        
    def pstreams(self):
        self.__git(self.__path('pstreams'))

    def cppitertools(self):
        self.__git(self.__path('cppitertools'))
        i = self.__path('cppitertools')
        cmd = "cd {d} && git checkout d4f79321842dd584f799a7d51d3e066a2cdb7cac".format(d=shellquote(i.local_dir))
        Utils.run(cmd)

    def cppprogutils(self):
        self.__git(self.__path('cppprogutils'))

    def catch(self):
        self.__git(self.__path('catch'))


def ubuntu(self):
        pkgs = """libbz2-dev python2.7-dev cmake libpcre3-dev zlib1g-dev libgcrypt11-dev libicu-dev
python doxygen doxygen-gui auctex xindy graphviz libcurl4-openssl-dev""".split()

def startSrc():
    if not os.path.isdir("src/"):
        os.mkdir("src/")
    if not os.path.isfile("src/main.cpp"):
        cmd = "./scripts/genHelloWorldCpp.sh src/main.cpp"
        Utils.run(cmd)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-compfile', type=str, nargs=1)
    parser.add_argument('-libs', type=str, help="The libraries to install")
    parser.add_argument('-printLibs', action = "store_true", help="Print Available Libs")
    parser.add_argument('-CC', type=str, nargs=1)
    parser.add_argument('-CXX', type=str, nargs=1)
    parser.add_argument('-instRPackageName',type=str, nargs=1)
    parser.add_argument('-instRPackageSource',type=str, nargs=1)
    parser.add_argument('-addBashCompletion', dest = 'addBashCompletion', action = 'store_true' )
    return parser.parse_args()

def main():
    args = parse_args()
    s = Setup(args)
    if(args.instRPackageName):
        s.installRPackageName(args.instRPackageName[0])
        return 0
    if(args.instRPackageSource):
        s.installRPackageSource(args.instRPackageSource[0])
        return 0
    if args.printLibs:
        s.printAvailableSetUps()
    elif args.addBashCompletion:
        if(os.path.isdir("./bashCompletes")):
            cmd = "cat bashCompletes/* >> ~/.bash_completion"
            Utils.run(cmd)
    else:
        s.setup()

main()
