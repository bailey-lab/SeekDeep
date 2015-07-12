#!/usr/bin/env python

# by purcaro@gmail.com

import subprocess, sys, os, argparse
from collections import namedtuple
import shutil
sys.path.append("scripts/pyUtils")
sys.path.append("scripts/setUpScripts")
from utils import Utils
from genFuncs import genHelper 
from color_text import ColorText as CT


BuildPaths = namedtuple("BuildPaths", 'url build_dir build_sub_dir local_dir')

def shellquote(s):
    # from http://stackoverflow.com/a/35857
    return "'" + s.replace("'", "'\\''") + "'"

def isMac():
    return Utils.isMac()

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
        self.paths["dlib"] = self.__dlib()
        self.paths["libsvm"] = self.__libsvm()
        self.paths["mongoc"] = self.__mongoc()
        self.paths["mongocxx"] = self.__mongocxx()
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
    
    def __mongoc(self):
        url = "https://github.com/mongodb/mongo-c-driver"
        name = "mongoc"
        return self.__package_dirs(url, name)
    
    def __mongocxx(self):
        url = "https://github.com/mongodb/mongo-cxx-driver"
        name = "mongocxx"
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
        #url = "http://sourceforge.net/projects/boost/files/boost/1.57.0/boost_1_57_0.tar.gz"
        url = "http://sourceforge.net/projects/boost/files/boost/1.58.0/boost_1_58_0.tar.gz"
        return self.__package_dirs(url, "boost")

    def __armadillo(self):
        #url = "http://freefr.dl.sourceforge.net/project/arma/armadillo-4.000.2.tar.gz"
        url = "http://freefr.dl.sourceforge.net/project/arma/armadillo-5.100.2.tar.gz"
        return self.__package_dirs(url, "armadillo")

    def __mlpack(self):
        url = "http://www.mlpack.org/files/mlpack-1.0.8.tar.gz"
        return self.__package_dirs(url, "mlpack")

    def __liblinear(self):
        url = "http://www.csie.ntu.edu.tw/~cjlin/liblinear/oldfiles/liblinear-1.94.tar.gz"
        return self.__package_dirs(url, "liblinear")

    def __cppcms(self):
        url = "http://freefr.dl.sourceforge.net/project/cppcms/cppcms/1.0.5/cppcms-1.0.5.tar.bz2"
        return self.__package_dirs(url, "cppcms") 

    def __mathgl(self):
        url = "http://freefr.dl.sourceforge.net/project/mathgl/mathgl/mathgl%202.2.1/mathgl-2.2.1.tar.gz"
        return self.__package_dirs(url, "mathgl")
    
    def __dlib(self):
        url = "http://freefr.dl.sourceforge.net/project/dclib/dlib/v18.7/dlib-18.7.tar.bz2"
        return self.__package_dirs(url, "dlib")
    
    def __libsvm(self):
        url = "http://www.csie.ntu.edu.tw/~cjlin/libsvm/oldfiles/libsvm-3.18.tar.gz"
        return self.__package_dirs(url, "libsvm")

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
        name = "seqServer"
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
        self.extDirLoc = "" # the location where the libraries will be installed
        #if no compile file set up and assume external is next to setup.py
        if not args.compfile:
            self.extDirLoc = os.path.abspath(os.path.join(os.path.dirname(__file__), "external"))
        else:
            self.extDirLoc = os.path.abspath(self.parseForExtPath(args.compfile[0]))
        self.paths = Paths(self.extDirLoc) # path object to hold the paths for install
        self.args = args # command line arguments parsed by argument parser
        self.setUps = {} # all available set ups
        self.setUpsNeeded = [] # the setups that need to be done
        self.installed = [] # the setups that able to install
        self.failedInstall = [] # the setups that failed
        self.CC = "" # the c compilier being used
        self.CXX = "" # the c++ compilier being used
        self.njhProjects = ["bibcpp", "bibcppdev", "bibseq", "bibseqdev", "seekdeep", "seekdeepdev", "njhrinside", "seqserver"]
        self.__initSetUpFuncs()
        self.__processArgs()
        
    def setup(self):
        if self.args.forceUpdate:
            for set in self.setUpsNeeded:
                if not set in self.setUps.keys():
                    print CT.boldBlack( "Unrecognized option ") + CT.boldRed(set)
                else:
                    self.rmDirsForLib(set)
                        
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
                       "pstreams": self.pstreams,
                       "dlib": self.dlib,
                       "libsvm": self.libsvm,
                       "mongoc": self.mongoc,
                       "mongocxx": self.mongocxx
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
        # if no compfile need to determine compiler, will default to env CC and CXX
        else:
            self.CC = genHelper.determineCC(self.args)
            self.CXX = genHelper.determineCXX(self.args)
        if "clang" in self.CXX:
            self.args.clang = True
        else:
            self.args.clang = False

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
    
    def rmDirsForLibs(self,libs):
        for l in libs:
            self.rmDirsForLib(l)
    
    def rmDirsForLib(self,lib):
        if lib not in self.setUps:
            print CT.boldBlack( "Unrecognized lib: ") + CT.boldRed(lib)
        else:
            p = self.__path(lib)
            if p.build_dir:
                print "Removing " + CT.boldBlack(p.build_dir)
                Utils.rm_rf(p.build_dir)
            if p.local_dir:
                print "Removing " + CT.boldBlack(p.local_dir)
                Utils.rm_rf(p.local_dir)
    

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
            except Exception as inst:
                print type(inst)
                print inst 
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
        retCores = Utils.num_cores()
        if self.args.numCores:
            if not self.args.numCores > retCores:
                retCores = self.args.numCores
        else:
            if retCores > 8:
                retCores  = retCores/2
            if 1 != retCores:
                retCores -= 1
        return retCores

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


    def __buildNjhProject(self,libPaths):
        cmd = """
        python ./configure.py -CC {CC} -CXX {CXX} -externalLibDir {external} -prefix {localTop} 
        && python ./setup.py -compfile compfile.mk 
        && make -j {num_cores} && make install""".format(localTop=shellquote(self.paths.install_dir),
                                                          num_cores=self.num_cores(), CC=self.CC, CXX=self.CXX,
                                                           external=self.extDirLoc)
        cmd = " ".join(cmd.split())
        self.__buildFromGit(libPaths, cmd)
        
    def __buildNjhProjectTag(self,libPaths,tagName):
        cmd = """
        python ./configure.py -CC {CC} -CXX {CXX} -externalLibDir {external} -prefix {localTop} 
        && python ./setup.py -compfile compfile.mk 
        && make -j {num_cores} && make install""".format(localTop=shellquote(self.paths.install_dir),
                                                          num_cores=self.num_cores(), CC=self.CC, CXX=self.CXX,
                                                           external=self.extDirLoc)
        cmd = " ".join(cmd.split())
        self.__buildFromGitTag(libPaths, cmd,tagName)
        
    
    def __buildNjhProjectBranch(self,libPaths,branchName):
        cmd = """
        python ./configure.py -CC {CC} -CXX {CXX} -externalLibDir {external} -prefix {localTop} 
        && python ./setup.py -compfile compfile.mk 
        && make -j {num_cores} && make install""".format(localTop=shellquote(self.paths.install_dir),
                                                          num_cores=self.num_cores(), CC=self.CC, CXX=self.CXX,
                                                           external=self.extDirLoc)
        cmd = " ".join(cmd.split())
        self.__buildFromGitTag(libPaths, cmd,branchName)
        
    def updateNjhProject(self,lib):
        cmd = """
        python ./configure.py -CC {CC} -CXX {CXX} -externalLibDir {external} -prefix {localTop} 
        && python ./setup.py -compfile compfile.mk
        && make clean 
        && make -j {num_cores} && make install""".format(localTop=shellquote(self.paths.install_dir),
                                                          num_cores=self.num_cores(), CC=self.CC, CXX=self.CXX,
                                                           external=self.extDirLoc)
        cmd = " ".join(cmd.split())
        libPaths = self.__path(lib.lower())
        self.__buildFromGit(libPaths, cmd)
    
    def updateNjhProjects(self, libs):
        for l in libs:
            libLower = l.lower()
            if libLower in self.njhProjects:
                self.updateNjhProject(libLower)
    
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
            
    def __buildFromGitBranch(self, i, cmd, branchName):
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
            cCmd = "git clone -b "+ branchName + " {url} {d}".format(url=i.url, d=i.build_dir)
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
    
    def __buildFromGitTag(self, i, cmd, tagName):
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
            tagCmd = "git checkout {tag}".format(tag=tagName)
            try:
                print self.paths.ext_build
                Utils.run_in_dir(cCmd, self.paths.ext_build)
                Utils.run_in_dir(tagCmd, i.build_dir)
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
        print "start"
        i = self.__path("boost")
        #boostLibs = "date_time,filesystem,iostreams,math,program_options,random,regex,serialization,signals,system,test,thread,log"
        boostLibs = "filesystem,iostreams,system"
        if isMac():
            print "here"
            setUpDir = os.path.dirname(os.path.abspath(__file__))
            gccJamLoc =  os.path.join(setUpDir, "scripts/etc/boost/gcc.jam")
            gccJamOutLoc = os.path.abspath("{build_dir}/tools/build/src/tools/gcc.jam".format(build_dir = i.build_sub_dir))
            print gccJamLoc
            print gccJamOutLoc
            installNameToolCmd  = """ 
            && install_name_tool -change $(otool -L {local_dir}/lib/libboost_filesystem.dylib | egrep -o "\\S.*libboost_system.dylib") {local_dir}/lib/libboost_system.dylib {local_dir}/lib/libboost_filesystem.dylib
            && install_name_tool -id {local_dir}/lib/libboost_filesystem.dylib {local_dir}/lib/libboost_filesystem.dylib
            && install_name_tool -id {local_dir}/lib/libboost_iostreams.dylib {local_dir}/lib/libboost_iostreams.dylib
            && install_name_tool -id {local_dir}/lib/libboost_system.dylib {local_dir}/lib/libboost_system.dylib
            """.format(local_dir=i.local_dir)
        if self.args.clang:
             if isMac():
                cmd = """./bootstrap.sh --with-toolset=clang --prefix={local_dir} --with-libraries=""" + boostLibs + """
                  &&  ./b2  toolset=clang cxxflags=\"-stdlib=libc++ -std=c++14\" linkflags=\"-stdlib=libc++\" -j {num_cores} install 
                  &&  install_name_tool -change libboost_system.dylib {local_dir}/lib/libboost_system.dylib {local_dir}/lib/libboost_filesystem.dylib
                  """
                  #&&  install_name_tool -change libboost_system.dylib {local_dir}/lib/libboost_system.dylib {local_dir}/lib/libboost_thread.dylib
                #cmd = """wget https://github.com/boostorg/atomic/commit/6bb71fdd.diff && wget https://github.com/boostorg/atomic/commit/e4bde20f.diff&&  wget https://gist.githubusercontent.com/philacs/375303205d5f8918e700/raw/d6ded52c3a927b6558984d22efe0a5cf9e59cd8c/0005-Boost.S11n-include-missing-algorithm.patch&&  patch -p2 -i 6bb71fdd.diff&&  patch -p2 -i e4bde20f.diff&&  patch -p1 -i 0005-Boost.S11n-include-missing-algorithm.patch&&  echo "using clang;  " >> tools/build/v2/user-config.jam&&  ./bootstrap.sh --with-toolset=clang --prefix={local_dir} --with-libraries=""" + boostLibs + """  &&  ./b2  -d 2 toolset=clang cxxflags=\"-stdlib=libc++\" linkflags=\"-stdlib=libc++\" -j {num_cores} install &&  install_name_tool -change libboost_system.dylib {local_dir}/lib/libboost_system.dylib {local_dir}/lib/libboost_thread.dylib&&  install_name_tool -change libboost_system.dylib {local_dir}/lib/libboost_system.dylib {local_dir}/lib/libboost_filesystem.dylib""".format(local_dir=shellquote(i.local_dir).replace(' ', '\ '), num_cores=self.num_cores())
             else:
                cmd = """./bootstrap.sh --with-toolset=clang --prefix={local_dir}  --with-libraries=""" + boostLibs + """ &&  ./b2 toolset=clang cxxflags=\"-std=c++14\" -j {num_cores} install""".format(local_dir=shellquote(i.local_dir).replace(' ', '\ '), num_cores=self.num_cores())
                #cmd = """wget https://github.com/boostorg/atomic/commit/6bb71fdd.diff && wget https://github.com/boostorg/atomic/commit/e4bde20f.diff&&  wget https://gist.githubusercontent.com/philacs/375303205d5f8918e700/raw/d6ded52c3a927b6558984d22efe0a5cf9e59cd8c/0005-Boost.S11n-include-missing-algorithm.patch&&  patch -p2 -i 6bb71fdd.diff&&  patch -p2 -i e4bde20f.diff&&  patch -p1 -i 0005-Boost.S11n-include-missing-algorithm.patch&&  echo "using clang;  " >> tools/build/v2/user-config.jam&&  ./bootstrap.sh --with-toolset=clang --prefix={local_dir}  --with-libraries=""" + boostLibs + """ &&  ./b2  -d 2 toolset=clang -j {num_cores} install""".format(local_dir=shellquote(i.local_dir).replace(' ', '\ '), num_cores=self.num_cores())
        elif self.CXX == "g++-4.8":
            if isMac():
                cmd = "cp " + gccJamLoc + "  " + gccJamOutLoc + """ && echo "using gcc : 4.8 : g++-4.8 : <linker-type>darwin ;" >> project-config.jam && ./bootstrap.sh --with-toolset=gcc --prefix={local_dir} --with-libraries=""" + boostLibs + """
                 && ./b2 --toolset=gcc-4.8 -j {num_cores} install 
                 """ + installNameToolCmd
            else:
                cmd = """echo "using gcc : 4.8 : g++-4.8;" >> project-config.jam && ./bootstrap.sh --with-toolset=gcc --prefix={local_dir} --with-libraries=""" + boostLibs + """
                 && ./b2 --toolset=gcc-4.8 -j {num_cores} install 
                 """
        elif self.CXX == "g++-4.9":
            if isMac():
                cmd = "cp " + gccJamLoc + "  " + gccJamOutLoc + """ && echo "using gcc : 4.9 : g++-4.9 : <linker-type>darwin ;" >> project-config.jam && ./bootstrap.sh --with-toolset=gcc --prefix={local_dir} --with-libraries=""" + boostLibs + """
                 && ./b2 --toolset=gcc-4.9 -j {num_cores} install 
                 """ + installNameToolCmd
            else:
                cmd = """echo "using gcc : 4.9 : g++-4.9;" >> project-config.jam && CC={CC} CXX={CXX} ./bootstrap.sh --with-toolset=gcc --prefix={local_dir} --with-libraries=""" + boostLibs + """
                 && CC={CC} CXX={CXX}  ./b2 --toolset=gcc-4.9 -j {num_cores} install 
                 """
        elif self.CXX == "g++-5":
            if isMac():
                cmd = """echo "using gcc : 5 : g++-5;" >> project-config.jam && CC={CC} CXX={CXX} ./bootstrap.sh --with-toolset=gcc --prefix={local_dir} --with-libraries=""" + boostLibs + """
                 && CC={CC} CXX={CXX} ./b2 --toolset=gcc-5 -j {num_cores} install 
                 """ + installNameToolCmd
            else:
                cmd = """echo "using gcc : 5 : g++-5;" >> project-config.jam && CC={CC} CXX={CXX} ./bootstrap.sh --with-toolset=gcc --prefix={local_dir} --with-libraries=""" + boostLibs + """
                 && CC={CC} CXX={CXX}  ./b2 --toolset=gcc-5 -j {num_cores} install 
                 """
        elif self.CXX == "g++":
            if isMac():
                cmd = "cp " + gccJamLoc + "  " + gccJamOutLoc + """ && echo "using gcc : 4.9 : g++ : <linker-type>darwin ;" >> project-config.jam && ./bootstrap.sh --with-toolset=gcc --prefix={local_dir} --with-libraries=""" + boostLibs + """
                 && ./b2 --toolset=gcc-4.9 -j {num_cores} install 
                 """ + installNameToolCmd
            else:
                cmd = """echo "using gcc : 4.9 : g++;" >> project-config.jam && ./bootstrap.sh --with-toolset=gcc --prefix={local_dir} --with-libraries=""" + boostLibs + """
                 && ./b2 --toolset=gcc-4.9 -j {num_cores} install 
                 """
        
        cmd = " ".join(cmd.split())
        cmd = cmd.format(local_dir=shellquote(i.local_dir).replace(' ', '\ '), num_cores=self.num_cores(), CC=self.CC, CXX=self.CXX)
        print cmd
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
        cmd = """git checkout v2.4.0 && mkdir -p build && cd build && CC={CC} CXX={CXX} cmake -DCMAKE_INSTALL_PREFIX:PATH={local_dir} .. 
        && make -j {num_cores} install""".format(
            local_dir=shellquote(i.local_dir), num_cores=self.num_cores(), CC=self.CC, CXX=self.CXX)
        cmd = " ".join(cmd.split())
        self.__buildFromGit(i, cmd)

    def bibcpp(self):
        i = self.__path('bibcpp')
        branch = "release/2"
        version = "2"
        self.__buildNjhProject(i)
        #self.__buildNjhProjectTag(i, version)
        #self.__buildNjhProjectBranch(i, branch)
    
    def bibcppDev(self):
        i = self.__path('bibcppdev')
        self.__buildNjhProject(i)

    def bibseq(self):
        i = self.__path('bibseq')
        branch = "release/2"
        version = "2"
        self.__buildNjhProject(i)
        #self.__buildNjhProjectTag(i, version)
        #self.__buildNjhProjectBranch(i, branch)
        
    def bibseqDev(self):
        i = self.__path('bibseqdev')
        self.__buildNjhProject(i)
        
    def SeekDeep(self):
        i = self.__path('seekdeep')
        branch = "release/2"
        version = "2"
        self.__buildNjhProject(i)
        #self.__buildNjhProjectTag(i, version)
        #self.__buildNjhProjectBranch(i, branch)
    
    
    def SeekDeepDev(self):
        i = self.__path('seekdeepdev')
        self.__buildNjhProject(i)
        
    def seqserver(self):
        i = self.__path('seqserver')
        branch = "develop"
        version = "2"
        #self.__buildNjhProject(i)
        #self.__buildNjhProjectTag(i, version)
        self.__buildNjhProjectBranch(i, branch)
        
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
        
    def mongoc(self):
        i = self.__path('mongoc')
        cmd = """sed -i.bak s/git:/http:/g .gitmodules && CC={CC} CXX={CXX} ./autogen.sh --prefix={local_dir}
        && make -j {num_cores}  && make install""".format(
            local_dir=shellquote(i.local_dir), num_cores=self.num_cores(),CC=self.CC, CXX=self.CXX)
        cmd = " ".join(cmd.split())
        branchName = "1.2.0-dev"
        self.__buildFromGitBranch(i, cmd, branchName)
        
    def mongocxx(self):
        i = self.__path('mongocxx')
        cmd = """cd build && CC={CC} CXX={CXX} PKG_CONFIG_PATH={ext_dir}/local/mongoc/lib/pkgconfig:$PKG_CONFIG_PATH cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX={local_dir} .. 
        && make -j {num_cores} && make install""".format(
            local_dir=i.local_dir, num_cores=self.num_cores(),CC=self.CC, CXX=self.CXX, ext_dir=self.extDirLoc)
        cmd = " ".join(cmd.split())
        branchName = "master"
        self.__buildFromGitBranch(i, cmd, branchName)
    
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
        cmd = """
perl -p -i -e 's/if\(check_probability_model/if\(1 || check_probability_model/' linear.cpp &&
make &&
mkdir -p {local_dir} &&
cp predict train {local_dir} &&
make lib &&
cp linear.h liblinear.so.1 README {local_dir} &&
ln -s {local_dir}/liblinear.so.1 {local_dir}/liblinear.so
""".format(local_dir=shellquote(i.local_dir))
        cmd = " ".join(cmd.split("\n"))
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
    
    def dlib(self):
        i = self.__path('dlib')
        cmd = """
mkdir {local_dir} &&
cp -a * {local_dir}/
""".format(local_dir=shellquote(i.local_dir), num_cores=self.num_cores())
        cmd = " ".join(cmd.split('\n'))
        self.__build(i, cmd)
        
    def libsvm(self):
        i = self.__path('libsvm')
        cmd = "make && make lib && mkdir -p {local_dir} && cp -a * {local_dir}".format(
            local_dir=shellquote(i.local_dir))
        self.__build(i, cmd)

    def cppprogutils(self):
        self.__git(self.__path('cppprogutils'))

    def catch(self):
        self.__git(self.__path('catch'))


def ubuntu(self):
        pkgs = """libbz2-dev python2.7-dev cmake libpcre3-dev zlib1g-dev libgcrypt11-dev libicu-dev
python doxygen doxygen-gui auctex xindy graphviz libcurl4-openssl-dev""".split()



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-compfile', type=str, nargs=1)
    parser.add_argument('-libs', type=str, help="The libraries to install")
    parser.add_argument('-printLibs', action = "store_true", help="Print Available Libs")
    parser.add_argument('-forceUpdate', action = "store_true", help="Remove already installed libs and re-install")
    parser.add_argument('-updateNjhProjects', type = str, help="Remove already installed libs and re-install")
    parser.add_argument('-CC', type=str, nargs=1)
    parser.add_argument('-CXX', type=str, nargs=1)
    parser.add_argument('-instRPackageName',type=str, nargs=1)
    parser.add_argument('-instRPackageSource',type=str, nargs=1) 
    parser.add_argument('-addBashCompletion', dest = 'addBashCompletion', action = 'store_true')
    parser.add_argument('-numCores', type=str)
    return parser.parse_args()

def main():
    args = parse_args()
    s = Setup(args)
    ccWhich = Utils.which(s.CC)
    cxxWhich = Utils.which(s.CXX)
    cmakeWhich = Utils.which("cmake")
    if not ccWhich or not cxxWhich or not cmakeWhich:
        if not ccWhich:
            print CT.boldRed("Could not find c compiler " + CT.purple + s.CC)
            if args.compfile:
                print "Change CC in " + args.compfile
            else:
                print "Can supply another c compiler by using -CC [option] or by defining bash environmental CC "
            print ""
        if not cxxWhich:
            print CT.boldRed("Could not find c++ compiler " + CT.purple + s.CXX)
            if args.compfile:
                print "Change CXX in " + args.compfile
            else:
                print "Can supply another c++ compiler by using -CXX [option] or by defining bash environmental CXX "
            print ""
        if not cmakeWhich:
            print CT.boldRed("Could not find " + CT.purple + "cmake")
            if Utils.isMac():
                print "If you have brew install via, brew update && brew install cmake, otherwise you can follow instructions from http://www.cmake.org/install/"
            else:
                print "On ubuntu to install latest cmake do the following"
                print "sudo add-apt-repository ppa:george-edison55/cmake-3.x"
                print "sudo apt-get update"
                print "sudo apt-get install cmake"
        return 1
        
    
    if(args.instRPackageName):
        s.installRPackageName(args.instRPackageName[0])
        return 0
    if(args.instRPackageSource):
        s.installRPackageSource(args.instRPackageSource[0])
        return 0
    if args.updateNjhProjects:
        projectsSplit = args.updateNjhProjects.split(",")
        s.updateNjhProjects(projectsSplit)
    if args.printLibs:
        s.printAvailableSetUps()
    elif args.addBashCompletion:
        if(os.path.isdir("./bashCompletes")):
            cmd = "cat bashCompletes/* >> ~/.bash_completion"
            Utils.run(cmd)
        if(os.path.isdir("./bash_completion.d")):
            cmd = "cat bash_completion.d/* >> ~/.bash_completion"
            Utils.run(cmd)
    else:
        s.setup()

main()
