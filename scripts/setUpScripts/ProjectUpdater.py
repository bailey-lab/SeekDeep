#!/usr/bin/env python

import shutil, os, argparse, sys, stat, time, platform
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), "pyUtils"))
from utils import Utils
from color_text import ColorText as CT
from shutil import ignore_patterns


class ProjectUpdater():
    def __init__(self, externalLoc, numCores, setupFrom, projectDir = "./"):
        self.numCores = numCores
        #make sure not to take all cores or more than the amount of cores available
        if not self.numCores or self.numCores > Utils.num_cores():
            self.numCores = Utils.num_cores()
        if self.numCores == Utils.num_cores():
            self.numCores -= 1
        self.externalLoc = os.path.abspath(externalLoc)
        self.projectDir = os.path.abspath(projectDir)
        self.setupFrom = os.path.abspath(setupFrom)
    
    def remakeBibseqProject(self, subProjectDir):
        makeCmd = ""
        cmd =   """
                if [ -f .git ]; then git pull; fi
                && if [ -f .gitmodules ]; then git submodule init && git submodule update; fi
                && ./configure.py -externalLibDir {EXT_LOC} 
                && ./setup.py --compfile compfile.mk --outMakefile makefile-common.mk --overWrite  
                && make clean && make -j {NUMCORES} 
                && if [ -d $(make printInstallDir) ]; then rm -rf $(make printInstallDir); fi
                && mkdir -p $(make printInstallDir)
                && make -j {NUMCORES} install""".format(EXT_LOC = self.externalLoc, NUMCORES = self.numCores)
        cmd = cmd + makeCmd
        cmd = " ".join(cmd.split())
        Utils.run_in_dir(cmd, subProjectDir)
    
    def remakeBibseqProjects(self, dirs):
        """
            Remake dirs in the project directory
        """
        installed = []
        failed = []
        nonExistent = []
        for dir in dirs:
            dir = os.path.join(self.projectDir, dir)
            if os.path.exists(dir):
                try:
                    self.remakeBibseqProject(dir)
                    installed.append(dir)
                except Exception as inst:
                    print inst
                    failed.append(dir)
            else:
                nonExistent.append(dir)
        #installed
        for d in installed:
            print CT.green + "Updated " + d + CT.reset
        #failed
        for d in failed:
            print CT.red + "Failed to update " + d + CT.reset
        #didn't exist
        for d in nonExistent:
            print CT.cyan + os.path.basename(d) + " doesn't exist in " + self.projectDir + ", skipped" + CT.reset
            
    def copyDir(self, src, dist, overWrite):
        if(os.path.exists(dist)):
            if(overWrite):
                shutil.rmtree(dist)
            else:
                raise Exception( "Error, directory " + str(dist) + " already exist, set overwrite to overwrite dir")
        shutil.copytree(src, dist, ignore=ignore_patterns('*.pyc', '*~', '.*'))
    
    def copySetUp(self, distDir, overWrite=False):
        if not os.path.exists(distDir):
            print CT.boldGreen(distDir) + CT.boldRed(" doesn't exist skipping")
            return
        if os.path.exists(os.path.join(distDir, "setup.py")) and not overWrite:
            raise Exception("Error, file " + os.path.join(distDir, "setup.py") + " already exist, use --overWrite overwrite setup.py")
        shutil.copy(os.path.join(self.setupFrom, "setup.py"), os.path.join(distDir, "setup.py"))    
        self.copyDir(os.path.join(self.setupFrom ,"scripts"), os.path.join(distDir,"scripts"), overWrite)
    
    def copySetUpToSub(self, distDir, overWrite=False):
        distDir = os.path.join(self.projectDir, distDir)
        self.copySetUp(distDir, overWrite)
    
    def copySetUpToMain(self, overWrite=False):
        distDir = self.projectDir
        self.copySetUp(distDir, overWrite)
        
    def copyMakefile(self, distDir, overWrite=False):
        distDir = os.path.join(self.projectDir, distDir)
        if not os.path.exists(distDir):
            print CT.boldGreen(distDir) + CT.boldRed(" doesn't exist skipping")
            return
        if os.path.exists(os.path.join(distDir, "Makefile")) and not overWrite:
            raise Exception("Error, file " + os.path.join(distDir, "Makefile") + " already exist, use --overWrite overwrite Makefile")
        shutil.copy(os.path.join(self.setupFrom, "scripts/cppMakefiles/Makefile"), os.path.join(distDir, "Makefile"))
    
    def reatchHeadAndPull(self, subModuleDir, branch):
        cmds = "cd " + os.path.join(self.projectDir, subModuleDir) + " && git checkout " + branch + " && git pull"
        Utils.run(cmds)
    
    def reatchHeadCommitPush(self, subModuleDir, branch, commitMessage):
        try:
            self.reatchHeadAndPull(subModuleDir, branch)
        except Exception as excp:
            print "Failed to reatch head and pull for " + subModuleDir + ":" + branch
            raise(excp)
        cmds = "cd " + os.path.join(self.projectDir, subModuleDir) + " && git commit -a -m \"" + commitMessage + "\" && git push"
        try:
            Utils.run(cmds)
        except Exception as pushExcep:
            print "Failed to commit and push for " + subModuleDir + ":" + branch
            raise(pushExcep)
        
   