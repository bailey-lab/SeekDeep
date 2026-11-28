#!/usr/bin/env python3

import shutil, os, argparse, sys, stat
import os, subprocess, sys

projectDir = os.path.dirname(os.path.abspath(__file__))
setupToolsDir = os.path.join(projectDir, "cppSetupTools")
if not os.path.exists(os.path.join(setupToolsDir, "setup.py")):
    if not os.path.exists(os.path.join(projectDir, ".git")):
        sys.exit("cppSetupTools is missing; it is a git submodule, so get this project with: git clone --recursive <url>")
    subprocess.run(["git", "submodule", "update", "--init", "cppSetupTools"], cwd=projectDir, check=True)
sys.path.append(os.path.join(setupToolsDir, "scripts", "pyUtils"))
sys.path.append(os.path.join(setupToolsDir, "scripts", "setUpScripts"))
from utils import Utils
from genFuncs import genHelper
def main():
    name = "SeekDeep"
    libs = "seqServer:develop"
    args = genHelper.parseNjhConfigureArgs()
    cmd = genHelper.mkConfigCmd(name, libs, sys.argv)
    Utils.run(cmd)

main()
