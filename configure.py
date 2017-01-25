#!/usr/bin/env python

import shutil, os, argparse, sys, stat
sys.path.append("scripts/pyUtils")
sys.path.append("scripts/setUpScripts")
from utils import Utils
from genFuncs import genHelper
def main():
    name = "SeekDeep"
    libs = "seqServer:release/v1.5.0"
    args = genHelper.parseNjhConfigureArgs()
    cmd = genHelper.mkConfigCmd(name, libs, sys.argv)
    Utils.run(cmd)
    
main()
