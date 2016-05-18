#!/usr/bin/env python

import shutil, os, argparse, sys
sys.path.append(os.path.join(os.path.dirname(__file__), "pyUtils"))
from utils import Utils
from color_text import ColorText as CT


def genPyHello(outFileName):
    with open(outFileName, "w") as f:
        f.write("""#!/usr/bin/env python
import shutil, os, argparse, sys, stat

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--name', type=str)
    return parser.parse_args()

def main():
    args = parse_args()
    if(args.name):
        print "Hello " + args.name + "!"
    else:
        print "Hello World!"
    
if __name__ == "__main__":
    main()

""")
        

def genCppHello(outFileName):
    with open(outFileName, "w") as f:
        f.write("""#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>
#include <cstdint>
#include <cstdio>
#include <cstddef>
#include <utility>

int main(int argc, char* argv[])
{
    std::cout << \"Hello World!\" << std::endl;
    return 0;
}""")
        
def parse_args_genHello():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outFilename',type=str, required = True)
    parser.add_argument('--overWrite', action = "store_true", help = "Over write file if it exists")
    parser.add_argument('-p','--python', action = "store_true", help = "Generate a python hello world script, default is a cpp hello world")
    return parser.parse_args()

def genHellos(outFilename, overWrite = False, python = False):
    if python:
        if not outFilename.endswith(".py"):
            outFilename = outFilename + ".py"
        if os.path.exists(outFilename) and not overWrite:
            raise Exception(CT.boldRed("File ") + CT.boldBlack(outFilename) + CT.boldRed(" already exists, use --overWrite to over write it"))
        genPyHello(outFilename)
        #rwx to user and group, r-x to everyone
        fd = os.open( outFilename, os.O_RDONLY )
        os.fchmod( fd, 0775)
        os.close( fd )
        print (CT.boldGreen("Now run"))
        print ("./" + outFilename + " ")
        print (CT.boldGreen("or"))
        print ("./" + outFilename + " --name Nick")
    else:
        if not outFilename.endswith(".cpp"):
            outFilename = outFilename + ".cpp"
        if os.path.exists(outFilename) and not overWrite:
            raise Exception("File " + outFilename + " already exists, use --overWrite to over write it")
        genCppHello(outFilename)
        print (CT.boldGreen("Now run"))
        print ("g++ -std=c++11 " + outFilename + " -o hello #-std=c++11 needed for cstdint include in libstdc++")
        print ("./hello")
        print (CT.boldGreen("or run"))
        #mac
        if Utils.isMac():
            print ("clang++ " + outFilename + " -o hello")
        else:
            print ("clang++-3.6 -std=c++11 " + outFilename + " -o hello #-std=c++11 needed for cstdint include in libstdc++")
        print ("./hello")
if __name__ == '__main__':
    args = parse_args_genHello()
    genHellos(args.outFilename, args.overWrite, args.python)

