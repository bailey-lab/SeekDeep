#!/usr/bin/env python

import shutil, os, argparse, sys


def genPyHello(outFileName):
    with open(outFileName, "w") as f:
        f.write("#!/usr/bin/env python\n")
        f.write("\n")
        f.write("import shutil, os, argparse, sys, stat\n")
        f.write("\n")
        f.write("def parse_args():\n")
        f.write("    parser = argparse.ArgumentParser()\n")
        f.write("    parser.add_argument('-name', type=str, nargs=1)\n")
        f.write("    return parser.parse_args()\n")
        f.write("\n")
        f.write("def main():\n")
        f.write("    args = parse_args()\n")
        f.write("    if(args.name):\n")
        f.write("        print \"Hello \" + args.name[0] + \"!\"\n")
        f.write("    else:\n")
        f.write("        print \"Hello World!\"\n")
        f.write("    \n")
        f.write("main()\n")
        
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-outFilename',type=str, required = True)
    return parser.parse_args()

def main():
    args = parse_args()
    genPyHello(args.outFilename[0])
    #give execute rights to everyone 
    fd = os.open( args.outFilename[0], os.O_RDONLY )
    os.fchmod( fd, 0775)
    os.close( fd )

    
main()

