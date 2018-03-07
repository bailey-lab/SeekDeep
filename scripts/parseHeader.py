#!/usr/bin/env python3

import shutil, os, argparse, sys, stat
import CppHeaderParser
sys.path.append(os.path.join(os.path.dirname(__file__), "pyUtils"))
from color_text import ColorText as CT

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str, required = True)
    return parser.parse_args()

def main():
    args = parse_args()
    try:
        cppHeader = CppHeaderParser.CppHeader(args.file)
    except CppHeaderParser.CppParseError as e:
        print(e)
        sys.exit(1)
    
    print(CT.boldBlack("Class public methods"))
    for k in list(cppHeader.classes.keys()):
        print(CT.boldBlack(k))
        for i in range(len(cppHeader.classes[k]["methods"]["public"])):
            print("\t",cppHeader.classes[k]["methods"]["public"][i]["name"])
    print("")
    print(CT.boldBlack("Includes"))        
    for include in cppHeader.includes:
        if "//" in include:
            print("\t" + include[0:include.find("//")].strip())
        elif "/*" in include:
            print("\t" + include[0:include.find("/*")].strip())
        else:
            print("\t" + include)
    
    print(("\n" + CT.boldBlack("Free functions are:")))
    for func in cppHeader.functions:
        print(("\t%s"%func["name"]))
main()
