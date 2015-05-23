#!/usr/bin/env python

from __future__ import print_function
import sys, os, argparse
from sets import Set
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), "pyUtils"))
from color_text import ColorText as CT

def warning(*objs):
    # from http://stackoverflow.com/a/14981125
    print("\nWARNING: ", *objs, file=sys.stderr)

objToHeader = {
               "std::string": "#include<string>",
               "std::basic_string": "#include<string>",
               "std::copy" : "#include<algorithm>",
               "std::find" : "#include<algorithm>",
               "std::sort" : "#include<algorithm>",
               "std::copy" : "#include<algorithm>",
               "std::find" : "#include<algorithm>",
               "std::sort" : "#include<algorithm>",
               "std::array" : "#include<array>",
               "std::duration" : "#include<chrono>",
               "std::time_point" : "#include<chrono>",
               "std::sqrt" : "#include<cmath>",
               "std::pow" : "#include<cmath>",
               "std::fstream" : "#include<fstream>",
               "std::ifstream" : "#include<fstream>",
               "std::ofstream" : "#include<fstream>",
               "std::future" : "#include<future>",
               "std::promise" : "#include<future>",
               "std::thread" : "#include<thread>",
               "std::mutex" : "#include<mutex>",
               "std::lock" : "#include<mutex>",
               "std::shared_timed_mutex" : "#include<shared_mutex>",
               "std::condition_variable" : "#include<condition_variable>",
               "std::istream" : "#include<iostream>",
               "std::ostream" : "#include<iostream>",
               "std::cin" : "#include<iostream>",
               "std::cout" : "#include<iostream>",
               "std::istrstream" : "#include<sstream>",
               "std::ostrstream" : "#include<sstream>",
               "std::stringstream" : "#include<sstream>",
               "std::map" : "#include<map>",
               "std::multimap" : "#include<map>",
               "std::unordered_map" : "#include<unordered_map>",
               "std::unordered_multimap" : "#include<unordered_map>",
               "std::unique_ptr" : "#include<memory>",
               "std::shared_ptr" : "#include<memory>",
               "std::allocator" : "#include<memory>",
               "std::default_random_engine" : "#include<random>",
               "std::normal_distribution" : "#include<random>",
               "std::random_device" : "#include<random>",
               "std::std::mt19937_64" : "#include<random>",
               "std::std::mt19937" : "#include<random>",
               "std::regex" : "#include<regex>",
               "std::smatch" : "#include<regex>",
               "std::set" : "#include<set>",
               "std::multiset" : "#include<set>",
               "std::unordered_set" : "#include <unordered_set>",
               "std::unordered_multiset" : "#include<unordered_set>",
               "std::pair" : "#include <utility>",
               "std::move" : "#include<utility>",
               "std::swap" : "#include<utility>",
               "std::vector" : "#include<vector>"
               }

def parseForHeaders(args):
    headers = set()
    unknowns = []

    for arg in args:
        if arg in objToHeader:
            headers.add(objToHeader[arg])
        elif ("std::" + arg) in objToHeader:
            headers.add(objToHeader["std::" + arg])
        else:
            unknowns.append(arg)
    print(CT.boldRed("Unknowns:"))
    for u in unknowns:
        print(CT.boldBlack("unknown which header file contains"), CT.boldRed(u))

    if headers:
        print(CT.boldGreen("found headers:"))
        for h in sorted(headers):
            print(h)

def printBashMatches(args):
    if not args: # show all available function names
        for k in [x[5:] for x in sorted(objToHeader.keys())]:
            print(k)
        return

    # we are completing a partial function name
    partial = args[0]
    def matchForBashAutocomplete(x):
        return x.startswith("std::" + partial)
    matches = filter(matchForBashAutocomplete, objToHeader.keys())
    for x in sorted(matches):
        print(x[5:]) # skip "std::"

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bashAutoComplete', action="store_true", default=False)
    parser.add_argument('args', type=str, nargs='*')
    return parser.parse_args()

def main():
    args = parse_args()

    #warning(','.join(sys.argv))

    if args.bashAutoComplete:
        printBashMatches(args.args)
        sys.exit(1)

    parseForHeaders(args.args)

if __name__ == "__main__":
    main()
