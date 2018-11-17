#!/usr/bin/env python2

import shutil, os, argparse, sys, stat,time


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-name', type=str, required = True)
    parser.add_argument('-author', type=str, required = True)
    parser.add_argument('-prepath', type=str, required = True)
    return parser.parse_args()


def fileInfoHeader(headerName, author):
    
    return """
//  {headerName}
//
//  Created by {author} on {date}.
//  Copyright (c) {year} {author}. All rights reserved.
//
""".format(headerName=headerName, author=author,year=time.strftime("%Y"),date=time.strftime("%Y/%m/%d"))

def startHeader(headerName, author):
    
    return """#pragma once
//
""" + fileInfoHeader(headerName, author)

def startCpp(nameStub, author):
    return fileInfoHeader(nameStub + ".cpp", author) + """
    
#include "{name}.hpp"
    
    """.format(name = nameStub)



def main():
    args = parse_args()
    name = args.name
    prepath = args.prepath
    if not prepath.endswith("/"):
        prepath += "/"
    if os.path.exists(name) or os.path.exists(name + ".h"):
        print "Error, " + name + " already exists"
        exit(1)
    #create main dir
    os.mkdir(name)
    #create main header to include the ring
    with open(name + ".h", "w") as f:
        mainHeaderOut = startHeader( name + ".h", args.author) + """
#include "{prepath}{name}/{name}SetUp.hpp"
#include "{prepath}{name}/{name}Runner.hpp"
""".format(name=args.name,prepath = prepath)
        f.write(mainHeaderOut)
    #create setUp header
    with open(os.path.join(name,name + "SetUp.hpp"), "w") as f:
        defaultHeader = startHeader(name + "SetUp.hpp", args.author)
        defaultHeader += """
#include <njhseq.h>
#include <njhseq/programUtils/seqSetUp.hpp>
#include <njhcpp.h>
namespace njhseq {{

class {name}SetUp : public seqSetUp {{

 public:
    using seqSetUp::seqSetUp;
}};
}} // namespace njhseq
""".format(name =name)
        f.write(defaultHeader)
    #create setUp cpp
    with open(os.path.join(name,name + "SetUp.cpp"), "w") as f:
        infoHeader = startCpp(name + "SetUp", args.author)
        infoHeader +="""
namespace njhseq {

} // namespace njhseq
"""
        f.write(infoHeader)
    #create runner header
    with open(os.path.join(name,name + "Runner.hpp"), "w") as f:
        infoHeader = startHeader(name + "Runner.hpp", args.author)
        infoHeader +="""
#include "{name}SetUp.hpp"

namespace njhseq {{

class {name}Runner : public njh::progutils::programRunner {{
 public:
  {name}Runner();
  
  static int placeHolder(MapStrStr inputCommands);

}};
}} // namespace njhseq
""".format(name = name)
        f.write(infoHeader)
    #create runner cpp
    with open(os.path.join(name,name + "Runner.cpp"), "w") as f:
        infoHeader = startCpp(name + "Runner", args.author)
        infoHeader +="""
namespace njhseq {{

{name}Runner::{name}Runner()
    : njh::progutils::programRunner({{addFunc("placeHolder", placeHolder, false)}},
                    "{name}") {{}}
                    
int {name}Runner::placeHolder(MapStrStr inputCommands) {{
  {name}SetUp setUp(inputCommands);
  setUp.finishSetUp(std::cout);
  return 0;
}}
                    
}} // namespace njhseq
""".format(name = name)
        f.write(infoHeader)
    
    
main()
