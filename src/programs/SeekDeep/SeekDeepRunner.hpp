#pragma once
//
//  SeekDeepRunner.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/24/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//

//#include <cppprogutils/programRunner.hpp>
#include <bibcpp.h>
namespace bibseq {

class SeekDeepRunner : public bib::progutils::programRunner {

 public:
  SeekDeepRunner();


  static int extractor(MapStrStr inputCommands);
  static int qluster(MapStrStr inputCommands);
  static int processClusters(MapStrStr inputCommands);
  static int makeSampleDirectories(MapStrStr inputCommands);
};
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "SeekDeepRunner.cpp"
#endif
