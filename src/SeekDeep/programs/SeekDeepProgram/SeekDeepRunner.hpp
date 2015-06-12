#pragma once
//
//  SeekDeepRunner.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/24/14.
//  Copyright (c) 2014 Nicholas Hathaway. All rights reserved.
//



#include <bibcpp.h>
#include "SeekDeep/programs/SeekDeepProgram/SeekDeepSetUp.hpp"

namespace bibseq {

class SeekDeepRunner : public bib::progutils::programRunner {

 public:
  SeekDeepRunner();


  static int sffExtractor(MapStrStr inputCommands);
  static int extractor(MapStrStr inputCommands);
  static int qluster(MapStrStr inputCommands);
  static int processClusters(MapStrStr inputCommands);
  static int makeSampleDirectories(MapStrStr inputCommands);
};
}  // namespace bibseq

#ifndef NOT_HEADER_ONLY
#include "SeekDeepRunner.cpp"
#endif
