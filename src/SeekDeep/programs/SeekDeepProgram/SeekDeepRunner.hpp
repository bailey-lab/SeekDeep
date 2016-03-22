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

class SeekDeepRunner : public bib::progutils::oneRing {

 public:
  SeekDeepRunner();
  static int extractor(const bib::progutils::CmdArgs & inputCommands);
  static int qluster(const bib::progutils::CmdArgs & inputCommands);
  static int processClusters(const bib::progutils::CmdArgs & inputCommands);
  static int makeSampleDirectories(const bib::progutils::CmdArgs & inputCommands);
};
}  // namespace bibseq

#ifndef NOT_HEADER_ONLY
#include "SeekDeepRunner.cpp"
#endif
