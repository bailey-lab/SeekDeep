#pragma once
//

//  SeekDeepUtilsRunner.hpp
//
//  Created by Nick Hathaway on 2015/06/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

#include "SeekDeepUtilsSetUp.hpp"
#include "SeekDeep/server.h"

namespace bibseq {

class SeekDeepUtilsRunner : public bib::progutils::programRunner {
 public:
  SeekDeepUtilsRunner();
  
  static int dryRunQaulityFiltering(const bib::progutils::CmdArgs & inputCommands);
	static int runMultipleCommands(const bib::progutils::CmdArgs & inputCommands);

};

} // namespace bibseq
