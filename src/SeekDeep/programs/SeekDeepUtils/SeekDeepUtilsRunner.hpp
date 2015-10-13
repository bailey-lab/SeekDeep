#pragma once
//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012, 2015 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
//
// This file is part of SeekDeep.
//
// SeekDeep is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// SeekDeep is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with SeekDeep.  If not, see <http://www.gnu.org/licenses/>.
//
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
  
  static int dryRunQaulityFiltering(MapStrStr inputCommands);
	static int runMultipleCommands(MapStrStr inputCommands);

};

} // namespace bibseq
