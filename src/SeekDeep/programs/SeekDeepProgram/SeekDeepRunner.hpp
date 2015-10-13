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
//  SeekDeepRunner.hpp
//  SeekDeep
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
