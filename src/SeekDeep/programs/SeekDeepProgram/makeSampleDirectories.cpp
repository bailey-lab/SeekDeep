//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012-2016 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
//  main.cpp
//  SeekDeep
//
//  Created by Nicholas Hathaway on 8/11/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//
#include "SeekDeepRunner.hpp"


namespace bibseq {

int SeekDeepRunner::makeSampleDirectories(const bib::progutils::CmdArgs & inputCommands) {
	SeekDeepSetUp setUp(inputCommands);
	std::string sampleNameFilename = "";
	bool separatedDirs = false;
	setUp.setOption(separatedDirs, "--separatedDirs",
			"Create a separate directory for each index");
	setUp.setUpMakeSampleDirectories(sampleNameFilename);
	setUpSampleDirs(sampleNameFilename, setUp.pars_.directoryName_, separatedDirs);
	setUp.startARunLog(setUp.pars_.directoryName_);
	setUp.logRunTime(std::cout);
	return 0;
}

}  // namespace bib
