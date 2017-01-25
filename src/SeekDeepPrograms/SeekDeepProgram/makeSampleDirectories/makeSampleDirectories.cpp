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
#include "SeekDeepPrograms/SeekDeepProgram/SeekDeepRunner.hpp"


namespace bibseq {

int SeekDeepRunner::makeSampleDirectories(const bib::progutils::CmdArgs & inputCommands) {
	SeekDeepSetUp setUp(inputCommands);
	makeSampleDirectoriesPars pars;

	setUp.setUpMakeSampleDirectories(pars);
	setUpSampleDirs(pars.sampleNameFilename.string(), setUp.pars_.directoryName_, pars.separatedDirs);
	auto inputInfoDirPath = bib::files::make_path(setUp.pars_.directoryName_, "inputInfo");
	bib::files::makeDir(bib::files::MkdirPar(inputInfoDirPath.string()));
	bfs::copy(pars.sampleNameFilename, bib::files::make_path(inputInfoDirPath, "inputSampleNames.tab.txt"));
	setUp.startARunLog(setUp.pars_.directoryName_);
	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}
	return 0;
}

}  // namespace bibseq
