

//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012-2019 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include "SeekDeepPrograms/SeekDeepProgram/SeekDeepSetUp.hpp"
#include <njhcpp/bashUtils.h>

namespace njhseq {

void SeekDeepSetUp::setUpMakeSampleDirectories(
		makeSampleDirectoriesPars & pars) {
	description_ = "Make the directory tree processClusters reads from a table of sample names, one directory per sample with a subdirectory per replicate, "
			"plus locationByIndex/ with a file per index listing where each of its replicates went";
	examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --file names.tab.txt --dout clustering");
	processVerbose();
	setOption(pars.separatedDirs, "--separatedDirs",
				"Create a separate directory for each index, with the sample directories under it");
	setOption(pars.sampleNameFilename, "--file", "Tab (or whitespace) delimited sample names file, at least 3 columns: 1) index (or input file) name, 2) sample name, 3+) the replicate (e.g. MID) names for that sample, "
			"e.g. a row of 1 090-00 MID01 MID02 is sample 090-00 with replicates MID01 and MID02 from index 1, rows starting with # are skipped", true);
	setOption(pars_.directoryName_, "--dout", "Main Out Directory Name", true);
	setOption(pars_.overWriteDir_, "--overWriteDir",
			"If the directory already exists over write it");
	if (!failed_ && !needsHelp()) {
		std::string newDirectoryName = "./"
				+ njh::replaceString(njh::replaceString(pars_.directoryName_, "./", ""), "TODAY",
						getCurrentDate()) + "/";
		pars_.directoryName_ = njh::files::makeDir("./",
				njh::files::MkdirPar(njh::replaceString(pars_.directoryName_, "./", ""),
						pars_.overWriteDir_)).string();
	}
	finishSetUp(std::cout);
}

}  // namespace njhseq

