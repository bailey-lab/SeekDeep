


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
#include "SeekDeepPrograms/SeekDeepProgram/SeekDeepSetUp.hpp"
#include <bibcpp/bashUtils.h>

namespace bibseq {

void SeekDeepSetUp::setUpMakeSampleDirectories(
		makeSampleDirectoriesPars & pars) {
	if (needsHelp()) {
		std::stringstream tempStream;
		tempStream << "makeSampleDirectoires" << std::endl;
		tempStream << "Set up a directory tree for processClusters" << std::endl;
		tempStream << "Commands, order not necessary, flags are case insensitive"
				<< std::endl;
		tempStream << "Required commands" << std::endl;
		tempStream << "-file [option], name of the file of sample names to read in"
				<< std::endl;
		tempStream << "-dout [option], name of the main directory to create"
				<< std::endl;
		tempStream << "File should be tab delimited and a few examples are below"
				<< std::endl;
		tempStream << "File should have at least three columns" << std::endl;
		tempStream << "Where first column is the name of the index or sff file "
				"used, second column is the sample names, and all following "
				"columns are the MIDs for that samples " << std::endl;
		std::cout << cleanOut(tempStream.str(), width_, indent_);
		tempStream.str(std::string());
		std::cout << "Example with two replicates and two separate master indexes"
				<< std::endl;
		std::cout << "1\t090-00\tMID01\tMID02" << std::endl;
		std::cout << "1\t090-24\tMID03\tMID04" << std::endl;
		std::cout << "1\t090-48\tMID05\tMID06" << std::endl;
		std::cout << "1\t090-72\tMID07\tMID08" << std::endl;
		std::cout << "2\t095-00\tMID01\tMID02" << std::endl;
		std::cout << "2\t095-24\tMID03\tMID04" << std::endl;
		std::cout << "2\t095-48\tMID05\tMID06" << std::endl;
		std::cout << "2\t095-72\tMID07\tMID08" << std::endl;
		std::cout << std::endl;
		std::cout << "Example with one replicate" << std::endl;
		std::cout << "1\t090-00\tMID01" << std::endl;
		std::cout << "1\t090-24\tMID02" << std::endl;
		std::cout << "2\t095-00\tMID01" << std::endl;
		std::cout << "2\t095-24\tMID02" << std::endl;
		std::cout << std::endl;
		std::cout << "Example with with mix amount of replicates" << std::endl;
		std::cout << "1\t090-00\tMID01" << std::endl;
		std::cout << "1\t090-24\tMID02" << std::endl;
		std::cout << "2\t095-00\tMID01\tMID02" << std::endl;
		std::cout << "2\t095-24\tMID03\tMID04" << std::endl;

		std::cout << "examples, SeekDeep makesampleDirectories -file "
				"names.tab.txt -dout clustering" << std::endl;
		exit(0);
	}
	processVerbose();
	setOption(pars.separatedDirs, "--separatedDirs",
				"Create a separate directory for each index");
	setOption(pars.sampleNameFilename, "--file", "Sample Names Filename", true);
	setOption(pars_.directoryName_, "-dout", "Main Out Directory Name", true);
	setOption(pars_.overWriteDir_, "--overWriteDir",
			"If the directory already exists over write it");
	if (!failed_) {
		std::string newDirectoryName = "./"
				+ bib::replaceString(bib::replaceString(pars_.directoryName_, "./", ""), "TODAY",
						getCurrentDate()) + "/";
		pars_.directoryName_ = bib::files::makeDir("./",
				bib::files::MkdirPar(bib::replaceString(pars_.directoryName_, "./", ""),
						pars_.overWriteDir_));
	}
	finishSetUp(std::cout);
}

}  // namespace bibseq

