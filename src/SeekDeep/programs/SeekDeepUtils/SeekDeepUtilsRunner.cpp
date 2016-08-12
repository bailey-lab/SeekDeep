
//  SeekDeepUtilsRunner.cpp
//
//  Created by Nick Hathaway on 2015/06/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//
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
    
#include "SeekDeepUtilsRunner.hpp"

    
namespace bibseq {

SeekDeepUtilsRunner::SeekDeepUtilsRunner()
    : bib::progutils::programRunner(
    		{
				addFunc("dryRunQaulityFiltering", dryRunQaulityFiltering, false),
				addFunc("popClusteringViewer", popClusteringViewer, false),
				addFunc("genProjectConfig", genProjectConfig, false),
				addFunc("runMultipleCommands", runMultipleCommands, false)
    		},
                    "SeekDeepUtils") {}
                    
int SeekDeepUtilsRunner::dryRunQaulityFiltering(const bib::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	setUp.processDefaultReader(true);
	std::string qualWindow;
	uint32_t qualityWindowLength;
	uint32_t qualityWindowStep;
	uint32_t qualityWindowThres;
  if (setUp.setOption(qualWindow, "-qualWindow", "SlidingQualityWindow")) {
    seqUtil::processQualityWindowString(qualWindow, qualityWindowLength,
                                        qualityWindowStep, qualityWindowThres);
  } else {
    qualityWindowLength = 50;
    qualityWindowStep = 5;
    qualityWindowThres = 25;
  }
  uint32_t qualCheck = 30;
  setUp.setOption(qualCheck, "-qualCheck", "Qual Check Level");
  double qualCheckCutOff = 0.90;
  setUp.setOption(qualCheckCutOff, "-qualCheckCutOff", "Cut Off for fraction of bases above qual check"	);
  bool plot = false;
  setUp.setOption(plot, "--plot", "Plot");
	setUp.finishSetUp(std::cout);
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	std::unordered_map<uint32_t, uint32_t> qualWindowCounts;
	uint32_t count = 0;
	std::vector<double> qualChecks;
	uint32_t failedQualCheck = 0;
	uint32_t failedQualWindow = 0;
	readObject read;
	while(reader.readNextRead(read)){
		read.setBaseCountOnQualCheck(qualCheck);
		std::cout << "Currently on " << count << "\r";
		std::cout.flush();
		qualChecks.emplace_back(read.fractionAboveQualCheck_);
		if(read.fractionAboveQualCheck_ < qualCheckCutOff){
			++failedQualCheck;
		}
		if(!seqUtil::checkQualityWindow(qualityWindowLength, qualityWindowThres,
        qualityWindowStep, read.seqBase_.qual_)){
			++failedQualWindow;
		}
		++count;
	}
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << bib::bashCT::boldBlack("Sliding Quality Window") << std::endl;
	std::cout << "Window Size: "  << qualityWindowLength << ", Window Step: " << qualityWindowStep << ", Window Thresdhold: " << qualityWindowThres << std::endl;
	std::cout << "FailedWindow: " << getPercentageString(failedQualWindow, count) << std::endl;
	std::cout << std::endl;
	std::cout << bib::bashCT::boldBlack("Quality Fraction above cut off") << std::endl;
	std::cout << "Q"<< qualCheck << ">" << qualCheckCutOff << std::endl;
	std::cout << "FailedCheck: " << getPercentageString(failedQualCheck, count) << std::endl;
	auto qualCheckStats = getStatsOnVec(qualChecks);
	table qualChecksTab(qualCheckStats, VecStr{"stat", "value"});
	qualChecksTab.outPutContentOrganized(std::cout);
	return 0;
}



int SeekDeepUtilsRunner::runMultipleCommands(
		const bib::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	std::string filename = "";
	std::string logFile = "";
	uint32_t numThreads = 1;
	bool raw = false;
	setUp.setOption(numThreads, "--numThreads", "Number of threads to use");
	setUp.setOption(logFile, "--logFile",
			"Name of a file to log the output of the commands");
	setUp.setOption(filename, "--cmdFile",
			"Name of the file, first line is command with REPLACETHIS, the next lines are the cmd to run with that line replacing REPLACETHIS",
			true);
	setUp.setOption(raw, "--raw,-r",
			"If if the file is simply just a list of commands to run");
	setUp.processVerbose();
	setUp.processWritingOptions();
	setUp.processDebug();
	if (setUp.needsHelp()) {
		std::cout
				<< "Input cmdFile should start with CMD: and then command and "
						"each line after that should be some sort of replacement to run the commnd"
				<< std::endl;
		std::cout << "Example" << std::endl;
		std::cout << "CMD:echo hello REPLACETHIS1 from REPLACETHIS2" << std::endl;
		std::cout << "REPLACETHIS1:nick,jon,mike" << std::endl;
		std::cout << "REPLACETHIS2:world,everyone" << std::endl;
		setUp.printFlags(std::cout);
		exit(1);
	}
	if (setUp.commands_.hasFlagCaseInsenNoDash("--gen")) {
		std::cout << "CMD:echo hello REPLACETHIS1 from REPLACETHIS2" << std::endl;
		std::cout << "REPLACETHIS1:nick,jon,mike" << std::endl;
		std::cout << "REPLACETHIS2:world,everyone" << std::endl;
		exit(1);
	}
	setUp.finishSetUp(std::cout);
	if (logFile == "") {
		logFile = bfs::path(filename).filename().replace_extension("").string() + "Log.json";
	}
	std::ofstream outFile;
	openTextFile(outFile, logFile, ".json", setUp.pars_.ioOptions_.out_);
	std::ifstream inFile(filename);
	if (!inFile) {
		std::stringstream ss;
		ss << bib::bashCT::bold << "Error in opening "
				<< bib::bashCT::boldRed(filename) << std::endl;
		throw std::runtime_error { ss.str() };
	}
	std::string cmd;
	VecStr cmds;
	std::string line;
	if (raw) {
		while (bib::files::crossPlatGetline(inFile, line)) {
			if ("" == line || allWhiteSpaceStr(line) || bib::beginsWith(line, "#")) {
				continue;
			}
			cmds.emplace_back(line);
		}
	} else {
		std::map<std::string, VecStr> replacements;
		while (bib::files::crossPlatGetline(inFile, line)) {
			if ("" == line || allWhiteSpaceStr(line) || bib::beginsWith(line, "#")) {
				continue;
			}
			auto colonPos = line.find(":");

			if (std::string::npos == colonPos) {
				std::stringstream ss;
				ss << "Error in processing line: " << bib::bashCT::boldRed(line) << std::endl;
				ss << "Need at least one colon" << std::endl;
				throw std::runtime_error{ss.str()};
			}
			VecStr toks{line.substr(0,colonPos), line.substr(colonPos + 1)};
			if (toks.front() == "CMD") {
				cmd = toks.back();
			} else {
				replacements[toks.front()] = tokenizeString(toks.back(), ",");
			}
		}
		for (const auto & r : replacements) {
			if (cmds.empty()) {
				for (const auto & subR : r.second) {
					cmds.emplace_back(bib::replaceString(cmd, r.first, subR));
				}
			} else {
				VecStr newCmds;
				for (const auto & c : cmds) {
					for (const auto & subR : r.second) {
						newCmds.emplace_back(bib::replaceString(c, r.first, subR));
					}
				}
				cmds = newCmds;
			}
		}
	}

	if (setUp.pars_.verbose_) {
		printVector(cmds, "\n");
	}

	auto allRunOutputs = bib::sys::runCmdsThreaded(cmds, numThreads,
			setUp.pars_.verbose_, setUp.pars_.debug_);
	Json::Value allLog;

	allLog["totalTime"] = setUp.timer_.totalTime();
	allLog["cmdsfile"] = bib::json::toJson(bfs::absolute(filename));
	auto & cmdsLog = allLog["cmdsLog"];
	for (const auto & out : allRunOutputs) {
		cmdsLog.append(out.toJson());
	}
	outFile << cmdsLog << std::endl;
	return 0;
}
} // namespace bibseq
