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
//  SeekDeepUtilsRunner.cpp
//
//  Created by Nick Hathaway on 2015/06/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

    
#include "SeekDeepUtilsRunner.hpp"

    
namespace bibseq {

SeekDeepUtilsRunner::SeekDeepUtilsRunner()
    : bib::progutils::programRunner({addFunc("dryRunQaulityFiltering", dryRunQaulityFiltering, false),
																		 addFunc("popClusteringViewer", popClusteringViewer, false)},
                    "SeekDeepUtils") {}
                    
int SeekDeepUtilsRunner::dryRunQaulityFiltering(MapStrStr inputCommands){
	seqSetUp setUp(inputCommands);
	setUp.processDefaultReader(true);
	std::string qualWindow;
	int32_t qualityWindowLength;
	int32_t qualityWindowStep;
	int32_t qualityWindowThres;
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
	readObjectIO reader;

	std::unordered_map<uint32_t, uint32_t> qualWindowCounts;
	uint32_t count = 0;
	uint32_t tenPer = 2000;
	std::vector<double> qualChecks;
	uint32_t failedQualCheck = 0;
	uint32_t failedQualWindow = 0;
	if(setUp.ioOptions_.inFormat_ != "fastq"){
		reader.read(setUp.ioOptions_);
		readVec::allSetQualCheck(reader.reads, qualCheck);
		for(const auto & read : reader.reads){
			if(count % tenPer == 0){
				std::cout << "Currently on " << count << " of " << len(reader.reads) << "\r";
				std::cout.flush();
			}
			qualChecks.emplace_back(read.fractionAboveQualCheck_);
			if(read.fractionAboveQualCheck_ < qualCheckCutOff){
				++failedQualCheck;
			}
			if(seqUtil::checkQualityWindow(qualityWindowLength, qualityWindowThres,
	        qualityWindowStep, read.seqBase_.qual_)){
				++failedQualWindow;
			}
			++count;
		}
	}else{
		readObject read;
		std::ifstream inFile(setUp.ioOptions_.firstName_);
		while(reader.readNextFastqStream(inFile, reader.SangerQualOffset, read, false)){
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
                    
} // namespace bibseq
