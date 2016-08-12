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


std::map<std::string, double> processCustomCutOffs(const bfs::path & customCutOffsFnp){
	table customCutOffsTab;
	std::map<std::string, double> ret;
	if ("" != customCutOffsFnp.string()) {
		customCutOffsTab = table(customCutOffsFnp.string(), "whitespace", true);
		if (!bib::in(std::string("sample"), customCutOffsTab.columnNames_)
				|| !bib::in(std::string("cutOff"), customCutOffsTab.columnNames_)) {
			std::stringstream ss;
			ss << "Error in loading custom cut off file, "
					<<  customCutOffsFnp << ", missing sample or cutOff\n"
					<< vectorToString(customCutOffsTab.columnNames_, ",");
			throw std::runtime_error { ss.str() };
		}
		for (const auto & rowPos : iter::range(customCutOffsTab.content_.size())) {
			ret[customCutOffsTab.content_[rowPos][customCutOffsTab.getColPos(
					"sample")]] =
					bib::lexical_cast<double>(
							customCutOffsTab.content_[rowPos][customCutOffsTab.getColPos(
									"cutOff")]);
		}
	}
	return ret;
}




int SeekDeepRunner::processClusters(const bib::progutils::CmdArgs & inputCommands) {
	// parameters
	SeekDeepSetUp setUp(inputCommands);
	processClustersPars pars;
	setUp.setUpMultipleSampleCluster(pars);
	// start a run log
	setUp.startARunLog(setUp.pars_.directoryName_);
	// parameters file
	setUp.writeParametersFile(setUp.pars_.directoryName_ + "parametersUsed.txt", false,
			false);
	//process custom cut offs
	std::map<std::string, double> customCutOffsMap = processCustomCutOffs(pars.customCutOffs);
	//read in the files in the corresponding sample directories

	auto analysisFiles = bib::files::listAllFiles(pars.masterDir, true,
			{ std::regex { "^" + setUp.pars_.ioOptions_.firstName_ + "$" } }, 3);

	std::set<std::string> samplesDirsSet;
	for (const auto & af : analysisFiles) {
		auto fileToks = bib::tokenizeString(bfs::relative(af.first, pars.masterDir).string(), "/");
		if (3 != fileToks.size()) {
			std::stringstream ss;
			ss << "File path should be three levels deep, not " << fileToks.size()
					<< " for " << bfs::relative(af.first, pars.masterDir).string() << std::endl;
			throw std::runtime_error { ss.str() };
		}
		samplesDirsSet.insert(fileToks[0]);
	}
	VecStr samplesDirs(samplesDirsSet.begin(), samplesDirsSet.end());
	VecStr specificFiles;

	for (const auto& fileIter : analysisFiles) {
		specificFiles.push_back(fileIter.first.string());
	}
	if(setUp.pars_.verbose_){
		std::cout << "Reading from" << std::endl;
		for (const auto& sfIter : specificFiles) {
			std::cout << sfIter << std::endl;
		}
	}
	// reading in reads
	uint64_t maxSize = 0;
	// reading expected sequences to compare to
	bool checkingExpected = setUp.pars_.refIoOptions_.firstName_ != "";
	std::vector<readObject> expectedSeqs;
	if (checkingExpected) {
		expectedSeqs = SeqInput::getReferenceSeq(setUp.pars_.refIoOptions_, maxSize);
	}
	// get max size for aligner
	for (const auto& sf : specificFiles) {
		SeqIOOptions inOpts(sf, setUp.pars_.ioOptions_.inFormat_, true);
		SeqInput reader(inOpts);
		reader.openIn();
		seqInfo seq;
		while(reader.readNextRead(seq)){
			readVec::getMaxLength(seq, maxSize);
		}
	}
	// create aligner class object
	aligner alignerObj(maxSize, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
			KmerMaps(setUp.pars_.colOpts_.kmerOpts_.kLength_),
			setUp.pars_.qScorePars_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
	// create collapserObj used for clustering
	collapser collapserObj(setUp.pars_.colOpts_);
	collapserObj.opts_.kmerOpts_.checkKmers_ = false;
	// output info about the read In reads

	collapse::SampleCollapseCollection sampColl(setUp.pars_.ioOptions_, pars.masterDir,
			setUp.pars_.directoryName_, PopNamesInfo(pars.experimentName,samplesDirs),  pars.clusterCutOff);

	if("" != pars.groupingsFile){
		sampColl.addGroupMetaData(pars.groupingsFile);
	}

	for (const auto & samp : samplesDirs) {
		if(setUp.pars_.verbose_){
			std::cout << "Starting: " << samp << std::endl;
		}
		sampColl.setUpSample(samp, alignerObj, collapserObj, setUp.pars_.chiOpts_);
		sampColl.clusterSample(samp, alignerObj, collapserObj, pars.iteratorMap);
		//exclude clusters that don't have the necessary replicate number
		//defaults to the number of input replicates if none supplied
		if (0 != pars.runsRequired) {
			sampColl.sampleCollapses_[samp]->excludeBySampNum(pars.runsRequired, false);
		} else {
			sampColl.sampleCollapses_[samp]->excludeBySampNum(
					sampColl.sampleCollapses_[samp]->input_.info_.infos_.size(), false);
		}
		if (!expectedSeqs.empty()) {
			sampColl.sampleCollapses_[samp]->collapsed_.checkAgainstExpected(
					expectedSeqs, alignerObj, false, false);
		}
		sampColl.dumpSample(samp);
		if(setUp.pars_.verbose_){
			std::cout << "Ending: " << samp << std::endl;
		}
	}

	if (pars.investigateChimeras) {
		sampColl.investigateChimeras(pars.chiCutOff, alignerObj);
	}

	//exclude
	for(const auto & sampleName : samplesDirs){
		sampColl.setUpSampleFromPrevious(sampleName);
		if (!pars.keepChimeras) {
			//now exclude all marked chimeras, currently this will also remark chimeras unnecessarily
			sampColl.sampleCollapses_[sampleName]->excludeChimeras(false, pars.chiCutOff);
		}
		if (bib::in(sampleName, customCutOffsMap)) {
			if (setUp.pars_.debug_) {
				std::cout << "Custom Cut off for " << sampleName << " : "
						<< customCutOffsMap[sampleName] << std::endl;
			}
			sampColl.sampleCollapses_[sampleName]->excludeFraction(customCutOffsMap[sampleName],
					true);
		} else {
			sampColl.sampleCollapses_[sampleName]->excludeFraction(pars.fracCutoff, true);
		}

		std::string sortBy = "fraction";
		sampColl.sampleCollapses_[sampleName]->renameClusters(sortBy);
		sampColl.dumpSample(sampleName);
	}
	if(setUp.pars_.verbose_){
		std::cout << bib::bashCT::boldGreen("Pop Clustering") << std::endl;
	}
	sampColl.doPopulationClustering(sampColl.createPopInput(),
			alignerObj, collapserObj, pars.popIteratorMap);

	if ("" != pars.previousPopFilename) {
		sampColl.renamePopWithSeqs(getSeqs<readObject>(pars.previousPopFilename));
	}

	if (!expectedSeqs.empty()) {
		sampColl.comparePopToRefSeqs(expectedSeqs, alignerObj);
	}

	sampColl.printSampleCollapseInfo(
			bib::files::join(sampColl.masterOutputDir_.string(),
					"selectedClustersInfo.tab.txt"));

	sampColl.symlinkInSampleFinals();
	sampColl.outputRepAgreementInfo();
	sampColl.dumpPopulation();

	if("" != pars.groupingsFile){
		sampColl.createGroupInfoFiles();
	}
	alignerObj.processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_,
			setUp.pars_.verbose_);
	setUp.rLog_ << alignerObj.numberOfAlingmentsDone_ << "\n";
	if (setUp.pars_.verbose_) {
		std::cout << alignerObj.numberOfAlingmentsDone_ << std::endl;
		setUp.logRunTime(std::cout);
	}
	return 0;
}

}  // namespace bibseq
