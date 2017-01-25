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
	setUp.writeParametersFile(setUp.pars_.directoryName_ + "parametersUsed.txt",
			false, false);
	//write clustering parameters
	auto parsDir = bib::files::makeDir(setUp.pars_.directoryName_, bib::files::MkdirPar("pars"));
	std::ofstream parsOutFile;
	openTextFile(parsOutFile, OutOptions(bib::files::make_path(parsDir, "pars.tab.txt").string()));
	pars.iteratorMap.writePars(parsOutFile);
	std::ofstream popParsOutFile;
	openTextFile(popParsOutFile, OutOptions(bib::files::make_path(parsDir, "popPars.tab.txt").string()));
	pars.popIteratorMap.writePars(popParsOutFile);


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
			setUp.pars_.directoryName_,
			PopNamesInfo(pars.experimentName, samplesDirs),
			pars.clusterCutOff);

	if("" != pars.groupingsFile){
		sampColl.addGroupMetaData(pars.groupingsFile);
	}

	{
		bib::concurrent::LockableQueue<std::string> sampleQueue(samplesDirs);
		bibseq::concurrent::AlignerPool alnPool(alignerObj, pars.numThreads);
		alnPool.initAligners();
		alnPool.outAlnDir_ = setUp.pars_.outAlnInfoDirName_;

		auto setupClusterSamples = [&sampleQueue, &alnPool,&collapserObj,&pars, &setUp,&expectedSeqs,&sampColl](){
			std::string samp = "";
			auto currentAligner = alnPool.popAligner();
			while(sampleQueue.getVal(samp)){
				if(setUp.pars_.verbose_){
					std::cout << "Starting: " << samp << std::endl;
				}

				sampColl.setUpSample(samp, *currentAligner, collapserObj, setUp.pars_.chiOpts_);
				sampColl.clusterSample(samp, *currentAligner, collapserObj, pars.iteratorMap);

				//exclude clusters that don't have the necessary replicate number
				//defaults to the number of input replicates if none supplied
				if (0 != pars.runsRequired) {
					sampColl.sampleCollapses_.at(samp)->excludeBySampNum(pars.runsRequired, true);
				} else {
					sampColl.sampleCollapses_.at(samp)->excludeBySampNum(
							sampColl.sampleCollapses_.at(samp)->input_.info_.infos_.size(), true);
				}
				std::string sortBy = "fraction";
				sampColl.sampleCollapses_.at(samp)->renameClusters(sortBy);


				if (!expectedSeqs.empty()) {
					sampColl.sampleCollapses_.at(samp)->excluded_.checkAgainstExpected(
							expectedSeqs, *currentAligner, false);
					sampColl.sampleCollapses_.at(samp)->collapsed_.checkAgainstExpected(
							expectedSeqs, *currentAligner, false);
					if(setUp.pars_.debug_){
						std::cout << "sample: " << samp << std::endl;
					}
					for(const auto & clus : sampColl.sampleCollapses_.at(samp)->collapsed_.clusters_){
						if(setUp.pars_.debug_){
							std::cout << clus.seqBase_.name_ << " : " << clus.expectsString << std::endl;
						}
						if("" ==  clus.expectsString ){
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ": Error, expects string is blank" << std::endl;
							ss << clus.seqBase_.name_ << std::endl;
							throw std::runtime_error{ss.str()};
						}
					}
				}

				sampColl.dumpSample(samp);

				if(setUp.pars_.verbose_){
					std::cout << "Ending: " << samp << std::endl;
				}
			}
		};
		std::vector<std::thread> threads;
		for(uint32_t t = 0; t < pars.numThreads; ++t){
			threads.emplace_back(std::thread(setupClusterSamples));
		}
		for(auto & t : threads){
			t.join();
		}
	}

	//read in the dump alignment cache
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);

	if (pars.investigateChimeras) {
		sampColl.investigateChimeras(pars.chiCutOff, alignerObj);
	}

	//exclude
	for(const auto & sampleName : samplesDirs){

		sampColl.setUpSampleFromPrevious(sampleName);
		if(setUp.pars_.debug_){
			std::cout << "sample: " << sampleName << std::endl;
		}
		for(const auto & clus : sampColl.sampleCollapses_.at(sampleName)->collapsed_.clusters_){
			if(setUp.pars_.debug_){
				std::cout << clus.seqBase_.name_ << " : " << clus.expectsString << std::endl;
			}
		}
		if (!pars.keepChimeras) {
			//now exclude all marked chimeras, currently this will also remark chimeras unnecessarily
			sampColl.sampleCollapses_.at(sampleName)->excludeChimeras(false, pars.chiCutOff);
		}
		if (bib::in(sampleName, customCutOffsMap)) {
			if (setUp.pars_.debug_) {
				std::cout << "Custom Cut off for " << sampleName << " : "
						<< customCutOffsMap[sampleName] << std::endl;
			}
			sampColl.sampleCollapses_.at(sampleName)->excludeFraction(customCutOffsMap[sampleName],
					true);
		} else {
			sampColl.sampleCollapses_.at(sampleName)->excludeFraction(pars.fracCutoff, true);
		}

		std::string sortBy = "fraction";
		sampColl.sampleCollapses_.at(sampleName)->renameClusters(sortBy);
		sampColl.dumpSample(sampleName);
	}
	if(setUp.pars_.verbose_){
		std::cout << bib::bashCT::boldGreen("Pop Clustering") << std::endl;
	}
	if(!pars.noPopulation){
		sampColl.doPopulationClustering(sampColl.createPopInput(),
				alignerObj, collapserObj, pars.popIteratorMap);
	}
	if(setUp.pars_.verbose_){
		std::cout << bib::bashCT::boldRed("Done Pop Clustering") << std::endl;
	}
	if ("" != pars.previousPopFilename && !pars.noPopulation) {
		sampColl.renamePopWithSeqs(getSeqs<readObject>(pars.previousPopFilename), pars.previousPopErrors);
	}

	if (!expectedSeqs.empty()) {
		sampColl.comparePopToRefSeqs(expectedSeqs, alignerObj);
	}

	sampColl.printSampleCollapseInfo(
			bib::files::join(sampColl.masterOutputDir_.string(),
					"selectedClustersInfo.tab.txt"));

	sampColl.symlinkInSampleFinals();
	sampColl.outputRepAgreementInfo();
	if(!pars.noPopulation){
		sampColl.dumpPopulation();
	}
	if("" != pars.groupingsFile){
		sampColl.createGroupInfoFiles();
	}

	sampColl.createCoreJsonFile();

	//collect extraction dirs
	std::set<bfs::path> extractionDirs;
	for(const auto & file : analysisFiles){
		auto metaDataJsonFnp = bib::files::make_path(file.first.parent_path(), "metaData.json");
		if(bfs::exists(metaDataJsonFnp)){
			auto metaJson = bib::json::parseFile(metaDataJsonFnp.string());
			if(metaJson.isMember("extractionDir")){
				extractionDirs.emplace(metaJson["extractionDir"].asString());
			}
		}
	}
	if(setUp.pars_.verbose_){
		std::cout << "Extraction Dirs" << std::endl;
		std::cout << bib::conToStr(extractionDirs, "\n") << std::endl;
	}
	table profileTab;
	table statsTab;
	for(const auto & extractDir : extractionDirs){
		auto profileFnp = bib::files::make_path(extractDir, "extractionProfile.tab.txt");
		auto statsFnp = bib::files::make_path(extractDir, "extractionStats.tab.txt");
		if(bfs::exists(profileFnp)){
			table currentProfileTab(profileFnp.string(), "\t", true);
			auto minLenPos = getPositionsOfTargetStartsWith(currentProfileTab.columnNames_, "len<");
			auto maxLenPos = getPositionsOfTargetStartsWith(currentProfileTab.columnNames_, "len>");
			currentProfileTab.columnNames_[minLenPos.front()] = "minlen";
			currentProfileTab.columnNames_[maxLenPos.front()] = "maxlen";
			auto oldColumnNames = currentProfileTab.columnNames_;
			currentProfileTab.addColumn(VecStr{extractDir.filename().string()}, "extractionDir");
			currentProfileTab = currentProfileTab.getColumns(catenateVectors(VecStr{"extractionDir"}, oldColumnNames));

			if(profileTab.empty()){
				profileTab = currentProfileTab;
			}else{
				profileTab.rbind(currentProfileTab, false);
			}
		}
		if(bfs::exists(statsFnp)){
			table curentStatsTab(statsFnp.string(), "\t", true);
			auto oldColumnNames = curentStatsTab.columnNames_;
			curentStatsTab.addColumn(VecStr{extractDir.filename().string()}, "extractionDir");
			curentStatsTab = curentStatsTab.getColumns(catenateVectors(VecStr{"extractionDir"}, oldColumnNames));
			if(statsTab.empty()){
				statsTab = curentStatsTab;
			}else{
				statsTab.rbind(curentStatsTab, false);
			}
		}
	}

	auto extractionOutputDir = bib::files::make_path(setUp.pars_.directoryName_,
			"extractionInfo");
	bib::files::makeDirP(bib::files::MkdirPar(extractionOutputDir.string()));
	if (!profileTab.empty()) {
		profileTab.sortTable("extractionDir", false);
		auto profileTabOpts =
				TableIOOpts::genTabFileOut(
						bib::files::make_path(extractionOutputDir,
								"extractionProfile.tab.txt").string(), true);
		profileTabOpts.out_.overWriteFile_ = true;
		profileTab.outPutContents(profileTabOpts);
	}
	if (!statsTab.empty()) {
		auto statsTabOpts =
				TableIOOpts::genTabFileOut(
						bib::files::make_path(extractionOutputDir,
								"extractionStats.tab.txt").string(), true);
		statsTabOpts.out_.overWriteFile_ = true;
		statsTab.sortTable("extractionDir", false);
		statsTab.outPutContents(statsTabOpts);
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
