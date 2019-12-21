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
//  main.cpp
//  SeekDeep
//
//  Created by Nicholas Hathaway on 8/11/13.
//
#include "SeekDeepPrograms/SeekDeepProgram/SeekDeepRunner.hpp"



namespace njhseq {


std::unordered_map<std::string, double> processCustomCutOffs(const bfs::path & customCutOffsFnp, const VecStr & allSamples, double defaultFracCutOff){
	std::unordered_map<std::string, double> ret;
	if ("" != customCutOffsFnp.string()) {
		table customCutOffsTab(customCutOffsFnp.string(), "\t", true);
		customCutOffsTab.checkForColumnsThrow(VecStr{"sample", "cutOff"}, __PRETTY_FUNCTION__);
		for (const auto & rowPos : iter::range(customCutOffsTab.content_.size())) {
			ret[customCutOffsTab.content_[rowPos][customCutOffsTab.getColPos(
					"sample")]] =
					njh::lexical_cast<double>(
							customCutOffsTab.content_[rowPos][customCutOffsTab.getColPos(
									"cutOff")]);
		}
	}
	for(const auto & samp : allSamples){
		if(!njh::in(samp, ret)){
			ret[samp] = defaultFracCutOff;
		}
	}
	return ret;
}




int SeekDeepRunner::processClusters(const njh::progutils::CmdArgs & inputCommands) {
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
	auto parsDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("pars"));
	std::ofstream parsOutFile;
	openTextFile(parsOutFile, OutOptions(njh::files::make_path(parsDir, "pars.tab.txt")));
	pars.iteratorMap.writePars(parsOutFile);
	std::ofstream popParsOutFile;
	openTextFile(popParsOutFile, OutOptions(njh::files::make_path(parsDir, "popPars.tab.txt")));
	pars.popIteratorMap.writePars(popParsOutFile);

	std::unique_ptr<TranslatorByAlignment> translator;
	if("" != pars.transPars.gffFnp_){
		translator = std::make_unique<TranslatorByAlignment>(pars.transPars);
		if("" == pars.transPars.lzPars_.genomeFnp){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << " if supplying gff file, must also supply --genomeFnp"<< "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	table knownAminoAcidChanges;
	if("" != pars.knownAminoAcidChangesFnp){
		if("" == pars.transPars.gffFnp_){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "if supplying known amino acid positions than must also supply --gffFnp file"<< "\n";
			throw std::runtime_error{ss.str()};
		}
		knownAminoAcidChanges = table(pars.knownAminoAcidChangesFnp, "\t", true);
		njh::for_each(knownAminoAcidChanges.columnNames_, [](std::string & col){
			njh::strToLower(col);
		});
		knownAminoAcidChanges.setColNamePositions();
		knownAminoAcidChanges.checkForColumnsThrow(VecStr{"transcriptid", "aaposition"}, __PRETTY_FUNCTION__);
	}
	std::unordered_map<std::string, std::set<uint32_t>> knownMutationsLocationsMap;
	if(knownAminoAcidChanges.nRow() > 0){
		for(const auto & row : knownAminoAcidChanges){
			if(std::all_of(row.begin(), row.end(), [](const std::string & element){
				return "" ==element;
			})){
				continue;
			}
			knownMutationsLocationsMap[row[knownAminoAcidChanges.getColPos("transcriptid")]].emplace(njh::StrToNumConverter::stoToNum<uint32_t>(row[knownAminoAcidChanges.getColPos("aaposition")]));
		}
	}

	//population seqs;
	std::vector<seqInfo> globalPopSeqs;
	if("" != pars.popSeqsFnp){
		globalPopSeqs = SeqInput::getSeqVec<seqInfo>(SeqIOOptions::genFastaIn(pars.popSeqsFnp));
	}

	//read in the files in the corresponding sample directories
	auto analysisFiles = njh::files::listAllFiles(pars.masterDir, true,
			{ std::regex { "^" + setUp.pars_.ioOptions_.firstName_.string() + "$" } }, 2);

	std::set<std::string> samplesDirsSet;
	for (const auto & af : analysisFiles) {
		auto fileToks = njh::tokenizeString(bfs::relative(af.first, pars.masterDir).string(), "/");
		if (3 != fileToks.size()) {
			std::stringstream ss;
			ss << "File path should be three levels deep, not " << fileToks.size()
					<< " for " << bfs::relative(af.first, pars.masterDir).string() << std::endl;
			throw std::runtime_error { ss.str() };
		}
		if(njh::in(fileToks[0], pars.excludeSamples)){
			continue;
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
			pars.preFiltCutOffs);



	if("" != pars.groupingsFile){
		sampColl.addGroupMetaData(pars.groupingsFile);
	}

	//process custom cut offs
	std::unordered_map<std::string, double> customCutOffsMap = processCustomCutOffs(pars.customCutOffs, samplesDirs, pars.fracCutoff);
	std::unordered_map<std::string, double> customCutOffsMapPerRep = processCustomCutOffs(pars.customCutOffs, samplesDirs, pars.withinReplicateFracCutOff);

	{
		njh::concurrent::LockableQueue<std::string> sampleQueue(samplesDirs);
		njhseq::concurrent::AlignerPool alnPool(alignerObj, pars.numThreads);
		alnPool.initAligners();
		alnPool.outAlnDir_ = setUp.pars_.outAlnInfoDirName_;

		auto setupClusterSamples = [&sampleQueue, &alnPool,&collapserObj,&pars,&setUp,
																&expectedSeqs,&sampColl,&customCutOffsMap,
																&customCutOffsMapPerRep](){
			std::string samp = "";
			auto currentAligner = alnPool.popAligner();
			while(sampleQueue.getVal(samp)){
				if(setUp.pars_.verbose_){
					std::cout << "Starting: " << samp << std::endl;
				}

				sampColl.setUpSample(samp, *currentAligner, collapserObj, setUp.pars_.chiOpts_);

				sampColl.clusterSample(samp, *currentAligner, collapserObj, pars.iteratorMap);

				sampColl.sampleCollapses_.at(samp)->markChimeras(pars.chiCutOff);

				//exclude clusters that don't have the necessary replicate number
				//defaults to the number of input replicates if none supplied

				if (0 != pars.runsRequired) {
					sampColl.sampleCollapses_.at(samp)->excludeBySampNum(pars.runsRequired, true);

				} else {
					sampColl.sampleCollapses_.at(samp)->excludeBySampNum(sampColl.sampleCollapses_.at(samp)->input_.info_.infos_.size(), true);
				}

				if(pars.collapseLowFreqOneOffs){
					bool skipChimeras = true;
					sampColl.sampleCollapses_.at(samp)->excludeLowFreqOneOffs(true, pars.lowFreqMultiplier, *currentAligner,
							skipChimeras, customCutOffsMap.at(samp));
				}

				sampColl.sampleCollapses_.at(samp)->excludeFractionAnyRep(customCutOffsMapPerRep.at(samp), true);
				sampColl.sampleCollapses_.at(samp)->excludeFraction(customCutOffsMap.at(samp), true);

				if (!pars.keepChimeras) {
					//now exclude all marked chimeras
					sampColl.sampleCollapses_.at(samp)->excludeChimerasNoReMark(true);
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

//	if (pars.investigateChimeras) {
//		sampColl.investigateChimeras(pars.chiCutOff, alignerObj);
//	}




	if(setUp.pars_.verbose_){
		std::cout << njh::bashCT::boldGreen("Pop Clustering") << std::endl;
	}

	sampColl.doPopulationClustering(sampColl.createPopInput(),
			alignerObj, collapserObj, pars.popIteratorMap);

	if(pars.rescueExcludedChimericHaplotypes || pars.rescueExcludedOneOffLowFreqHaplotypes || pars.rescueExcludedLowFreqHaplotypes){
		//first gather major haplotypes
		std::set<std::string> majorHaps;
		std::set<std::string> majorHapsForChi;

		for(const auto & sampleName : sampColl.passingSamples_){
			if(njh::in(sampleName, pars.controlSamples)){
				continue;
			}
			sampColl.setUpSampleFromPrevious(sampleName);
			auto sampPtr = sampColl.sampleCollapses_.at(sampleName);
			for(uint32_t clusPos = 0;  clusPos < sampPtr->collapsed_.clusters_.size(); ++clusPos){
				if(sampPtr->collapsed_.clusters_[clusPos].seqBase_.frac_ >= pars.majorHaplotypeFracForRescue){
					majorHaps.emplace(sampColl.popCollapse_->collapsed_.clusters_[sampColl.popCollapse_->collapsed_.subClustersPositions_.at(sampPtr->collapsed_.clusters_[clusPos].getStubName(true))].seqBase_.name_);
					if(clusPos < 2){
						majorHapsForChi.emplace(sampColl.popCollapse_->collapsed_.clusters_[sampColl.popCollapse_->collapsed_.subClustersPositions_.at(sampPtr->collapsed_.clusters_[clusPos].getStubName(true))].seqBase_.name_);
					}
				}
			}
			sampColl.dumpSample(sampleName);
		}
		if(setUp.pars_.debug_){
			std::cout << "majorHaps: " << njh::conToStr(majorHaps, ",") << std::endl;
		}
		bool rescuedHaplotypes = false;
		for(const auto & sampleName : sampColl.passingSamples_){
			sampColl.setUpSampleFromPrevious(sampleName);
			auto sampPtr = sampColl.sampleCollapses_.at(sampleName);
			std::vector<uint32_t> toBeRescued;
			//iterator over haplotypes, determine if they should be considered for rescue, if they should be then check to see if they match a major haplotype
			for(const auto excludedPos : iter::range(sampPtr->excluded_.clusters_.size())){
				const auto & excluded = sampPtr->excluded_.clusters_[excludedPos];
				if(excluded.nameHasMetaData()){
					MetaDataInName excludedMeta(excluded.seqBase_.name_);
					std::set<std::string> otherExcludedCriteria;
					bool chimeriaExcludedRescue = false;
					bool oneOffExcludedRescue = false;
					bool lowFreqExcludedRescue = false;

					for(const auto & excMeta : excludedMeta.meta_){
						if(njh::beginsWith(excMeta.first, "Exclude") ){
							bool other = true;
							if(pars.rescueExcludedChimericHaplotypes && "ExcludeIsChimeric" == excMeta.first){
								chimeriaExcludedRescue = true;
								other = false;
							}
							if(pars.rescueExcludedOneOffLowFreqHaplotypes && "ExcludeFailedLowFreqOneOff" == excMeta.first){
								oneOffExcludedRescue = true;
								other = false;
							}
							if(pars.rescueExcludedLowFreqHaplotypes && "ExcludeFailedFracCutOff" == excMeta.first){
								lowFreqExcludedRescue = true;
								other = false;
							}
							if(other){
								otherExcludedCriteria.emplace(excMeta.first);
							}
						}
					}
					//check if it should be considered for rescue
					//std::cout << excluded.seqBase_.name_ << " consider for rescue: " << njh::colorBool((chimeriaExcludedRescue || oneOffExcludedRescue) && otherExcludedCriteria.empty()) << std::endl;
					if((chimeriaExcludedRescue || oneOffExcludedRescue || lowFreqExcludedRescue) && otherExcludedCriteria.empty()){
						//see if it matches a major haplotype
						bool rescue = false;
						if(chimeriaExcludedRescue){
							for(const auto & popHap : sampColl.popCollapse_->collapsed_.clusters_){
								if(popHap.seqBase_.seq_ == excluded.seqBase_.seq_ &&
										popHap.seqBase_.cnt_ > excluded.seqBase_.cnt_ &&
										njh::in(popHap.seqBase_.name_, majorHapsForChi)){
									rescue = true;
									break;
								}
							}
						}else{
							for(const auto & popHap : sampColl.popCollapse_->collapsed_.clusters_){
								if(popHap.seqBase_.seq_ == excluded.seqBase_.seq_ &&
										popHap.seqBase_.cnt_ > excluded.seqBase_.cnt_ &&
										njh::in(popHap.seqBase_.name_, majorHaps)){
									rescue = true;
									break;
								}
							}
						}
						if(rescue){
							toBeRescued.emplace_back(excludedPos);
						}
					}
				}
			}
			if(!toBeRescued.empty()){
				rescuedHaplotypes = true;
				std::sort(toBeRescued.rbegin(), toBeRescued.rend());
				for(const auto toRescue : toBeRescued){
					MetaDataInName excludedMeta(sampPtr->excluded_.clusters_[toRescue].seqBase_.name_);
					excludedMeta.addMeta("rescue", "TRUE");
					excludedMeta.resetMetaInName(sampPtr->excluded_.clusters_[toRescue].seqBase_.name_);
					//unmarking so as not to mess up chimera numbers
					sampPtr->excluded_.clusters_[toRescue].seqBase_.unmarkAsChimeric();
					for (auto & subRead : sampPtr->excluded_.clusters_[toRescue].reads_) {
						subRead->seqBase_.unmarkAsChimeric();
					}
					sampPtr->collapsed_.clusters_.emplace_back(sampPtr->excluded_.clusters_[toRescue]);
					sampPtr->excluded_.clusters_.erase(sampPtr->excluded_.clusters_.begin() + toRescue);
				}
				sampPtr->updateAfterExclustion();
				sampPtr->renameClusters("fraction");
			}
			sampColl.dumpSample(sampleName);
		}
		if(rescuedHaplotypes){
			//if excluded run pop clustering again
			sampColl.doPopulationClustering(sampColl.createPopInput(),
					alignerObj, collapserObj, pars.popIteratorMap);
		}
	}// end resuce operations for chimeria and low freq haplotypes



	if(pars.removeCommonlyLowFreqHaplotypes_){
		while(sampColl.excludeCommonlyLowFreqHaps(pars.lowFreqHaplotypeFracCutOff_)){
			//if excluded run pop clustering again
			sampColl.doPopulationClustering(sampColl.createPopInput(),
					alignerObj, collapserObj, pars.popIteratorMap);
		}
	}

	if(pars.removeOneSampOnlyOneOffHaps){
		if(sampColl.excludeOneSampOnlyOneOffHaps(pars.oneSampOnlyOneOffHapsFrac, alignerObj)){
			//if excluded run pop clustering again
			sampColl.doPopulationClustering(sampColl.createPopInput(),
					alignerObj, collapserObj, pars.popIteratorMap);
		}
	}

	if(pars.removeOneSampOnlyHaps){
		if(sampColl.excludeOneSampOnlyHaps(pars.oneSampOnlyHapsFrac)){
			//if excluded run pop clustering again
			sampColl.doPopulationClustering(sampColl.createPopInput(),
					alignerObj, collapserObj, pars.popIteratorMap);
		}
	}

	if(pars.rescueMatchingExpected && !expectedSeqs.empty()){
		bool rescuedHaplotypes = false;
		for(const auto & sampleName : sampColl.passingSamples_){
			sampColl.setUpSampleFromPrevious(sampleName);
			auto sampPtr = sampColl.sampleCollapses_.at(sampleName);
			std::vector<uint32_t> toBeRescued;
			//iterator over haplotypes, determine if they should be considered for rescue, if they should be then check to see if they match a major haplotype
			for(const auto excludedPos : iter::range(sampPtr->excluded_.clusters_.size())){
				const auto & excluded = sampPtr->excluded_.clusters_[excludedPos];
				if(excluded.nameHasMetaData()){
					MetaDataInName excludedMeta(excluded.seqBase_.name_);
					std::set<std::string> otherExcludedCriteria;
					bool chimeriaExcludedRescue = false;
					bool oneOffExcludedRescue = false;
					bool commonlyLowExcludedRescue = false;
					bool lowFreqExcludedRescue = false;

					for(const auto & excMeta : excludedMeta.meta_){
						if(njh::beginsWith(excMeta.first, "Exclude") ){
							bool other = true;
							if("ExcludeIsChimeric" == excMeta.first){
								chimeriaExcludedRescue = true;
								other = false;
							}
							if("ExcludeFailedLowFreqOneOff" == excMeta.first){
								oneOffExcludedRescue = true;
								other = false;
							}
							if("ExcludeCommonlyLowFreq" == excMeta.first){
								commonlyLowExcludedRescue = true;
								other = false;
							}
							if("ExcludeFailedFracCutOff" == excMeta.first){
								lowFreqExcludedRescue = true;
								other = false;
							}
							if(other){
								otherExcludedCriteria.emplace(excMeta.first);
							}
						}
					}
					//check if it should be considered for resuce
					if((chimeriaExcludedRescue || oneOffExcludedRescue || commonlyLowExcludedRescue || lowFreqExcludedRescue) && otherExcludedCriteria.empty()){
						//see if it matches a major haplotype
						bool rescue = false;
						for(const auto & expectedHap : expectedSeqs){
							if(expectedHap.seqBase_.seq_ == excluded.seqBase_.seq_){
								rescue = true;
								break;
							}
						}
						if(rescue){
							toBeRescued.emplace_back(excludedPos);
						}
					}
				}
			}
			if(!toBeRescued.empty()){
				rescuedHaplotypes = true;
				std::sort(toBeRescued.rbegin(), toBeRescued.rend());
				for(const auto toRescue : toBeRescued){
					MetaDataInName excludedMeta(sampPtr->excluded_.clusters_[toRescue].seqBase_.name_);
					excludedMeta.addMeta("rescue", "TRUE");
					excludedMeta.resetMetaInName(sampPtr->excluded_.clusters_[toRescue].seqBase_.name_);
					//unmarking so as not to mess up chimera numbers
					sampPtr->excluded_.clusters_[toRescue].seqBase_.unmarkAsChimeric();
					for (auto & subRead : sampPtr->excluded_.clusters_[toRescue].reads_) {
						subRead->seqBase_.unmarkAsChimeric();
					}
					sampPtr->collapsed_.clusters_.emplace_back(sampPtr->excluded_.clusters_[toRescue]);
					sampPtr->excluded_.clusters_.erase(sampPtr->excluded_.clusters_.begin() + toRescue);
				}
				sampPtr->updateAfterExclustion();
				sampPtr->renameClusters("fraction");
			}
			sampColl.dumpSample(sampleName);
		}
		if(rescuedHaplotypes){
			//if excluded run pop clustering again
			sampColl.doPopulationClustering(sampColl.createPopInput(),
					alignerObj, collapserObj, pars.popIteratorMap);
		}
	} //end resue of matching expected



	if(setUp.pars_.verbose_){
		std::cout << njh::bashCT::boldRed("Done Pop Clustering") << std::endl;
	}
	//if ("" != pars.previousPopFilename && !pars.noPopulation) {
	if ("" != pars.previousPopFilename) {
		sampColl.renamePopWithSeqs(getSeqs<readObject>(pars.previousPopFilename), pars.previousPopErrors);
	}

	if (!expectedSeqs.empty()) {
		sampColl.comparePopToRefSeqs(expectedSeqs, alignerObj);
	}


	std::vector<seqInfo> popSeqsPerSamp;
	std::vector<seqInfo> outPopSeqsPerSamp;

	std::unordered_map<std::string, uint32_t> sampCountsForPopHaps;
	uint32_t totalPopCount = 0;
	std::set<std::string> samplesCount;
	//if(!pars.noPopulation){
	popSeqsPerSamp = sampColl.genOutPopSeqsPerSample();
	outPopSeqsPerSamp = popSeqsPerSamp;
	for(const auto & popClus : sampColl.popCollapse_->collapsed_.clusters_){
		sampCountsForPopHaps[popClus.seqBase_.name_] = popClus.sampleClusters().size();
		totalPopCount += popClus.sampleClusters().size();
	}
	for(const auto & seq : popSeqsPerSamp){
		MetaDataInName seqMeta(seq.name_);
		samplesCount.emplace(seqMeta.getMeta("sample"));
	}

	//if(!pars.noPopulation){



	std::map<std::string, std::map<std::string, MetaDataInName>> knownAAMeta;
	//       seqName               transcript   amino acid positions and amino acid
	std::map<std::string, std::map<std::string, std::string>> knownAATyped;
	std::map<std::string, std::map<std::string, std::string>> fullAATyped;

	if("" != pars.transPars.gffFnp_){

		auto variantInfoDir =  njh::files::make_path(sampColl.masterOutputDir_, "variantInfo");
		njh::files::makeDir(njh::files::MkdirPar{variantInfoDir});
		translator->pars_.keepTemporaryFiles_ = true;
		translator->pars_.workingDirtory_ = variantInfoDir;
		auto popSeqsOpts = setUp.pars_.ioOptions_;
		SeqIOOptions popOutOpts(
				njh::files::make_path(variantInfoDir,
						"PopSeqs" + setUp.pars_.ioOptions_.getOutExtension()),
						setUp.pars_.ioOptions_.outFormat_);
		SeqOutput::write(sampColl.popCollapse_->collapsed_.clusters_, popOutOpts);

		popSeqsOpts.firstName_ = popOutOpts.out_.outName();
		auto translatedRes = translator->run(popSeqsOpts, sampCountsForPopHaps, pars.variantCallerRunPars);
		SeqOutput transwriter(SeqIOOptions::genFastaOut(njh::files::make_path(variantInfoDir, "translatedInput.fasta")));
		for(const auto & seqName : translatedRes.translations_){
			for(const auto & transcript : seqName.second){
				transwriter.openWrite(transcript.second.translation_);
			}
		}
		SeqInput popReader(popSeqsOpts);
		auto popSeqs = popReader.readAllReads<seqInfo>();
		std::unordered_map<std::string, uint32_t> popSeqsPosition;
		for(const auto & popPos : iter::range(popSeqs.size())){
			popSeqsPosition[popSeqs[popPos].name_] = popPos;
		}
		OutputStream popBedLocs(njh::files::make_path(variantInfoDir, "PopSeqs.bed"));
		for(const auto & seqLocs : translatedRes.seqAlns_){
			for(const auto & loc : seqLocs.second){
				popBedLocs << loc.gRegion_.genBedRecordCore().toDelimStrWithExtra() << std::endl;
			}
		}

		for(const auto & pop : popSeqs){
			if(!njh::in(pop.name_, translatedRes.seqAlns_)){
				popBedLocs << "*"
						<< "\t" << "*"
						<< "\t" << "*"
						<< "\t" << pop.name_
						<< "\t" << "*"
						<< "\t" << "*" << std::endl;
			}
		}

		{
			//protein
			for(auto & varPerTrans : translatedRes.proteinVariants_){
				auto snpsPositions = getVectorOfMapKeys(varPerTrans.second.snpsFinal);
				njh::sort(snpsPositions);
				OutputStream snpTabOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidVariable.tab.txt"))));
				snpTabOut << "transcript\tposition(1-based)\trefAA\tAA\tcount\tfraction\talleleDepth\tsamples" << std::endl;
				for(const auto & snpPos : snpsPositions){
					for(const auto & aa : varPerTrans.second.allBases[snpPos]){
						snpTabOut << varPerTrans.first
								<< "\t" << snpPos + 1
								<< "\t" << translatedRes.proteinForTranscript_[varPerTrans.first][snpPos]
								<< "\t" << aa.first
								<< "\t" << aa.second
								<< "\t" << aa.second/static_cast<double>(totalPopCount)
								<< "\t" << totalPopCount
								<< "\t" << samplesCount.size() << std::endl;
					}
				}
				std::set<uint32_t> knownMutationsLocations;
				OutputStream allAATabOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidsAll.tab.txt"))));
				allAATabOut << "transcript\tposition(1-based)\trefAA\tAA\tcount\tfraction\talleleDepth\tsamples" << std::endl;
				for(const auto & snpPos : varPerTrans.second.allBases){
					if(njh::in(snpPos.first + 1, knownMutationsLocationsMap[varPerTrans.first])){
						knownMutationsLocations.emplace(snpPos.first);
					}
					for(const auto & aa : snpPos.second){
						allAATabOut << varPerTrans.first
								<< "\t" << snpPos.first + 1
								<< "\t" << translatedRes.proteinForTranscript_[varPerTrans.first][snpPos.first]
								<< "\t" << aa.first
								<< "\t" << aa.second
								<< "\t" << aa.second/static_cast<double>(totalPopCount)
								<< "\t" << totalPopCount
								<< "\t" << samplesCount.size() << std::endl;
					}
				}
				if(!varPerTrans.second.variablePositons_.empty()){
					uint32_t variableStart = vectorMinimum(std::vector<uint32_t>(varPerTrans.second.variablePositons_.begin(), varPerTrans.second.variablePositons_.end()));
					uint32_t variableStop = vectorMaximum(std::vector<uint32_t>(varPerTrans.second.variablePositons_.begin(), varPerTrans.second.variablePositons_.end()));

					GenomicRegion variableRegion;
					variableRegion.chrom_ = varPerTrans.first;
					variableRegion.start_ = variableStart +1;
					variableRegion.end_ = variableStop;

					OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_variableRegion.bed"))));
					bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
				}
				std::set<uint32_t> allLocations(knownMutationsLocations.begin(), knownMutationsLocations.end());
				for(const auto & variablePos : varPerTrans.second.snpsFinal){
					allLocations.emplace(variablePos.first);
				}
				std::map<std::string, MetaDataInName> aaMeta;


				for(auto & seqName : translatedRes.translations_){
					if(njh::in(varPerTrans.first, seqName.second)){
						for(const auto & variablePos : allLocations){
							//std::cout << __FILE__ << " " << __LINE__ << std::endl;
							auto aa = seqName.second[varPerTrans.first].queryAlnTranslation_.seq_[getAlnPosForRealPos(seqName.second[varPerTrans.first].refAlnTranslation_.seq_, variablePos)];
							aaMeta[seqName.first].addMeta(estd::to_string(variablePos), aa, false);
							//std::cout << __FILE__ << " " << __LINE__ << std::endl;
						}
					}
					for(const auto & knownLoc : knownMutationsLocations){

						//std::cout << __FILE__ << " " << __LINE__ << std::endl;
						auto aa = seqName.second[varPerTrans.first].queryAlnTranslation_.seq_[getAlnPosForRealPos(seqName.second[varPerTrans.first].refAlnTranslation_.seq_, knownLoc)];
						knownAAMeta[seqName.first.substr(0, seqName.first.rfind("_f"))][varPerTrans.first].addMeta(estd::to_string(knownLoc), aa, false);
						//std::cout << __FILE__ << " " << __LINE__ << std::endl;
					}
				}

				for (auto & seqName : translatedRes.translations_) {
					if (njh::in(varPerTrans.first, seqName.second)) {
						VecStr allAAPosCoded;
						std::string popName = seqName.first.substr(0, seqName.first.rfind("_f"));
						std::string transcript = varPerTrans.first;
						for (const auto & loc : allLocations) {
							//std::cout << __FILE__ << " " << __LINE__ << std::endl;
							auto aa =seqName.second[varPerTrans.first].queryAlnTranslation_.seq_[getAlnPosForRealPos(seqName.second[varPerTrans.first].refAlnTranslation_.seq_,loc)];
							allAAPosCoded.emplace_back(njh::pasteAsStr(loc + 1, "-", aa));
							//std::cout << __FILE__ << " " << __LINE__ << std::endl;
						}
						if(!allAAPosCoded.empty()){
							fullAATyped[popName][transcript] = njh::conToStr(allAAPosCoded, ":");
						}else{
							fullAATyped[popName][transcript] = "NONE";
						}
					}
				}
				if(!knownMutationsLocations.empty()){
					OutputStream knownTabOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidKnownMutations.tab.txt"))));
					knownTabOut << "transcript\tposition(1-based)\trefAA\tAA\tcount\tfraction\talleleDepth\tsamples" << std::endl;
					for(const auto & snpPos : knownMutationsLocations){
						for(const auto & aa : varPerTrans.second.allBases[snpPos]){
							knownTabOut << varPerTrans.first
									<< "\t" << snpPos + 1
									<< "\t" << translatedRes.proteinForTranscript_[varPerTrans.first][snpPos]
									<< "\t" << aa.first
									<< "\t" << aa.second
									<< "\t" << aa.second/static_cast<double>(totalPopCount)
									<< "\t" << totalPopCount
									<< "\t" << samplesCount.size() << std::endl;
						}
					}
					for(const auto & pop : knownAAMeta){
						std::string popName = pop.first;
						std::string transcript = varPerTrans.first;
						std::vector<std::string> aaPos;
						for(const auto & variablePos : knownMutationsLocations){
							aaPos.emplace_back(estd::to_string(variablePos + 1) + "-" + pop.second.at(transcript).getMeta(estd::to_string(variablePos)));
						}
						knownAATyped[popName][transcript] = njh::conToStr(aaPos, ":");
					}
				}




				OutputStream outPopHapAminos(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-popHapToAmino.tab.txt")));
				outPopHapAminos << "h_PopUID" ;
				VecStr aminoPositionsHeader;
				for(const auto & variablePos : allLocations){
					outPopHapAminos << "\t" << varPerTrans.first << "-aa" << variablePos + 1;
					aminoPositionsHeader.emplace_back(njh::pasteAsStr(varPerTrans.first, "-aa",  variablePos + 1));
				}
				outPopHapAminos << std::endl;
				std::unordered_map<std::string, std::string> popHapAminoTyped;
				std::unordered_map<std::string, VecStr> popHapAminoTypedRow;

				for(const auto & pop : aaMeta){
					outPopHapAminos << pop.first.substr(0, pop.first.rfind("_f")) ;
					std::string typed = varPerTrans.first + "=";
					std::vector<std::string> aaPos;
					std::vector<std::string> aa;

					for(const auto & variablePos : allLocations){
						outPopHapAminos << "\t" << pop.second.getMeta(estd::to_string(variablePos));
						aaPos.emplace_back(estd::to_string(variablePos + 1) + "-" + pop.second.getMeta(estd::to_string(variablePos)));
						aa.emplace_back(pop.second.getMeta(estd::to_string(variablePos)));

					}
					if (!aaPos.empty()) {
						typed += njh::conToStr(aaPos, ":");
					} else {
						typed += "NONE";
					}
					//typed +=";";

					popHapAminoTyped[pop.first.substr(0, pop.first.rfind("_f"))] = typed;
					popHapAminoTypedRow[pop.first.substr(0, pop.first.rfind("_f"))] = aa;
					outPopHapAminos << std::endl;
				}
				auto popMetaTable = seqsToMetaTable(popSeqsPerSamp);
				popMetaTable.deleteColumn("seq");
				popMetaTable.deleteColumn("count");
				addOtherVec(popMetaTable.columnNames_, aminoPositionsHeader);
				for(auto & row : popMetaTable){
					MetaDataInName::removeMetaDataInName(row[popMetaTable.getColPos("name")]);
					auto popName = row[popMetaTable.getColPos("PopUID")];
					VecStr aminos;
					if(njh::in(popName, popHapAminoTypedRow)){
						aminos = popHapAminoTypedRow[popName];
					}else{
						//pop uid was untypable, didn't align
						aminos = VecStr{aminoPositionsHeader.size(), std::string("NA")};
					}
					addOtherVec(row, aminos);
				}

				popMetaTable.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(variantInfoDir, varPerTrans.first + "-popSeqsWithMetaAndVariableAAInfoTable.tab.txt")));

//				for(auto & popHapSamp : outPopSeqsPerSamp){
//					MetaDataInName meta(popHapSamp.name_);
//					auto typed = popHapAminoTyped[meta.getMeta("PopUID")].substr(popHapAminoTyped[meta.getMeta("PopUID")].find("=") + 1);
//					if("" ==typed){
//						typed = "NONE";
//					}
//					meta.addMeta(varPerTrans.first, typed);
//					meta.resetMetaInName(popHapSamp.name_);
//				}
			}
		}

		{
			//snps
			OutputStream outSnpDepthPerSample(njh::files::make_path(variantInfoDir, njh::pasteAsStr("snpDepthPerSample.tab.txt")));
			outSnpDepthPerSample << "AnalysisName\tsample\tchromosome\tposition\trefBAse\tbase\treadDepth" << std::endl;

			for( auto & varPerChrom : translatedRes.seqVariants_){
				auto snpsPositions = getVectorOfMapKeys(varPerChrom.second.snpsFinal);
				njh::sort(snpsPositions);
				OutputStream snpTabOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-SNPs.tab.txt"))));
				snpTabOut << "chromosome\tposition(0-based)\trefBase\tbase\tcount\tfraction\talleleDepth\tsamples" << std::endl;
				for(const auto & snpPos : snpsPositions){
					for(const auto & base : varPerChrom.second.allBases[snpPos]){
						snpTabOut << varPerChrom.first
								<< "\t" << snpPos
								<< "\t" << translatedRes.baseForPosition_[varPerChrom.first][snpPos]
								<< "\t" << base.first
								<< "\t" << base.second
								<< "\t" << base.second/static_cast<double>(totalPopCount)
								<< "\t" << totalPopCount
								<< "\t" << samplesCount.size() << std::endl;
					}
				}
				OutputStream allBasesTabOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-allBases.tab.txt"))));
				allBasesTabOut << "chromosome\tposition(0-based)\trefBase\tbase\tcount\tfraction\talleleDepth\tsamples" << std::endl;
				for(const auto & snpPos : varPerChrom.second.allBases){
					for(const auto & base : snpPos.second){
						allBasesTabOut << varPerChrom.first
								<< "\t" << snpPos.first
								<< "\t" << translatedRes.baseForPosition_[varPerChrom.first][snpPos.first]
								<< "\t" << base.first
								<< "\t" << base.second
								<< "\t" << base.second/static_cast<double>(totalPopCount)
								<< "\t" << totalPopCount
								<< "\t" << samplesCount.size() << std::endl;
					}
				}

				if(!varPerChrom.second.variablePositons_.empty()){
					uint32_t variableStart = vectorMinimum(std::vector<uint32_t>(varPerChrom.second.variablePositons_.begin(), varPerChrom.second.variablePositons_.end()));
					uint32_t variableStop = vectorMaximum(std::vector<uint32_t>(varPerChrom.second.variablePositons_.begin(), varPerChrom.second.variablePositons_.end()));

					GenomicRegion variableRegion;
					variableRegion.chrom_ = varPerChrom.first;
					variableRegion.start_ = variableStart;
					variableRegion.end_ = variableStop;

					OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-chromosome_variableRegion.bed"))));
					bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
				}

				std::map<std::string, MetaDataInName> snpMeta;
				for(auto & seqName : translatedRes.seqAlns_){
					for(const auto & variablePos : varPerChrom.second.snpsFinal){
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
//						std::cout << "seqName.second.front()->gRegion_.start_: " << seqName.second.front()->gRegion_.start_ << std::endl;
//						std::cout << "variablePos: " << variablePos.first << std::endl;
						if(variablePos.first < seqName.second.front().gRegion_.start_ || variablePos.first >= seqName.second.front().gRegion_.end_){
							snpMeta[seqName.first].addMeta(estd::to_string(variablePos.first), "X", false);
						} else {
							//std::cout << __FILE__ << " " << __LINE__ << std::endl;
							auto aa = seqName.second.front().querySeq_.seq_[getAlnPosForRealPos(seqName.second.front().alnRefSeq_.seq_, variablePos.first - seqName.second.front().gRegion_.start_)];
							snpMeta[seqName.first].addMeta(estd::to_string(variablePos.first), aa, false);
							//std::cout << __FILE__ << " " << __LINE__ << std::endl;
						}
					}
				}

				OutputStream outPopHapAminos(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-popHapToSNPs.tab.txt")));
				outPopHapAminos << "h_PopUID" ;
				VecStr aminoPositionsHeader;
				for(const auto & variablePos : varPerChrom.second.snpsFinal){
					outPopHapAminos << "\t" << varPerChrom.first << "-" << variablePos.first;
					aminoPositionsHeader.emplace_back(njh::pasteAsStr(varPerChrom.first, "-",  variablePos.first));
				}
				outPopHapAminos << std::endl;
				std::unordered_map<std::string, std::string> popHapAminoTyped;
				std::unordered_map<std::string, VecStr> popHapAminoTypedRow;

				for(const auto & pop : snpMeta){
					outPopHapAminos << pop.first.substr(0, pop.first.rfind("_f")) ;
					std::string typed = varPerChrom.first + "=";
					std::vector<std::string> aaPos;
					std::vector<std::string> aa;

					for(const auto & variablePos : varPerChrom.second.snpsFinal){
						outPopHapAminos << "\t" << pop.second.getMeta(estd::to_string(variablePos.first));
						aaPos.emplace_back(estd::to_string(variablePos.first) + "-" + pop.second.getMeta(estd::to_string(variablePos.first)));
						aa.emplace_back(pop.second.getMeta(estd::to_string(variablePos.first)));
					}
					if (!aaPos.empty()) {
						typed += njh::conToStr(aaPos, ":");
					} else {
						typed += "NONE";
					}
					//typed +=";";
					popHapAminoTyped[pop.first.substr(0, pop.first.rfind("_f"))] = typed;
					popHapAminoTypedRow[pop.first.substr(0, pop.first.rfind("_f"))] = aa;
					outPopHapAminos << std::endl;
				}


				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
				auto popMetaTable = seqsToMetaTable(popSeqsPerSamp);
				popMetaTable.deleteColumn("seq");
				popMetaTable.deleteColumn("count");
				addOtherVec(popMetaTable.columnNames_, aminoPositionsHeader);

				std::map<std::string, std::map<std::string, std::map<uint32_t, std::map<std::string, double>>>> snpsPerSampleDepth;

				for(auto & row : popMetaTable){
					MetaDataInName::removeMetaDataInName(row[popMetaTable.getColPos("name")]);
					auto popName = row[popMetaTable.getColPos("PopUID")];
					VecStr aminos;
					if(njh::in(popName, popHapAminoTypedRow)){
						aminos = popHapAminoTypedRow[popName];
						uint32_t snpPosition = 0;
						for(const auto & snp : aminos){
							if("X" == snp){
								continue;
							}
							auto sample = row[popMetaTable.getColPos("sample")];
							auto chrom = aminoPositionsHeader[snpPosition].substr(0, aminoPositionsHeader[snpPosition].rfind("-"));
							auto position = njh::StrToNumConverter::stoToNum<uint32_t>(aminoPositionsHeader[snpPosition].substr(aminoPositionsHeader[snpPosition].rfind("-") + 1));
							auto depth = njh::StrToNumConverter::stoToNum<uint32_t>(row[popMetaTable.getColPos("readCount")]);

							snpsPerSampleDepth[sample][chrom][position][snp] += depth;
							++snpPosition;
						}
					}else{
						//pop uid was untypable, didn't align
						aminos = VecStr{aminoPositionsHeader.size(), std::string("NA")};
					}
					addOtherVec(row, aminos);
				}
				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
				popMetaTable.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(variantInfoDir, varPerChrom.first + "-popSeqsWithMetaAndVariableSNPInfoTable.tab.txt")));
				for(const auto & sample : snpsPerSampleDepth){
					for(const auto & chrom : sample.second){
						for(const auto & position : chrom.second){
							for(const auto & snp : position.second){
								outSnpDepthPerSample
										<< "\t" << pars.experimentName
								    << "\t" << sample.first
										<< "\t" << chrom.first
										<< "\t" << position.first
										<< "\t" << translatedRes.baseForPosition_[varPerChrom.first][position.first]
										<< "\t" << snp.first
										<< "\t" << snp.second << std::endl;
							}
						}
					}
				}
//				for(auto & popHapSamp : outPopSeqsPerSamp){
////					std::cout << __FILE__ << " " << __LINE__ << std::endl;
//					MetaDataInName meta(popHapSamp.name_);
//					auto typed = popHapAminoTyped[meta.getMeta("PopUID")].substr(popHapAminoTyped[meta.getMeta("PopUID")].find("=") + 1);
//					if("" ==typed){
//						typed = "NONE";
//					}
//					meta.addMeta(varPerChrom.first, typed);
//					meta.resetMetaInName(popHapSamp.name_);
//				}
			}
		}
	}

	for(auto & clus : sampColl.popCollapse_->collapsed_.clusters_){
		auto popName = clus.seqBase_.name_.substr(0, clus.seqBase_.name_.rfind("_f"));
		std::string typed = "";
		for(const auto & tran : fullAATyped[popName]){
			if("" != typed){
				typed += ";";
			}
			typed += tran.first + "=" + tran.second;
//			std::cout << tran.second << std::endl;
			clus.meta_.addMeta("h_AATyped", typed);
		}
	}
	sampColl.printSampleCollapseInfo(
			njh::files::make_path(sampColl.masterOutputDir_,
					"selectedClustersInfo.tab.txt.gz"));

	if(pars.writeOutAllInfoFile){
		sampColl.printAllSubClusterInfo(
					njh::files::make_path(sampColl.masterOutputDir_,
							"allClustersInfo.tab.txt.gz"));
	}
	sampColl.symlinkInSampleFinals();
	sampColl.outputRepAgreementInfo();

	table hapIdTab = sampColl.genHapIdTable();
	hapIdTab.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(sampColl.masterOutputDir_,
			"hapIdTable.tab.txt.gz"), true));
	sampColl.dumpPopulation();

	SeqOutput::write(outPopSeqsPerSamp,SeqIOOptions(njh::files::make_path(sampColl.masterOutputDir_, "population", "popSeqsWithMetaWtihSampleName"), setUp.pars_.ioOptions_.outFormat_));

	auto popMetaTable = seqsToMetaTable(outPopSeqsPerSamp);
	for(auto & row : popMetaTable){
		MetaDataInName::removeMetaDataInName(row[popMetaTable.getColPos("name")]);
	}

	popMetaTable.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(sampColl.masterOutputDir_, "population", "popSeqsWithMetaWtihSampleNameTable.tab.txt.gz")));


	if("" != pars.groupingsFile){
		if(pars.noWriteGroupInfoFiles){
			//still write the meta info for use with downstream analysis
			auto groupsTopDir = njh::files::make_path(sampColl.masterOutputDir_, "groups");
			njh::files::makeDir(njh::files::MkdirPar(groupsTopDir, true));
			OutputStream groupMetaJsonFile(njh::files::make_path(groupsTopDir, "groupMetaData.json"));
			groupMetaJsonFile << sampColl.groupMetaData_->toJson() << std::endl;
		} else {
			sampColl.createGroupInfoFiles();
		}
	}

	sampColl.createCoreJsonFile();

	//collect extraction dirs
	std::set<bfs::path> extractionDirs;
	for(const auto & file : analysisFiles){
		auto fileToks = njh::tokenizeString(bfs::relative(file.first, pars.masterDir).string(), "/");
		if(njh::in(fileToks[0], pars.excludeSamples)){
			continue;
		}
		auto metaDataJsonFnp = njh::files::make_path(file.first.parent_path(), "metaData.json");
		if(bfs::exists(metaDataJsonFnp)){
			auto metaJson = njh::json::parseFile(metaDataJsonFnp.string());
			if(metaJson.isMember("extractionDir")){
				extractionDirs.emplace(metaJson["extractionDir"].asString());
			}
		}
	}
	if(setUp.pars_.verbose_){
		std::cout << "Extraction Dirs" << std::endl;
		std::cout << njh::conToStr(extractionDirs, "\n") << std::endl;
	}
	table profileTab;
	table statsTab;
	for(const auto & extractDir : extractionDirs){
		auto profileFnp = njh::files::make_path(extractDir, "extractionProfile.tab.txt");
		auto statsFnp = njh::files::make_path(extractDir, "extractionStats.tab.txt");
		if(bfs::exists(profileFnp)){
			table currentProfileTab(profileFnp.string(), "\t", true);
			auto oldColumnNames = currentProfileTab.columnNames_;
			currentProfileTab.addColumn(VecStr{extractDir.filename().string()}, "extractionDir");
			currentProfileTab = currentProfileTab.getColumns(concatVecs(VecStr{"extractionDir"}, oldColumnNames));
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
			curentStatsTab = curentStatsTab.getColumns(concatVecs(VecStr{"extractionDir"}, oldColumnNames));
			if(statsTab.empty()){
				statsTab = curentStatsTab;
			}else{
				statsTab.rbind(curentStatsTab, false);
			}
		}
	}

	auto extractionOutputDir = njh::files::make_path(setUp.pars_.directoryName_,
			"extractionInfo");
	njh::files::makeDirP(njh::files::MkdirPar(extractionOutputDir.string()));
	if (!profileTab.empty()) {
		profileTab.sortTable("extractionDir", false);
		auto profileTabOpts =
				TableIOOpts::genTabFileOut(
						njh::files::make_path(extractionOutputDir,
								"extractionProfile.tab.txt").string(), true);
		profileTabOpts.out_.overWriteFile_ = true;
		profileTab.outPutContents(profileTabOpts);
	}
	if (!statsTab.empty()) {
		auto statsTabOpts =
				TableIOOpts::genTabFileOut(
						njh::files::make_path(extractionOutputDir,
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

}  // namespace njhseq
