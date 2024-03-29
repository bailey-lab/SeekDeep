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
#include <njhseq/objects/Meta/MetaUtils.hpp>
#include <njhseq/objects/dataContainers/tables/TableReader.hpp>

namespace njhseq {







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
	pars.experimentNames.samples_ = std::set<std::string>{samplesDirs.begin(), samplesDirs.end()};
	collapse::SampleCollapseCollection sampColl(setUp.pars_.ioOptions_, pars.masterDir,
			setUp.pars_.directoryName_,
			pars.experimentNames,
			pars.preFiltCutOffs);
	sampColl.keepSampleInfoInMemory_ = pars.keepSampleInfoInMemory_;


	if("" != pars.groupingsFile){
		sampColl.addGroupMetaData(pars.groupingsFile);
	}

	//process custom cut offs
	std::unordered_map<std::string, double> customCutOffsMap = collapse::SampleCollapseCollection::processCustomCutOffs(pars.customCutOffs, samplesDirs, pars.fracCutoff);
	std::unordered_map<std::string, double> customCutOffsMapPerRep = collapse::SampleCollapseCollection::processCustomCutOffs(pars.customCutOffs, samplesDirs, pars.withinReplicateFracCutOff);

	{
		njh::concurrent::LockableQueue<std::string> sampleQueue(samplesDirs);
		njhseq::concurrent::AlignerPool alnPool(alignerObj, pars.numThreads);
		alnPool.initAligners();
		alnPool.outAlnDir_ = setUp.pars_.outAlnInfoDirName_;

		std::function<void()> setupClusterSamples = [&sampleQueue, &alnPool,&collapserObj,&pars,&setUp,
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
					sampColl.sampleCollapses_.at(samp)->excluded_.checkAgainstExpected(expectedSeqs, *currentAligner, false);
					sampColl.sampleCollapses_.at(samp)->collapsed_.checkAgainstExpected(expectedSeqs, *currentAligner, false);
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

				if(!sampColl.keepSampleInfoInMemory_){
					sampColl.dumpSample(samp);
				}
				if(setUp.pars_.verbose_){
					std::cout << "Ending: " << samp << std::endl;
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(setupClusterSamples, pars.numThreads);
	}

	//read in the dump alignment cache
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);


	if(setUp.pars_.verbose_){
		std::cout << njh::bashCT::boldGreen("Pop Clustering") << std::endl;
	}

	//first population clustering
	sampColl.doPopulationClustering(sampColl.createPopInput(), alignerObj, collapserObj, pars.popIteratorMap);


	if(pars.rescuePars_.performResuce()){
		sampColl.conductResuceOperations(pars.rescuePars_, alignerObj, collapserObj, pars.popIteratorMap);
	}
	sampColl.performLowLevelFilters(pars.lowLevelPopFiltPars_, alignerObj, collapserObj, pars.popIteratorMap);

	if(pars.rescueMatchingExpected && !expectedSeqs.empty()){
		sampColl.rescueMatchingSeqs(expectedSeqs, alignerObj, collapserObj, pars.popIteratorMap);
	}

	if(setUp.pars_.verbose_){
		std::cout << njh::bashCT::boldRed("Done Pop Clustering") << std::endl;
	}
	//if ("" != pars.previousPopFilename && !pars.noPopulation) {
	if ("" != pars.previousPopFilename) {
		auto previousPopSeqsRaw = getSeqs<readObject>(pars.previousPopFilename);
		//collapse indentical seqs
		std::vector<readObject> previousPopSeqs;
		std::vector<std::set<std::string>> allNamesForPreviousPops;
		for(const auto & seq : previousPopSeqsRaw){
			bool matchPrevious = false;
			for(const auto  other : iter::enumerate(previousPopSeqs)){
				if(other.element.seqBase_.seq_ == seq.seqBase_.seq_){
					matchPrevious = true;
					allNamesForPreviousPops[other.first].emplace(seq.seqBase_.name_);
					break;
				}
			}
			if(!matchPrevious){
				previousPopSeqs.emplace_back(seq);
				allNamesForPreviousPops.emplace_back(std::set<std::string>{seq.seqBase_.name_});
			}
		}
		//rename
		for(const auto  other : iter::enumerate(previousPopSeqs)){
			if(allNamesForPreviousPops[other.index].size() > 1){
				other.element.seqBase_.name_ = njh::conToStr(allNamesForPreviousPops[other.index], ":");
			}
		}
		sampColl.renamePopWithSeqs(previousPopSeqs, pars.previousPopErrors);
	}

	if (!expectedSeqs.empty()) {
		sampColl.comparePopToRefSeqs(expectedSeqs, alignerObj);
	}

  alignerObj.processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_,
                                  setUp.pars_.verbose_);
  setUp.rLog_ << alignerObj.numberOfAlingmentsDone_ << "\n";
  if (setUp.pars_.verbose_) {
    std::cout << alignerObj.numberOfAlingmentsDone_ << std::endl;
    setUp.logRunTime(std::cout);
  }
	std::vector<seqInfo> popSeqsPerSamp;
  std::vector<seqInfo> outPopSeqsPerSamp = sampColl.genOutPopSeqsPerSample();



  auto varCallDirPath = njh::files::make_path(sampColl.masterOutputDir_,  "variantCalling");
  if(!pars.collapseVarCallPars.transPars.lzPars_.genomeFnp.empty()){
    //if genome set call variants against genome
    pars.collapseVarCallPars.identifier = pars.experimentNames.populationName_;
    pars.collapseVarCallPars.outputDirectory = varCallDirPath;
    collapseAndCallVariants(pars.collapseVarCallPars, outPopSeqsPerSamp);
  }

//  for(const auto & popClus : sampColl.popCollapse_->collapsed_.clusters_){
//		auto samples = njh::getVecOfMapKeys(popClus.sampleClusters());
//		sampCountsForPopHaps[popClus.seqBase_.name_].insert(samples.begin(), samples.end());
//		totalPopCount += popClus.sampleClusters().size();
//	}
//	for(const auto & seq : popSeqsPerSamp){
//		MetaDataInName seqMeta(seq.name_);
//		samplesCount.emplace(seqMeta.getMeta("sample"));
//	}




//	std::map<std::string, std::map<std::string, MetaDataInName>> knownAAMeta;
//	//       seqName               transcript   amino acid positions and amino acid
//	std::map<std::string, std::map<std::string, std::string>> knownAATyped;
//	std::map<std::string, std::map<std::string, std::string>> fullAATyped;
//
//	if("" != pars.transPars.gffFnp_){
//
//		auto variantInfoDir =  njh::files::make_path(sampColl.masterOutputDir_, "variantInfo");
//		njh::files::makeDir(njh::files::MkdirPar{variantInfoDir});
//		translator->pars_.keepTemporaryFiles_ = true;
//		translator->pars_.workingDirtory_ = variantInfoDir;
//		auto popSeqsOpts = setUp.pars_.ioOptions_;
//		SeqIOOptions popOutOpts(
//				njh::files::make_path(variantInfoDir,
//						"PopSeqs" + setUp.pars_.ioOptions_.getOutExtension()),
//						setUp.pars_.ioOptions_.outFormat_);
//		SeqOutput::write(sampColl.popCollapse_->collapsed_.clusters_, popOutOpts);
//
//		popSeqsOpts.firstName_ = popOutOpts.out_.outName();
//		auto translatedRes = translator->run(popSeqsOpts, sampCountsForPopHaps, pars.variantCallerRunPars);
//		SeqOutput transwriter(SeqIOOptions::genFastaOut(njh::files::make_path(variantInfoDir, "translatedInput.fasta")));
//		for(const auto & seqName : translatedRes.translations_){
//			for(const auto & transcript : seqName.second){
//				transwriter.openWrite(transcript.second.translation_);
//			}
//		}
//		SeqInput popReader(popSeqsOpts);
//		auto popSeqs = popReader.readAllReads<seqInfo>();
//		std::unordered_map<std::string, uint32_t> popSeqsPosition;
//		for(const auto popPos : iter::range(popSeqs.size())){
//			popSeqsPosition[popSeqs[popPos].name_] = popPos;
//		}
//		OutputStream popBedLocs(njh::files::make_path(variantInfoDir, "PopSeqs.bed"));
//		for(const auto & seqLocs : translatedRes.seqAlns_){
//			for(const auto & loc : seqLocs.second){
//				popBedLocs << loc.gRegion_.genBedRecordCore().toDelimStrWithExtra() << std::endl;
//			}
//		}
//
//		for(const auto & pop : popSeqs){
//			if(!njh::in(pop.name_, translatedRes.seqAlns_)){
//				totalPopCount -= sampCountsForPopHaps[pop.name_].size();
//				popBedLocs << "*"
//						<< "\t" << "*"
//						<< "\t" << "*"
//						<< "\t" << pop.name_
//						<< "\t" << "*"
//						<< "\t" << "*" << std::endl;
//			}
//		}
//
//		{
//			//protein
//			for(auto & varPerTrans : translatedRes.proteinVariants_){
//				varPerTrans.second.writeOutSNPsFinalInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidVariable.tab.txt")), varPerTrans.first, true	);
//				varPerTrans.second.writeOutSNPsAllInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidsAll.tab.txt")), varPerTrans.first, true	);
//
//				std::set<uint32_t> knownMutationsLocations;
//				for(const auto & snpPos : varPerTrans.second.allBases){
//					if(njh::in(snpPos.first + 1, knownMutationsLocationsMap[varPerTrans.first])){
//						knownMutationsLocations.emplace(snpPos.first);
//					}
//				}
//				if(!varPerTrans.second.variablePositons_.empty()){
//					GenomicRegion variableRegion = varPerTrans.second.getVariableRegion();
//					variableRegion.start_ = variableRegion.start_ +1;
//					OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_variableRegion.bed"))));
//					bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
//				}
//				std::set<uint32_t> allLocations(knownMutationsLocations.begin(), knownMutationsLocations.end());
//				for(const auto & variablePos : varPerTrans.second.snpsFinal){
//					allLocations.emplace(variablePos.first);
//				}
//				std::map<std::string, MetaDataInName> aaMeta;
//
//
//				for(auto & seqName : translatedRes.translations_){
//					if(njh::in(varPerTrans.first, seqName.second)){
//						for(const auto & variablePos : allLocations){
//							//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//							auto aa = seqName.second[varPerTrans.first].queryAlnTranslation_.seq_[getAlnPosForRealPos(seqName.second[varPerTrans.first].refAlnTranslation_.seq_, variablePos)];
//							aaMeta[seqName.first].addMeta(estd::to_string(variablePos), aa, false);
//							//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//						}
//					}
//					for(const auto & knownLoc : knownMutationsLocations){
//
//						//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//						auto aa = seqName.second[varPerTrans.first].queryAlnTranslation_.seq_[getAlnPosForRealPos(seqName.second[varPerTrans.first].refAlnTranslation_.seq_, knownLoc)];
//						knownAAMeta[seqName.first.substr(0, seqName.first.rfind("_f"))][varPerTrans.first].addMeta(estd::to_string(knownLoc), aa, false);
//						//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//					}
//				}
//
//				for (auto & seqName : translatedRes.translations_) {
//					if (njh::in(varPerTrans.first, seqName.second)) {
//						VecStr allAAPosCoded;
//						std::string popName = seqName.first.substr(0, seqName.first.rfind("_f"));
//						std::string transcript = varPerTrans.first;
//						for (const auto & loc : allLocations) {
//							//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//							auto aa =seqName.second[varPerTrans.first].queryAlnTranslation_.seq_[getAlnPosForRealPos(seqName.second[varPerTrans.first].refAlnTranslation_.seq_,loc)];
//							allAAPosCoded.emplace_back(njh::pasteAsStr(loc + 1, "-", aa));
//							//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//						}
//						if(!allAAPosCoded.empty()){
//							fullAATyped[popName][transcript] = njh::conToStr(allAAPosCoded, ":");
//						}else{
//							fullAATyped[popName][transcript] = "NONE";
//						}
//					}
//				}
//				if(!knownMutationsLocations.empty()){
//					varPerTrans.second.writeOutSNPsInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidKnownMutations.tab.txt")), varPerTrans.first,knownMutationsLocations, true	);
//					for(const auto & pop : knownAAMeta){
//						std::string popName = pop.first;
//						std::string transcript = varPerTrans.first;
//						std::vector<std::string> aaPos;
//						for(const auto & variablePos : knownMutationsLocations){
//							aaPos.emplace_back(estd::to_string(variablePos + 1) + "-" + pop.second.at(transcript).getMeta(estd::to_string(variablePos)));
//						}
//						knownAATyped[popName][transcript] = njh::conToStr(aaPos, ":");
//					}
//				}
//				OutputStream outPopHapAminos(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-popHapToAmino.tab.txt")));
//				outPopHapAminos << "h_PopUID" ;
//				VecStr aminoPositionsHeader;
//				for(const auto & variablePos : allLocations){
//					outPopHapAminos << "\t" << varPerTrans.first << "-aa" << variablePos + 1;
//					aminoPositionsHeader.emplace_back(njh::pasteAsStr(varPerTrans.first, "-aa",  variablePos + 1));
//				}
//				outPopHapAminos << std::endl;
//				std::unordered_map<std::string, std::string> popHapAminoTyped;
//				std::unordered_map<std::string, VecStr> popHapAminoTypedRow;
//
//				for(const auto & pop : aaMeta){
//					outPopHapAminos << pop.first.substr(0, pop.first.rfind("_f")) ;
//					std::string typed = varPerTrans.first + "=";
//					std::vector<std::string> aaPos;
//					std::vector<std::string> aa;
//
//					for(const auto & variablePos : allLocations){
//						outPopHapAminos << "\t" << pop.second.getMeta(estd::to_string(variablePos));
//						aaPos.emplace_back(estd::to_string(variablePos + 1) + "-" + pop.second.getMeta(estd::to_string(variablePos)));
//						aa.emplace_back(pop.second.getMeta(estd::to_string(variablePos)));
//
//					}
//					if (!aaPos.empty()) {
//						typed += njh::conToStr(aaPos, ":");
//					} else {
//						typed += "NONE";
//					}
//					//typed +=";";
//
//					popHapAminoTyped[pop.first.substr(0, pop.first.rfind("_f"))] = typed;
//					popHapAminoTypedRow[pop.first.substr(0, pop.first.rfind("_f"))] = aa;
//					outPopHapAminos << std::endl;
//				}
//				auto popMetaTable = seqsToMetaTable(popSeqsPerSamp);
//				popMetaTable.deleteColumn("seq");
//				popMetaTable.deleteColumn("count");
//				addOtherVec(popMetaTable.columnNames_, aminoPositionsHeader);
//				for(auto & row : popMetaTable){
//					MetaDataInName::removeMetaDataInName(row[popMetaTable.getColPos("name")]);
//					auto popName = row[popMetaTable.getColPos("PopUID")];
//					VecStr aminos;
//					if(njh::in(popName, popHapAminoTypedRow)){
//						aminos = popHapAminoTypedRow[popName];
//					}else{
//						//pop uid was untypable, didn't align
//						aminos = VecStr{aminoPositionsHeader.size(), std::string("NA")};
//					}
//					addOtherVec(row, aminos);
//				}
//
//				popMetaTable.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(variantInfoDir, varPerTrans.first + "-popSeqsWithMetaAndVariableAAInfoTable.tab.txt")));
//
////				for(auto & popHapSamp : outPopSeqsPerSamp){
////					MetaDataInName meta(popHapSamp.name_);
////					auto typed = popHapAminoTyped[meta.getMeta("PopUID")].substr(popHapAminoTyped[meta.getMeta("PopUID")].find("=") + 1);
////					if("" ==typed){
////						typed = "NONE";
////					}
////					meta.addMeta(varPerTrans.first, typed);
////					meta.resetMetaInName(popHapSamp.name_);
////				}
//			}
//		}
//
//		{
//			//snps
//			OutputStream outSnpDepthPerSample(njh::files::make_path(variantInfoDir, njh::pasteAsStr("snpDepthPerSample.tab.txt")));
//			outSnpDepthPerSample << "AnalysisName\tsample\tchromosome\tposition\trefBase\tbase\treadDepth" ;
//			VecStr metaLevels;
//			if(nullptr != sampColl.groupMetaData_){
//				metaLevels = getVectorOfMapKeys(sampColl.groupMetaData_->groupData_);
//				for(const auto & meta : metaLevels){
//					outSnpDepthPerSample << "\t" << meta;
//				}
//			}
//			outSnpDepthPerSample << std::endl;
//
//			for( auto & varPerChrom : translatedRes.seqVariants_){
//				varPerChrom.second.writeOutSNPsFinalInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-SNPs.tab.txt")), varPerChrom.first, false);
//				varPerChrom.second.writeOutSNPsAllInfo(  njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-allBases.tab.txt")), varPerChrom.first, false);
//				if(!varPerChrom.second.variablePositons_.empty()){
//					GenomicRegion variableRegion = varPerChrom.second.getVariableRegion();
//					OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-chromosome_variableRegion.bed"))));
//					bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
//				}
//
//				std::map<std::string, MetaDataInName> snpMeta;
//				for(auto & seqName : translatedRes.seqAlns_){
//					for(const auto & variablePos : varPerChrom.second.snpsFinal){
////						std::cout << __FILE__ << " " << __LINE__ << std::endl;
////						std::cout << "seqName.second.front()->gRegion_.start_: " << seqName.second.front().gRegion_.start_ << std::endl;
////						std::cout << "variablePos: " << variablePos.first << std::endl;
//						if(variablePos.first < seqName.second.front().gRegion_.start_ || variablePos.first >= seqName.second.front().gRegion_.end_){
//							snpMeta[seqName.first].addMeta(estd::to_string(variablePos.first), "X", false);
//						} else {
////							std::cout << __FILE__ << " " << __LINE__ << std::endl;
//							auto aa = seqName.second.front().alnQuerySeq_.seq_[getAlnPosForRealPos(seqName.second.front().alnRefSeq_.seq_, variablePos.first - seqName.second.front().gRegion_.start_)];
//							snpMeta[seqName.first].addMeta(estd::to_string(variablePos.first), aa, false);
////							std::cout << __FILE__ << " " << __LINE__ << std::endl;
//						}
//					}
//				}
////				std::cout << __FILE__ << " " << __LINE__ << std::endl;
//				OutputStream outPopHapAminos(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-popHapToSNPs.tab.txt")));
//				outPopHapAminos << "h_PopUID" ;
//				VecStr aminoPositionsHeader;
//				for(const auto & variablePos : varPerChrom.second.snpsFinal){
//					outPopHapAminos << "\t" << varPerChrom.first << "-" << variablePos.first;
//					aminoPositionsHeader.emplace_back(njh::pasteAsStr(varPerChrom.first, "-",  variablePos.first));
//				}
//				outPopHapAminos << std::endl;
//				std::unordered_map<std::string, std::string> popHapAminoTyped;
//				std::unordered_map<std::string, VecStr> popHapAminoTypedRow;
//
//				for(const auto & pop : snpMeta){
//					outPopHapAminos << pop.first.substr(0, pop.first.rfind("_f")) ;
//					std::string typed = varPerChrom.first + "=";
//					std::vector<std::string> aaPos;
//					std::vector<std::string> aa;
//
//					for(const auto & variablePos : varPerChrom.second.snpsFinal){
//						outPopHapAminos << "\t" << pop.second.getMeta(estd::to_string(variablePos.first));
//						aaPos.emplace_back(estd::to_string(variablePos.first) + "-" + pop.second.getMeta(estd::to_string(variablePos.first)));
//						aa.emplace_back(pop.second.getMeta(estd::to_string(variablePos.first)));
//					}
//					if (!aaPos.empty()) {
//						typed += njh::conToStr(aaPos, ":");
//					} else {
//						typed += "NONE";
//					}
//					//typed +=";";
//					popHapAminoTyped[pop.first.substr(0, pop.first.rfind("_f"))] = typed;
//					popHapAminoTypedRow[pop.first.substr(0, pop.first.rfind("_f"))] = aa;
//					outPopHapAminos << std::endl;
//				}
//
//
////				std::cout << __FILE__ << " " << __LINE__ << std::endl;
//				auto popMetaTable = seqsToMetaTable(popSeqsPerSamp);
//				popMetaTable.deleteColumn("seq");
//				popMetaTable.deleteColumn("count");
//				addOtherVec(popMetaTable.columnNames_, aminoPositionsHeader);
//
//				std::map<std::string, std::map<std::string, std::map<uint32_t, std::map<std::string, double>>>> snpsPerSampleDepth;
//
//				for(auto & row : popMetaTable){
//					MetaDataInName::removeMetaDataInName(row[popMetaTable.getColPos("name")]);
//					auto popName = row[popMetaTable.getColPos("PopUID")];
//					VecStr aminos;
//					if(njh::in(popName, popHapAminoTypedRow)){
//						aminos = popHapAminoTypedRow[popName];
//						uint32_t snpPosition = 0;
//						for(const auto & snp : aminos){
//							if("X" == snp){
//								continue;
//							}
//							auto sample = row[popMetaTable.getColPos("sample")];
//							auto chrom = aminoPositionsHeader[snpPosition].substr(0, aminoPositionsHeader[snpPosition].rfind("-"));
////							std::cout << __FILE__ << " " << __LINE__ << std::endl;
////							std::cout << "snpPosition: " << snpPosition << std::endl;
////							std::cout << "aminoPositionsHeader[snpPosition]: " <<  aminoPositionsHeader[snpPosition] << std::endl;
////							std::cout << "aminoPositionsHeader[snpPosition].substr(aminoPositionsHeader[snpPosition].rfind(\"-\") + 1): " << aminoPositionsHeader[snpPosition].substr(aminoPositionsHeader[snpPosition].rfind("-") + 1)<< std::endl;
////							std::cout << "row[popMetaTable.getColPos(\"readCount\")]: " << row[popMetaTable.getColPos("readCount")] << std::endl;
//							//std::cout << aminoPositionsHeader[snpPosition].substr(aminoPositionsHeader[snpPosition].rfind("-") + 1) << std::endl;
//							auto position = njh::StrToNumConverter::stoToNum<uint32_t>(aminoPositionsHeader[snpPosition].substr(aminoPositionsHeader[snpPosition].rfind("-") + 1));
//							double depth = njh::StrToNumConverter::stoToNum<double>(row[popMetaTable.getColPos("readCount")]);
////							std::cout << __FILE__ << " " << __LINE__ << std::endl;
//							snpsPerSampleDepth[sample][chrom][position][snp] += depth;
//							++snpPosition;
//						}
//					} else {
//						//pop uid was untypable, didn't align
//						aminos = VecStr{aminoPositionsHeader.size(), std::string("NA")};
//					}
//					addOtherVec(row, aminos);
//				}
//				//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//				popMetaTable.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(variantInfoDir, varPerChrom.first + "-popSeqsWithMetaAndVariableSNPInfoTable.tab.txt")));
//
//
//				for(const auto & sample : snpsPerSampleDepth){
//					for(const auto & chrom : sample.second){
//						for(const auto & position : chrom.second){
//							for(const auto & snp : position.second){
//								outSnpDepthPerSample
//										<< pars.experimentNames.populationName_
//								    << "\t" << sample.first
//										<< "\t" << chrom.first
//										<< "\t" << position.first
//										<< "\t" << varPerChrom.second.getBaseForGenomicRegion(position.first) //translatedRes.baseForPosition_[varPerChrom.first][position.first]
//										<< "\t" << snp.first
//										<< "\t" << snp.second ;
//
//								if (nullptr != sampColl.groupMetaData_) {
//									auto metaForSample = sampColl.groupMetaData_->getMetaForSample(sample.first, metaLevels);
//									for (const auto & meta : metaLevels) {
//										outSnpDepthPerSample << "\t" << metaForSample.getMeta(meta);
//									}
//								}
//								outSnpDepthPerSample << std::endl;
//							}
//						}
//					}
//				}
////				for(auto & popHapSamp : outPopSeqsPerSamp){
//////					std::cout << __FILE__ << " " << __LINE__ << std::endl;
////					MetaDataInName meta(popHapSamp.name_);
////					auto typed = popHapAminoTyped[meta.getMeta("PopUID")].substr(popHapAminoTyped[meta.getMeta("PopUID")].find("=") + 1);
////					if("" ==typed){
////						typed = "NONE";
////					}
////					meta.addMeta(varPerChrom.first, typed);
////					meta.resetMetaInName(popHapSamp.name_);
////				}
//			}
//		}
//	}

  std::map<std::string, std::string> fullAATyped;
  //if typing file exists, read it in and set in map
  auto seqTypingFnp = njh::files::make_path(varCallDirPath, "variantCalls/seqsAATyped.tab.txt.gz");
  if(bfs::exists(seqTypingFnp)){
//    std::cout << __FILE__ << " " << __PRETTY_FUNCTION__  << " "  << __LINE__ << std::endl;
    std::unordered_map<std::string, std::string> nameKey;
    auto seqNamesFnp = njh::files::make_path(varCallDirPath, "uniqueSeqs_meta.tab.txt.gz");
    {
      VecStr row;
      TableReader seqNamesTab(TableIOOpts::genTabFileIn(seqNamesFnp, true));
      while(seqNamesTab.getNextRow(row)){
//        std::cout << "CollapsedName: " << row[seqNamesTab.header_.getColPos("CollapsedName")] << std::endl;
//        std::cout << "PopUID: " << row[seqNamesTab.header_.getColPos("PopUID")] << std::endl;

        nameKey[row[seqNamesTab.header_.getColPos("CollapsedName")]] = row[seqNamesTab.header_.getColPos("PopUID")];
      }
    }
    {
      VecStr row;
      TableReader seqTypingTab(TableIOOpts::genTabFileIn(seqTypingFnp, true));
      while(seqTypingTab.getNextRow(row)){
        auto collapsedName = row[seqTypingTab.header_.getColPos("name")];
        auto typed = row[seqTypingTab.header_.getColPos("fullTyped")];
        auto popName = nameKey[collapsedName];
//        std::cout << "collapsedName: " << collapsedName << std::endl;
//        std::cout << "typed: " << typed << std::endl;
//        std::cout << "popName: " << popName << std::endl;
        fullAATyped[popName] = typed;
      }
    }
  }

	for(auto & clus : sampColl.popCollapse_->collapsed_.clusters_){
		auto popName = clus.seqBase_.name_.substr(0, clus.seqBase_.name_.rfind("_f"));
		std::string typed = fullAATyped[popName];
// 		std::cout << "clus.seqBase_.name_: " << clus.seqBase_.name_ << std::endl;
//		std::cout << "clus.meta_.toJson(): " << clus.meta_.toJson() << std::endl;
//		std::cout << "popName: " << popName << std::endl;
//    std::cout << njh::conToStr(njh::getVecOfMapKeys(fullAATyped), ",") << std::endl;
		clus.meta_.addMeta("h_AATyped", typed);
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
  auto outPopSeqsPerSampIoOpts = SeqIOOptions(njh::files::make_path(sampColl.masterOutputDir_, "population", "popSeqsWithMetaWtihSampleName"), setUp.pars_.ioOptions_.outFormat_);
  SeqOutput::write(outPopSeqsPerSamp,outPopSeqsPerSampIoOpts);

//	auto popMetaTable = seqsToMetaTable(outPopSeqsPerSamp);
//	for(auto & row : popMetaTable){
//		MetaDataInName::removeMetaDataInName(row[popMetaTable.getColPos("name")]);
//	}
//
//	popMetaTable.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(sampColl.masterOutputDir_, "population", "popSeqsWithMetaWtihSampleNameTable.tab.txt.gz")));


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




	return 0;
}

}  // namespace njhseq
