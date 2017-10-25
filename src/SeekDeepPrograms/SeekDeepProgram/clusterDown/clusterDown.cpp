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



int SeekDeepRunner::qluster(const bib::progutils::CmdArgs & inputCommands) {
	SeekDeepSetUp setUp(inputCommands);
	// parameters
	clusterDownPars pars;

	setUp.setUpClusterDown(pars);
	// make the runLog, this is what is seen on the terminal screen at run time
	setUp.startARunLog(setUp.pars_.directoryName_);
	// parameter file
	setUp.writeParametersFile(setUp.pars_.directoryName_ + "parametersUsed.txt", false,
			true);
	//add some meta data about file and analysis paths so latter trace back an happen
	Json::Value metaData;
	auto analysisDirPath = bib::files::bfs::canonical(setUp.pars_.directoryName_);
	metaData["analysisDirPath"] = bib::json::toJson(analysisDirPath.string());
	auto fullPathToInput = bib::files::bfs::canonical(setUp.pars_.ioOptions_.firstName_);
	metaData["inputFile"] = bib::json::toJson(fullPathToInput);
	metaData["extractionDir"] = bib::json::toJson(fullPathToInput.parent_path());
	if("" != setUp.pars_.ioOptions_.secondName_){
		metaData["inputFile2"] = bib::json::toJson(bib::files::bfs::canonical(setUp.pars_.ioOptions_.secondName_));
	}
	// print out the parameters read in
	if (setUp.pars_.colOpts_.nucCompBinOpts_.useNucComp_ && setUp.pars_.verbose_) {
		std::cout << "Nucleotide Composition Binning cut offs" << std::endl;
		printVector(setUp.pars_.colOpts_.nucCompBinOpts_.diffCutOffVec_, ", ", std::cout);
	}
	setUp.rLog_.setCurrentLapName("initialSetUp");
	setUp.rLog_.logCurrentTime("Reading In Sequences");

	//write out clustering parameters

	std::string parDir = bib::files::makeDir(setUp.pars_.directoryName_, bib::files::MkdirPar("pars")).string();
	std::ofstream parsOutFile;
	openTextFile(parsOutFile, OutOptions(bib::files::join(parDir, "pars.txt")));
	pars.iteratorMap.writePars(parsOutFile);
	if("" != pars.binParameters){
		std::ofstream binParsOutFile;
		openTextFile(binParsOutFile, OutOptions(bib::files::join(parDir, "binPars.txt")));
		pars.binIteratorMap.writePars(binParsOutFile);
	}

	// read in the sequences
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto reads = reader.readAllReads<readObject>();
	setUp.rLog_.logCurrentTime("Various filtering and little modifications");
	auto splitOnSize = readVecSplitter::splitVectorBellowLength(reads,
			pars.smallReadSize);
	reads = splitOnSize.first;
	if (!splitOnSize.second.empty()) {
		SeqOutput::write(splitOnSize.second,SeqIOOptions(setUp.pars_.directoryName_ + "smallReads",
				setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
	}
	if (setUp.pars_.colOpts_.iTOpts_.removeLowQualityBases_) {
		readVec::allRemoveLowQualityBases(reads, setUp.pars_.colOpts_.iTOpts_.lowQualityBaseTrim_);
	}
	if (setUp.pars_.colOpts_.iTOpts_.adjustHomopolyerRuns_) {
		readVec::allAdjustHomopolymerRunsQualities(reads);
	}
	bool containsCompReads = false;
	int compCount = 0;
	readVec::getCountOfReadNameContaining(reads, "_Comp", compCount);
	if (compCount > 0) {
		containsCompReads = true;
	}

	readVecSorter::sortReadVector(reads, pars.sortBy);
	// get the count of reads read in and the max length so far
	int counter = readVec::getTotalReadCount(reads);
	uint64_t maxSize = 0;
	readVec::getMaxLength(reads, maxSize);
	std::vector<readObject> refSequences;
	if (setUp.pars_.refIoOptions_.firstName_ != "") {
		refSequences = SeqInput::getReferenceSeq(setUp.pars_.refIoOptions_, maxSize);
	}
	maxSize = maxSize * 2;
	// calculate the runCutoff if necessary
	processRunCutoff(setUp.pars_.colOpts_.kmerOpts_.runCutOff_, setUp.pars_.colOpts_.kmerOpts_.runCutOffString_,
			readVec::getTotalReadCount(reads));
	if (setUp.pars_.verbose_ && !pars.onPerId) {
		std::cout << "Kmer Low Frequency Error Cut off Is: " << setUp.pars_.colOpts_.kmerOpts_.runCutOff_
				<< std::endl;
	}

	// read in the paramteres from the parameters file
	setUp.rLog_ << "Parameters used" << "\n";
	pars.iteratorMap.writePars(setUp.rLog_.runLogFile_);
	if (setUp.pars_.verbose_) {
		std::cout << "Reading clusters from " << setUp.pars_.ioOptions_.firstName_ << " "
				<< ("" == setUp.pars_.ioOptions_.secondName_ ? "": setUp.pars_.ioOptions_.secondName_.string()) << std::endl;
		std::cout << "Read in " << counter << " reads" << std::endl;
	}
	setUp.rLog_ << "Reading clusters from " << setUp.pars_.ioOptions_.firstName_ << " "
			<< setUp.pars_.ioOptions_.secondName_ << "\n";
	setUp.rLog_ << "Read in " << counter << " reads" << "\n";
	setUp.rLog_.logCurrentTime("Collapsing to unique sequences");
	// create cluster vector
	std::vector<identicalCluster> identicalClusters;
	std::vector<cluster> clusters;
	if (setUp.pars_.ioOptions_.processed_) {
		clusters = baseCluster::convertVectorToClusterVector<cluster>(reads);
	} else {
		identicalClusters = clusterCollapser::collapseIdenticalReads(reads,
				pars.qualRep);
		for (const auto & read : identicalClusters) {
			clusters.push_back(cluster(read.seqBase_));
		}
		//clusters = baseCluster::convertVectorToClusterVector<cluster>(identicalClusters);
	}

	if (setUp.pars_.verbose_) {
		std::cout << "Unique clusters numbers: " << clusters.size() << std::endl;
	}
	setUp.rLog_ << "Unique clusters numbers: " << clusters.size() << "\n";
	std::sort(clusters.begin(), clusters.end());
	//readVecSorter::sortReadVector(clusters, sortBy);
	setUp.rLog_.logCurrentTime("Indexing kmers");
	KmerMaps kMaps = indexKmers(clusters, setUp.pars_.colOpts_.kmerOpts_.kLength_, setUp.pars_.colOpts_.kmerOpts_.runCutOff_,
			setUp.pars_.colOpts_.kmerOpts_.kmersByPosition_, setUp.pars_.expandKmerPos_, setUp.pars_.expandKmerSize_);
	setUp.rLog_.logCurrentTime("Creating aligner");
	// create aligner class object
	aligner alignerObj(maxSize,
			setUp.pars_.gapInfo_, setUp.pars_.scoring_, kMaps,
			setUp.pars_.qScorePars_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
	if (setUp.pars_.verbose_ && !pars.onPerId) {
		std::cout << bib::bashCT::bold << "Primary Qual: "
				<< alignerObj.qScorePars_.primaryQual_ << std::endl;
		std::cout << "Secondary Qual: " << alignerObj.qScorePars_.secondaryQual_
				<< bib::bashCT::reset << std::endl;
	}
	if (setUp.pars_.debug_) {
		std::ofstream scoreArrayFile;
		openTextFile(scoreArrayFile, "scoreArrayTest.tab.txt", ".txt", true, false);
		std::ofstream scoreMapFile;
		openTextFile(scoreMapFile, "scoreMapFile.tab.txt", ".txt", true, false);
		for (const auto & row : iter::range(len(alignerObj.parts_.scoring_.mat_))) {
			std::vector<int32_t> currentRow;
			for (const auto & col : iter::range(
					len(alignerObj.parts_.scoring_.mat_[row]))) {
				currentRow.emplace_back(alignerObj.parts_.scoring_.mat_[row][col]);
			}
			printVector(currentRow, "\t", scoreArrayFile);
		}
		std::cout << alignerObj.parts_.gapScores_.toJson() << std::endl;
	}
	setUp.rLog_.logCurrentTime("Reading in previous alignments");
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
	setUp.rLog_.logCurrentTime("Removing singlets");
	collapser collapserObj = collapser(setUp.pars_.colOpts_);

	uint32_t singletonNum = 0;
	std::vector<cluster> singletons;
	if (!pars.startWithSingles) {
		if(setUp.pars_.verbose_){
			std::cout << "Removing Singlets for initial analysis" << std::endl;
		}
		for (auto& clus : clusters) {
			if (clus.seqBase_.cnt_ <= 1.01) {
				clus.remove = true;
			}
		}
		clusters = readVecSplitter::splitVectorOnRemoveAdd(clusters, singletons,
				singletonNum, "none", false);
		for (auto& clus : singletons) {
			clus.remove = false;
		}
		if(setUp.pars_.verbose_){
			std::cout << "Removed " << singletons.size() << " singlets" << std::endl;
		}
	}
	setUp.rLog_.logCurrentTime("Running initial clustering");
	//run clustering
	pars.snapShotsOpts_.snapShotsDirName_ = "firstSnaps";
	collapserObj.runFullClustering(clusters, pars.intialParameters,
			pars.binIteratorMap, alignerObj, setUp.pars_.directoryName_,
			setUp.pars_.ioOptions_, setUp.pars_.refIoOptions_, pars.snapShotsOpts_);
	//run again with singlets if needed
	if (!pars.startWithSingles && !pars.leaveOutSinglets) {
		setUp.rLog_.logCurrentTime("Running singlet clustering");
		addOtherVec(clusters, singletons);
		pars.snapShotsOpts_.snapShotsDirName_ = "secondSnaps";
		collapserObj.runFullClustering(clusters, pars.iteratorMap,
				pars.binIteratorMap, alignerObj, setUp.pars_.directoryName_,
				setUp.pars_.ioOptions_, setUp.pars_.refIoOptions_, pars.snapShotsOpts_);
	}

	//remove reads if they are made up of reads only in one direction
	if (containsCompReads && pars.useCompPerCutOff) {
		for (auto& clus : clusters) {
			if (containsCompReads) {
				uint32_t currentCompAmount = 0;
				for(const auto & read : clus.reads_){
					if(bib::containsSubString(read->seqBase_.name_, "_Comp")){
						currentCompAmount+= read->seqBase_.cnt_;
					}
				}
				double currentCompPer = currentCompAmount / clus.seqBase_.cnt_;
				if (currentCompPer > pars.compPerCutOff && clus.seqBase_.cnt_ != 1) {
					clus.remove = true;
				}
			}
		}
		auto splitClus = readVecSplitter::splitVectorOnRemove(clusters);
		if (!splitClus.second.empty()) {
			clusters = splitClus.first;
			clusterVec::allSetFractionClusters(clusters);
			SeqOutput::write(splitClus.second,SeqIOOptions(setUp.pars_.directoryName_ + "allOneDirection",
					setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
		}
	}

	readVecSorter::sortReadVector(clusters, "totalCount");
	std::string seqName = bfs::basename(setUp.pars_.ioOptions_.firstName_);
	readVecSorter::sort(clusters);
	renameReadNames(clusters, seqName, true, false, false);


	if (setUp.pars_.chiOpts_.checkChimeras_) {
		setUp.rLog_.logCurrentTime("Checking chimeras");
		std::ofstream chimerasInfoFile;
		openTextFile(chimerasInfoFile,
				setUp.pars_.directoryName_ + "chimeraNumberInfo.txt", ".txt", false, false);
		chimerasInfoFile << "#chimericClusters\t#chimericReads" << std::endl;
		setUp.pars_.chiOpts_.chiOverlap_.largeBaseIndel_ = .99;

//		collapserObj.opts_.verboseOpts_.verbose_ = true;
//		collapserObj.opts_.verboseOpts_.debug_ = true;
		auto chiInfoTab = collapserObj.markChimeras(clusters, alignerObj,
				setUp.pars_.chiOpts_);

		chiInfoTab.outPutContents(
				TableIOOpts(
						OutOptions(setUp.pars_.directoryName_ + "chiParentsInfo.txt",
								".txt"), "\t", true));
		/*clusterCollapser::markChimerasAdvanced(
		 clusters, alignerObj,pars.parFreqs, 1, setUp.pars_.local_, chiOverlap,
		 overLapSizeCutoff, setUp.pars_.weightHomopolymers_, chiCount, allowableError);*/

		int clusterCount = 0;
		int readCount = 0;
		readVec::getCountOfReadNameContaining(clusters, "CHI", clusterCount);
		readVec::getReadCountOfReadNameContaining(clusters, "CHI", readCount);
		chimerasInfoFile << getPercentageString(clusterCount, clusters.size())
				<< "\t"
				<< getPercentageString(readCount, readVec::getTotalReadCount(clusters))
				<< std::endl;
		if (setUp.pars_.verbose_) {
			std::cout << "Marked " << clusterCount << " as chimeric" << std::endl;
		}
	}

	if (pars.collapsingTandems) {
		setUp.rLog_.logCurrentTime("Collapsing tandems");
		SeqOutput::write(clusters,
				SeqIOOptions(
						setUp.pars_.directoryName_
								+ setUp.pars_.ioOptions_.out_.outFilename_.string() + "_befroeTanCol",
						setUp.pars_.ioOptions_.outFormat_, setUp.pars_.ioOptions_.out_));
		std::cout << "Collapsing on tandem repeat gaps" << std::endl;
		std::cout << "Starting with " << clusters.size() << " clusters"
				<< std::endl;
		clusterCollapser::collapseTandems(clusters, alignerObj,
				setUp.pars_.colOpts_.kmerOpts_.runCutOff_,
				setUp.pars_.colOpts_.kmerOpts_.kLength_,
				setUp.pars_.colOpts_.kmerOpts_.kmersByPosition_, setUp.pars_.chiOpts_.parentFreqs_,
				setUp.pars_.local_, true);
		clusters = readVecSplitter::splitVectorOnRemove(clusters).first;
		std::cout << "Collapsed down to " << clusters.size() << " clusters"
				<< std::endl;
	}
	setUp.rLog_.logCurrentTime("Writing outputs");
	if (setUp.pars_.refIoOptions_.firstName_ == "") {
		profiler::getFractionInfoCluster(clusters, setUp.pars_.directoryName_,
				"outputInfo");
	} else {
		profiler::getFractionInfoCluster(clusters, setUp.pars_.directoryName_,
				"outputInfo", setUp.pars_.refIoOptions_.firstName_.string(), alignerObj,
				setUp.pars_.local_);
	}
	std::ofstream startingInfo;
	openTextFile(startingInfo, setUp.pars_.directoryName_ + "startingInfo.txt", ".txt",
			false, false);
	startingInfo << "cluster\tstartingClusterName\tstartingSize" << std::endl;
	for (const auto& clus : clusters) {
		VecStr toks = tokenizeString(clus.firstReadName_, "_t");
		startingInfo << clus.seqBase_.name_ << "\t" << clus.firstReadName_ << "\t"
				<< toks.back() << std::endl;
	}
	if (pars.additionalOut) {
		std::string additionalOutDir = findAdditonalOutLocation(
				pars.additionalOutLocationFile, setUp.pars_.ioOptions_.firstName_.string());
		if (additionalOutDir == "") {
			std::cerr << bib::bashCT::red << bib::bashCT::bold;
			std::cerr << "No additional out directory found for: "
					<< setUp.pars_.ioOptions_.firstName_ << std::endl;
			std::cerr << bib::bashCT::reset;
		} else {
			SeqOutput::write(clusters, SeqIOOptions(additionalOutDir + setUp.pars_.ioOptions_.out_.outFilename_.string(),
					setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
			std::ofstream metaDataFile;
			openTextFile(metaDataFile, additionalOutDir + "/" + "metaData", ".json",
					setUp.pars_.ioOptions_.out_);
			metaDataFile << metaData;
		}
	}

	SeqOutput::write(clusters,
			SeqIOOptions(
					setUp.pars_.directoryName_ + setUp.pars_.ioOptions_.out_.outFilename_.string(),
					setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
	if(pars.writeOutFinalInternalSnps){
		setUp.rLog_.logCurrentTime("Calling internal snps");
		std::string snpDir = bib::files::makeDir(setUp.pars_.directoryName_,
				bib::files::MkdirPar("internalSnpInfo", false)).string();
		for (const auto & readPos : iter::range(clusters.size())) {
			std::unordered_map<uint32_t,
					std::unordered_map<char, std::vector<baseReadObject>>>mismatches;

			for (const auto & subReadPos : iter::range(
							clusters[readPos].reads_.size())) {
				const auto & subRead = clusters[readPos].reads_[subReadPos];
				alignerObj.alignCacheGlobal(clusters[readPos], subRead);
				//count gaps and mismatches and get identity
				alignerObj.profilePrimerAlignment(clusters[readPos], subRead);
				for (const auto & m : alignerObj.comp_.distances_.mismatches_) {
					mismatches[m.second.refBasePos][m.second.seqBase].emplace_back(
							subRead->seqBase_);
				}
			}
			table misTab {VecStr {"refPos", "refBase", "seqBase", "freq",
					"fraction", "seqs", "clusterName"}};
			for (const auto & m : mismatches) {
				for (const auto & seqM : m.second) {
					double totalCount = 0;
					for(const auto & read : seqM.second) {
						totalCount += read.seqBase_.cnt_;
					}
					misTab.content_.emplace_back(
							toVecStr(m.first, clusters[readPos].seqBase_.seq_[m.first],
									seqM.first, totalCount,
									totalCount / clusters[readPos].seqBase_.cnt_,
									vectorToString(readVec::getNames(seqM.second), ","),
									clusters[readPos].seqBase_.name_)
					);
				}
			}
			misTab.sortTable("seqBase", false);
			misTab.sortTable("refPos", false);
			misTab.outPutContents(
					TableIOOpts(OutOptions(snpDir + clusters[readPos].seqBase_.name_,
							".tab.txt"), "\t", misTab.hasHeader_));
		}
	}
	if (pars.createMinTree) {
		setUp.rLog_.logCurrentTime("Creating minimum spanning trees");
		std::string minTreeDirname = bib::files::makeDir(setUp.pars_.directoryName_,
				bib::files::MkdirPar("minTree", false)).string();
		auto clusSplit = readVecSplitter::splitVectorOnReadFraction(clusters,
				0.005);
		std::vector<readObject> tempReads;
		for (const auto & clus : clusSplit.first) {
			tempReads.emplace_back(readObject(clus.seqBase_));
		}
	  std::unordered_map<std::string, std::unique_ptr<aligner>> aligners;
	  std::mutex alignerLock;
	  uint32_t numThreads = 2;
		auto graph = genReadComparisonGraph(tempReads, alignerObj, aligners, alignerLock,
				numThreads);
		std::vector<std::string> popNames;
		for (const auto & n : graph.nodes_) {
			if (n->on_) {
				popNames.emplace_back(n->name_);
			}
		}
		auto nameColors = getColorsForNames(popNames);
		comparison maxEvents = graph.setMinimumEventConnections();
		if(setUp.pars_.debug_){
			std::cout << maxEvents.toJson() << std::endl;
		}
		auto treeData = graph.toD3Json(bib::color("#000000"), nameColors);
		std::ofstream outJson(minTreeDirname + "tree.json");
		std::ofstream outHtml(minTreeDirname + "tree.html");
		outJson << treeData;
		auto outHtmlStr = genHtmlStrForPsuedoMintree("tree.json",
				"http://bib8.umassmed.edu/~hathawan/js/psuedoMinTreeWithIndels.js");
		outHtml << outHtmlStr;
	}



	if (pars.writeOutInitalSeqs) {
		std::string clusterDirectoryName = bib::files::makeDir(setUp.pars_.directoryName_,
				bib::files::MkdirPar("clusters", false)).string();
		/*
		 * 	clusterVec::allWriteClustersInDir(clusters, clusterDirectoryName,
		 setUp.pars_.ioOptions_);
		 */
		SeqIOOptions clustersIoOpts(
				bib::files::make_path(clusterDirectoryName, "initialClusters"),
				setUp.pars_.ioOptions_.outFormat_);
		SeqOutput subClusterWriter(clustersIoOpts);
		subClusterWriter.openOut();
		for( const auto & clus : clusters){
			MetaDataInName clusMeta;
			clusMeta.addMeta("clusterName", clus.seqBase_.name_);
			for( const auto & subClus : clus.reads_){
				seqInfo subClusCopy = subClus->seqBase_;
				clusMeta.resetMetaInName(subClusCopy.name_, subClusCopy.name_.find("_t"));
				subClusterWriter.write(subClusCopy);
			}
		}
	}



	if (!setUp.pars_.ioOptions_.processed_) {
		std::ofstream compStats;
		if (containsCompReads) {
			openTextFile(compStats, setUp.pars_.directoryName_ + "compStats.tab.txt",
					".txt", false, false);
			compStats << "cluster\tcompAmount" << std::endl;
		}
		std::string allInputReadsDir;
		if(pars.writeOutInitalSeqs){
			allInputReadsDir = bib::files::makeDir(setUp.pars_.directoryName_,
					bib::files::MkdirPar("allInputReadsForEachCluster", false)).string();
		}
		SeqIOOptions clustersIoOpts(
				bib::files::make_path(allInputReadsDir, "allInitialReads"),
				setUp.pars_.ioOptions_.outFormat_);
		SeqOutput subClusterWriter(clustersIoOpts);
		if(pars.writeOutInitalSeqs){
			subClusterWriter.openOut();
		}

		for (const auto& clus : clusters) {
			MetaDataInName clusMeta;
			clusMeta.addMeta("clusterName", clus.seqBase_.name_);
			uint32_t currentCompAmount = 0;
			for (const auto & seq : clus.reads_) {
				auto input = readVec::getReadByName(identicalClusters, seq->seqBase_.name_);
				for( auto & inputRead : input.reads_){
					if (bib::containsSubString(inputRead->seqBase_.name_, "_Comp")) {
						currentCompAmount += inputRead->seqBase_.cnt_;
					}
					if(pars.writeOutInitalSeqs){
						seqInfo subClusCopy = inputRead->seqBase_;
						clusMeta.resetMetaInName(subClusCopy.name_);
						subClusterWriter.write(subClusCopy);
					}
				}
			}
			if (containsCompReads) {
				compStats << clus.seqBase_.name_ << "\t"
						<< getPercentageString(currentCompAmount, clus.seqBase_.cnt_)
						<< std::endl;
			}
		}
	}

	if (!pars.startWithSingles && pars.leaveOutSinglets) {
		SeqOutput::write(singletons, SeqIOOptions(setUp.pars_.directoryName_ + "singletons",
				setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
		std::ofstream singletonsInfoFile;
		openTextFile(singletonsInfoFile, setUp.pars_.directoryName_ + "singletonsInfo",
				".tab.txt", false, true);
		singletonsInfoFile << "readCnt\treadFrac\n";
		singletonsInfoFile << singletons.size() << "\t"
				<< singletons.size() / static_cast<double>(reads.size())
				<< std::endl;
	}

	if (setUp.pars_.writingOutAlnInfo_) {
		setUp.rLog_.logCurrentTime("Writing previous alignments");
		alignerObj.alnHolder_.write(setUp.pars_.outAlnInfoDirName_);
	}
	//log number of alignments done
	setUp.rLog_ << "Number of Alignments Done: "
			<< alignerObj.numberOfAlingmentsDone_ << "\n";
	if (setUp.pars_.verbose_) {
		std::cout << "Number of Alignments Done: "
				<< alignerObj.numberOfAlingmentsDone_ << std::endl;
		//log time
		setUp.logRunTime(std::cout);
	}

	return 0;
}



}  // namespace bib
