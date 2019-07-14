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



int SeekDeepRunner::clusterDown(const njh::progutils::CmdArgs & inputCommands) {
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
	auto analysisDirPath = njh::files::bfs::canonical(setUp.pars_.directoryName_);
	metaData["analysisDirPath"] = njh::json::toJson(analysisDirPath.string());
	auto fullPathToInput = njh::files::bfs::canonical(setUp.pars_.ioOptions_.firstName_);
	metaData["inputFile"] = njh::json::toJson(fullPathToInput);
	metaData["extractionDir"] = njh::json::toJson(fullPathToInput.parent_path());
	if("" != setUp.pars_.ioOptions_.secondName_){
		metaData["inputFile2"] = njh::json::toJson(njh::files::bfs::canonical(setUp.pars_.ioOptions_.secondName_));
	}
	// print out the parameters read in
	if (setUp.pars_.colOpts_.nucCompBinOpts_.useNucComp_ && setUp.pars_.verbose_) {
		std::cout << "Nucleotide Composition Binning cut offs" << std::endl;
		printVector(setUp.pars_.colOpts_.nucCompBinOpts_.diffCutOffVec_, ", ", std::cout);
	}
	setUp.rLog_.setCurrentLapName("initialSetUp");
	setUp.rLog_.logCurrentTime("Reading In Sequences");

	//write out clustering parameters

	std::string parDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("pars")).string();
	std::ofstream parsOutFile;
	openTextFile(parsOutFile, OutOptions(njh::files::join(parDir, "pars.txt")));
	pars.iteratorMap.writePars(parsOutFile);
	if("" != pars.binParameters){
		std::ofstream binParsOutFile;
		openTextFile(binParsOutFile, OutOptions(njh::files::join(parDir, "binPars.txt")));
		pars.binIteratorMap.writePars(binParsOutFile);
	}
	bool containsCompReads = false;
	int compCount = 0;
	// read in the sequences
	setUp.rLog_.logCurrentTime("Various filtering and little modifications");
	auto inputOpts = setUp.pars_.ioOptions_;

	//
	std::unordered_map<std::string, uint32_t> sampleNumberCounts;
	if(!pars.dontFilterToMostCommonIlluminaSampleNumber_){
		uint32_t totalInputCount = 0;
		{
			SeqInput counterIo(inputOpts);
			counterIo.openIn();
			seqInfo seq;
			while(counterIo.readNextRead(seq)){
				IlluminaNameFormatDecoder decoder(seq.name_, pars.IlluminaSampleRegPatStr_, pars.IlluminaSampleNumberPos_);
				if(0 == decoder.match_.size()){
					decoder = IlluminaNameFormatDecoder(seq.name_, pars.BackUpIlluminaSampleRegPatStr_, pars.BackUpIlluminaSampleNumberPos_);
				}
				++sampleNumberCounts[decoder.getSampleNumber()];
				++totalInputCount;
				if(totalInputCount > pars.useCutOff){
					break;
				}
			}
		}
	}
	VecStr sampleNames;
	if(!pars.dontFilterToMostCommonIlluminaSampleNumber_){
		sampleNames = getVectorOfMapKeys(sampleNumberCounts);
		njh::sort(sampleNames,[&sampleNumberCounts](const std::string & k1, const std::string & k2){
			if(sampleNumberCounts[k1] == sampleNumberCounts[k2]){
				return k1 > k2;
			}else{
				return sampleNumberCounts[k1] > sampleNumberCounts[k2];
			}
		});
	}


	//downsample input file to save on memory usage
	bfs::path downsampledFnp;
	if(!pars.useAllInput){
		//for limiting large number of input sequences
		uint32_t totalInputCount = 0;
		{
			SeqInput counterIo(inputOpts);
			counterIo.openIn();
			seqInfo seq;
			while(counterIo.readNextRead(seq)){
				if(!pars.dontFilterToMostCommonIlluminaSampleNumber_){
					IlluminaNameFormatDecoder decoder(seq.name_, pars.IlluminaSampleRegPatStr_, pars.IlluminaSampleNumberPos_);
					if(0 == decoder.match_.size()){
						decoder = IlluminaNameFormatDecoder(seq.name_, pars.BackUpIlluminaSampleRegPatStr_, pars.BackUpIlluminaSampleNumberPos_);
					}
					if(decoder.getSampleNumber() != sampleNames.front()){
						continue;
					}
				}
				++totalInputCount;
			}
		}
		if(setUp.pars_.verbose_){
			std::cout << "totalInputCount: " << totalInputCount << std::endl;
		}
		if(totalInputCount > pars.useCutOff){

			std::vector<uint32_t> randomSel;
			{
				njh::randomGenerator rGen;
				std::vector<uint32_t> allIndexes(totalInputCount);
				njh::iota<uint32_t>(allIndexes, 0);
				std::shuffle( std::begin(allIndexes), std::end(allIndexes), rGen.mtGen_) ;
				randomSel = std::vector<uint32_t>{ std::begin(allIndexes), std::begin(allIndexes)+pars.useCutOff} ;
			}
			SeqIOOptions outOpts(njh::files::make_path(setUp.pars_.directoryName_, "downsampledFile"), SeqIOOptions::getOutFormat(inputOpts.inFormat_));
			SeqOutput writer(outOpts);
			writer.openOut();

			SeqInput reader(inputOpts);
			reader.openIn();
			uint32_t seqCount = 0;
			seqInfo seq;

			while(reader.readNextRead(seq)){
				if(!pars.dontFilterToMostCommonIlluminaSampleNumber_){
					IlluminaNameFormatDecoder decoder(seq.name_, pars.IlluminaSampleRegPatStr_, pars.IlluminaSampleNumberPos_);
					if(0 == decoder.match_.size()){
						decoder = IlluminaNameFormatDecoder(seq.name_, pars.BackUpIlluminaSampleRegPatStr_, pars.BackUpIlluminaSampleNumberPos_);
					}
					if(decoder.getSampleNumber() != sampleNames.front()){
						continue;
					}
				}
				if(njh::in(seqCount, randomSel)){
					writer.write(seq);
				}
				++seqCount;
			}
			downsampledFnp = outOpts.out_.outFilename_.string() + outOpts.getOutExtension();
			inputOpts.firstName_ = downsampledFnp;
		}
	}

	bfs::path tempOutFnp;
	if(setUp.pars_.ioOptions_.inFormat_ == SeqIOOptions::inFormats::FASTAGZ ||
		 setUp.pars_.ioOptions_.inFormat_ == SeqIOOptions::inFormats::FASTQGZ){
		if (setUp.pars_.ioOptions_.inFormat_ == SeqIOOptions::inFormats::FASTAGZ) {
			inputOpts.inFormat_ = SeqIOOptions::inFormats::FASTA;
			tempOutFnp =njh::files::make_path(setUp.pars_.directoryName_, "tempinput.fasta");
		} else if (setUp.pars_.ioOptions_.inFormat_ == SeqIOOptions::inFormats::FASTQGZ) {
			inputOpts.inFormat_ = SeqIOOptions::inFormats::FASTQ;
			tempOutFnp =njh::files::make_path(setUp.pars_.directoryName_, "tempinput.fastq");
		}
		OutputStream tempOut(tempOutFnp);
		InputStream in(inputOpts.firstName_);
		std::string line;
		while(njh::files::crossPlatGetline(in, line)){
			tempOut << line << "\n";
		}
		inputOpts.firstName_ = tempOutFnp;
	}



	SeqInput reader(inputOpts);
	reader.openIn();
	std::vector<cluster> clusters;
	uint32_t counter = 0;
	std::unordered_map<uint32_t, std::vector<uint64_t>> filepositions;
	std::unordered_map<std::string, uint32_t> clusterNameToFilePosKey;
	{
		SeqIOOptions smallSeqOpts(setUp.pars_.directoryName_ + "smallReads",
						setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_);
		SeqOutput smallWriter(smallSeqOpts);
		seqInfo seq;
		uint64_t fPos = reader.tellgPri();
		while(reader.readNextRead(seq)){
			if(!pars.dontFilterToMostCommonIlluminaSampleNumber_){
				IlluminaNameFormatDecoder decoder(seq.name_, pars.IlluminaSampleRegPatStr_, pars.IlluminaSampleNumberPos_);
				if(0 == decoder.match_.size()){
					decoder = IlluminaNameFormatDecoder(seq.name_, pars.BackUpIlluminaSampleRegPatStr_, pars.BackUpIlluminaSampleNumberPos_);
				}
				if(decoder.getSampleNumber() != sampleNames.front()){
					fPos = reader.tellgPri();
					continue;
				}
			}
			++counter;
			if(len(seq) <= pars.smallReadSize){
				smallWriter.openWrite(seq);
			}else{
				if (setUp.pars_.colOpts_.iTOpts_.removeLowQualityBases_) {
					seq.removeLowQualityBases(setUp.pars_.colOpts_.iTOpts_.lowQualityBaseTrim_);
				}
				if (setUp.pars_.colOpts_.iTOpts_.adjustHomopolyerRuns_) {
					seq.adjustHomopolymerRunQualities();
				}
				if(std::string::npos != seq.name_.find("_Comp")){
					++compCount;
				}
			  readVec::handelLowerCaseBases(seq, setUp.pars_.ioOptions_.lowerCaseBases_);
			  if (setUp.pars_.ioOptions_.removeGaps_) {
			    seq.removeGaps();
			  }
			  bool found = false;
			  for(const auto & clusPos : iter::range(clusters.size())){
			  	if(seq.seq_ == clusters[clusPos].seqBase_.seq_){
			  		clusters[clusPos].seqBase_.cnt_ += seq.cnt_;
			  		clusters[clusPos].firstReadCount_ += seq.cnt_;
			  		found = true;
			  		filepositions[clusPos].emplace_back(fPos);
			  		break;
			  	}
			  }
			  if(!found){
		  		filepositions[clusters.size()].emplace_back(fPos);
			  	clusters.emplace_back(seq);
			  }
			}
			fPos = reader.tellgPri();
		}

		reader.reOpenIn();
		//now calculate the qualities if fastq
		for(const auto & fPositions : filepositions){
			if(setUp.pars_.ioOptions_.inFormat_ != SeqIOOptions::inFormats::FASTA && setUp.pars_.ioOptions_.inFormat_ != SeqIOOptions::inFormats::FASTAGZ){
				//first iterator over the files positions for the seq and read in their qualities
				std::vector<std::vector<uint32_t>> qualities(clusters[fPositions.first].seqBase_.qual_.size());
				for(const auto seqPos : fPositions.second){
					reader.seekgPri(seqPos);
					reader.readNextRead(seq);
					if (setUp.pars_.colOpts_.iTOpts_.removeLowQualityBases_) {
						seq.removeLowQualityBases(setUp.pars_.colOpts_.iTOpts_.lowQualityBaseTrim_);
					}
					if (setUp.pars_.colOpts_.iTOpts_.adjustHomopolyerRuns_) {
						seq.adjustHomopolymerRunQualities();
					}
					if(std::string::npos != seq.name_.find("_Comp")){
						++compCount;
					}
				  readVec::handelLowerCaseBases(seq, setUp.pars_.ioOptions_.lowerCaseBases_);
				  if (setUp.pars_.ioOptions_.removeGaps_) {
				    seq.removeGaps();
				  }

					for(const auto i : iter::range(seq.qual_.size())){
						qualities[i].emplace_back(seq.qual_[i]);
					}
				}

				//calculate qualities
		    if (pars.qualRep == "worst") {
		    	clusters[fPositions.first].seqBase_.qual_.clear();
				  for (const auto i : iter::range(clusters[fPositions.first].seqBase_.seq_.size())) {
				  	clusters[fPositions.first].seqBase_.qual_.push_back(vectorMinimum(qualities[i]));
				  }
		    } else if (pars.qualRep == "median") {
		    	clusters[fPositions.first].seqBase_.qual_.clear();
				  for (const auto i : iter::range(clusters[fPositions.first].seqBase_.seq_.size())) {
				  	clusters[fPositions.first].seqBase_.qual_.push_back(vectorMedianRef(qualities[i]));
				  }
		    } else if (pars.qualRep == "average") {
		    	clusters[fPositions.first].seqBase_.qual_.clear();
				  for (const auto i : iter::range(clusters[fPositions.first].seqBase_.seq_.size())) {
				  	clusters[fPositions.first].seqBase_.qual_.push_back(vectorMean(qualities[i]));
				  }
		    } else if (pars.qualRep == "bestQual") {
		    	clusters[fPositions.first].seqBase_.qual_.clear();
				  for (const auto i : iter::range(clusters[fPositions.first].seqBase_.seq_.size())) {
				  	clusters[fPositions.first].seqBase_.qual_.push_back(vectorMaximum(qualities[i]));
				  }
		    } else {
		    	std::stringstream ss;
		      ss << __PRETTY_FUNCTION__ << ' '<< "Unrecognized qualRep: " << pars.qualRep << std::endl;
		      ss << "Needs to be median, average, bestQual, or worst"
		                << std::endl;
		      throw std::runtime_error{ss.str()};
		    }
			}
			clusters[fPositions.first].averageErrorRate = clusters[fPositions.first].getAverageErrorRate();
			clusters[fPositions.first].updateName();
			clusters[fPositions.first].reads_.front()->averageErrorRate = clusters[fPositions.first].getAverageErrorRate();
			clusters[fPositions.first].reads_.front()->seqBase_ = clusters[fPositions.first].seqBase_;
			clusterNameToFilePosKey[clusters[fPositions.first].seqBase_.name_] = fPositions.first;
		}
	}
	if(!pars.countIlluminaSampleNumbers_ && !pars.writeOutInitalSeqs){
		filepositions.clear();
		clusterNameToFilePosKey.clear();
	}

	if (compCount > 0) {
		containsCompReads = true;
	}

	uint64_t maxSize = 0;
	readVec::getMaxLength(clusters, maxSize);

	if (setUp.pars_.verbose_) {
		std::cout << "Reading clusters from " << setUp.pars_.ioOptions_.firstName_ << " "
				<< ("" == setUp.pars_.ioOptions_.secondName_ ? "": setUp.pars_.ioOptions_.secondName_.string()) << std::endl;
		std::cout << "Read in " << counter << " reads" << std::endl;
	}
	setUp.rLog_ << "Reading clusters from " << setUp.pars_.ioOptions_.firstName_ << " "
			<< setUp.pars_.ioOptions_.secondName_ << "\n";
	setUp.rLog_ << "Read in " << counter << " reads" << "\n";
	setUp.rLog_.logCurrentTime("Collapsing to unique sequences");



	if (setUp.pars_.verbose_) {
		std::cout << "Unique clusters numbers: " << clusters.size() << std::endl;
	}
	setUp.rLog_ << "Unique clusters numbers: " << clusters.size() << "\n";
	std::sort(clusters.begin(), clusters.end());

	std::vector<readObject> refSequences;
	if (setUp.pars_.refIoOptions_.firstName_ != "") {
		refSequences = SeqInput::getReferenceSeq(setUp.pars_.refIoOptions_, maxSize);
	}
	maxSize = maxSize * 2;
	// calculate the runCutoff if necessary
	processRunCutoff(setUp.pars_.colOpts_.kmerOpts_.runCutOff_, setUp.pars_.colOpts_.kmerOpts_.runCutOffString_,
			readVec::getTotalReadCount(clusters));
	if (setUp.pars_.verbose_ && !pars.onPerId) {
		std::cout << "Kmer Low Frequency Error Cut off Is: " << setUp.pars_.colOpts_.kmerOpts_.runCutOff_
				<< std::endl;
	}

	// read in the parameters from the parameters file
	setUp.rLog_ << "Parameters used" << "\n";
	pars.iteratorMap.writePars(setUp.rLog_.runLogFile_);



	//readVecSorter::sortReadVector(clusters, sortBy);
	setUp.rLog_.logCurrentTime("Indexing kmers");
	KmerMaps kMaps = indexKmers(clusters,
			setUp.pars_.colOpts_.kmerOpts_.kLength_,
			setUp.pars_.colOpts_.kmerOpts_.runCutOff_,
			setUp.pars_.colOpts_.kmerOpts_.kmersByPosition_,
			setUp.pars_.expandKmerPos_,
			setUp.pars_.expandKmerSize_);
	setUp.rLog_.logCurrentTime("Creating aligner");
	// create aligner class object
	aligner alignerObj(maxSize,
			setUp.pars_.gapInfo_, setUp.pars_.scoring_, kMaps,
			setUp.pars_.qScorePars_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
	if (setUp.pars_.verbose_ && !pars.onPerId) {
		std::cout << njh::bashCT::bold << "Primary Qual: "
				<< alignerObj.qScorePars_.primaryQual_ << std::endl;
		std::cout << "Secondary Qual: " << alignerObj.qScorePars_.secondaryQual_
				<< njh::bashCT::reset << std::endl;
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
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
	if(setUp.pars_.verbose_){
		uint32_t alignmentsReadIn = 0;
		for(const auto & alnHold : alignerObj.alnHolder_.globalHolder_){
			alignmentsReadIn += alnHold.second.infos_.size();
		}
		for(const auto & alnHold : alignerObj.alnHolder_.localHolder_){
			alignmentsReadIn += alnHold.second.infos_.size();
		}
		std::cout << "Global Alignments Holders: " << alignerObj.alnHolder_.globalHolder_.size() << std::endl;
		std::cout << "Local Alignments Holders: " << alignerObj.alnHolder_.localHolder_.size() << std::endl;
		std::cout << "Read in: " << alignmentsReadIn << "alignments" << std::endl;
	}
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
					if(njh::containsSubString(read->seqBase_.name_, "_Comp")){
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
		auto fnp = setUp.pars_.ioOptions_.firstName_.filename().string();
		if(njh::endsWith(fnp, ".gz")){
			fnp = fnp.substr(0, fnp.rfind(".gz"));
		}
		std::string additionalOutDir = findAdditonalOutLocation(
				pars.additionalOutLocationFile, fnp);
		if (additionalOutDir == "") {
			std::cerr << njh::bashCT::red << njh::bashCT::bold;
			std::cerr << "No additional out directory found for: "
					<< setUp.pars_.ioOptions_.firstName_ << std::endl;
			std::cerr << njh::bashCT::reset;
			std::cerr << processFileNameForID(setUp.pars_.ioOptions_.firstName_.string())<< std::endl;
			std::cerr << "Options:" << std::endl;
			table inTab(pars.additionalOutLocationFile, "\t");
		  MapStrStr additionalOutNames;
		  for (const auto& fIter : inTab.content_) {
		    additionalOutNames[makeIDNameComparable(fIter[0])] = fIter[1];
		  }
		  std::cerr << njh::conToStr(njh::getVecOfMapKeys(additionalOutNames)) << std::endl;
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
		std::string snpDir = njh::files::makeDir(setUp.pars_.directoryName_,
				njh::files::MkdirPar("internalSnpInfo", false)).string();
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
		std::string minTreeDirname = njh::files::makeDir(setUp.pars_.directoryName_,
				njh::files::MkdirPar("minTree", false)).string();
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
		auto treeData = graph.toD3Json(njh::color("#000000"), nameColors);
		std::ofstream outJson(minTreeDirname + "tree.json");
		std::ofstream outHtml(minTreeDirname + "tree.html");
		outJson << treeData;
		auto outHtmlStr = genHtmlStrForPsuedoMintree("tree.json",
				"http://njh8.umassmed.edu/~hathawan/js/psuedoMinTreeWithIndels.js");
		outHtml << outHtmlStr;
	}

	if(pars.countIlluminaSampleNumbers_){
		OutputStream sampleCountsOut(njh::files::make_path(setUp.pars_.directoryName_, "illuminaSampleNumbersCounts.tab.txt"));
		SeqInput reader(inputOpts);
		reader.openIn();
		seqInfo inSeq;
		sampleCountsOut << "SeqName\tsampleNumber\tcount\tfrac" << std::endl;
		for (const auto& clus : clusters) {
			std::unordered_map<std::string, uint32_t> sampleNumberCounts;
			double total = 0;
			for (const auto & seq : clus.reads_) {

				for (const auto seqPos : filepositions[clusterNameToFilePosKey[seq->seqBase_.name_]]) {
					reader.seekgPri(seqPos);
					reader.readNextRead(inSeq);
					IlluminaNameFormatDecoder nameDecoded(inSeq.name_, pars.IlluminaSampleRegPatStr_, pars.IlluminaSampleNumberPos_);
					if(0 == nameDecoded.match_.size()){
						nameDecoded = IlluminaNameFormatDecoder(inSeq.name_, pars.BackUpIlluminaSampleRegPatStr_, pars.BackUpIlluminaSampleNumberPos_);
					}
					++sampleNumberCounts[nameDecoded.getSampleNumber()];
					++total;
				}
			}

			VecStr sampleNumbersNames = getVectorOfMapKeys(sampleNumberCounts);
			njh::sort(sampleNumbersNames,[&sampleNumberCounts](const std::string &k1, const std::string &k2){
				if(sampleNumberCounts[k1] == sampleNumberCounts[k2]){
					return k1 > k2;
				}else{
					return sampleNumberCounts[k1] > sampleNumberCounts[k2];
				}
			});
			for(const auto & k : sampleNumbersNames){
				sampleCountsOut << clus.seqBase_.name_
						<< "\t" << k
						<< "\t" << sampleNumberCounts[k]
						<< "\t" << sampleNumberCounts[k]/total<< std::endl;
			}
		}
	}

	if (pars.writeOutInitalSeqs) {
		std::string clusterDirectoryName = njh::files::makeDir(setUp.pars_.directoryName_,
				njh::files::MkdirPar("clusters", false)).string();
		/*
		 * 	clusterVec::allWriteClustersInDir(clusters, clusterDirectoryName,
		 setUp.pars_.ioOptions_);
		 */
		SeqIOOptions clustersIoOpts(
				njh::files::make_path(clusterDirectoryName, "initialClusters"),
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



	if (pars.writeOutInitalSeqs) {
		std::ofstream compStats;
		if (containsCompReads) {
			openTextFile(compStats, setUp.pars_.directoryName_ + "compStats.tab.txt",
					".txt", false, false);
			compStats << "cluster\tcompAmount" << std::endl;
		}
		std::string allInputReadsDir;
		allInputReadsDir = njh::files::makeDir(setUp.pars_.directoryName_,
				njh::files::MkdirPar("allInputReadsForEachCluster", false)).string();
		SeqIOOptions clustersIoOpts(
				njh::files::make_path(allInputReadsDir, "allInitialReads"),
				setUp.pars_.ioOptions_.outFormat_);
		SeqOutput subClusterWriter(clustersIoOpts);
		subClusterWriter.openOut();
		SeqInput reader(inputOpts);
		reader.openIn();
		seqInfo inSeq;
		for (const auto& clus : clusters) {
			MetaDataInName clusMeta;
			clusMeta.addMeta("clusterName", clus.seqBase_.name_);
			uint32_t currentCompAmount = 0;
			for (const auto & seq : clus.reads_) {
				for (const auto seqPos : filepositions[clusterNameToFilePosKey[seq->seqBase_.name_]]) {
					reader.seekgPri(seqPos);
					reader.readNextRead(inSeq);
					if (setUp.pars_.colOpts_.iTOpts_.removeLowQualityBases_) {
						inSeq.removeLowQualityBases(setUp.pars_.colOpts_.iTOpts_.lowQualityBaseTrim_);
					}
					if (setUp.pars_.colOpts_.iTOpts_.adjustHomopolyerRuns_) {
						inSeq.adjustHomopolymerRunQualities();
					}
					if(std::string::npos != inSeq.name_.find("_Comp")){
						++compCount;
					}
				  readVec::handelLowerCaseBases(inSeq, setUp.pars_.ioOptions_.lowerCaseBases_);
				  if (setUp.pars_.ioOptions_.removeGaps_) {
				    inSeq.removeGaps();
				  }
					if (njh::containsSubString(inSeq.name_, "_Comp")) {
						currentCompAmount += inSeq.cnt_;
					}
					if (pars.writeOutInitalSeqs) {
						clusMeta.resetMetaInName(inSeq.name_);
						subClusterWriter.write(inSeq);
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
				<< singletons.size() / static_cast<double>(counter)
				<< std::endl;
	}

	if (setUp.pars_.writingOutAlnInfo_) {
		setUp.rLog_.logCurrentTime("Writing previous alignments");
		alignerObj.alnHolder_.write(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);
	}
	if(bfs::exists(tempOutFnp)){
		bfs::remove(tempOutFnp);
	}
	if(bfs::exists(downsampledFnp) && !pars.keepDownSampledFile){
		bfs::remove(downsampledFnp);
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



}  // namespace njh
