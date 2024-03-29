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
//  main.cpp
//  SeekDeep
//
//  Created by Nicholas Hathaway on 8/11/13.
//

#include "SeekDeepPrograms/SeekDeepProgram/SeekDeepRunner.hpp"
namespace njhseq {


int SeekDeepRunner::extractorPairedEnd(const njh::progutils::CmdArgs & inputCommands) {
	SeekDeepSetUp setUp(inputCommands);
	ExtractorPairedEndPars pars;
	//std::cout << "pars.corePars_.sampleName: " << pars.corePars_.sampleName << std::endl;
	setUp.setUpExtractorPairedEnd(pars);
	//std::cout << "pars.corePars_.sampleName: " << pars.corePars_.sampleName << std::endl;
	// run log
	setUp.startARunLog(setUp.pars_.directoryName_);
	// parameter file
	setUp.writeParametersFile(setUp.pars_.directoryName_ + "parametersUsed.txt", false, false);
	PrimersAndMids ids(pars.corePars_.primIdsPars.idFile_);


	ids.checkIfMIdsOrPrimersReadInThrow(__PRETTY_FUNCTION__);
	if(setUp.pars_.verbose_){
		if(ids.getMids().size()> 0){
			std::cout << "Found: " << ids.getMids().size() << " MIDs to de-multiplex on" << std::endl;
		}
		if(ids.getTargets().size() > 0){
			std::cout << "Found: " << ids.getTargets().size() << " target primer pairs to de-multiplex on" << std::endl;
		}
	}

	//init mids
	//init primers
	//add in any length cuts if any
	//add in ref sequences if any
	ids.initAllAddLenCutsRefs(pars.corePars_.primIdsPars);
	//add in overlap status
	if(!pars.defaultStatuses_.empty()){
		ids.addOverLapStatuses(pars.defaultStatuses_);
	}else{
		ids.addOverLapStatuses(pars.corePars_.primIdsPars.overlapStatusFnp_);
	}

	if (pars.corePars_.noPrimers_) {
		if (nullptr == ids.pDeterminator_
				|| ids.pDeterminator_->primers_.size() != 1) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__
					<< ", error if setting --noPrimers then there should be just target listed after the target/gene header in id file" << "\n"
					<< "the forward and reverse primer sequences columns will be ignored but a name must still appear"
					<< "\n";
			throw std::runtime_error { ss.str() };
		}
	}

	//default checks for Ns and quality
	ReadCheckerOnSeqContaining nChecker("N", pars.corePars_.numberOfNs, true);
	ReadCheckerQualCheck qualChecker(pars.corePars_.qPars_.qualCheck_, pars.corePars_.qPars_.qualCheckCutOff_, true);


	// make some directories for outputs
	bfs::path unfilteredReadsDir = njh::files::makeDir(setUp.pars_.directoryName_,
			njh::files::MkdirPar("unfilteredReads", false));
	bfs::path unfilteredByBarcodesDir = njh::files::makeDir(unfilteredReadsDir,
			njh::files::MkdirPar("byBarcodes", false));
	bfs::path unfilteredByPrimersDir = njh::files::makeDir(unfilteredReadsDir,
			njh::files::MkdirPar("byPrimers", false));
	bfs::path unfilteredByPairsProcessedDir = njh::files::makeDir(unfilteredReadsDir,
			njh::files::MkdirPar("pairsProcessed", false));
	bfs::path filteredOffDir = njh::files::makeDir(setUp.pars_.directoryName_,
			njh::files::MkdirPar("filteredOff", false));
	bfs::path overHansDir;
	if (pars.pairProcessorParams_.writeOverHangs_) {
		overHansDir = njh::files::makeDir(setUp.pars_.directoryName_,
				njh::files::MkdirPar("overhangs", false));
	}
	bfs::path badDir = njh::files::makeDir(filteredOffDir,
			njh::files::MkdirPar("bad", false));
	bfs::path unrecognizedPrimerDir = njh::files::makeDir(filteredOffDir,
			njh::files::MkdirPar("unrecognizedPrimer", false));
	bfs::path contaminationDir = "";


	PairedRead seq;

	// read in reads and remove lower case bases indicating tech low quality like
	// tags and such

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto smallOpts = setUp.pars_.ioOptions_;
	smallOpts.out_.outFilename_ = njh::files::make_path(badDir, "smallFragments").string();
	SeqIO smallFragMentOut(smallOpts);
	//smallFragMentOut.openOut();
	auto startsWtihBadQualOpts = setUp.pars_.ioOptions_;
	startsWtihBadQualOpts.out_.outFilename_ = njh::files::make_path(badDir , "startsWtihBadQual").string();
	SeqOutput startsWtihBadQualOut(startsWtihBadQualOpts);

	uint32_t smallFragmentCount = 0;
	uint64_t maxReadSize = 0;

	MultiSeqIO readerOuts;

	std::map<std::string, std::pair<uint32_t, uint32_t>> counts;
	std::unordered_map<std::string, uint32_t>  failBarCodeCounts;
	std::unordered_map<std::string, uint32_t>  failBarCodeCountsPossibleContamination;

	if (ids.containsMids() && pars.corePars_.primIdsPars.mPars_.allowableErrors_> 0) {
		if (setUp.pars_.debug_) {
			std::cout << "Allowing " << pars.corePars_.primIdsPars.mPars_.allowableErrors_ << " errors in barcode" << std::endl;
		}
	}
	if (ids.containsMids()) {
		for (const auto & mid : ids.mDeterminator_->mids_) {
			auto midOpts = setUp.pars_.ioOptions_;
			midOpts.out_.outFilename_ = njh::files::make_path(unfilteredByBarcodesDir, mid.first).string();
			if (setUp.pars_.debug_) {
				std::cout << "Inserting: " << mid.first << std::endl;
			}
			readerOuts.addReader(mid.first, midOpts);
		}
	} else {
		auto midOpts = setUp.pars_.ioOptions_;
		midOpts.out_.outFilename_ = njh::files::make_path(unfilteredByBarcodesDir, "all").string();
		if (setUp.pars_.debug_) {
			std::cout << "Inserting: " << "all" << std::endl;
		}
		readerOuts.addReader("all", midOpts);
	}

	if (ids.containsMids()) {
		auto failureCases = MidDeterminator::ProcessedRes::getProcessedCaseNames();
		for(const auto & failureCase : failureCases){
			std::string unRecName = "unrecognizedBarcode_" + failureCase;
			auto midOpts = setUp.pars_.ioOptions_;
			midOpts.out_.outFilename_ = njh::files::make_path(badDir, unRecName).string();
			if (setUp.pars_.debug_){
				std::cout << "Inserting: " << unRecName << std::endl;
			}
			readerOuts.addReader(unRecName, midOpts);
			//std::cout << "midOpts.outFormat_: " << SeqIOOptions::getOutFormat(midOpts.outFormat_) << std::endl;
		}
	}


	if(setUp.pars_.verbose_){
		std::cout << njh::bashCT::boldGreen("Extracting on MIDs") << std::endl;
	}
	if(setUp.pars_.verbose_){
		std::cout << "Reading in reads:" << std::endl;
	}
	uint32_t count = 0;
	uint32_t readsNotMatchedToBarcode = 0;
	uint32_t unrecognizedPrimers = 0;
	uint32_t mismatchedPrimers = 0;
	uint32_t mismatchedPrimersDimers = 0;
	uint32_t contamination = 0;
	uint32_t qualityFilters = 0;
	uint32_t used = 0;

	std::string seqName = bfs::basename(setUp.pars_.ioOptions_.firstName_);
	seqName = seqName.substr(0,seqName.find("_"));


	while (reader.readNextRead(seq)) {
//		std::cout << seq.seqBase_.name_ << std::endl;
//		bool print = false;
//		if("M02551:63:000000000-D3YB2:1:1102:10373:15072 1:N:0:1" == seq.seqBase_.name_){
//			print = true;
//		}
		//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		++count;
		if (setUp.pars_.verbose_ && count % 50 == 0) {
			std::cout << "\r" << count ;
			std::cout.flush();
		}
		readVec::handelLowerCaseBases(seq, setUp.pars_.ioOptions_.lowerCaseBases_);

		if (len(seq) < pars.corePars_.smallFragmentCutoff) {
			smallFragMentOut.openWrite(seq);
			++smallFragmentCount;
			continue;
		}
		readVec::getMaxLength(seq.seqBase_, maxReadSize);
		readVec::getMaxLength(seq.mateSeqBase_, maxReadSize);

		if (ids.containsMids()) {
			auto searchRes = ids.mDeterminator_->searchPairedEndRead(seq);

			auto processRes = ids.mDeterminator_->processSearchPairedEndRead(seq, searchRes);
//			if(MidDeterminator::ProcessedRes::PROCESSED_CASE::MISMATCHING_MIDS == processRes.case_){
//				std::cout << seq.seqBase_.name_ << std::endl;
//				std::cout << "Forward: ";
//				for(const auto & f : searchRes.forward_){
//					std::cout << f.toJson() << std::endl << std::endl;
//				}
//				std::cout << "Reverse: ";
//				for(const auto & f : searchRes.reverse_){
//					std::cout << f.toJson() << std::endl << std::endl;
//				}
//				exit(1);
//			}
			if(MidDeterminator::ProcessedRes::PROCESSED_CASE::MATCH == processRes.case_){
				if (processRes.rcomplement_) {
					++counts[processRes.midName_].second;
				} else {
					++counts[processRes.midName_].first;
				}
				readerOuts.openWrite(processRes.midName_, seq);
			}else{
				std::string unRecName = "unrecognizedBarcode_" + MidDeterminator::ProcessedRes::getProcessedCaseName(processRes.case_);
				++readsNotMatchedToBarcode;
				++failBarCodeCounts[MidDeterminator::ProcessedRes::getProcessedCaseName(processRes.case_)];
				readerOuts.openWrite(unRecName, seq);
			}
		} else {
			++counts["all"].first;
			readerOuts.openWrite("all", seq);
		}
	}
	if(ids.containsMids()){
		OutOptions barcodeCountOpts(njh::files::make_path(setUp.pars_.directoryName_, "midCounts.tab.txt"));
		OutputStream barcodeCountOut(barcodeCountOpts);
		barcodeCountOut << "inputName\tMID\tforwardCount\tforwardCountPerc\treverseCount\treverseCountPerc\ttotal\tfraction" << std::endl;
		std::set<std::string> barKeysSet{"all"};
		if(ids.containsMids()){
			auto inputMNames = ids.getMids();
			barKeysSet = std::set<std::string>(inputMNames.begin(), inputMNames.end());
		}
		for(const auto & count : counts){
			barKeysSet.emplace(count.first);
		}
		for(const auto & countkey : barKeysSet){
			double total = counts[countkey].first + counts[countkey].second;
			if(total > 0){
				barcodeCountOut
						<< seqName
						<< "\t" << countkey
						<< "\t" << counts[countkey].first
						<< "\t" << 100 * (counts[countkey].first/total)
						<< "\t" << counts[countkey].second
						<< "\t" << 100 * (counts[countkey].second/total)
						<< "\t" << total
						<< "\t" << total/count << std::endl;
			} else {
				barcodeCountOut << seqName
						<< "\t" << countkey
						<< "\t" << "0"
						<< "\t" << "0"
						<< "\t" << "0"
						<< "\t" << "0"
						<< "\t" << "0"
						<< "\t" << "0" << std::endl;
			}
		}
	}
	if (setUp.pars_.verbose_) {
		std::cout << std::endl;
	}
	//close mid outs;
	readerOuts.closeOutAll();


	std::ofstream renameKeyFile;
	if (pars.corePars_.rename) {
		openTextFile(renameKeyFile, setUp.pars_.directoryName_ + "renameKey.tab.txt", ".tab.txt", false, false);
		renameKeyFile << "originalName\tnewName\n";
	}

	// create aligner for primer identification
	if(setUp.pars_.debug_){
		std::cout << "setUp.pars_.generalMatch_: " << setUp.pars_.generalMatch_ << std::endl;
		std::cout << "setUp.pars_.generalMismatch_: " << setUp.pars_.generalMismatch_ << std::endl;
	}


	auto scoreMatrix = substituteMatrix::createDegenScoreMatrixLessN(
			setUp.pars_.generalMatch_, setUp.pars_.generalMismatch_);

//	auto scoreMatrix = substituteMatrix::createDegenScoreMatrixNoNInRef(
//			setUp.pars_.generalMatch_, setUp.pars_.generalMismatch_);

	gapScoringParameters gapPars(setUp.pars_.gapInfo_);
	KmerMaps emptyMaps;
	bool countEndGaps = false;
	//to avoid allocating an extremely large aligner matrix;

	aligner alignObj(maxReadSize, gapPars, scoreMatrix, emptyMaps,
			setUp.pars_.qScorePars_, countEndGaps, false);

	alignObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
	ExtractionStator stats(count, readsNotMatchedToBarcode, 0, smallFragmentCount);
	std::map<std::string, uint32_t> allPrimerCounts;
	std::map<std::string, uint32_t> matchingPrimerCounts;
	std::vector<std::string> expectedSamples;
	if(ids.containsMids()){
		expectedSamples = getVectorOfMapKeys(ids.mDeterminator_->mids_);
	}else{
		expectedSamples = {"all"};
	}
	auto barcodeFiles = njh::files::listAllFiles(unfilteredByBarcodesDir, false, VecStr { });
	ReadPairsOrganizer prOrg(expectedSamples);
	prOrg.processFiles(barcodeFiles);
	auto readsByPairs = prOrg.processReadPairs();
	auto barcodeReadsKeys = getVectorOfMapKeys(readsByPairs);
	njh::sort(barcodeReadsKeys);
	if(setUp.pars_.verbose_){
		std::cout << "The following mids were found with reads" << std::endl;
		printVector(barcodeReadsKeys);
	}

	std::unordered_map<std::string, std::set<std::string>> primersInMids;

	//now post processing
	PairedReadProcessor pairProcessor(pars.pairProcessorParams_);
	auto alnGapPars = gapScoringParameters(
			10,
			1,
			0,0,
			0,0);

	PairedReadProcessor::ProcessedResultsCounts mismatchedPrimerPairProcessCounts;

	auto pairedProcessingScoring = substituteMatrix::createScoreMatrix(2, -2, true, true, true);
	aligner processingPairsAligner(maxReadSize, alnGapPars, pairedProcessingScoring, false);
	processingPairsAligner.qScorePars_.qualThresWindow_ = 0;
	std::unordered_map<std::string, std::vector<uint32_t>> lengthsPerStitchedTarget;
	std::unordered_map<std::string, std::pair<std::string, PairedReadProcessor::ProcessedResultsCounts>> resultsPerMidTarPair;


	std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, std::vector<uint32_t>>>>> failedPrimer;

	for (const auto & barcodeReadsKey : barcodeReadsKeys) {
		const auto & barcodeReadPairs = readsByPairs[barcodeReadsKey];
		std::string barcodeName = barcodeReadsKey;
		if (setUp.pars_.verbose_) {
			if (ids.containsMids()) {
				std::cout
						<< njh::bashCT::boldGreen("Filtering on barcode: " + barcodeName)
						<< std::endl;
			} else {
				std::cout << njh::bashCT::boldGreen("Filtering") << std::endl;
			}
		}

		//create outputs
		MultiSeqIO midReaderOuts;
		auto unrecogPrimerOutOpts = setUp.pars_.ioOptions_;
		unrecogPrimerOutOpts.out_.outFilename_ = njh::files::make_path(
				unrecognizedPrimerDir, barcodeName).string();
		midReaderOuts.addReader("unrecognized", unrecogPrimerOutOpts);
		for (const auto & primerName : getVectorOfMapKeys(
				ids.pDeterminator_->primers_)) {
			//determine full name
			std::string fullname = primerName;
			if (ids.containsMids()) {
				fullname += barcodeName;
			} else if (pars.corePars_.sampleName != "") {
				fullname += pars.corePars_.sampleName;
			}
			//std::cout << "fullname: " << fullname << std::endl;
			//bad out
			auto badDirOutOpts = setUp.pars_.ioOptions_;
			badDirOutOpts.out_.outFilename_ = njh::files::make_path(badDir, fullname).string();
			midReaderOuts.addReader(fullname + "bad", badDirOutOpts);
			//good out
			auto goodDirOutOpts = setUp.pars_.ioOptions_;
			goodDirOutOpts.out_.outFilename_ = njh::files::make_path(unfilteredByPrimersDir,  fullname);
			midReaderOuts.addReader(fullname + "good", goodDirOutOpts);
		}

		uint32_t barcodeCount = 1;
		njh::ProgressBar pbar(counts[barcodeName].first + counts[barcodeName].second);
		pbar.progColors_ = pbar.RdYlGn_;

		SeqIOOptions barcodePairsReaderOpts = SeqIOOptions::genPairedIn(
				barcodeReadPairs.first.front(), barcodeReadPairs.second.front());
		SeqInput barcodePairsReader(barcodePairsReaderOpts);
		barcodePairsReader.openIn();
		bool primerCheckComplement = pars.corePars_.pDetPars.checkComplement_;
		//turn on auto determination for dual barcoded system where the barcodes are the same since direction cannot be determined by the MID barcode
		if("all" != barcodeName && ids.containsMids() && ids.mDeterminator_->mids_.at(barcodeName).forSameAsRev_){
			primerCheckComplement = true;
		}
		while(barcodePairsReader.readNextRead(seq)){
			//std::cout << barcodeCount << std::endl;
			if (setUp.pars_.verbose_) {
				pbar.outputProgAdd(std::cout, 1, true);
			}
			++barcodeCount;
			//find primers
			//forward
			std::string forwardPrimerName = "";
			bool foundInReverse = false;
			if(pars.corePars_.noPrimers_){
				forwardPrimerName = ids.pDeterminator_->primers_.begin()->first;
			}else{
				forwardPrimerName = ids.pDeterminator_->determineForwardPrimer(seq.seqBase_, pars.corePars_.pDetPars, alignObj);
				if ("unrecognized" ==  forwardPrimerName && primerCheckComplement) {
					forwardPrimerName = ids.pDeterminator_->determineForwardPrimer(seq.mateSeqBase_, pars.corePars_.pDetPars, alignObj);
					if (seq.mateSeqBase_.on_) {
						foundInReverse = true;
					}
				}
			}



			//reverse primer
			std::string reversePrimerName = "";
//			std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//			std::cout << "\tbarcodeName: " << barcodeName << std::endl;

			if(pars.corePars_.noPrimers_){
				reversePrimerName = ids.pDeterminator_->primers_.begin()->first;
			}else{
				if (!foundInReverse) {
					reversePrimerName = ids.pDeterminator_->determineWithReversePrimer(seq.mateSeqBase_, pars.corePars_.pDetPars, alignObj);
				} else {
					reversePrimerName = ids.pDeterminator_->determineWithReversePrimer(seq.seqBase_,     pars.corePars_.pDetPars, alignObj);
				}
			}



			std::string fullname = "";
			if(forwardPrimerName != reversePrimerName){
				fullname = forwardPrimerName + "-" + reversePrimerName;
			}else{
				fullname = forwardPrimerName;
			}
			if (ids.containsMids()) {
				fullname += barcodeName;
			} else if ("" != pars.corePars_.sampleName) {
				fullname += pars.corePars_.sampleName;
			}


			if("unrecognized" == forwardPrimerName ||
				 "unrecognized" == reversePrimerName){
				//check for unrecognized primers
				stats.increaseFailedForward(barcodeName, seq.seqBase_.name_);
				midReaderOuts.openWrite("unrecognized", seq);
				//++allPrimerCounts[fullname];
				++unrecognizedPrimers;
				seq.mateSeqBase_.reverseComplementRead(false, true);
				seq.mateRComplemented_ = true;
				auto mismatchedPairRes = pairProcessor.processPairedEnd(seq, mismatchedPrimerPairProcessCounts, processingPairsAligner);
				if(mismatchedPairRes.status_ != PairedReadProcessor::ReadPairOverLapStatus::NONE &&
						mismatchedPairRes.status_ != PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP){
					readVec::handelLowerCaseBases(mismatchedPairRes.combinedSeq_, "remove");
					failedPrimer[barcodeName][forwardPrimerName][reversePrimerName][PairedReadProcessor::getOverlapStatusStr(mismatchedPairRes.status_)].emplace_back(len(*mismatchedPairRes.combinedSeq_));
				}else{
					failedPrimer[barcodeName][forwardPrimerName][reversePrimerName][PairedReadProcessor::getOverlapStatusStr(mismatchedPairRes.status_)].emplace_back(0);
				}
			}else if(forwardPrimerName != reversePrimerName){
				//check for primer mismatch
				//stats.increasePrimerMismatch(barcodeName, seq.seqBase_.name_);
				stats.increaseCounts(barcodeName, seq.seqBase_.name_, ExtractionStator::extractCase::MISMATCHPRIMERS);
				auto badDirOutOpts = setUp.pars_.ioOptions_;
				badDirOutOpts.out_.outFilename_ =
						njh::files::make_path(badDir, fullname).string();
				if(pars.corePars_.keepFilteredOff){
					if(!midReaderOuts.containsReader(fullname + "bad")){
						midReaderOuts.addReader(fullname + "bad", badDirOutOpts);
					}
					midReaderOuts.openWrite(fullname + "bad", seq);
				}
				//++allPrimerCounts[fullname];
				seq.mateSeqBase_.reverseComplementRead(false, true);
				seq.mateRComplemented_ = true;
				auto mismatchedPairRes = pairProcessor.processPairedEnd(seq, mismatchedPrimerPairProcessCounts, processingPairsAligner);
				if(mismatchedPairRes.status_ != PairedReadProcessor::ReadPairOverLapStatus::NONE &&
						mismatchedPairRes.status_ != PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP){
					readVec::handelLowerCaseBases(mismatchedPairRes.combinedSeq_, "remove");
					failedPrimer[barcodeName][forwardPrimerName][reversePrimerName][PairedReadProcessor::getOverlapStatusStr(mismatchedPairRes.status_)].emplace_back(len(*mismatchedPairRes.combinedSeq_));
				}else{
					failedPrimer[barcodeName][forwardPrimerName][reversePrimerName][PairedReadProcessor::getOverlapStatusStr(mismatchedPairRes.status_)].emplace_back(0);
				}
				if(PairedReadProcessor::ReadPairOverLapStatus::R1BEGINSINR2 == mismatchedPairRes.status_ &&
						len(*mismatchedPairRes.combinedSeq_)<= pars.pairProcessorParams_.primerDimmerSize_ ){
					++mismatchedPrimersDimers;
				}else{
					++mismatchedPrimers;
				}
			}else{
				primersInMids[barcodeName].emplace(forwardPrimerName);
				//primer match
				if(foundInReverse){
					seq.seqBase_.name_.append("_Comp");
					seq.mateSeqBase_.name_.append("_Comp");
					auto tempMate = seq.seqBase_;
					seq.seqBase_ = seq.mateSeqBase_;
					seq.mateSeqBase_ = tempMate;

				}
				stats.increaseCounts(fullname, seq.seqBase_.name_,
						ExtractionStator::extractCase::GOOD);
				midReaderOuts.openWrite(fullname + "good", seq);
				++allPrimerCounts[fullname];
				++matchingPrimerCounts[fullname];
			}
		}
	}

//	{
//		auto primerCountsKeys = getVectorOfMapKeys(allPrimerCounts);
//		njh::sort(primerCountsKeys);
//		OutOptions allPrimerCountsOpts(njh::files::make_path(setUp.pars_.directoryName_, "allPrimerCounts.tab.txt"));
//		OutputStream allPrimerCountsOut(allPrimerCountsOpts);
//		allPrimerCountsOut << "inputName\tname\tcount" << std::endl;
//		for (const auto & good : allPrimerCounts) {
//			allPrimerCountsOut << seqName << '\t'
//					<< good.first << "\t" << good.second << std::endl;
//		}
//	}
	{
		OutOptions allFailedPrimerCountsOpts(njh::files::make_path(setUp.pars_.directoryName_, "allFailedPrimerCounts.tab.txt"));
		OutputStream allFailedPrimerCountsOut(allFailedPrimerCountsOpts);
		allFailedPrimerCountsOut << "inputName\tMIDName\tforwardPrimer\treversePrimer\toverlapStatus\tcount\tMIDtotal\tavgLength" << "\n";
		for(const auto & samp : failedPrimer){
			for(const auto & forw : samp.second){
				for(const auto & rev : forw.second){
					for(const auto & status : rev.second){
						std::string meanStr = "";
						if(status.first != "NONE" && status.first != "NOOVERLAP"){
							meanStr = estd::to_string(vectorMean(status.second));
						}
						allFailedPrimerCountsOut
						<< seqName
								<< "\t" << samp.first
								<< "\t" << forw.first
								<< "\t" << rev.first
								<< "\t" << status.first
								<< "\t" << status.second.size()
								<< "\t" << counts[samp.first].first + counts[samp.first].second
								<< "\t" << meanStr << std::endl;
					}
				}
			}
		}
		failedPrimer.clear();
	}

	//key1 = target, key2 = overlap status, value = count of that status
	std::unordered_map<std::string, std::unordered_map<PairedReadProcessor::ReadPairOverLapStatus, uint32_t>> overlapStatusCounts;

	for(const auto & extractedMid : primersInMids){
		for(const auto & extractedPrimer : extractedMid.second){
			std::string name = extractedPrimer + extractedMid.first;
			if(!ids.containsMids() && "" != pars.corePars_.sampleName){
				name = extractedPrimer + pars.corePars_.sampleName;
			}else if(!ids.containsMids() && "all" == extractedMid.first){
				name = extractedPrimer;
			}

			SeqIOOptions currentPairOpts;
			if(setUp.pars_.ioOptions_.inFormat_ == SeqIOOptions::inFormats::FASTQPAIREDGZ){
				currentPairOpts = SeqIOOptions::genPairedInGz(
						njh::files::make_path(unfilteredByPrimersDir, name + "_R1.fastq.gz"),
						njh::files::make_path(unfilteredByPrimersDir, name + "_R2.fastq.gz"));
			}else{
				currentPairOpts = SeqIOOptions::genPairedIn(
						njh::files::make_path(unfilteredByPrimersDir, name + "_R1.fastq"),
						njh::files::make_path(unfilteredByPrimersDir, name + "_R2.fastq"));
			}
			if(pars.corePars_.primIdsPars.noOverlapProcessForNoOverlapStatusTargets_ && njh::in(PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_)) {
				PairedReadProcessor::ProcessedResultsCounts currentProcessResults;
				currentProcessResults.notCombinedOpts = std::make_shared<SeqIOOptions>(currentPairOpts);
				resultsPerMidTarPair[name] = std::make_pair(extractedPrimer, currentProcessResults);
				continue;
			}
			OutOptions currentOutOpts(njh::files::make_path(unfilteredByPairsProcessedDir, name));
			PairedReadProcessor::ProcessorOutWriters processWriter;
			if(SeqIOOptions::inFormats::FASTQPAIREDGZ == setUp.pars_.ioOptions_.inFormat_){
				processWriter.overhangsWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOutGz(njh::files::make_path(overHansDir, name + "_overhangs")));
				if(njh::in(PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_)){
					processWriter.notCombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genPairedOutGz(njh::files::make_path(unfilteredByPairsProcessedDir, name)));
					processWriter.r1BeginsInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOutGz(njh::files::make_path(badDir, name + "_r1BeginsInR2")));
					processWriter.r1EndsInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOutGz(njh::files::make_path(badDir, name + "_r1EndsInR2")));
					processWriter.r1AllInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOutGz(njh::files::make_path(badDir, name + "_r1AllInR2")));
					processWriter.r2AllInR1CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOutGz(njh::files::make_path(badDir, name + "_r2AllInR1")));
				} else {
					//not combined
					processWriter.notCombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genPairedOutGz(njh::files::make_path(badDir, name + "_notCombined")));
					if(njh::in(PairedReadProcessor::ReadPairOverLapStatus::R1BEGINSINR2, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_) ||
						 njh::in(PairedReadProcessor::ReadPairOverLapStatus::PERFECTOVERLAP, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_)){
						processWriter.r1AllInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOutGz(njh::files::make_path(unfilteredByPairsProcessedDir, name + "_r1AllInR2")));
						processWriter.r2AllInR1CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOutGz(njh::files::make_path(unfilteredByPairsProcessedDir, name + "_r2AllInR1")));
					}else{
						processWriter.r1AllInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOutGz(njh::files::make_path(badDir, name + "_r1AllInR2")));
						processWriter.r2AllInR1CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOutGz(njh::files::make_path(badDir, name + "_r2AllInR1")));
					}
					if(njh::in(PairedReadProcessor::ReadPairOverLapStatus::R1BEGINSINR2, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_)){
						processWriter.r1BeginsInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOutGz(njh::files::make_path(unfilteredByPairsProcessedDir, name + "_r1BeginsInR2")));
					}else{
						processWriter.r1BeginsInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOutGz(njh::files::make_path(badDir, name + "_r1BeginsInR2")));
					}
					if(njh::in(PairedReadProcessor::ReadPairOverLapStatus::R1ENDSINR2, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_)){
						processWriter.r1EndsInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOutGz(njh::files::make_path(unfilteredByPairsProcessedDir, name + "_r1EndsInR2")));
					}else{
						processWriter.r1EndsInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOutGz(njh::files::make_path(badDir, name + "_r1EndsInR2")));
					}
					if(njh::in(PairedReadProcessor::ReadPairOverLapStatus::PERFECTOVERLAP, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_)){
						processWriter.perfectOverlapCombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOutGz(njh::files::make_path(unfilteredByPairsProcessedDir, name + "_perfectOverlap")));
					}else{
						processWriter.perfectOverlapCombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOutGz(njh::files::make_path(badDir, name + "_perfectOverlap")));
					}
				}
			}else{
				processWriter.overhangsWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(njh::files::make_path(overHansDir, name + "_overhangs")));
				if(njh::in(PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_)){
					processWriter.notCombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genPairedOut(njh::files::make_path(unfilteredByPairsProcessedDir, name)));
					processWriter.r1BeginsInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(njh::files::make_path(badDir, name + "_r1BeginsInR2")));
					processWriter.r1EndsInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(njh::files::make_path(badDir, name + "_r1EndsInR2")));
					processWriter.r1AllInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(njh::files::make_path(badDir, name + "_r1AllInR2")));
					processWriter.r2AllInR1CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(njh::files::make_path(badDir, name + "_r2AllInR1")));
				} else {
					//not combined
					processWriter.notCombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genPairedOut(njh::files::make_path(badDir, name + "_notCombined")));
					if(njh::in(PairedReadProcessor::ReadPairOverLapStatus::R1BEGINSINR2, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_) ||
						 njh::in(PairedReadProcessor::ReadPairOverLapStatus::PERFECTOVERLAP, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_)){
						processWriter.r1AllInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(njh::files::make_path(unfilteredByPairsProcessedDir, name + "_r1AllInR2")));
						processWriter.r2AllInR1CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(njh::files::make_path(unfilteredByPairsProcessedDir, name + "_r2AllInR1")));
					}else{
						processWriter.r1AllInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(njh::files::make_path(badDir, name + "_r1AllInR2")));
						processWriter.r2AllInR1CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(njh::files::make_path(badDir, name + "_r2AllInR1")));
					}
					if(njh::in(PairedReadProcessor::ReadPairOverLapStatus::R1BEGINSINR2, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_)){
						processWriter.r1BeginsInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(njh::files::make_path(unfilteredByPairsProcessedDir, name + "_r1BeginsInR2")));
					}else{
						processWriter.r1BeginsInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(njh::files::make_path(badDir, name + "_r1BeginsInR2")));
					}
					if(njh::in(PairedReadProcessor::ReadPairOverLapStatus::R1ENDSINR2, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_)){
						processWriter.r1EndsInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(njh::files::make_path(unfilteredByPairsProcessedDir, name + "_r1EndsInR2")));
					}else{
						processWriter.r1EndsInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(njh::files::make_path(badDir, name + "_r1EndsInR2")));
					}
					if(njh::in(PairedReadProcessor::ReadPairOverLapStatus::PERFECTOVERLAP, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_)){
						processWriter.perfectOverlapCombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(njh::files::make_path(unfilteredByPairsProcessedDir, name + "_perfectOverlap")));
					}else{
						processWriter.perfectOverlapCombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(njh::files::make_path(badDir, name + "_perfectOverlap")));
					}
				}
			}
			currentPairOpts.revComplMate_ = true;
			SeqInput currentReader(currentPairOpts);
			currentReader.openIn();
			if(setUp.pars_.verbose_){
				std::cout << "Pair Processing " << name << std::endl;
			}
			auto currentProcessResults = pairProcessor.processPairedEnd(currentReader, processWriter, processingPairsAligner);
			if(setUp.pars_.verbose_){
				std::cout << "Done Pair Processing for " << name << std::endl;
			}

			if(njh::in(PairedReadProcessor::ReadPairOverLapStatus::AUTO, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_)){
				uint32_t highestCount = 0;
				PairedReadProcessor::ReadPairOverLapStatus autoDetStat = PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP;

				overlapStatusCounts[extractedPrimer][PairedReadProcessor::ReadPairOverLapStatus::R1ENDSINR2] += currentProcessResults.r1EndsInR2Combined;
				overlapStatusCounts[extractedPrimer][PairedReadProcessor::ReadPairOverLapStatus::R1BEGINSINR2] += currentProcessResults.r1EndsInR2Combined;

				if(currentProcessResults.r1EndsInR2Combined > highestCount){
					highestCount = currentProcessResults.r1EndsInR2Combined;
					autoDetStat = PairedReadProcessor::ReadPairOverLapStatus::R1ENDSINR2;
				}
				if(currentProcessResults.r1BeginsInR2CombinedAboveCutOff > highestCount){
					highestCount = currentProcessResults.r1BeginsInR2CombinedAboveCutOff;
					autoDetStat = PairedReadProcessor::ReadPairOverLapStatus::R1BEGINSINR2;
				}
				njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_ = {autoDetStat};
			}


			resultsPerMidTarPair[name] = std::make_pair(extractedPrimer,currentProcessResults);
		}
	}

	//determine AUTO overlap statuses
	//this is somewhat dangerous because if even one case of r1endsinr2 or r1beginsinr2 happens that will be the auto detected status but 99% could have not stitched
	//but i don't want to count nooverlap status because if it's a low quality run it could mean much more don't stitch
	for(auto & tar : ids.targets_){
		if(njh::in(PairedReadProcessor::ReadPairOverLapStatus::AUTO,tar.second.overlapStatuses_)){
			uint32_t highestCount = 0;
			tar.second.overlapStatuses_ = {PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP};
			for(const auto & count : overlapStatusCounts[tar.first]){
				if(count.second > highestCount){
					highestCount = count.second;
					tar.second.overlapStatuses_ = {count.first};
				}
			}
		}
	}



	OutputStream outPairProcessOut(njh::files::make_path(setUp.pars_.directoryName_, "processPairsCounts.tab.txt"));
	outPairProcessOut << "inputName\tName\ttarget\ttotal\tnotCombined\tperfectOverlapCombined\tr1BeginsInR2Combined\tr1BeginsInR2CombinedAboveCutOff\tr1EndsInR2Combined\tr1AllInR2Combined\tr2AllInR1Combined" << std::endl;;
	auto processedPairsKeys = njh::getVecOfMapKeys(resultsPerMidTarPair);
	njh::sort(processedPairsKeys);

	for(const auto & processedResultsKey : processedPairsKeys){
		const auto & processedResults = resultsPerMidTarPair[processedResultsKey];
		outPairProcessOut << seqName
				<< "\t" << processedResultsKey
				<< "\t" << processedResults.first
				<< "\t" << processedResults.second.total
				<< "\t" << getPercentageString(processedResults.second.overlapFail, processedResults.second.total)
				<< "\t" << getPercentageString(processedResults.second.perfectOverlapCombined, processedResults.second.total)
				<< "\t" << getPercentageString(processedResults.second.r1BeginsInR2Combined, processedResults.second.total)
				<< "\t" << getPercentageString(processedResults.second.r1BeginsInR2CombinedAboveCutOff, processedResults.second.total)
				<< "\t" << getPercentageString(processedResults.second.r1EndsInR2Combined, processedResults.second.total)
				<< "\t" << getPercentageString(processedResults.second.r1AllInR2Combined, processedResults.second.total)
				<< "\t" << getPercentageString(processedResults.second.r2AllInR1Combined, processedResults.second.total)
				<< std::endl;
	}

	VecStr lengthNeeded;
	for(const auto & t : ids.targets_){
		if(nullptr == t.second.lenCuts_){
			lengthNeeded.emplace_back(t.first);
		}
	}

	std::unordered_map<std::string, SeqIOOptions> tempOuts;
	std::unordered_map<std::string, std::vector<uint32_t>> readLengthsPerTarget;
	std::unordered_map<std::string, uint32_t> possibleContaminationCounts;
	std::unordered_map<std::string, uint32_t> failedPairProcessingOverlap;
	uint32_t failedPairProcessingOverlapTotal = 0;

	std::unordered_map<std::string, uint32_t> failedPairProcessingUnexpectedStatus;
	uint32_t failedPairProcessingUnexpectedStatusTotal = 0;

	for(const auto & extractedMid : primersInMids){
		for(const auto & extractedPrimer : extractedMid.second){
			std::string name = extractedPrimer + extractedMid.first;
			if(!ids.containsMids() && "" != pars.corePars_.sampleName){
				name = extractedPrimer + pars.corePars_.sampleName;
			}else if(!ids.containsMids() && "all" == extractedMid.first){
				name = extractedPrimer;
			}
			if(njh::in(PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_)){
				if(setUp.pars_.verbose_ && setUp.pars_.debug_){
					std::cout << extractedPrimer << " " << "NOOVERLAP" << std::endl;
				}
				PairedRead filteringSeq;
				if(nullptr == resultsPerMidTarPair[name].second.notCombinedOpts){
					if(setUp.pars_.verbose_){
						std::cout << "No reads for " << name << std::endl;
					}
					continue;
				}
				//add failed pair processing
				uint32_t failedPairProcessingUnexpectedStatusForName = resultsPerMidTarPair[name].second.total - resultsPerMidTarPair[name].second.overlapFail;
				failedPairProcessingUnexpectedStatus[name] = failedPairProcessingUnexpectedStatusForName;
				failedPairProcessingUnexpectedStatusTotal += failedPairProcessingUnexpectedStatusForName;

				//get seq options for expected pair status
				SeqIOOptions filterSeqOpts = *resultsPerMidTarPair[name].second.notCombinedOpts;
				SeqInput processedReader(filterSeqOpts);

				SeqIOOptions contaminationOutOpts;
				SeqIOOptions tempOutOpts;

				if (SeqIOOptions::inFormats::FASTQPAIREDGZ == setUp.pars_.ioOptions_.inFormat_) {
					contaminationOutOpts = SeqIOOptions::genPairedOutGz(njh::files::make_path(unfilteredByPairsProcessedDir, "possibleContamination_" + name));
					tempOutOpts = SeqIOOptions::genPairedOutGz(njh::files::make_path(unfilteredByPairsProcessedDir, "temp_" + name));
				}else{
					contaminationOutOpts = SeqIOOptions::genPairedOut(njh::files::make_path(unfilteredByPairsProcessedDir, "possibleContamination_" + name));
					tempOutOpts = SeqIOOptions::genPairedOut(njh::files::make_path(unfilteredByPairsProcessedDir, "temp_" + name));
				}
				SeqOutput contaminationWriter(contaminationOutOpts);
				SeqOutput tempWriter(tempOutOpts);
				processedReader.openIn();
				tempWriter.openOut();
				tempOuts[name] = SeqIOOptions::genPairedIn(tempWriter.getPrimaryOutFnp(),
						tempWriter.getSecondaryOutFnp());
				while(processedReader.readNextRead(filteringSeq)
						|| filteringSeq.seqBase_.seq_.size() <= pars.corePars_.primIdsPars.compKmerLen_
						|| filteringSeq.mateSeqBase_.seq_.size() <= pars.corePars_.primIdsPars.compKmerLen_){
					bool pass = false;
					if(pairProcessor.params_.r1Trim_ > 0 && pairProcessor.params_.r1Trim_ < len(filteringSeq.seqBase_)){
						readVecTrimmer::trimOffEndBases(filteringSeq.seqBase_, pairProcessor.params_.r1Trim_);
					}
					if(pairProcessor.params_.r2Trim_ > 0 && pairProcessor.params_.r2Trim_ < len(filteringSeq.mateSeqBase_)){
						readVecTrimmer::trimOffEndBases(filteringSeq.mateSeqBase_, pairProcessor.params_.r2Trim_);
					}

					if(ids.targets_.at(extractedPrimer).refs_.empty()){
						pass = true;
					}else{
						auto firstMateCopy = filteringSeq.seqBase_;
						auto secodnMateCopy = filteringSeq.mateSeqBase_;
						seqUtil::removeLowerCase(firstMateCopy.seq_,firstMateCopy.qual_);
						seqUtil::removeLowerCase(secodnMateCopy.seq_,secodnMateCopy.qual_);
						kmerInfo seqKInfo(firstMateCopy.seq_, pars.corePars_.primIdsPars.compKmerLen_, false);
						kmerInfo mateSeqInfo(secodnMateCopy.seq_, pars.corePars_.primIdsPars.compKmerLen_, true);

						for(const auto & refKinfo : ids.targets_.at(extractedPrimer).refKInfos_){
							if(refKinfo.compareKmers(seqKInfo).second >= pars.corePars_.primIdsPars.compKmerSimCutOff_ &&
							   refKinfo.compareKmersRevComp(mateSeqInfo).second >= pars.corePars_.primIdsPars.compKmerSimCutOff_){
								pass = true;
								break;
							}
						}
					}
					if(pass){
						tempWriter.write(filteringSeq);
					}else{
						++contamination;
						contaminationWriter.openWrite(filteringSeq);
						++possibleContaminationCounts[name];
					}
				}
			} else {
				//add failed pair processing
				uint32_t failedPairProcessingForNameOverlapFail = resultsPerMidTarPair[name].second.overlapFail;
				failedPairProcessingOverlap[name] = failedPairProcessingForNameOverlapFail;
				failedPairProcessingOverlapTotal += failedPairProcessingForNameOverlapFail;

				uint32_t failedPairProcessingForNameUnexpectedStatus = resultsPerMidTarPair[name].second.total - resultsPerMidTarPair[name].second.overlapFail;


				std::vector<SeqIOOptions> allFilterSeqOpts;
				if(njh::in(PairedReadProcessor::ReadPairOverLapStatus::R1BEGINSINR2, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_) ||
					 njh::in(PairedReadProcessor::ReadPairOverLapStatus::PERFECTOVERLAP, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_)){
					failedPairProcessingForNameUnexpectedStatus -= resultsPerMidTarPair[name].second.r1AllInR2Combined;
					if(resultsPerMidTarPair[name].second.r1AllInR2Combined > 0){
						allFilterSeqOpts.emplace_back(*resultsPerMidTarPair[name].second.r1AllInR2CombinedOpts);
					}
					failedPairProcessingForNameUnexpectedStatus -= resultsPerMidTarPair[name].second.r2AllInR1Combined;
					if(resultsPerMidTarPair[name].second.r2AllInR1Combined > 0){
						allFilterSeqOpts.emplace_back(*resultsPerMidTarPair[name].second.r2AllInR1CombinedOpts);
					}
				}

				if(njh::in(PairedReadProcessor::ReadPairOverLapStatus::R1BEGINSINR2, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_)){
					failedPairProcessingForNameUnexpectedStatus -= resultsPerMidTarPair[name].second.r1BeginsInR2Combined;
					if(resultsPerMidTarPair[name].second.r1BeginsInR2Combined > 0){
						allFilterSeqOpts.emplace_back(*resultsPerMidTarPair[name].second.r1BeginsInR2CombinedOpts);
					}
				}
				if(njh::in(PairedReadProcessor::ReadPairOverLapStatus::R1ENDSINR2, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_)){
					failedPairProcessingForNameUnexpectedStatus -= resultsPerMidTarPair[name].second.r1EndsInR2Combined;
					if(resultsPerMidTarPair[name].second.r1EndsInR2Combined > 0){
						allFilterSeqOpts.emplace_back(*resultsPerMidTarPair[name].second.r1EndsInR2CombinedOpts);
					}
				}
				if(njh::in(PairedReadProcessor::ReadPairOverLapStatus::PERFECTOVERLAP, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_) ){
					failedPairProcessingForNameUnexpectedStatus -= resultsPerMidTarPair[name].second.perfectOverlapCombined;
					if(resultsPerMidTarPair[name].second.perfectOverlapCombined > 0){
						allFilterSeqOpts.emplace_back(*resultsPerMidTarPair[name].second.perfectOverlapCombinedOpts);
					}
				}
				failedPairProcessingUnexpectedStatus[name] = failedPairProcessingForNameUnexpectedStatus;
				failedPairProcessingUnexpectedStatusTotal += failedPairProcessingForNameUnexpectedStatus;

				SeqIOOptions contaminationOutOpts;
				SeqIOOptions tempOutOpts;

				if (SeqIOOptions::inFormats::FASTQPAIREDGZ == setUp.pars_.ioOptions_.inFormat_) {
					contaminationOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(unfilteredByPairsProcessedDir, "possibleContamination_" + name + ".fastq.gz"));
					tempOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(unfilteredByPairsProcessedDir, "temp_" + name + ".fastq.gz"));
				}else{
					contaminationOutOpts = SeqIOOptions::genFastqOut(njh::files::make_path(unfilteredByPairsProcessedDir, "possibleContamination_" + name + ".fastq"));
					tempOutOpts = SeqIOOptions::genFastqOut(njh::files::make_path(unfilteredByPairsProcessedDir, "temp_" + name + ".fastq"));\
				}
				SeqOutput tempWriter(tempOutOpts);
				if(!allFilterSeqOpts.empty()){
					tempWriter.openOut();
					tempOuts[name] = SeqIOOptions::genFastqIn(tempWriter.getPrimaryOutFnp());
				}
				SeqOutput contaminationWriter(contaminationOutOpts);
				for(const auto & filterSeqOpts : allFilterSeqOpts){
					seqInfo filteringSeq;
					SeqInput processedReader(filterSeqOpts);
					processedReader.openIn();

					while(processedReader.readNextRead(filteringSeq)){
						bool pass = false;
						if(ids.targets_.at(extractedPrimer).refs_.empty() || filteringSeq.seq_.size() <= pars.corePars_.primIdsPars.compKmerLen_){
							pass = true;
						}else{
							kmerInfo seqKInfo(filteringSeq.seq_, pars.corePars_.primIdsPars.compKmerLen_, false);
							for(const auto & refKinfo : ids.targets_.at(extractedPrimer).refKInfos_){
								if(refKinfo.compareKmers(seqKInfo).second >= pars.corePars_.primIdsPars.compKmerSimCutOff_){
									pass = true;
									break;
								}
							}
						}
						if(pass){
							if(njh::in(extractedPrimer, lengthNeeded)){
								readLengthsPerTarget[extractedPrimer].emplace_back(len(filteringSeq));
							}
							tempWriter.write(filteringSeq);
						}else{
							++contamination;
							contaminationWriter.openWrite(filteringSeq);
							++possibleContaminationCounts[name];
						}
					}
				}
			}
		}
	}

	for(const auto & tar : lengthNeeded){
		//will only be in here if any reads pass
		if(njh::in(tar, readLengthsPerTarget)){
			auto medianlength = vectorMedianRef(njh::mapAt(readLengthsPerTarget,tar));
			njh::mapAt(ids.targets_, tar).addLenCutOff(medianlength - (medianlength * .20), medianlength + (medianlength * .20));
		}
	}
	//log read lengths used as cut offs
	OutOptions readLengthOpts(njh::files::make_path(setUp.pars_.directoryName_, "readLengthsUsed.tab.txt"));
	OutputStream readLengthOut(readLengthOpts);
	readLengthOut << "target\tminlen\tmaxlen" << std::endl;
	for(const auto & tar : ids.targets_){
		if(nullptr != tar.second.lenCuts_){
			readLengthOut << tar.first
					<< "\t" << tar.second.lenCuts_->minLenChecker_.minLen_
					<< "\t" << tar.second.lenCuts_->maxLenChecker_.maxLen_
					<< std::endl;
		}
	}

	std::unordered_map<std::string, uint32_t> badQual;
	std::unordered_map<std::string, uint32_t> badNs;
	std::unordered_map<std::string, uint32_t> badMaxLen;
	std::unordered_map<std::string, uint32_t> badMinLen;
	std::unordered_map<std::string, uint32_t> goodFinal;

	std::set<std::string> allNames;

	for(const auto & extractedMid : primersInMids){
		for(const auto & extractedPrimer : extractedMid.second){
			std::string name = extractedPrimer + extractedMid.first;
			if(!ids.containsMids() && "" != pars.corePars_.sampleName){
				name = extractedPrimer + pars.corePars_.sampleName;
			}else if(!ids.containsMids() && "all" == extractedMid.first){
				name = extractedPrimer;
			}
			if(njh::in(PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP, njh::mapAt(ids.targets_, extractedPrimer).overlapStatuses_)){
				if(setUp.pars_.verbose_&& setUp.pars_.debug_){
					std::cout << extractedPrimer << " " << "NOOVERLAP" << std::endl;
				}
				if(!njh::in(name, tempOuts)){
					continue;
				}
				allNames.emplace(name);
				PairedRead filteringSeq;
				SeqInput tempReader(tempOuts[name]);
				tempReader.openIn();
				SeqIOOptions finalSeqOut = SeqIOOptions::genPairedOut(njh::files::make_path(setUp.pars_.directoryName_, name));
				SeqIOOptions badSeqOut =  SeqIOOptions::genPairedOut(njh::files::make_path(badDir, name));

				if(SeqIOOptions::inFormats::FASTQPAIREDGZ == setUp.pars_.ioOptions_.inFormat_){
					finalSeqOut = SeqIOOptions::genPairedOutGz(njh::files::make_path(setUp.pars_.directoryName_, name));
					badSeqOut =  SeqIOOptions::genPairedOutGz(njh::files::make_path(badDir, name));
				}else{
					finalSeqOut = SeqIOOptions::genPairedOut(njh::files::make_path(setUp.pars_.directoryName_, name));
					badSeqOut =  SeqIOOptions::genPairedOut(njh::files::make_path(badDir, name));
				}
				SeqOutput finalWriter(finalSeqOut);
				SeqOutput badWriter(badSeqOut);
				while(tempReader.readNextRead(filteringSeq)){
					bool bad = false;
					if(!nChecker.checkRead(filteringSeq)){
						++badNs[name];
						bad = true;
					}else if(!qualChecker.checkRead(filteringSeq)){
						++badQual[name];
						bad = true;
					}else{
						++used;
						++goodFinal[name];
						finalWriter.openWrite(filteringSeq);
					}
					if(bad){
						++qualityFilters;
						if(pars.corePars_.keepFilteredOff){
							badWriter.openWrite(filteringSeq);
						}
					}
				}
			} else {
				if(!njh::in(name, tempOuts)){
					//no reads
					continue;
				}
				allNames.emplace(name);
				seqInfo filteringSeq;
				SeqInput tempReader(tempOuts[name]);
				tempReader.openIn();
				SeqIOOptions finalSeqOut;
				SeqIOOptions badSeqOUt;
				if(SeqIOOptions::inFormats::FASTQPAIREDGZ == setUp.pars_.ioOptions_.inFormat_){
					finalSeqOut = SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, name));
					badSeqOUt =   SeqIOOptions::genFastqOutGz(njh::files::make_path(badDir, name));
				}else{
					finalSeqOut = SeqIOOptions::genFastqOut(njh::files::make_path(setUp.pars_.directoryName_, name));
					badSeqOUt =   SeqIOOptions::genFastqOut(njh::files::make_path(badDir, name));
				}
				SeqOutput finalWriter(finalSeqOut);
				SeqOutput badWriter(badSeqOUt);
				while(tempReader.readNextRead(filteringSeq)){
					bool bad = false;
					if(!nChecker.checkRead(filteringSeq)){
						bad = true;
						++badNs[name];
					}else if(!qualChecker.checkRead(filteringSeq)){
						bad = true;
						++badQual[name];
					}else if(!ids.targets_.at(extractedPrimer).lenCuts_->minLenChecker_.checkRead(filteringSeq)){
						bad = true;
						++badMinLen[name];
					}else if(!ids.targets_.at(extractedPrimer).lenCuts_->maxLenChecker_.checkRead(filteringSeq)){
						bad = true;
						++badMaxLen[name];
					}else{
						++used;
						++goodFinal[name];
						finalWriter.openWrite(filteringSeq);
					}
					if(bad){
						++qualityFilters;
						if(pars.corePars_.keepFilteredOff){
							badWriter.openWrite(filteringSeq);
						}
					}
				}
			}
		}
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

	OutOptions extractionStatsOpts(njh::files::make_path(setUp.pars_.directoryName_, "extractionStats.tab.txt"));
	OutputStream extractionStatsOut(extractionStatsOpts);
	extractionStatsOut << "inputName\tTotal\tfailedBarcode\tfailedPrimers\tmismatchedPrimers\tmismatchedPrimersDimers\tfailedOverlapPairProcessing\tfailedStatusPairProcessing\tpossibleContamination\tfilteredOff\tused" << std::endl;
	extractionStatsOut << seqName
			<< "\t" << count
			<< "\t" << getPercentageString(readsNotMatchedToBarcode, count)
			<< "\t" << getPercentageString(unrecognizedPrimers, count)
			<< "\t" << getPercentageString(mismatchedPrimers, count)
			<< "\t" << getPercentageString(mismatchedPrimersDimers, count)
			<< "\t" << getPercentageString(failedPairProcessingOverlapTotal, count)
			<< "\t" << getPercentageString(failedPairProcessingUnexpectedStatusTotal, count)
			<< "\t" << getPercentageString(contamination, count)
			<< "\t" << getPercentageString(qualityFilters, count)
			<< "\t" << getPercentageString(used, count) << std::endl;

	OutOptions extractionProfileOpts(njh::files::make_path(setUp.pars_.directoryName_, "extractionProfile.tab.txt"));
	OutputStream extractionProfileOut(extractionProfileOpts);
	extractionProfileOut << "inputName\tname\ttotalMatching\tgood\tbad\tfailedMinLen\tfailedMaxLen\tfailedQuality\tfailedNs\tfailedPossibleContamination\tfailedOverlapPairProcessing\tfailedStatusPairProcessing";

	extractionProfileOut << std::endl;
	for(const auto & name : allNames){
		uint32_t totalBad = badQual[name] + badNs[name] + badMinLen[name] + badMaxLen[name] + possibleContaminationCounts[name];
		extractionProfileOut << seqName
				<< "\t" << name
				<< "\t" << matchingPrimerCounts[name]
				<< "\t" <<  getPercentageString(goodFinal[name], matchingPrimerCounts[name])
				<< "\t" <<  getPercentageString(totalBad, matchingPrimerCounts[name])
				<< "\t" <<  getPercentageString(badMinLen[name], totalBad)
				<< "\t" <<  getPercentageString(badMaxLen[name], totalBad)
				<< "\t" <<  getPercentageString(badQual[name], totalBad)
				<< "\t" <<  getPercentageString(badNs[name], totalBad)
				<< "\t" <<  getPercentageString(possibleContaminationCounts[name], totalBad)
				<< "\t" <<  getPercentageString(failedPairProcessingOverlap[name], matchingPrimerCounts[name])
				<< "\t" <<  getPercentageString(failedPairProcessingUnexpectedStatus[name], matchingPrimerCounts[name])
				<< std::endl;
	}

	auto writeOutUnrecCounts = [&seqName](const SeqIOOptions & opts, const bfs::path & outFilename){
		if(opts.inExists()){
			auto outTable = getSeqPortionCounts(opts, 0, 20);
			OutputStream countsOuts(outFilename);
			uint32_t count = 0;
			if(opts.isPairedIn()){
				countsOuts << "inputName\tr1Seq\tr2Seq\tcount\tfractionOfUnrecognizedReads" << std::endl;
			}else{
				countsOuts << "inputName\tseq\tcount\tfractionOfUnrecognizedReads" << std::endl;
			}
			while(count < 10 && count < outTable.content_.size()){
				countsOuts << seqName
						<< "\t" << njh::conToStr(outTable.content_[count], "\t") << std::endl;
				++count;
			}
		}
	};
	if(SeqIOOptions::inFormats::FASTQPAIREDGZ == setUp.pars_.ioOptions_.inFormat_){
		auto unRecBarR1Opts = SeqIOOptions::genFastqIn(njh::files::make_path(setUp.pars_.directoryName_, "filteredOff/bad/unrecognizedBarcode_NO_MATCHING_R1.fastq.gz"));
		auto unRecBarR2Opts = SeqIOOptions::genFastqIn(njh::files::make_path(setUp.pars_.directoryName_, "filteredOff/bad/unrecognizedBarcode_NO_MATCHING_R2.fastq.gz"));
		auto unRecBarPaired = SeqIOOptions::genPairedIn(
				njh::files::make_path(setUp.pars_.directoryName_, "filteredOff/bad/unrecognizedBarcode_NO_MATCHING_R1.fastq.gz"),
				njh::files::make_path(setUp.pars_.directoryName_, "filteredOff/bad/unrecognizedBarcode_NO_MATCHING_R2.fastq.gz"));

		writeOutUnrecCounts(unRecBarR1Opts, njh::files::make_path(setUp.pars_.directoryName_, "top_mostCommonR1Starts_for_unrecognizedBarcodes.tab.txt"));
		writeOutUnrecCounts(unRecBarR2Opts, njh::files::make_path(setUp.pars_.directoryName_, "top_mostCommonR2Starts_for_unrecognizedBarcodes.tab.txt"));
		writeOutUnrecCounts(unRecBarPaired, njh::files::make_path(setUp.pars_.directoryName_, "top_mostCommonR1AndR2Starts_for_unrecognizedBarcodes.tab.txt"));

	}else{
		auto unRecBarR1Opts = SeqIOOptions::genFastqIn(njh::files::make_path(setUp.pars_.directoryName_, "filteredOff/bad/unrecognizedBarcode_NO_MATCHING_R1.fastq"));
		auto unRecBarR2Opts = SeqIOOptions::genFastqIn(njh::files::make_path(setUp.pars_.directoryName_, "filteredOff/bad/unrecognizedBarcode_NO_MATCHING_R2.fastq"));
		auto unRecBarPaired = SeqIOOptions::genPairedIn(
				njh::files::make_path(setUp.pars_.directoryName_, "filteredOff/bad/unrecognizedBarcode_NO_MATCHING_R1.fastq"),
				njh::files::make_path(setUp.pars_.directoryName_, "filteredOff/bad/unrecognizedBarcode_NO_MATCHING_R2.fastq"));

		writeOutUnrecCounts(unRecBarR1Opts, njh::files::make_path(setUp.pars_.directoryName_, "top_mostCommonR1Starts_for_unrecognizedBarcodes.tab.txt"));
		writeOutUnrecCounts(unRecBarR2Opts, njh::files::make_path(setUp.pars_.directoryName_, "top_mostCommonR2Starts_for_unrecognizedBarcodes.tab.txt"));
		writeOutUnrecCounts(unRecBarPaired, njh::files::make_path(setUp.pars_.directoryName_, "top_mostCommonR1AndR2Starts_for_unrecognizedBarcodes.tab.txt"));

	}
	if(!pars.corePars_.keepUnfilteredReads){
		njh::files::rmDirForce(unfilteredReadsDir);
	}
	if(!pars.corePars_.keepFilteredOff){
		njh::files::rmDirForce(filteredOffDir);
	}
	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}

	return 0;
}



}  // namespace njh
