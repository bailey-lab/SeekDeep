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


int SeekDeepRunner::extractorPairedEnd(const bib::progutils::CmdArgs & inputCommands) {
	SeekDeepSetUp setUp(inputCommands);
	ExtractorPairedEndPars pars;

	setUp.setUpExtractorPairedEnd(pars);

	// run log
	setUp.startARunLog(setUp.pars_.directoryName_);
	// parameter file
	setUp.writeParametersFile(setUp.pars_.directoryName_ + "parametersUsed.txt", false,
			false);

	PrimersAndMids ids(pars.idFilename);
	if(0 == ids.getMids().size() && 0 == ids.getTargets().size()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", failed to read either targets or primers from " << pars.idFilename << "\n";
		throw std::runtime_error{ss.str()};
	}
	if(setUp.pars_.verbose_){
		if(ids.getMids().size()> 0){
			std::cout << "Found: " << ids.getMids().size() <<" MIDs to de-multiplex on" << std::endl;
		}
		if(ids.getTargets().size() > 0){
			std::cout << "Found: " << ids.getTargets().size() <<" target primer pairs to de-multiplex on" << std::endl;
		}
	}
	//init mids
	if(ids.containsMids()){
		ids.initMidDeterminator();
		ids.mDeterminator_->setAllowableMismatches(pars.barcodeErrors);
		ids.mDeterminator_->setMidEndsRevComp(pars.midEndsRevComp);
	}
	//init primers
	if(ids.containsTargets()){
		ids.initPrimerDeterminator();
	}
	//add in any length cuts if any
	if("" != pars.lenCutOffFilename_){
		ids.addLenCutOffs(pars.lenCutOffFilename_);
	}
	//add in ref sequences if any
	if("" != pars.comparisonSeqFnp_){
		ids.addRefSeqs(pars.comparisonSeqFnp_);
		ids.setRefSeqsKInfos(pars.compKmerLen, pars.compKmerSimCutOff);
	}

	//add in overlap status
	ids.addOverLapStatuses(pars.overlapStatusFnp_);

	//default checks for Ns and quality
	ReadCheckerOnSeqContaining nChecker("N", pars.numberOfNs, true);
	ReadCheckerQualCheck qualChecker(pars.qPars_.qualCheck_, pars.qPars_.qualCheckCutOff_, true);


	// make some directories for outputs
	bfs::path unfilteredReadsDir = bib::files::makeDir(setUp.pars_.directoryName_,
			bib::files::MkdirPar("unfilteredReads", false));
	bfs::path unfilteredByBarcodesDir = bib::files::makeDir(unfilteredReadsDir,
			bib::files::MkdirPar("byBarcodes", false));
	bfs::path unfilteredByPrimersDir = bib::files::makeDir(unfilteredReadsDir,
			bib::files::MkdirPar("byPrimers", false));
	bfs::path unfilteredByPairsProcessedDir = bib::files::makeDir(unfilteredReadsDir,
			bib::files::MkdirPar("pairsProcessed", false));
	bfs::path filteredOffDir = bib::files::makeDir(setUp.pars_.directoryName_,
			bib::files::MkdirPar("filteredOff", false));
	bfs::path overHansDir;
	if (pars.pairProcessorParams_.writeOverHangs_) {
		overHansDir = bib::files::makeDir(setUp.pars_.directoryName_,
				bib::files::MkdirPar("overhangs", false));
	}
	bfs::path badDir = bib::files::makeDir(filteredOffDir,
			bib::files::MkdirPar("bad", false));
	bfs::path unrecognizedPrimerDir = bib::files::makeDir(filteredOffDir,
			bib::files::MkdirPar("unrecognizedPrimer", false));
	bfs::path contaminationDir = "";


	PairedRead seq;

	// read in reads and remove lower case bases indicating tech low quality like
	// tags and such

	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto smallOpts = setUp.pars_.ioOptions_;
	smallOpts.out_.outFilename_ = bib::files::make_path(badDir, "smallFragments").string();
	SeqIO smallFragMentOut(smallOpts);
	smallFragMentOut.openOut();
	auto startsWtihBadQualOpts = setUp.pars_.ioOptions_;
	startsWtihBadQualOpts.out_.outFilename_ = bib::files::make_path(badDir , "startsWtihBadQual").string();
	SeqOutput startsWtihBadQualOut(startsWtihBadQualOpts);

	uint32_t smallFragmentCount = 0;
	uint64_t maxReadSize = 0;

	MultiSeqIO readerOuts;

	std::map<std::string, std::pair<uint32_t, uint32_t>> counts;
	std::unordered_map<std::string, uint32_t>  failBarCodeCounts;
	std::unordered_map<std::string, uint32_t>  failBarCodeCountsPossibleContamination;

	if (ids.containsMids() && pars.barcodeErrors > 0) {
		if (setUp.pars_.debug_) {
			std::cout << "Allowing " << pars.barcodeErrors << " errors in barcode"
					<< std::endl;
		}
	}
	if (ids.containsMids()) {
		for (const auto & mid : ids.mDeterminator_->mids_) {
			auto midOpts = setUp.pars_.ioOptions_;
			midOpts.out_.outFilename_ = bib::files::make_path(unfilteredByBarcodesDir, mid.first).string();
			if (setUp.pars_.debug_) {
				std::cout << "Inserting: " << mid.first << std::endl;
			}
			readerOuts.addReader(mid.first, midOpts);
		}
	} else {
		auto midOpts = setUp.pars_.ioOptions_;
		midOpts.out_.outFilename_ = bib::files::make_path(unfilteredByBarcodesDir, "all").string();
		if (setUp.pars_.debug_) {
			std::cout << "Inserting: " << "all" << std::endl;
		}
		readerOuts.addReader("all", midOpts);
	}

	if (ids.containsMids()) {
		auto failureCases = MidDeterminator::midPos::getFailureCaseNames();
		for(const auto & failureCase : failureCases){
			std::string unRecName = "unrecognizedBarcode_" + failureCase;
			auto midOpts = setUp.pars_.ioOptions_;
			midOpts.out_.outFilename_ = bib::files::make_path(badDir, unRecName).string();
			if (setUp.pars_.debug_){
				std::cout << "Inserting: " << unRecName << std::endl;
			}
			readerOuts.addReader(unRecName, midOpts);
		}
	}


	if(setUp.pars_.verbose_){
		std::cout << bib::bashCT::boldGreen("Extracting on MIDs") << std::endl;
	}
	if(setUp.pars_.verbose_){
		std::cout << "Reading in reads:" << std::endl;
	}
	uint32_t count = 0;
	uint32_t readsNotMatchedToBarcode = 0;
	uint32_t unrecognizedPrimers = 0;
	uint32_t contamination = 0;
	uint32_t qualityFilters = 0;
	uint32_t used = 0;


	while (reader.readNextRead(seq)) {
		++count;
		if (setUp.pars_.verbose_ && count % 50 == 0) {
			std::cout << "\r" << count ;
			std::cout.flush();
		}

		readVec::handelLowerCaseBases(seq, setUp.pars_.ioOptions_.lowerCaseBases_);

		if (len(seq) < pars.smallFragmentCutoff) {
			smallFragMentOut.write(seq);
			++smallFragmentCount;
			continue;
		}
		readVec::getMaxLength(seq.seqBase_, maxReadSize);
		readVec::getMaxLength(seq.mateSeqBase_, maxReadSize);

		std::pair<MidDeterminator::midPos, MidDeterminator::midPos> currentMid;
		if (ids.containsMids()) {
			currentMid = ids.mDeterminator_->fullDetermine(seq, pars.mDetPars);
		} else {
			currentMid = {MidDeterminator::midPos("all", 0, 0, 0), MidDeterminator::midPos("all", 0, 0, 0)};
		}
		if (seq.seqBase_.name_.find("_Comp") != std::string::npos) {
			++counts[currentMid.first.midName_].second;
		} else {
			++counts[currentMid.first.midName_].first;
		}
		if (!currentMid.first) {
			std::string unRecName = "unrecognizedBarcode_" + MidDeterminator::midPos::getFailureCaseName(currentMid.first.fCase_);
			++readsNotMatchedToBarcode;
			MidDeterminator::increaseFailedBarcodeCounts(currentMid.first, failBarCodeCounts);
			readerOuts.openWrite(unRecName, seq);
		}else{
			readerOuts.openWrite(currentMid.first.midName_, seq);
		}
	}

	OutOptions barcodeCountOpts(bib::files::make_path(setUp.pars_.directoryName_, "midCounts.tab.txt"));
	OutputStream barcodeCountOut(barcodeCountOpts);
	barcodeCountOut << "MID\tforwardCount\tforwardCountPerc\treverseCount\treverseCountPerc\ttotal" << std::endl;
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
			barcodeCountOut << countkey
					<< "\t" << counts[countkey].first
					<< "\t" << 100 * (counts[countkey].first/total)
					<< "\t" << counts[countkey].second
					<< "\t" << 100 * (counts[countkey].second/total)
					<< "\t" << total << std::endl;
		}else{
			barcodeCountOut << countkey
					<< "\t" << "0"
					<< "\t" << "0"
					<< "\t" << "0"
					<< "\t" << "0"
					<< "\t" << "0" << std::endl;
		}
	}
	if (setUp.pars_.verbose_) {
		std::cout << std::endl;
	}
	//close mid outs;
	readerOuts.closeOutAll();


	std::ofstream renameKeyFile;
	if (pars.rename) {
		openTextFile(renameKeyFile, setUp.pars_.directoryName_ + "renameKey.tab.txt", ".tab.txt", false, false);
		renameKeyFile << "originalName\tnewName\n";
	}

	// create aligner for primer identification
	auto scoreMatrix = substituteMatrix::createDegenScoreMatrixNoNInRef(
			setUp.pars_.generalMatch_, setUp.pars_.generalMismatch_);
	gapScoringParameters gapPars(setUp.pars_.gapInfo_);
	KmerMaps emptyMaps;
	bool countEndGaps = false;
	//to avoid allocating an extremely large aligner matrix;



	aligner alignObj = aligner(maxReadSize, gapPars, scoreMatrix, emptyMaps,
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
	auto barcodeFiles = bib::files::listAllFiles(unfilteredByBarcodesDir, false,
			VecStr { });
	ReadPairsOrganizer prOrg(expectedSamples);
	prOrg.processFiles(barcodeFiles);
	auto readsByPairs = prOrg.processReadPairs();
	auto barcodeReadsKeys = getVectorOfMapKeys(readsByPairs);
	bib::sort(barcodeReadsKeys);
	if(setUp.pars_.verbose_){
		std::cout << "The following mids were found with reads" << std::endl;
		printVector(barcodeReadsKeys);
	}

	std::unordered_map<std::string, std::set<std::string>> primersInMids;

	for (const auto & barcodeReadsKey : barcodeReadsKeys) {
		const auto & barcodeReadPairs = readsByPairs[barcodeReadsKey];
		std::string barcodeName = barcodeReadsKey;
		if (setUp.pars_.verbose_) {
			if (ids.containsMids()) {
				std::cout
						<< bib::bashCT::boldGreen("Filtering on barcode: " + barcodeName)
						<< std::endl;
			} else {
				std::cout << bib::bashCT::boldGreen("Filtering") << std::endl;
			}
		}

		//create outputs
		MultiSeqIO midReaderOuts;
		auto unrecogPrimerOutOpts = setUp.pars_.ioOptions_;
		unrecogPrimerOutOpts.out_.outFilename_ = bib::files::make_path(
				unrecognizedPrimerDir, barcodeName).string();
		midReaderOuts.addReader("unrecognized", unrecogPrimerOutOpts);
		for (const auto & primerName : getVectorOfMapKeys(
				ids.pDeterminator_->primers_)) {
			//determine full name
			std::string fullname = primerName;
			if (ids.containsMids()) {
				fullname += barcodeName;
			} else if (pars.sampleName != "") {
				fullname += pars.sampleName;
			}
			//bad out
			auto badDirOutOpts = setUp.pars_.ioOptions_;
			badDirOutOpts.out_.outFilename_ =
					bib::files::make_path(badDir, fullname).string();
			midReaderOuts.addReader(fullname + "bad", badDirOutOpts);
			//good out
			auto goodDirOutOpts = setUp.pars_.ioOptions_;

			goodDirOutOpts.out_.outFilename_ = bib::files::make_path(unfilteredByPrimersDir,  fullname);
			midReaderOuts.addReader(fullname + "good", goodDirOutOpts);
		}

		uint32_t barcodeCount = 1;
		bib::ProgressBar pbar(counts[barcodeName].first + counts[barcodeName].second);
		pbar.progColors_ = pbar.RdYlGn_;

		SeqIOOptions barcodePairsReaderOpts = SeqIOOptions::genPairedIn(
				barcodeReadPairs.first.front(), barcodeReadPairs.second.front());
		SeqInput barcodePairsReader(barcodePairsReaderOpts);
		barcodePairsReader.openIn();
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
			if (ids.containsMids()) {
				forwardPrimerName = ids.pDeterminator_->determineForwardPrimer(seq.seqBase_, pars.primerWithinStart_, alignObj,
						pars.primerErrors, !pars.primerToUpperCase);
				if ("unrecognized" ==  forwardPrimerName && pars.mDetPars.checkComplement_) {
					forwardPrimerName = ids.pDeterminator_->determineWithReversePrimer(seq.seqBase_, pars.primerWithinStart_, alignObj,
							pars.primerErrors, !pars.primerToUpperCase);
					if (seq.seqBase_.on_) {
						foundInReverse = true;
					}
				}
			} else {
				forwardPrimerName = ids.pDeterminator_->determineForwardPrimer(seq.seqBase_, pars.primerWithinStart_, alignObj,
						pars.primerErrors, !pars.primerToUpperCase);
				if ("unrecognized" ==  forwardPrimerName && pars.mDetPars.checkComplement_) {
					forwardPrimerName = ids.pDeterminator_->determineForwardPrimer(seq.mateSeqBase_, pars.primerWithinStart_,
							alignObj, pars.primerErrors, !pars.primerToUpperCase);
					if (seq.seqBase_.on_) {
						foundInReverse = true;
					}
				}
			}

			//reverse primer
			std::string reversePrimerName = "";
			if (!foundInReverse) {
				reversePrimerName = ids.pDeterminator_->determineWithReversePrimer(seq.mateSeqBase_,
						pars.primerWithinStart_,
						alignObj,
						pars.primerErrors,
						!pars.primerToUpperCase);
			} else {
				reversePrimerName = ids.pDeterminator_->determineWithReversePrimer(seq.seqBase_,
						pars.primerWithinStart_,
						alignObj,
						pars.primerErrors,
						!pars.primerToUpperCase);
			}

			std::string fullname = "";
			if(forwardPrimerName != reversePrimerName){
				fullname = forwardPrimerName + "-" + reversePrimerName;
			}else{
				fullname = forwardPrimerName;
			}
			if (ids.containsMids()) {
				fullname += barcodeName;
			}// else if (pars.sampleName != "") {
//				fullname += pars.sampleName;
//			}

			if("unrecognized" == forwardPrimerName ||
				 "unrecognized" == reversePrimerName){
				//check for unrecognized primers
				stats.increaseFailedForward(barcodeName, seq.seqBase_.name_);
				midReaderOuts.openWrite("unrecognized", seq);
				++allPrimerCounts[fullname];
				++unrecognizedPrimers;
			}else if(forwardPrimerName != reversePrimerName){
				//check for primer mismatch
				//stats.increasePrimerMismatch(barcodeName, seq.seqBase_.name_);
				stats.increaseCounts(barcodeName, seq.seqBase_.name_, ExtractionStator::extractCase::MISMATCHPRIMERS);
				auto badDirOutOpts = setUp.pars_.ioOptions_;
				badDirOutOpts.out_.outFilename_ =
						bib::files::make_path(badDir, fullname).string();
				if(!midReaderOuts.containsReader(fullname + "bad")){
					midReaderOuts.addReader(fullname + "bad", badDirOutOpts);
				}
				midReaderOuts.openWrite(fullname + "bad", seq);
				++allPrimerCounts[fullname];
				++unrecognizedPrimers;
			}else{
				primersInMids[barcodeName].emplace(forwardPrimerName);
				//primer match
				stats.increaseCounts(fullname, seq.seqBase_.name_,
						ExtractionStator::extractCase::GOOD);
				midReaderOuts.openWrite(fullname + "good", seq);
				++allPrimerCounts[fullname];
				++matchingPrimerCounts[fullname];
			}
		}
	}

	auto primerCountsKeys = getVectorOfMapKeys(allPrimerCounts);
	bib::sort(primerCountsKeys);
	OutOptions allPrimerCountsOpts(bib::files::make_path(setUp.pars_.directoryName_, "allPrimerCounts.tab.txt"));
	OutputStream allPrimerCountsOut(allPrimerCountsOpts);

	OutOptions goodCountsOpts(
			bib::files::make_path(setUp.pars_.directoryName_,
					"matchingPrimerCounts.tab.txt"));
	OutputStream goodCountsOut(goodCountsOpts);
	goodCountsOut << "name\tcount" << std::endl;
	for (const auto & good : allPrimerCounts) {
		allPrimerCountsOut << good.first << "\t" << good.second << std::endl;
		if (bib::in(good.first, matchingPrimerCounts)) {
			goodCountsOut << good.first << "\t" << good.second << std::endl;
		}
	}
	//now post processing
	PairedReadProcessor pairProcessor(pars.pairProcessorParams_);
	auto alnGapPars = gapScoringParameters(
			setUp.pars_.gapInfo_.gapOpen_,
			setUp.pars_.gapInfo_.gapExtend_,
			0,0,
			0,0);
	auto pairedProcessingScoring = substituteMatrix::createScoreMatrix(2, -2, true, true, true);
	aligner processingPairsAligner(maxReadSize, alnGapPars, pairedProcessingScoring, false);
	processingPairsAligner.qScorePars_.qualThresWindow_ = 0;
	std::unordered_map<std::string, std::vector<uint32_t>> lengthsPerStitchedTarget;
	std::unordered_map<std::string, PairedReadProcessor::ProcessedResults> resultsPerMidTarPair;
	for(const auto & extractedMid : primersInMids){
		for(const auto & extractedPrimer : extractedMid.second){
			std::string name = extractedPrimer + extractedMid.first;
			SeqIOOptions currentPairOpts;
			if(setUp.pars_.ioOptions_.inFormat_ == SeqIOOptions::inFormats::FASTQPAIREDGZ){
				currentPairOpts = SeqIOOptions::genPairedInGz(
						bib::files::make_path(unfilteredByPrimersDir, name + "_R1.fastq.gz"),
						bib::files::make_path(unfilteredByPrimersDir, name + "_R2.fastq.gz"));
			}else{
				currentPairOpts = SeqIOOptions::genPairedIn(
						bib::files::make_path(unfilteredByPrimersDir, name + "_R1.fastq"),
						bib::files::make_path(unfilteredByPrimersDir, name + "_R2.fastq"));
			}
			if(pars.noOverlapProcessForNoOverlapStatusTargets_ && PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP == bib::mapAt(ids.targets_, extractedPrimer).overlapStatus_){
				PairedReadProcessor::ProcessedResults currentProcessResults;
				currentProcessResults.notCombinedOpts = std::make_shared<SeqIOOptions>(currentPairOpts);
				resultsPerMidTarPair[name] = currentProcessResults;
				continue;
			}
			OutOptions currentOutOpts(bib::files::make_path(unfilteredByPairsProcessedDir, name));
			PairedReadProcessor::ProcessorOutWriters processWriter;
			processWriter.overhangsWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(bib::files::make_path(overHansDir, name + "_overhangs")));
			processWriter.perfectOverlapCombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(bib::files::make_path(unfilteredByPairsProcessedDir, name + "_perfectOverlap")));
			if(PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP == bib::mapAt(ids.targets_, extractedPrimer).overlapStatus_){
				processWriter.notCombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genPairedOut(bib::files::make_path(unfilteredByPairsProcessedDir, name)));
				processWriter.r1BeginsInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(bib::files::make_path(badDir, name + "_r1BeginsInR2.fastq")));
				processWriter.r1EndsInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(bib::files::make_path(badDir, name + "_r1EndsInR2.fastq")));
			}else if(PairedReadProcessor::ReadPairOverLapStatus::R1BEGINSINR2 == bib::mapAt(ids.targets_, extractedPrimer).overlapStatus_){
				processWriter.r1BeginsInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(bib::files::make_path(unfilteredByPairsProcessedDir, name + "")));
				processWriter.notCombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genPairedOut(bib::files::make_path(badDir, name + "_notCombined")));
				processWriter.r1EndsInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(bib::files::make_path(badDir, name + "_r1EndsInR2.fastq")));
			}else if(PairedReadProcessor::ReadPairOverLapStatus::R1ENDSINR2 == bib::mapAt(ids.targets_, extractedPrimer).overlapStatus_){
				processWriter.r1EndsInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(bib::files::make_path(unfilteredByPairsProcessedDir, name + "")));
				processWriter.notCombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genPairedOut(bib::files::make_path(badDir, name + "_notCombined")));
				processWriter.r1BeginsInR2CombinedWriter = std::make_unique<SeqOutput>(SeqIOOptions::genFastqOut(bib::files::make_path(badDir, name + "_r1BeginsInR2.fastq")));
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
			resultsPerMidTarPair[name] = currentProcessResults;
		}
	}

	OutputStream outPairProcessOut(bib::files::make_path(setUp.pars_.directoryName_, "processPairsCounts.tab.txt"));
	outPairProcessOut << "Name\ttotal\tnotCombined\tperfectOverlapCombined\tr1BeginsInR2Combined\tr1EndsInR2Combined" << std::endl;;
	for(const auto & n : resultsPerMidTarPair){
		outPairProcessOut << n.first
				<< "\t" << n.second.total
				<< "\t" << getPercentageString(n.second.overlapFail, n.second.total)
				<< "\t" << getPercentageString(n.second.perfectOverlapCombined, n.second.total)
				<< "\t" << getPercentageString(n.second.r1BeginsInR2Combined, n.second.total)
				<< "\t" << getPercentageString(n.second.r1EndsInR2Combined, n.second.total)
				<< std::endl;
	}

	bool checkingContamination = pars.comparisonSeqFnp_ != "";
	VecStr lengthNeeded;
	for(const auto & t : ids.targets_){
		if(nullptr == t.second.lenCuts_){
			lengthNeeded.emplace_back(t.first);
		}
	}
	std::unordered_map<std::string, SeqIOOptions> tempOuts;
	std::unordered_map<std::string, std::vector<uint32_t>> readLengthsPerTarget;
	std::unordered_map<std::string, uint32_t> possibleContaminationCounts;
	for(const auto & extractedMid : primersInMids){
		for(const auto & extractedPrimer : extractedMid.second){
			std::string name = extractedPrimer + extractedMid.first;
			if(PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP == bib::mapAt(ids.targets_, extractedPrimer).overlapStatus_){
				if(setUp.pars_.verbose_){
					std::cout << extractedPrimer << " " << "NOOVERLAP" << std::endl;
				}
				PairedRead filteringSeq;
				if(nullptr == resultsPerMidTarPair[name].notCombinedOpts){
					if(setUp.pars_.verbose_){
						std::cout << "No reads for " << name << std::endl;
					}
					continue;
				}
				SeqIOOptions filterSeqOpts = *resultsPerMidTarPair[name].notCombinedOpts;
				SeqInput processedReader(filterSeqOpts);
				auto contaminationOutOpts = SeqIOOptions::genPairedOut(bib::files::make_path(unfilteredByPairsProcessedDir, "possibleContamination_" + name));
				SeqOutput contaminationWriter(contaminationOutOpts);
				auto tempOutOpts = SeqIOOptions::genPairedOut(bib::files::make_path(unfilteredByPairsProcessedDir, "temp_" + name));

				SeqOutput tempWriter(tempOutOpts);
				processedReader.openIn();
				tempWriter.openOut();
				tempOuts[name] = SeqIOOptions::genPairedIn(tempWriter.getPrimaryOutFnp(),
						tempWriter.getSecondaryOutFnp());
				while(processedReader.readNextRead(filteringSeq)){
					if(checkingContamination){
						bool pass = false;
						if(ids.targets_.at(extractedPrimer).refs_.empty()){
							pass = true;
						}
						auto firstMateCopy = filteringSeq.seqBase_;
						auto secodnMateCopy = filteringSeq.mateSeqBase_;
						seqUtil::removeLowerCase(firstMateCopy.seq_,firstMateCopy.qual_);
						seqUtil::removeLowerCase(secodnMateCopy.seq_,secodnMateCopy.qual_);
						kmerInfo seqKInfo(firstMateCopy.seq_, pars.compKmerLen, false);
						kmerInfo mateSeqInfo(secodnMateCopy.seq_, pars.compKmerLen, true);

						for(const auto & refKinfo : ids.targets_.at(extractedPrimer).refKInfos_){
							if(refKinfo.compareKmers(seqKInfo).second >= pars.compKmerSimCutOff &&
							   refKinfo.compareKmersRevComp(mateSeqInfo).second >= pars.compKmerSimCutOff){
								pass = true;
								break;
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
				}
			}else if(PairedReadProcessor::ReadPairOverLapStatus::R1BEGINSINR2 == bib::mapAt(ids.targets_, extractedPrimer).overlapStatus_ ||
					PairedReadProcessor::ReadPairOverLapStatus::R1ENDSINR2 == bib::mapAt(ids.targets_, extractedPrimer).overlapStatus_){
				if(setUp.pars_.verbose_){
					std::cout << extractedPrimer << " " << "r1BeginsInR2" << std::endl;
				}
				seqInfo filteringSeq;
				SeqIOOptions filterSeqOpts;

				if(PairedReadProcessor::ReadPairOverLapStatus::R1BEGINSINR2 == bib::mapAt(ids.targets_, extractedPrimer).overlapStatus_){
					if(nullptr == resultsPerMidTarPair[name].r1BeginsInR2CombinedOpts){
						if(setUp.pars_.verbose_){
							std::cout << "No reads stitched for " << name << std::endl;
						}
						continue;
					}
					filterSeqOpts = *resultsPerMidTarPair[name].r1BeginsInR2CombinedOpts;
				}else{
					if(nullptr == resultsPerMidTarPair[name].r1EndsInR2CombinedOpts){
						if(setUp.pars_.verbose_){
							std::cout << "No reads stitched for " << name << std::endl;
						}
						continue;
					}
					filterSeqOpts = *resultsPerMidTarPair[name].r1EndsInR2CombinedOpts;
				}
				SeqInput processedReader(filterSeqOpts);
				auto contaminationOutOpts = SeqIOOptions::genFastqOut(bib::files::make_path(unfilteredByPairsProcessedDir, "possibleContamination_" + name + ".fastq"));
				SeqOutput contaminationWriter(contaminationOutOpts);
				auto tempOutOpts = SeqIOOptions::genFastqOut(bib::files::make_path(unfilteredByPairsProcessedDir, "temp_" + name + ".fastq"));

				SeqOutput tempWriter(tempOutOpts);
				processedReader.openIn();
				tempWriter.openOut();
				tempOuts[name] = SeqIOOptions::genFastqIn(tempWriter.getPrimaryOutFnp());
				while(processedReader.readNextRead(filteringSeq)){
					if(checkingContamination){
						bool pass = false;
						if(ids.targets_.at(extractedPrimer).refs_.empty()){
							pass = true;
						}
						kmerInfo seqKInfo(filteringSeq.seq_, pars.compKmerLen, false);
						for(const auto & refKinfo : ids.targets_.at(extractedPrimer).refKInfos_){
							if(refKinfo.compareKmers(seqKInfo).second >= pars.compKmerSimCutOff){
								pass = true;
								break;
							}
						}
						if(pass){
							if(bib::in(extractedPrimer, lengthNeeded)){
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
		auto medianlength = vectorMedianRef(readLengthsPerTarget.at(tar));
		ids.targets_.at(tar).addLenCutOff(medianlength - (medianlength * .10), medianlength + (medianlength * .10));
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

			if(PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP == bib::mapAt(ids.targets_, extractedPrimer).overlapStatus_){
				if(setUp.pars_.verbose_){
					std::cout << extractedPrimer << " " << "NOOVERLAP" << std::endl;
				}
				if(!bib::in(name, tempOuts)){
					continue;
				}
				allNames.emplace(name);
				PairedRead filteringSeq;
				SeqInput tempReader(tempOuts[name]);
				tempReader.openIn();
				auto finalSeqOut = SeqIOOptions::genPairedOut(bib::files::make_path(setUp.pars_.directoryName_, name));
				auto badSeqOut =  SeqIOOptions::genPairedOut(bib::files::make_path(badDir, name));
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
						badWriter.openWrite(filteringSeq);
					}
				}
			}else if(PairedReadProcessor::ReadPairOverLapStatus::R1BEGINSINR2 == bib::mapAt(ids.targets_, extractedPrimer).overlapStatus_ ||
					PairedReadProcessor::ReadPairOverLapStatus::R1ENDSINR2 == bib::mapAt(ids.targets_, extractedPrimer).overlapStatus_){
				if(!bib::in(name, tempOuts)){
					continue;
				}
				if(setUp.pars_.verbose_){
					std::cout << extractedPrimer << " " << "r1BeginsInR2" << std::endl;
				}
				allNames.emplace(name);
				seqInfo filteringSeq;
				SeqInput tempReader(tempOuts[name]);
				tempReader.openIn();
				auto finalSeqOut = SeqIOOptions::genFastqOut(bib::files::make_path(setUp.pars_.directoryName_, name));
				auto badSeqOUt =  SeqIOOptions::genFastqOut(bib::files::make_path(badDir, name));
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
						badWriter.openWrite(filteringSeq);
					}
				}
			}
		}
	}

	OutOptions extractionStatsOpts(bib::files::make_path(setUp.pars_.directoryName_, "extractionStats.tab.txt"));
	OutputStream extractionStatsOut(extractionStatsOpts);
	extractionStatsOut << "Total\tfailedBarcode\tfailedPrimers\tpossibleContamination\tfailedQualityLenFilters\tused" << std::endl;
	extractionStatsOut << count
			<< "\t" << getPercentageString(readsNotMatchedToBarcode, count)
			<< "\t" << getPercentageString(unrecognizedPrimers, count)
			<< "\t" << getPercentageString(contamination, count)
			<< "\t" << getPercentageString(qualityFilters, count)
			<< "\t" << getPercentageString(used, count) << std::endl;

	OutOptions extractionProfileOpts(bib::files::make_path(setUp.pars_.directoryName_, "extractionProfile.tab.txt"));
	OutputStream extractionProfileOut(extractionProfileOpts);
	extractionProfileOut << "name\ttotalMatching\tgood\tbad\tfailedMinLen\tfailedMaxLen\tfailedQuality\tfailedNs\tfailedPossibleContamination" << std::endl;
	for(const auto & name : allNames){
		uint32_t totalBad = badQual[name] + badNs[name] + badMinLen[name] + badMaxLen[name] + possibleContaminationCounts[name];
		extractionProfileOut << name
				<< "\t" << matchingPrimerCounts[name]
				<< "\t" <<  getPercentageString(goodFinal[name], matchingPrimerCounts[name])
				<< "\t" <<  getPercentageString(totalBad, matchingPrimerCounts[name])
				<< "\t" <<  getPercentageString(badMinLen[name], totalBad)
				<< "\t" <<  getPercentageString(badMaxLen[name], totalBad)
				<< "\t" <<  getPercentageString(badQual[name], totalBad)
				<< "\t" <<  getPercentageString(badNs[name], totalBad)
				<< "\t" <<  getPercentageString(possibleContaminationCounts[name], totalBad)
				<< std::endl;
	}

	{
		auto unRecBarOpts = SeqIOOptions::genFastqIn(bib::files::make_path(setUp.pars_.directoryName_, "filteredOff/bad/unrecognizedBarcode_NO_MATCHING_R1.fastq"));
		if(unRecBarOpts.inExists()){
			SeqIO reader(unRecBarOpts);
			reader.openIn();
			seqInfo read;
			strCounter counter;
			std::function<std::string(const seqInfo &, size_t, uint32_t)> getSubStr  = [](const seqInfo & read, size_t position, uint32_t size){
				std::string ret = "";
				if(position + size <= read.seq_.size()){
					ret = read.seq_.substr(position, size);
				}
				return ret;
			};
			while(reader.readNextRead(read)){
				counter.increaseCountByString(getSubStr(read, 0, 20));
			}
			counter.setFractions();
			table outTable(VecStr { "str", "count", "fraction" });
			for (const auto & str : counter.counts_) {
				if (str.second > 0 && str.first != "") {
					outTable.content_.emplace_back(
							toVecStr(str.first, str.second, counter.fractions_[str.first]));
				}
			}
			outTable.sortTable("count", true);
			OutputStream countsOuts(bib::files::make_path(setUp.pars_.directoryName_, "top_mostCommonR1Starts_for_unrecognizedBarcodes"));
			uint32_t count = 0;
			countsOuts << "seq\tcount\tpercentOfUnrecognizedReads" << std::endl;
			while(count < 10 && count < outTable.content_.size()){
				countsOuts << outTable.content_[count][0] << "\t" <<  outTable.content_[count][1] << "\t" << getPercentageString(bib::StrToNumConverter::stoToNum<uint32_t>(outTable.content_[count][1]), counter.getTotalCount()) << std::endl;
				++count;
			}
		}
	}
	{
		auto unRecBarOpts = SeqIOOptions::genFastqIn(bib::files::make_path(setUp.pars_.directoryName_, "filteredOff/bad/unrecognizedBarcode_NO_MATCHING_R2.fastq"));
		if(unRecBarOpts.inExists()){
			SeqIO reader(unRecBarOpts);
			reader.openIn();
			seqInfo read;
			strCounter counter;
			std::function<std::string(const seqInfo &, size_t, uint32_t)> getSubStr  = [](const seqInfo & read, size_t position, uint32_t size){
				std::string ret = "";
				if(position + size <= read.seq_.size()){
					ret = read.seq_.substr(position, size);
				}
				return ret;
			};
			while(reader.readNextRead(read)){
				counter.increaseCountByString(getSubStr(read, 0, 20));
			}
			counter.setFractions();
			table outTable(VecStr { "str", "count", "fraction" });
			for (const auto & str : counter.counts_) {
				if (str.second > 0 && str.first != "") {
					outTable.content_.emplace_back(
							toVecStr(str.first, str.second, counter.fractions_[str.first]));
				}
			}
			outTable.sortTable("count", true);
			OutputStream countsOuts(bib::files::make_path(setUp.pars_.directoryName_, "top_mostCommonR2Starts_for_unrecognizedBarcodes"));
			uint32_t count = 0;
			countsOuts << "seq\tcount\tpercentOfUnrecognizedReads" << std::endl;
			while(count < 10 && count < outTable.content_.size()){
				countsOuts << outTable.content_[count][0] << "\t" <<  outTable.content_[count][1] << "\t" << getPercentageString(bib::StrToNumConverter::stoToNum<uint32_t>(outTable.content_[count][1]), counter.getTotalCount()) << std::endl;
				++count;
			}
		}
	}

//		auto barcodeName = bfs::basename(f.first.string());
//		if ((counts[barcodeName].first + counts[barcodeName].second) == 0
//				&& ids.containsMids()) {
//			//no reads extracted for barcode so skip filtering step
//			continue;
//		}
//		if (setUp.pars_.verbose_) {
//			if (ids.containsMids()) {
//				std::cout
//						<< bib::bashCT::boldGreen("Filtering on barcode: " + barcodeName)
//						<< std::endl;
//			} else {
//				std::cout << bib::bashCT::boldGreen("Filtering") << std::endl;
//			}
//		}
//
//		auto barcodeOpts = setUp.pars_.ioOptions_;
//		barcodeOpts.firstName_ = f.first.string();
//		barcodeOpts.inFormat_ = SeqIOOptions::getInFormat(
//				bib::files::getExtension(f.first.string()));
//		SeqIO barcodeIn(barcodeOpts);
//		barcodeIn.openIn();
//
//		//create outputs
//		MultiSeqIO midReaderOuts;
//		auto unrecogPrimerOutOpts = setUp.pars_.ioOptions_;
//		unrecogPrimerOutOpts.out_.outFilename_ = bib::files::make_path(unrecognizedPrimerDir, barcodeName).string();
//		midReaderOuts.addReader("unrecognized", unrecogPrimerOutOpts);
//
//		for (const auto & primerName : getVectorOfMapKeys(ids.pDeterminator_->primers_)) {
//			std::string fullname = primerName;
//			if (ids.containsMids()) {
//				fullname += barcodeName;
//			} else if (pars.sampleName != "") {
//				fullname += pars.sampleName;
//			}
//			//bad out
//			auto badDirOutOpts = setUp.pars_.ioOptions_;
//			badDirOutOpts.out_.outFilename_ = bib::files::make_path(badDir, fullname).string();
//			midReaderOuts.addReader(fullname + "bad", badDirOutOpts);
//			//good out
//			auto goodDirOutOpts = setUp.pars_.ioOptions_;
//			goodDirOutOpts.out_.outFilename_ = setUp.pars_.directoryName_ + fullname;
//			midReaderOuts.addReader(fullname + "good", goodDirOutOpts);
//		}
//
//		uint32_t barcodeCount = 1;
//		bib::ProgressBar pbar(
//				counts[barcodeName].first + counts[barcodeName].second);
//		pbar.progColors_ = pbar.RdYlGn_;
//		std::string readFlows = "";
//		std::ifstream inFlowFile;
//
//
//
//		while (barcodeIn.readNextRead(seq)) {
//			//std::cout << barcodeCount << std::endl;
//			std::getline(inFlowFile, readFlows);
//			if(setUp.pars_.verbose_){
//				pbar.outputProgAdd(std::cout, 1, true);
//			}
//			++barcodeCount;
//			//filter on primers
//			//forward
//			std::string primerName = "";
//			bool foundInReverse = false;
//			if (pars.noForwardPrimer) {
//				primerName = primerTable.content_.front()[primerTable.getColPos(
//						"geneName")];
//			} else {
//				if (pars.multiplex) {
//
//					primerName = pDetermine.determineForwardPrimer(seq, 0, alignObj,
//							pars.primerErrors, !pars.primerToUpperCase);
//
//					if (primerName == "unrecognized" && pars.mDetPars.checkComplement_) {
//						primerName = pDetermine.determineWithReversePrimer(seq, 0,
//								alignObj, pars.primerErrors, !pars.primerToUpperCase);
//						if (seq.seqBase_.on_) {
//							foundInReverse = true;
//						}
//					}
//
//				} else {
//					uint32_t start = pars.mDetPars.variableStop_;
//					//std::cout << "Determining primer" << std::endl;
//
//					primerName = pDetermine.determineForwardPrimer(seq, start, alignObj,
//							pars.primerErrors, !pars.primerToUpperCase);
//					if (primerName == "unrecognized" && pars.mDetPars.checkComplement_) {
//						primerName = pDetermine.determineWithReversePrimer(seq, start,
//								alignObj, pars.primerErrors, !pars.primerToUpperCase);
//						if (seq.seqBase_.on_) {
//							foundInReverse = true;
//						}
//					}
//				}
//
//				if (!seq.seqBase_.on_) {
//					stats.increaseFailedForward(barcodeName, seq.seqBase_.name_);
//					midReaderOuts.openWrite("unrecognized", seq);
//					continue;
//				}
//			}
//			std::string fullname = primerName;
//			if (pars.multiplex) {
//				fullname += barcodeName;
//			} else if (pars.sampleName != "") {
//				fullname += pars.sampleName;
//			}
//
//
//			//min len
//			multipleLenCutOffsR1.at(primerName).minLenChecker_.checkRead(seq.seqBase_);
//			multipleLenCutOffsR2.at(primerName).minLenChecker_.checkRead(seq.mateSeqBase_);
//
//			if (!seq.seqBase_.on_) {
//				stats.increaseCounts(fullname, seq.seqBase_.name_,
//						ExtractionStator::extractCase::MINLENBAD);
//				midReaderOuts.openWrite(fullname + "bad", seq);
//				continue;
//			}
//
//			//reverse
//			if (!pars.noReversePrimer) {
//				if (foundInReverse) {
//					pDetermine.checkForForwardPrimerInRev(seq, primerName, alignObj,
//							pars.rPrimerErrors, !pars.reversePrimerToUpperCase,
//							pars.mDetPars.variableStop_, false);
//				} else {
//					pDetermine.checkForReversePrimer(seq, primerName, alignObj,
//							pars.rPrimerErrors, !pars.reversePrimerToUpperCase,
//							pars.mDetPars.variableStop_, false);
//				}
//				if (!seq.seqBase_.on_) {
//					stats.increaseCounts(fullname, seq.seqBase_.name_,
//							ExtractionStator::extractCase::BADREVERSE);
//					seq.seqBase_.name_.append("_badReverse");
//					midReaderOuts.openWrite(fullname + "bad", seq);
//					continue;
//				}
//			}
//
//			//min len again becuase the reverse primer search trims to the reverse primer so it could be short again
//			multipleLenCutOffsR1.at(primerName).minLenChecker_.checkRead(seq.seqBase_);
//			multipleLenCutOffsR2.at(primerName).minLenChecker_.checkRead(seq.mateSeqBase_);
//
//
//			if (!seq.seqBase_.on_) {
//				stats.increaseCounts(fullname, seq.seqBase_.name_,
//						ExtractionStator::extractCase::MINLENBAD);
//				midReaderOuts.openWrite(fullname + "bad", seq);
//				continue;
//			}
//
//			//if found in the reverse direction need to re-orient now
//			/**@todo handle finding in rev comp */
////			if (foundInReverse) {
////				seq.seqBase_.reverseComplementRead(true, true);
////			}
//
//
//
//
//			if (seq.seqBase_.on_) {
//				stats.increaseCounts(fullname, seq.seqBase_.name_,
//						ExtractionStator::extractCase::GOOD);
//				if (pars.rename) {
//					std::string oldName = bib::replaceString(seq.seqBase_.name_, "_Comp", "");
//					seq.seqBase_.name_ = fullname + "."
//							+ leftPadNumStr(goodCounts[fullname],
//									counts[barcodeName].first + counts[barcodeName].second);
//					if (bib::containsSubString(oldName, "_Comp")) {
//						seq.seqBase_.name_.append("_Comp");
//					}
//					renameKeyFile << oldName << "\t" << seq.seqBase_.name_ << "\n";
//				}
//				midReaderOuts.openWrite(fullname + "good", seq);
//				++goodCounts[fullname];
//			}
//		}
//		if(setUp.pars_.verbose_){
//			std::cout << std::endl;
//		}
//	}
//
//	std::ofstream profileLog;
//	openTextFile(profileLog, setUp.pars_.directoryName_ + "extractionProfile.tab.txt",
//			".txt", false, false);
//	profileLog
//			<< "name\ttotalReadsExtracted\t\n";
//	stats.outStatsPerName(profileLog, "\t");
//	std::ofstream extractionStatsFile;
//	openTextFile(extractionStatsFile,
//			setUp.pars_.directoryName_ + "extractionStats.tab.txt", ".txt",
//			setUp.pars_.ioOptions_.out_);
//	extractionStatsFile
//			<< "TotalReads\tReadsNotMatchedBarcodes\tReadsNotMatchedBarcodesPosContamination\tSmallFragments(len<"
//			<< pars.smallFragmentCutoff
//			<< ")\tfailedForwardPrimer\tfailedQualityFiltering\tused";
//	extractionStatsFile << "\tcontamination";
//	extractionStatsFile << std::endl;
//	stats.outTotalStats(extractionStatsFile, "\t");
//	std::ofstream failedForwardFile;
//	openTextFile(failedForwardFile, setUp.pars_.directoryName_ + "failedForward.tab.txt",
//			".txt", false, false);
//	failedForwardFile << "MidName\ttotalFailed\tfailedInFor\tfailedInRev"
//			<< std::endl;
//	stats.outFailedForwardStats(failedForwardFile, "\t");
//	if(pars.multiplex){
//		std::ofstream failedBarcodeFile;
//		openTextFile(failedBarcodeFile, setUp.pars_.directoryName_ + "failedBarcode.tab.txt",
//				".txt", false, false);
//		failedBarcodeFile << "Reason\tcount"
//				<< std::endl;
//		auto countKeys = getVectorOfMapKeys(failBarCodeCounts);
//		bib::sort(countKeys);
//		for(const auto & countKey : countKeys){
//			failedBarcodeFile << countKey << "\t" << getPercentageString(failBarCodeCounts.at(countKey), readsNotMatchedToBarcode)<< std::endl;
//		}
//	}
//
//	if (!setUp.pars_.debug_) {
//		bib::files::rmDirForce(unfilteredReadsDir);
//	}
//
//	if (setUp.pars_.writingOutAlnInfo_) {
//		setUp.rLog_ << "Number of alignments done" << "\n";
//		alignObj.alnHolder_.write(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);
//	}

	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}
	return 0;
}



}  // namespace bib
