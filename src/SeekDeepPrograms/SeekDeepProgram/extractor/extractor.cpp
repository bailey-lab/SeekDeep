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





int SeekDeepRunner::extractor(const njh::progutils::CmdArgs & inputCommands) {
	SeekDeepSetUp setUp(inputCommands);
	extractorPars pars;
	setUp.setUpExtractor(pars);

	uint32_t readsNotMatchedToBarcode = 0;
	uint32_t readsNotMatchedToBarcodePossContam = 0;

	// run log
	setUp.startARunLog(setUp.pars_.directoryName_);
	// parameter file
	setUp.writeParametersFile(setUp.pars_.directoryName_ + "parametersUsed.txt", false, false);

	// create Primers and MIDs
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
	// init
	ids.initAllAddLenCutsRefs(pars.corePars_.primIdsPars);
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
	// make some directories for outputs
	bfs::path unfilteredReadsDir = njh::files::makeDir(
			setUp.pars_.directoryName_,
			njh::files::MkdirPar("unfilteredReads", false));
	bfs::path unfilteredByBarcodesDir = njh::files::makeDir(unfilteredReadsDir,
			njh::files::MkdirPar("byBarcodes", false));
	bfs::path unfilteredByBarcodesFlowDir = njh::files::makeDir(
			unfilteredReadsDir, njh::files::MkdirPar("flowsByBarcodes", false));
	bfs::path unfilteredByPrimersDir = njh::files::makeDir(unfilteredReadsDir,
			njh::files::MkdirPar("byPrimers", false));
	bfs::path filteredOffDir = njh::files::makeDir(setUp.pars_.directoryName_,
			njh::files::MkdirPar("filteredOff", false));
	bfs::path badDir = njh::files::makeDir(filteredOffDir,
			njh::files::MkdirPar("bad", false));
	bfs::path unrecognizedPrimerDir = njh::files::makeDir(filteredOffDir,
			njh::files::MkdirPar("unrecognizedPrimer", false));
	bfs::path contaminationDir = "";
	if ("" != pars.corePars_.primIdsPars.comparisonSeqFnp_) {
		contaminationDir = njh::files::makeDir(filteredOffDir,
				njh::files::MkdirPar("contamination", false));
	}

	std::shared_ptr<readObject> seq = std::make_shared<readObject>();

	// read in reads and remove lower case bases indicating tech low quality like
	// tags and such
	if(setUp.pars_.verbose_){
		std::cout << "Reading in reads:" << std::endl;
	}
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto smallOpts = setUp.pars_.ioOptions_;
	smallOpts.out_.outFilename_ = njh::files::make_path(badDir,"smallFragments").string();
	SeqIO smallFragMentOut(smallOpts);
	smallFragMentOut.openOut();
	auto startsWtihBadQualOpts = setUp.pars_.ioOptions_;
	startsWtihBadQualOpts.out_.outFilename_ = njh::files::make_path(badDir,"startsWtihBadQual").string();
	SeqOutput startsWtihBadQualOut(startsWtihBadQualOpts);

	uint32_t smallFragmentCount = 0;
	// uint32_t startsWithBadQualCount = 0;
	uint64_t maxReadSize = 0;
	uint32_t count = 0;
	MultiSeqIO readerOuts;

	std::map<std::string, std::pair<uint32_t, uint32_t>> counts;
	std::unordered_map<std::string, uint32_t>  failBarCodeCounts;
	std::unordered_map<std::string, uint32_t>  failBarCodeCountsPossibleContamination;

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
			if(ids.screeningForPossibleContamination()){
				std::string unRecNamePosCon = "possible_contamination_unrecognizedBarcode_" + failureCase;
				auto midOpts = setUp.pars_.ioOptions_;
				midOpts.out_.outFilename_ = njh::files::make_path(contaminationDir, unRecNamePosCon).string();
				if (setUp.pars_.debug_){
					std::cout << "Inserting: " << unRecNamePosCon << std::endl;
				}
				readerOuts.addReader(unRecNamePosCon, midOpts);
			}
		}
	}
	ReadCheckerOnSeqContaining nChecker("N", pars.corePars_.numberOfNs, true);


	if(setUp.pars_.verbose_){
		std::cout << njh::bashCT::boldGreen("Extracting on MIDs") << std::endl;
	}

	std::vector<size_t> readLens;


	while (reader.readNextRead(seq)) {
		++count;
		if (setUp.pars_.verbose_ && count % 50 == 0) {
			std::cout << "\r" << count ;
			std::cout.flush();
		}
		readVec::handelLowerCaseBases(seq, setUp.pars_.ioOptions_.lowerCaseBases_);

		//possibly trim reads at low quality
		if(pars.trimAtQual){
			readVecTrimmer::trimAtFirstQualScore(seq->seqBase_, pars.trimAtQualCutOff);
			if(0 == len(*seq)){
				startsWtihBadQualOut.openWrite(seq);
				// ++startsWithBadQualCount;
				continue;
			}
		}

		if (len(*seq) < pars.corePars_.smallFragmentCutoff) {
			smallFragMentOut.write(seq);
			++smallFragmentCount;
			continue;
		}
		readVec::getMaxLength(seq, maxReadSize);


		if (ids.containsMids()) {
			auto searchRes = ids.mDeterminator_->searchRead(seq->seqBase_);
			auto processRes = ids.mDeterminator_->processSearchRead(seq->seqBase_, searchRes);
			if(MidDeterminator::ProcessedRes::PROCESSED_CASE::MATCH == processRes.case_){
				if (processRes.rcomplement_) {
					++counts[processRes.midName_].second;
				} else {
					++counts[processRes.midName_].first;
				}
				readLens.emplace_back(len(*seq));
				readerOuts.openWrite(processRes.midName_, seq);
			}else{
				std::string unRecName = "unrecognizedBarcode_" + MidDeterminator::ProcessedRes::getProcessedCaseName(processRes.case_);
				bool possibleContaimination = false;
				if(ids.screeningForPossibleContamination()){
					//this will check the read against all targets and their reverse complement so it will be a conservative estimate
					//of whether or not this is contamination, if the read is still on by the end then that it means it's not
					//considered possible contamination, could mark a lot seqs as contamination if not all seqs have comparison seqs
					kmerInfo seqKInfo(seq->seqBase_.seq_, pars.corePars_.primIdsPars.compKmerLen_, false);
					seq->seqBase_.on_ = false;
					for(const auto & tar : ids.targets_){
						for(const auto & refInfo : tar.second.refKInfos_){
							if(refInfo.compareKmers(seqKInfo).second >= pars.corePars_.primIdsPars.compKmerSimCutOff_){
								seq->seqBase_.on_ = true;
								break;
							}
						}
					}
					if(!seq->seqBase_.on_){
						possibleContaimination = true;
					}
				}
				if(possibleContaimination){
					unRecName = "possible_contamination_" + unRecName;
					++readsNotMatchedToBarcodePossContam;
					++failBarCodeCountsPossibleContamination[MidDeterminator::ProcessedRes::getProcessedCaseName(processRes.case_)];
				}else{
					++readsNotMatchedToBarcode;
					++failBarCodeCounts[MidDeterminator::ProcessedRes::getProcessedCaseName(processRes.case_)];
				}
				readerOuts.openWrite(unRecName, seq);
			}
		} else {
			++counts["all"].first;
			readLens.emplace_back(len(*seq));
			readerOuts.openWrite("all", seq);
		}
	}
	if (setUp.pars_.verbose_) {
		std::cout << std::endl;
	}
	//close mid outs;
	readerOuts.closeOutAll();

	//if no length was supplied, calculate a min and max length off of the median read length
	auto readLenMedian = vectorMedianRef(readLens);
	auto lenStep = readLenMedian * .20;
	if(std::numeric_limits<uint32_t>::max() == pars.minLen){
		if(lenStep > readLenMedian){
			pars.minLen = 0;
		}else{
			pars.minLen = ::round(readLenMedian - lenStep);
		}
	}
	if(std::numeric_limits<uint32_t>::max() == pars.maxLength){
		pars.maxLength = ::round(readLenMedian + lenStep);
	}
	ids.addDefaultLengthCutOffs(pars.minLen, pars.maxLength);

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



	// set up quality filtering
	std::unique_ptr<ReadChecker> qualChecker;
	if (pars.corePars_.qPars_.checkingQFrac_) {
		qualChecker = std::make_unique<ReadCheckerQualCheck>(pars.corePars_.qPars_.qualCheck_,
				pars.corePars_.qPars_.qualCheckCutOff_, true);
	} else {
		if (pars.qualWindowTrim) {
			qualChecker = std::make_unique<ReadCheckerOnQualityWindowTrim>(
					pars.corePars_.qPars_.qualityWindowLength_,
					pars.corePars_.qPars_.qualityWindowStep_,
					pars.corePars_.qPars_.qualityWindowThres_,
					pars.minLen, true);
		} else {
			qualChecker = std::make_unique<ReadCheckerOnQualityWindow>(
					pars.corePars_.qPars_.qualityWindowLength_,
					pars.corePars_.qPars_.qualityWindowStep_,
					pars.corePars_.qPars_.qualityWindowThres_, true);
		}
	}

	if (setUp.pars_.debug_) {
		table midCounts { VecStr { "MidName", "For", "Rev" } };
		for (const auto & mCount : counts) {
			midCounts.content_.emplace_back(
					toVecStr(mCount.first, mCount.second.first, mCount.second.second));
		}
		midCounts.outPutContentOrganized(std::cout);
	}

	std::ofstream renameKeyFile;
	if (pars.corePars_.rename) {
		openTextFile(renameKeyFile, setUp.pars_.directoryName_ + "renameKey.tab.txt",
				".tab.txt", false, false);
		renameKeyFile << "originalName\tnewName\n";
	}

	auto barcodeFiles = njh::files::listAllFiles(unfilteredByBarcodesDir, false, VecStr { });
	// creating aligner
	// create aligner for primer identification
	auto scoreMatrix = substituteMatrix::createDegenScoreMatrixNoNInRef(
			setUp.pars_.generalMatch_, setUp.pars_.generalMismatch_);
//	setUp.pars_.gapInfo_.gapLeftRefOpen_ = 5;
//	setUp.pars_.gapInfo_.gapLeftRefExtend_ = 1;
//	setUp.pars_.gapInfo_.gapRightRefOpen_ = 5;
//	setUp.pars_.gapInfo_.gapRightRefExtend_ = 1;
  setUp.pars_.gapInfo_.gapLeftQueryOpen_ = 0;
  setUp.pars_.gapInfo_.gapLeftQueryExtend_ = 0;

  setUp.pars_.gapInfo_.gapRightQueryOpen_ = 0;
  setUp.pars_.gapInfo_.gapRightQueryExtend_ = 0;

  setUp.pars_.gapInfo_.gapLeftRefOpen_ = 0;
  setUp.pars_.gapInfo_.gapLeftRefExtend_ = 0;

  setUp.pars_.gapInfo_.gapRightRefOpen_ = 0;
  setUp.pars_.gapInfo_.gapRightRefExtend_ = 0;

	gapScoringParameters gapPars(setUp.pars_.gapInfo_);
	KmerMaps emptyMaps;
	bool countEndGaps = false;
	//to avoid allocating an extremely large aligner matrix;
	if(maxReadSize > 1000){
		auto maxPrimerSize = ids.pDeterminator_->getMaxPrimerSize();
		if(setUp.pars_.debug_){
			std::cout << njh::bashCT::boldBlack("maxPrimerSize: ") << maxPrimerSize << std::endl;
		}
		maxReadSize =  maxPrimerSize * 4 + pars.corePars_.pDetPars.primerWithin_;
	}

	aligner alignObj(maxReadSize, gapPars, scoreMatrix, emptyMaps, setUp.pars_.qScorePars_, countEndGaps, false);
	alignObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
	bfs::path smallDir = "";
	if (pars.filterOffSmallReadCounts) {
		smallDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("smallReadCounts", false));
	}
	ExtractionStator stats(count, readsNotMatchedToBarcode,
			readsNotMatchedToBarcodePossContam, smallFragmentCount);
	std::map<std::string, uint32_t> goodCounts;

	for (const auto & barcodeFile : barcodeFiles) {
		std::string barcodeName = bfs::basename(barcodeFile.first.string());
		std::string extension =   njh::files::getExtension(barcodeFile.first.string());
		if(setUp.pars_.ioOptions_.isInGz()){
			barcodeName = bfs::basename(bfs::basename(barcodeFile.first.string()));
			extension = njh::files::getExtension(bfs::path(barcodeFile.first).replace_extension("")) + "." + njh::files::getExtension(bfs::path(barcodeFile.first));
			//std::cout << "extension: " << extension << std::endl;
		}

		if ((counts[barcodeName].first + counts[barcodeName].second) == 0) {
			//no reads extracted for barcode so skip filtering step
			continue;
		}
		if (setUp.pars_.verbose_) {
			if(ids.containsMids()){
				std::cout << "Starting Filtering On MID: " << barcodeName << std::endl;
			}else{
				std::cout << "Starting Filtering" << std::endl;
			}
		}
		if (pars.filterOffSmallReadCounts && (counts[barcodeName].first + counts[barcodeName].second) <= pars.smallExtractReadCount) {
			auto barcodeOpts = setUp.pars_.ioOptions_;
			barcodeOpts.firstName_ = barcodeFile.first.string();
			barcodeOpts.inFormat_ = SeqIOOptions::getInFormat(extension);
			barcodeOpts.out_.outFilename_ = njh::files::make_path(smallDir,  barcodeName).string();
			SeqIO barcodeIn(barcodeOpts);
			barcodeIn.openIn();
			readObject read;
			while (barcodeIn.readNextRead(read)) {
				barcodeIn.openWrite(read);
			}
			continue;
		}
		if (setUp.pars_.verbose_) {
			if (ids.containsMids()) {
				std::cout
						<< njh::bashCT::boldGreen("Filtering on barcode: " + barcodeName)
						<< std::endl;
			} else {
				std::cout << njh::bashCT::boldGreen("Filtering") << std::endl;
			}
		}

		auto barcodeOpts = setUp.pars_.ioOptions_;
		barcodeOpts.firstName_ = barcodeFile.first.string();
		barcodeOpts.inFormat_ = SeqIOOptions::getInFormat(extension);
		SeqIO barcodeIn(barcodeOpts);
		barcodeIn.openIn();

		//create outputs
		MultiSeqIO midReaderOuts;
		auto unrecogPrimerOutOpts = setUp.pars_.ioOptions_;
		unrecogPrimerOutOpts.out_.outFilename_ = njh::files::make_path(unrecognizedPrimerDir
				,barcodeName).string();
		midReaderOuts.addReader("unrecognized", unrecogPrimerOutOpts);

		for (const auto & primerName : getVectorOfMapKeys(ids.pDeterminator_->primers_)) {
			std::string fullname = primerName;
			if (ids.containsMids()) {
				fullname += barcodeName;
			} else if (pars.corePars_.sampleName != "") {
				fullname += pars.corePars_.sampleName;
			}
			//bad out
			auto badDirOutOpts = setUp.pars_.ioOptions_;
			badDirOutOpts.out_.outFilename_ = njh::files::make_path( badDir, fullname).string();
			midReaderOuts.addReader(fullname + "bad", badDirOutOpts);
			//good out
			auto goodDirOutOpts = setUp.pars_.ioOptions_;
			goodDirOutOpts.out_.outFilename_ = setUp.pars_.directoryName_ + fullname;
			midReaderOuts.addReader(fullname + "good", goodDirOutOpts);
			//contamination out
			if (ids.screeningForPossibleContamination()) {
				auto contamOutOpts = setUp.pars_.ioOptions_;
				contamOutOpts.out_.outFilename_ = njh::files::make_path(contaminationDir, fullname).string();
				midReaderOuts.addReader(fullname + "contamination", contamOutOpts);
			}
		}

		// uint32_t barcodeCount = 1;
		njh::ProgressBar pbar(
				counts[barcodeName].first + counts[barcodeName].second);
		pbar.progColors_ = pbar.RdYlGn_;

		while (barcodeIn.readNextRead(seq)) {
			if(setUp.pars_.verbose_){
				pbar.outputProgAdd(std::cout, 1, true);
			}
			// ++barcodeCount;
			//filter on primers
			//front primer determination
			std::string frontPrimerName = "unrecognized";
			std::string backPrimerName = "unrecognized";
			bool foundInReverse = false;
			std::string fullname = "";
			std::string targetName = "";
			if (pars.corePars_.noPrimers_) {
				frontPrimerName = ids.pDeterminator_->primers_.begin()->first;
				backPrimerName = ids.pDeterminator_->primers_.begin()->first;
				fullname = frontPrimerName;
				targetName = frontPrimerName;
				if (ids.containsMids()) {
					fullname += barcodeName;
				} else if (pars.corePars_.sampleName != "") {
					fullname += pars.corePars_.sampleName;
				}
			} else {
				bool primerCheckComplement = pars.corePars_.pDetPars.checkComplement_;
				//turn on auto determination for dual barcoded system where the barcodes are the same since direction cannot be determined by the MID barcode
				if("all" != barcodeName && ids.containsMids() && ids.mDeterminator_->mids_.at(barcodeName).forSameAsRev_){
					primerCheckComplement = true;
				}
				//front end primer
				frontPrimerName = ids.pDeterminator_->determineForwardPrimer(seq, pars.corePars_.pDetPars, alignObj);
				if (frontPrimerName == "unrecognized" && primerCheckComplement) {
					frontPrimerName = ids.pDeterminator_->determineWithReversePrimer(seq, pars.corePars_.pDetPars, alignObj);
					if (seq->seqBase_.on_) {
						foundInReverse = true;
					}
				}
				if ("unrecognized" == frontPrimerName) {
					stats.increaseFailedForward(barcodeName, seq->seqBase_.name_);
					midReaderOuts.openWrite("unrecognized", seq);
					continue;
				}

				//trim to max length
				if(pars.trimToMaxLength){
					readVecTrimmer::trimToMaxLength(seq, ids.targets_.at(frontPrimerName).lenCuts_->maxLenChecker_.maxLen_);
				}

				if(!pars.corePars_.noReversePrimer_){
					//back end primer
					seq->seqBase_.reverseComplementRead(true, true);
					if(foundInReverse){
						backPrimerName = ids.pDeterminator_->determineForwardPrimer(seq, pars.corePars_.backEndpDetPars, alignObj);
					}else{
						backPrimerName = ids.pDeterminator_->determineWithReversePrimer(seq, pars.corePars_.backEndpDetPars, alignObj);
						//if wasn't found in reverse, reverse back
						seq->seqBase_.reverseComplementRead(true, true);
					}
				}else{
					backPrimerName = frontPrimerName;
				}


				std::string primerName = "";
				if (frontPrimerName == backPrimerName) {
					primerName = frontPrimerName;
				} else {
					primerName = frontPrimerName + "-" + backPrimerName;
				}
				targetName = frontPrimerName;
				fullname = frontPrimerName;
				if (ids.containsMids()) {
					fullname += barcodeName;
				} else if (pars.corePars_.sampleName != "") {
					fullname += pars.corePars_.sampleName;
				}

				if (!seq->seqBase_.on_ || frontPrimerName != backPrimerName) {
					stats.increaseCounts(fullname, seq->seqBase_.name_,
							ExtractionStator::extractCase::BADREVERSE);
					if("unrecognized" == backPrimerName){
						seq->seqBase_.name_.append("_badReverse");
					}else{
						seq->seqBase_.name_.append("[backPrimer=" + backPrimerName + "]");
					}
					midReaderOuts.openWrite(fullname + "bad", seq);
					continue;
				}
			}

			//look for possible contamination
			if (!njh::mapAt(ids.targets_, targetName).refKInfos_.empty() ) {
				bool contamination = true;
				kmerInfo seqKInfo(seq->seqBase_.seq_, pars.corePars_.primIdsPars.compKmerLen_, false);
				for(const auto & refInfo : ids.targets_.at(targetName).refKInfos_){
					if(refInfo.compareKmers(seqKInfo).second >= pars.corePars_.primIdsPars.compKmerSimCutOff_){
						contamination = false;
						break;
					}
				}
				if(contamination){
					seq->seqBase_.on_ = false;
				}
				if (!seq->seqBase_.on_) {
					stats.increaseCounts(fullname, seq->seqBase_.name_,
							ExtractionStator::extractCase::CONTAMINATION);
					midReaderOuts.openWrite(fullname + "contamination", seq);
					continue;
				}
			}

			//min len
			ids.targets_.at(targetName).lenCuts_->minLenChecker_.checkRead(seq->seqBase_);

			if (!seq->seqBase_.on_) {
				stats.increaseCounts(fullname, seq->seqBase_.name_,
						ExtractionStator::extractCase::MINLENBAD);
				midReaderOuts.openWrite(fullname + "bad", seq);
				continue;
			}


			//contains n
			nChecker.checkRead(seq->seqBase_);
			if (!seq->seqBase_.on_) {
				stats.increaseCounts(fullname, seq->seqBase_.name_,
						ExtractionStator::extractCase::CONTAINSNS);
				midReaderOuts.openWrite(fullname + "bad", seq);
				continue;
			}

			//max len
			ids.targets_.at(targetName).lenCuts_->maxLenChecker_.checkRead(seq->seqBase_);
			if (!seq->seqBase_.on_) {
				stats.increaseCounts(fullname, seq->seqBase_.name_,
						ExtractionStator::extractCase::MAXLENBAD);
				midReaderOuts.openWrite(fullname + "bad", seq);
				continue;
			}
			//quality
			qualChecker->checkRead(seq->seqBase_);

			if (!seq->seqBase_.on_) {
				stats.increaseCounts(fullname, seq->seqBase_.name_,
						ExtractionStator::extractCase::QUALITYFAILED);
				midReaderOuts.openWrite(fullname + "bad", seq);
				continue;
			}

			if (seq->seqBase_.on_) {
				stats.increaseCounts(fullname, seq->seqBase_.name_,
						ExtractionStator::extractCase::GOOD);
				if (pars.corePars_.rename) {
					std::string oldName = njh::replaceString(seq->seqBase_.name_, "_Comp", "");
					seq->seqBase_.name_ = fullname + "."
							+ leftPadNumStr(goodCounts[fullname],
									counts[barcodeName].first + counts[barcodeName].second);
					if (njh::containsSubString(oldName, "_Comp")) {
						seq->seqBase_.name_.append("_Comp");
					}
					renameKeyFile << oldName << "\t" << seq->seqBase_.name_ << "\n";
				}
				midReaderOuts.openWrite(fullname + "good", seq);
				++goodCounts[fullname];
			}
		}
		if(setUp.pars_.verbose_){
			std::cout << std::endl;
		}
	}

	std::ofstream profileLog;
	openTextFile(profileLog, setUp.pars_.directoryName_ + "extractionProfile.tab.txt",
			".txt", false, false);
	profileLog
			<< "name\ttotalReadsExtracted\tgoodReadsExtracted\tforGood\trevGood\ttotalBadReads\t"
					"badBackEndPrimer\tcontainsNs\tfailedMinLen" << "\tfailedMaxLen";

	if (pars.corePars_.qPars_.checkingQFrac_) {
		profileLog << "\tq" + estd::to_string(pars.corePars_.qPars_.qualCheck_) + "<"
				<< pars.corePars_.qPars_.qualCheckCutOff_;
	} else {
		profileLog << "\tbadQualityWindow";
	}
	profileLog << "\tcontamination\n";
	stats.outStatsPerName(profileLog, "\t");
	std::ofstream extractionStatsFile;
	openTextFile(extractionStatsFile,
			setUp.pars_.directoryName_ + "extractionStats.tab.txt", ".txt",
			setUp.pars_.ioOptions_.out_);
	extractionStatsFile
			<< "TotalReads\tReadsNotMatchedBarcodes\tReadsNotMatchedBarcodesPosContamination\tSmallFragments(len<"
			<< pars.corePars_.smallFragmentCutoff
			<< ")\tfailedFrontPrimer\tfailedQualityFiltering\tused";
	extractionStatsFile << "\tcontamination";
	extractionStatsFile << std::endl;
	stats.outTotalStats(extractionStatsFile, "\t");
	std::ofstream failedForwardFile;
	openTextFile(failedForwardFile, setUp.pars_.directoryName_ + "failedForward.tab.txt",
			".txt", false, false);
	failedForwardFile << "MidName\ttotalFailed\tfailedInFor\tfailedInRev"
			<< std::endl;
	stats.outFailedForwardStats(failedForwardFile, "\t");
	if(ids.containsMids()){
		std::ofstream failedBarcodeFile;
		openTextFile(failedBarcodeFile, setUp.pars_.directoryName_ + "failedBarcode.tab.txt",
				".txt", false, false);
		failedBarcodeFile << "Reason\tcount"
				<< std::endl;
		auto countKeys = getVectorOfMapKeys(failBarCodeCounts);
		njh::sort(countKeys);
		for(const auto & countKey : countKeys){
			failedBarcodeFile << countKey << "\t" << getPercentageString(failBarCodeCounts.at(countKey), readsNotMatchedToBarcode)<< std::endl;
		}
		if(ids.screeningForPossibleContamination() && ! failBarCodeCountsPossibleContamination.empty()){
			std::ofstream failedBarcodePosContFile;
			openTextFile(failedBarcodePosContFile, setUp.pars_.directoryName_ + "failedBarcodePossibleContamination.tab.txt",
					".txt", false, false);
			failedBarcodePosContFile << "Reason\tcount"
					<< std::endl;
			auto countKeys = getVectorOfMapKeys(failBarCodeCountsPossibleContamination);
			njh::sort(countKeys);
			for(const auto & countKey : countKeys){
				failedBarcodePosContFile << countKey << "\t" << getPercentageString(failBarCodeCountsPossibleContamination.at(countKey), readsNotMatchedToBarcodePossContam)<< std::endl;
			}
		}
	}

	if(!pars.corePars_.keepUnfilteredReads){
		njh::files::rmDirForce(unfilteredReadsDir);
	}
	if (setUp.pars_.writingOutAlnInfo_) {
		setUp.rLog_ << "Number of alignments done" << "\n";
		alignObj.alnHolder_.write(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);
	}

	if(setUp.pars_.verbose_){
		setUp.logRunTime(std::cout);
	}
	return 0;
}



}  // namespace njh
