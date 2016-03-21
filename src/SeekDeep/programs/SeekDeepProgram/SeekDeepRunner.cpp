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
#include <bibcpp.h>
#include "SeekDeep/programs/SeekDeepUtils.h"

namespace bibseq {

SeekDeepRunner::SeekDeepRunner() :
		bib::progutils::oneRing( { addRing<SeekDeepUtilsRunner>() },
				{
						addFunc("extractor", extractor, false),
						addFunc("processClusters", processClusters,false),
						addFunc("qluster", qluster, false),
						addFunc("clusterDown",qluster, true),
						addFunc("makeSampleDirectories", makeSampleDirectories, false)
				}, "SeekDeep") {
}

int SeekDeepRunner::extractor(const bib::progutils::CmdArgs & inputCommands) {
	SeekDeepSetUp setUp(inputCommands);
	extractorPars pars;
	setUp.setUpExtractor(pars);

	uint32_t readsNotMatchedToBarcode = 0;
	if (pars.checkingQCheck) {
		std::cout << "Quality Check: " << pars.qualCheck << "\n";
		std::cout << "Quality Check Cut Off: " << pars.qualCheckCutOff << "\n";
		std::cout << "Q" << pars.qualCheck << ">" << pars.qualCheckCutOff << "\n";
	} else {
		std::cout << "Quality Window Length: " << pars.qualityWindowLength << "\n";
		std::cout << "Quality Window Step: " << pars.qualityWindowStep << "\n";
		std::cout << "Quality Window Threshold: " << pars.qualityWindowThres
				<< "\n";
	}
	// run log
	setUp.startARunLog(setUp.pars_.directoryName_);
	// parameter file
	setUp.writeParametersFile(setUp.pars_.directoryName_ + "parametersUsed.txt", false,
			false);
	readObject compareObject = readObject(seqInfo("Compare", pars.compareSeq));
	table primerTable = seqUtil::readPrimers(pars.idFilename, pars.idFileDelim,
			false);
	int midSize = 10;
	table mids = seqUtil::readBarcodes(pars.idFilename, pars.idFileDelim,
			midSize);
	std::unique_ptr<MidDeterminator> determinator;
	if (pars.multiplex) {
		determinator = std::make_unique<MidDeterminator>(mids);
	}
	PrimerDeterminator pDetermine(primerTable);

	std::string tcagLower = "tcag";
	std::string tcagUpper = "TCAG";
	auto checkForTcag = [&tcagLower, &tcagUpper](const std::string & seq)->bool {
		return std::equal(tcagLower.begin(), tcagLower.end(), seq.begin()) ||
		std::equal(tcagUpper.begin(), tcagUpper.end(), seq.begin());
	};
	// make some directories for outputs
	std::string unfilteredReadsDir = bib::files::makeDir(setUp.pars_.directoryName_,
			"unfilteredReads", false);
	std::string unfilteredByBarcodesDir = bib::files::makeDir(unfilteredReadsDir,
			"byBarcodes", false);
	std::string unfilteredByBarcodesFlowDir = bib::files::makeDir(unfilteredReadsDir,
			"flowsByBarcodes", false);
	std::string unfilteredByPrimersDir = bib::files::makeDir(unfilteredReadsDir,
			"byPrimers", false);
	std::string filteredOffDir = bib::files::makeDir(setUp.pars_.directoryName_,
			"filteredOff", false);
	std::string badDir = bib::files::makeDir(filteredOffDir, "bad", false);
	std::string unrecognizedPrimerDir = bib::files::makeDir(filteredOffDir,
			"unrecognizedPrimer", false);
	std::string contaminationDir = "";
	if (pars.screenForPossibleContamination) {
		contaminationDir = bib::files::makeDir(filteredOffDir, "contamination", false);
	}
	readObject read;
	// read in reads and remove lower case bases indicating tech low quality like
	// tags and such
	std::cout << "Reading in reads:" << std::endl;
	SeqIO reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto smallOpts = setUp.pars_.ioOptions_;
	smallOpts.out_.outFilename_ = badDir + "smallFragments";
	SeqIO smallFragMentOut(smallOpts);
	smallFragMentOut.openOut();
	uint32_t smallFragmentCount = 0;
	uint64_t maxReadSize = 0;
	uint32_t count = 0;
	MultiSeqIO readerOuts;
	std::map<std::string, std::pair<uint32_t, uint32_t>> counts;
	if (pars.multiplex && pars.barcodeErrors > 0) {
		if (setUp.pars_.debug_) {
			std::cout << "Allowing " << pars.barcodeErrors << " errors in barcode"
					<< std::endl;
		}
	}
	if (pars.multiplex) {
		for (const auto & mid : determinator->mids_) {
			auto midOpts = setUp.pars_.ioOptions_;
			midOpts.out_.outFilename_ = unfilteredByBarcodesDir + mid.first;
			if (setUp.pars_.debug_) {
				std::cout << "Inserting: " << mid.first << std::endl;
			}
			readerOuts.addReader(mid.first, midOpts);
			if (pars.mothurExtract || pars.pyroExtract) {
				auto flowExtractOutOpts = midOpts;
				flowExtractOutOpts.out_.outFilename_ = unfilteredByBarcodesFlowDir + mid.first + "_temp";
				flowExtractOutOpts.outFormat_ = SeqIOOptions::outFormats::FLOW;
				flowExtractOutOpts.out_.outExtention_ = ".dat";
				readerOuts.addReader(mid.first + "flow", flowExtractOutOpts);
			}
		}
	} else {
		auto midOpts = setUp.pars_.ioOptions_;
		midOpts.out_.outFilename_ = unfilteredByBarcodesDir + "all";
		if (setUp.pars_.debug_) {
			std::cout << "Inserting: " << "all" << std::endl;
		}
		readerOuts.addReader("all", midOpts);
		if (pars.mothurExtract || pars.pyroExtract) {
			auto flowExtractOutOpts = midOpts;
			flowExtractOutOpts.out_.outFilename_ = unfilteredByBarcodesFlowDir + "all" + "_temp";
			flowExtractOutOpts.outFormat_ = SeqIOOptions::outFormats::FLOW;
			flowExtractOutOpts.out_.outExtention_ = ".dat";
			readerOuts.addReader(std::string("all") + "flow", flowExtractOutOpts);
		}
	}

	if (pars.multiplex) {
		auto midOpts = setUp.pars_.ioOptions_;
		midOpts.out_.outFilename_ = badDir + "unrecognizedBarcode";
		if (setUp.pars_.debug_){
			std::cout << "Inserting: " << "unrecognized" << std::endl;
		}
		readerOuts.addReader("unrecognized", midOpts);
	}

	std::cout << bib::bashCT::boldGreen("Extracting on MIDs") << std::endl;
	while (reader.readNextRead(read)) {
		++count;
		if (count % 50 == 0) {
			std::cout << count << "\r";
			std::cout.flush();
		}

		if (pars.HMP) {
			subStrToUpper(read.seqBase_.seq_, 4, pars.primerLen);
		}
		if (pars.trimTcag) {
			if (checkForTcag(read.seqBase_.seq_)) {
				read.trimFront(4);
			}
		}
		readVec::handelLowerCaseBases(read, setUp.pars_.ioOptions_.lowerCaseBases_);

		if (len(read) < pars.smallFragmentCutoff) {
			smallFragMentOut.write(read);
			++smallFragmentCount;
			continue;
		}
		readVec::getMaxLength(read, maxReadSize);

		midPos currentMid;
		if (pars.multiplex) {
			if (pars.barcodeErrors > 0) {
				currentMid = determinator->fullDetermine(read, pars.variableStart,
						pars.variableStop, pars.checkComplement, pars.barcodesBothEnds,
						pars.barcodeErrors);
			} else {
				currentMid = determinator->fullDetermine(read, pars.variableStart,
						pars.variableStop, pars.checkComplement, pars.barcodesBothEnds);
			}
		} else {
			currentMid = midPos("all", 0, 0);
		}

		if (!currentMid) {
			++readsNotMatchedToBarcode;
		}

		if (read.seqBase_.name_.find("_Comp") != std::string::npos) {
			++counts[currentMid.midName_].second;
		} else {
			++counts[currentMid.midName_].first;
		}
		/**@todo need to reorient the reads here before outputing if that's needed*/
		readerOuts.openWrite(currentMid.midName_, read);
		if (currentMid.midName_ != "unrecognized" && (pars.mothurExtract || pars.pyroExtract)) {
			readerOuts.openWriteFlow(currentMid.midName_ + "flow", *reader.in_.lastSffRead_);
		}
	}
	std::cout << std::endl;
	//close mid outs;
	readerOuts.closeOutAll();

	if (setUp.pars_.debug_) {
		table midCounts { VecStr { "MidName", "For", "Rev" } };
		for (const auto & mCount : counts) {
			midCounts.content_.emplace_back(
					toVecStr(mCount.first, mCount.second.first, mCount.second.second));
		}
		midCounts.outPutContentOrganized(std::cout);
	}

	std::ofstream renameKeyFile;
	if (pars.rename) {
		openTextFile(renameKeyFile, setUp.pars_.directoryName_ + "renameKey.tab.txt",
				".tab.txt", false, false);
		renameKeyFile << "originalName\tnewName\n";
	}

	if (pars.noForwardPrimer) {
		if (primerTable.content_.size() > 1) {
			std::cerr
					<< "Error, if noForwardPrimer is turned on can only supply one gene name, curently have: "
					<< primerTable.content_.size() << std::endl;
			std::cerr << bib::conToStr(primerTable.getColumn("geneName"), ",")
					<< std::endl;
			exit(1);
		}
	}
	auto barcodeFiles = bib::files::listAllFiles(unfilteredByBarcodesDir, false,
			VecStr { });
	kmerInfo compareInfo;
	kmerInfo compareInfoRev;
	std::map<std::string, kmerInfo> compareInfos;
	std::map<std::string, kmerInfo> compareInfosRev;
	struct lenCutOffs {
		uint32_t minLen_;
		uint32_t maxLen_;
	};
	std::map<std::string, lenCutOffs> multipleLenCutOffs;
	if (pars.multipleTargets && pars.multipleLenCutOffFilename != "") {
		table lenCutTab = table(pars.multipleLenCutOffFilename, "whitespace", true);
		bib::for_each(lenCutTab.columnNames_,
				[](std::string & str) {stringToLower(str);});
		lenCutTab.setColNamePositions();
		if (!bib::in(std::string("target"), lenCutTab.columnNames_)
				|| !bib::in(std::string("minlen"), lenCutTab.columnNames_)
				|| !bib::in(std::string("maxlen"), lenCutTab.columnNames_)) {
			std::stringstream ss;
			ss << "need to have columns " << "target,minlen, and maxlen"
					<< " when reading in a table for multiple cut off lengths"
					<< std::endl;
			ss << "only have " << vectorToString(lenCutTab.columnNames_, ",")
					<< std::endl;
			throw std::runtime_error { bib::bashCT::boldRed(ss.str()) };
		}
		for (const auto & row : lenCutTab.content_) {
			multipleLenCutOffs[row[lenCutTab.getColPos("target")]] = {bib::lexical_cast<uint32_t>(row[lenCutTab.getColPos("minlen")]),bib::lexical_cast<uint32_t>(row[lenCutTab.getColPos("maxlen")])};
		}

	}
	if (pars.screenForPossibleContamination) {
		readVec::getMaxLength(compareObject, maxReadSize);
		compareInfo = kmerInfo(compareObject.seqBase_.seq_, pars.contaminationKLen,
				false);
		auto rev = compareObject;
		rev.seqBase_.reverseComplementRead(true, true);
		compareInfoRev = kmerInfo(rev.seqBase_.seq_, pars.contaminationKLen, false);
	}
	if (pars.multipleTargets && pars.screenForPossibleContamination) {
		std::cout << pars.compareSeqFilename << std::endl;
		SeqIOOptions conTamOpts;
		conTamOpts.inFormat_ = SeqIOOptions::getInFormat(bib::files::getExtension(pars.compareSeqFilename));
		conTamOpts.firstName_ = pars.compareSeqFilename;
		SeqInput readerCon(conTamOpts);
		readerCon.openIn();
		auto reads = readerCon.readAllReads<readObject>();
		for (const auto & read : reads) {
			readVec::getMaxLength(read, maxReadSize);
			compareInfos[read.seqBase_.name_] = kmerInfo(read.seqBase_.seq_,
					pars.contaminationKLen, false);
			auto rev = read;
			rev.seqBase_.reverseComplementRead(true, true);
			compareInfosRev[read.seqBase_.name_] = kmerInfo(rev.seqBase_.seq_,
					pars.contaminationKLen, false);
		}
	}

	// create aligner for primer identification
	auto scoreMatrix = substituteMatrix::createDegenScoreMatrixNoNInRef(
			setUp.pars_.generalMatch_, setUp.pars_.generalMismatch_);
	gapScoringParameters gapPars(setUp.pars_.gapInfo_);
	KmerMaps emptyMaps;
	bool countEndGaps = true;
	aligner alignObj = aligner(maxReadSize, gapPars, scoreMatrix, emptyMaps,
			setUp.pars_.qScorePars_, countEndGaps);
	alignObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);

	std::string smallDir = "";
	if (pars.filterOffSmallReadCounts) {
		smallDir = bib::files::makeDir(setUp.pars_.directoryName_, "smallReadCounts", false);
	}

	ExtractionStator stats = ExtractionStator(count, readsNotMatchedToBarcode,
			smallFragmentCount);
	std::map<std::string, uint32_t> goodCounts;

	for (const auto & f : barcodeFiles) {
		auto barcodeName = bib::files::getFileName(f.first.string());
		if ((counts[barcodeName].first + counts[barcodeName].second) == 0
				&& pars.multiplex) {
			//no reads extracted for barcode so skip filtering step
			continue;
		}

		if (pars.filterOffSmallReadCounts
				&& (counts[barcodeName].first + counts[barcodeName].second)
						<= pars.smallExtractReadCount) {
			auto barcodeOpts = setUp.pars_.ioOptions_;
			barcodeOpts.firstName_ = f.first.string();
			barcodeOpts.inFormat_ = SeqIOOptions::getInFormat(bib::files::getExtension(f.first.string()));
			barcodeOpts.out_.outFilename_ = smallDir + barcodeName;
			SeqIO barcodeIn(barcodeOpts);
			barcodeIn.openIn();
			readObject read;
			while (barcodeIn.readNextRead(read)) {
				barcodeIn.openWrite(read);
			}
			continue;
		}

		if (pars.multiplex) {
			std::cout
					<< bib::bashCT::boldGreen("Filtering on barcode: " + barcodeName)
					<< std::endl;
		} else {
			std::cout << bib::bashCT::boldGreen("Filtering") << std::endl;
		}
		auto barcodeOpts = setUp.pars_.ioOptions_;
		barcodeOpts.firstName_ = f.first.string();
		barcodeOpts.inFormat_ = SeqIOOptions::getInFormat(
				bib::files::getExtension(f.first.string()));
		SeqIO barcodeIn(barcodeOpts);
		barcodeIn.openIn();

		//create outputs
		MultiSeqIO midReaderOuts;
		auto unrecogPrimerOutOpts = setUp.pars_.ioOptions_;
		unrecogPrimerOutOpts.out_.outFilename_ = unrecognizedPrimerDir
				+ barcodeName;
		midReaderOuts.addReader("unrecognized", unrecogPrimerOutOpts);
		for (const auto & primerName : getVectorOfMapKeys(pDetermine.primers_)) {
			std::string fullname = primerName;
			if (pars.multiplex) {
				fullname += barcodeName;
			} else if (pars.sampleName != "") {
				fullname += pars.sampleName;
			}
			//bad out
			auto badDirOutOpts = setUp.pars_.ioOptions_;
			badDirOutOpts.out_.outFilename_ = badDir + fullname;
			midReaderOuts.addReader(fullname + "bad", badDirOutOpts);
			//good out
			auto goodDirOutOpts = setUp.pars_.ioOptions_;
			goodDirOutOpts.out_.outFilename_ = setUp.pars_.directoryName_ + fullname;
			midReaderOuts.addReader(fullname + "good", goodDirOutOpts);
			//contamination out
			if (pars.screenForPossibleContamination) {
				auto contamOutOpts = setUp.pars_.ioOptions_;
				contamOutOpts.out_.outFilename_ = contaminationDir + fullname;
				midReaderOuts.addReader(fullname + "contamination", contamOutOpts);
			}
			if (pars.mothurExtract || pars.pyroExtract) {
				auto flowExtractOutOpts = setUp.pars_.ioOptions_;
				flowExtractOutOpts.out_.outFilename_ = setUp.pars_.directoryName_ + fullname + "_temp";
				flowExtractOutOpts.outFormat_ = SeqIOOptions::outFormats::FLOW;
				flowExtractOutOpts.out_.outExtention_ = ".dat";
				midReaderOuts.addReader(fullname + "flow", flowExtractOutOpts);
			}
		}

		uint32_t barcodeCount = 1;
		bib::ProgressBar pbar(
				counts[barcodeName].first + counts[barcodeName].second);
		pbar.progColors_ = pbar.RdYlGn_;
		std::string readFlows = "";
		std::ifstream inFlowFile;
		if(pars.mothurExtract || pars.pyroExtract){
			inFlowFile.open(unfilteredByBarcodesFlowDir + barcodeName + "_temp.dat");
		}



		while (barcodeIn.readNextRead(read)) {
			std::getline(inFlowFile, readFlows);
			pbar.outputProg(std::cout, barcodeCount);
			++barcodeCount;
			//filter on primers
			//forward
			std::string primerName = "";
			bool foundInReverse = false;
			if (pars.noForwardPrimer) {
				primerName = primerTable.content_.front()[primerTable.getColPos(
						"geneName")];
			} else {
				if (pars.multiplex) {
					primerName = pDetermine.determineForwardPrimer(read, 0, alignObj,
							pars.fPrimerErrors, !pars.forwardPrimerToUpperCase,
							setUp.pars_.weightHomopolymers_);
					if (primerName == "unrecognized" && pars.checkComplement) {
						primerName = pDetermine.determineWithReversePrimer(read, 0,
								alignObj, pars.fPrimerErrors, !pars.forwardPrimerToUpperCase,
								setUp.pars_.weightHomopolymers_);
						if (read.seqBase_.on_) {
							foundInReverse = true;
						}
					}
				} else {
					uint32_t start = 0;
					if (pars.variableStart) {
						start = pars.variableStop;
					}
					primerName = pDetermine.determineForwardPrimer(read, start, alignObj,
							pars.fPrimerErrors, !pars.forwardPrimerToUpperCase,
							setUp.pars_.weightHomopolymers_);
					if (primerName == "unrecognized" && pars.checkComplement) {
						primerName = pDetermine.determineWithReversePrimer(read, start,
								alignObj, pars.fPrimerErrors, !pars.forwardPrimerToUpperCase,
								setUp.pars_.weightHomopolymers_);
						if (read.seqBase_.on_) {
							foundInReverse = true;
						}
					}
				}

				if (!read.seqBase_.on_) {
					stats.increaseFailedForward(barcodeName, read.seqBase_.name_);
					midReaderOuts.openWrite("unrecognized", read);
					continue;
				}
			}
			std::string fullname = primerName;
			if (pars.multiplex) {
				fullname += barcodeName;
			} else if (pars.sampleName != "") {
				fullname += pars.sampleName;
			}
			//look for possible contamination
			if (pars.screenForPossibleContamination) {
				kmerInfo compareInfo(compareObject.seqBase_.seq_,
						pars.contaminationKLen, false);
				if (pars.multipleTargets) {
					if (compareInfos.find(primerName) != compareInfos.end()) {
						if (foundInReverse) {
							readChecker::checkReadOnKmerComp(read.seqBase_,
									compareInfosRev[primerName], pars.contaminationKLen,
									pars.kmerCutOff, true);
						} else {
							readChecker::checkReadOnKmerComp(read.seqBase_,
									compareInfos[primerName], pars.contaminationKLen,
									pars.kmerCutOff, true);
						}

						if (!read.seqBase_.on_) {
							stats.increaseCounts(fullname, read.seqBase_.name_,
									ExtractionStator::extractCase::CONTAMINATION);
							midReaderOuts.openWrite(fullname + "contamination", read);
							continue;
						}
					} else {
						std::cerr
								<< "Error in screening for contamination, multiple targets turned on but no contamination found for "
								<< primerName << std::endl;
						std::cerr << "Options are: "
								<< vectorToString(getVectorOfMapKeys(compareInfos), ",")
								<< std::endl;
					}
				} else {
					if (foundInReverse) {
						readChecker::checkReadOnKmerComp(read.seqBase_, compareInfoRev,
								pars.contaminationKLen, pars.kmerCutOff, true);
					} else {
						readChecker::checkReadOnKmerComp(read.seqBase_, compareInfo,
								pars.contaminationKLen, pars.kmerCutOff, true);
					}
					if (!read.seqBase_.on_) {
						stats.increaseCounts(fullname, read.seqBase_.name_,
								ExtractionStator::extractCase::CONTAMINATION);
						midReaderOuts.openWrite(fullname + "contamination", read);
						continue;
					}
				}
			}
			//min len
			if (pars.multipleTargets) {
				if (multipleLenCutOffs.find(primerName) != multipleLenCutOffs.end()) {
					readChecker::checkReadLenAbove(read.seqBase_,
							multipleLenCutOffs[primerName].minLen_, true);
				} else {
					readChecker::checkReadLenAbove(read.seqBase_, pars.minLen, true);
				}
			} else {
				readChecker::checkReadLenAbove(read.seqBase_, pars.minLen, true);
			}
			if (!read.seqBase_.on_) {
				stats.increaseCounts(fullname, read.seqBase_.name_,
						ExtractionStator::extractCase::MINLENBAD);
				midReaderOuts.openWrite(fullname + "bad", read);
				continue;
			}

			//reverse
			if (!pars.noReversePrimer) {
				if (foundInReverse) {
					pDetermine.checkForForwardPrimerInRev(read, primerName, alignObj,
							pars.rPrimerErrors, !pars.reversePrimerToUpperCase,
							setUp.pars_.weightHomopolymers_, pars.variableStop, false);
				} else {
					pDetermine.checkForReversePrimer(read, primerName, alignObj,
							pars.rPrimerErrors, !pars.reversePrimerToUpperCase,
							setUp.pars_.weightHomopolymers_, pars.variableStop, false);
				}

				if (!read.seqBase_.on_) {
					stats.increaseCounts(fullname, read.seqBase_.name_,
							ExtractionStator::extractCase::BADREVERSE);
					read.seqBase_.name_.append("_badReverse");
					midReaderOuts.openWrite(fullname + "bad", read);
					continue;
				}
			}
			//min len again becuase the reverse primer search trims to the reverse primer so it could be short again
			if (pars.multipleTargets) {
				if (multipleLenCutOffs.find(primerName) != multipleLenCutOffs.end()) {
					readChecker::checkReadLenAbove(read.seqBase_,
							multipleLenCutOffs[primerName].minLen_, true);
				} else {
					readChecker::checkReadLenAbove(read.seqBase_, pars.minLen, true);
				}
			} else {
				readChecker::checkReadLenAbove(read.seqBase_, pars.minLen, true);
			}
			if (!read.seqBase_.on_) {
				stats.increaseCounts(fullname, read.seqBase_.name_,
						ExtractionStator::extractCase::MINLENBAD);
				midReaderOuts.openWrite(fullname + "bad", read);
				continue;
			}

			//if found in the reverse direction need to re-orient now
			if (foundInReverse) {
				read.seqBase_.reverseComplementRead(true, true);
			}

			//contains n
			readChecker::checkReadOnSeqContaining(read.seqBase_, "N", pars.numberOfNs,
					true);
			if (!read.seqBase_.on_) {
				stats.increaseCounts(fullname, read.seqBase_.name_,
						ExtractionStator::extractCase::CONTAINSNS);
				midReaderOuts.openWrite(fullname + "bad", read);
				continue;
			}

			//max len
			if (pars.multipleTargets) {
				if (multipleLenCutOffs.find(primerName) != multipleLenCutOffs.end()) {
					readChecker::checkReadLenBellow(read.seqBase_,
							multipleLenCutOffs[primerName].maxLen_, true);
				} else {
					readChecker::checkReadLenBellow(read.seqBase_, pars.maxLength, true);
				}
			} else {
				readChecker::checkReadLenBellow(read.seqBase_, pars.maxLength, true);
			}
			if (!read.seqBase_.on_) {
				stats.increaseCounts(fullname, read.seqBase_.name_,
						ExtractionStator::extractCase::MAXLENBAD);
				midReaderOuts.openWrite(fullname + "bad", read);
				continue;
			}

			//quality
			if (pars.checkingQCheck) {
				readChecker::checkReadQualCheck(read.seqBase_, pars.qualCheck,
						pars.qualCheckCutOff, true);
			} else {
				if (pars.qualWindowTrim) {
					readChecker::checkReadOnQualityWindowTrim(read.seqBase_,
							pars.qualityWindowLength, pars.qualityWindowStep,
							pars.qualityWindowThres, pars.minLen, true);
				} else {
					readChecker::checkReadOnQualityWindow(read.seqBase_,
							pars.qualityWindowLength, pars.qualityWindowStep,
							pars.qualityWindowThres, true);
				}
			}

			if (!read.seqBase_.on_) {
				stats.increaseCounts(fullname, read.seqBase_.name_,
						ExtractionStator::extractCase::QUALITYFAILED);
				midReaderOuts.openWrite(fullname + "bad", read);
				continue;
			}

			if (read.seqBase_.on_) {
				stats.increaseCounts(fullname, read.seqBase_.name_,
						ExtractionStator::extractCase::GOOD);
				if (pars.rename) {
					std::string oldName = replaceString(read.seqBase_.name_, "_Comp", "");
					read.seqBase_.name_ = fullname + "."
							+ leftPadNumStr(goodCounts[fullname],
									counts[barcodeName].first + counts[barcodeName].second);
					if (containsSubString(oldName, "_Comp")) {
						read.seqBase_.name_.append("_Comp");
					}
					renameKeyFile << oldName << "\t" << read.seqBase_.name_ << "\n";
				}
				midReaderOuts.openWrite(fullname + "good", read);
				if(pars.mothurExtract || pars.pyroExtract){
					midReaderOuts.openWrite(fullname + "flow", readFlows);
				}
				++goodCounts[fullname];
			}
		}
		std::cout << std::endl;
	}
	if(pars.mothurExtract || pars.pyroExtract){
		if(pars.mothurExtract){
			for(const auto & name : goodCounts){
				std::ofstream mothurOut;
				openTextFile(mothurOut,setUp.pars_.directoryName_ + name.first,".flow",false, true);
				uint32_t mothurExtractFlowNum = 800;
				if (pars.maxFlowCutoff == 720) {
					mothurExtractFlowNum = 800;
				} else if (pars.maxFlowCutoff <= 400) {
					mothurExtractFlowNum = 400;
				}
				mothurOut << mothurExtractFlowNum << std::endl;
				std::ifstream flowFile(setUp.pars_.directoryName_ + name.first + "_temp.dat");
				for(std::string line; std::getline(flowFile,line);){
					mothurOut << line << std::endl;
				}
			}
		}
		if(pars.pyroExtract){
			for(const auto & name : goodCounts){
				std::ofstream ampNoiseOut;
				openTextFile(ampNoiseOut,setUp.pars_.directoryName_ + name.first,".dat",false, true);
				ampNoiseOut << name.second << " " << pars.maxFlowCutoff << std::endl;
				std::ifstream flowFile(setUp.pars_.directoryName_ + name.first + "_temp.dat");
				for(std::string line; std::getline(flowFile,line);){
					ampNoiseOut << line << std::endl;
				}
			}
		}
		for(const auto & name : goodCounts){
			bib::files::bfs::remove(setUp.pars_.directoryName_ + name.first + "_temp.dat");
		}
	}
	std::ofstream profileLog;
	openTextFile(profileLog, setUp.pars_.directoryName_ + "extractionProfile.tab.txt",
			".txt", false, false);

	profileLog
			<< "name\ttotalReadsExtracted\tgoodReadsExtracted\tforGood\trevGood\ttotalBadReads\t"
					"badReverse\tcontainsNs\tlen<" << pars.minLen << "\tlen>"
			<< pars.maxLength;
	if (pars.checkingQCheck) {
		profileLog << "\tq" + estd::to_string(pars.qualCheck) + "<"
				<< pars.qualCheckCutOff;
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
			<< "TotalReads\tReadsNotMatchedBarcodes\tSmallFragments(len<"
			<< pars.smallFragmentCutoff
			<< ")\tfailedForwardPrimer\tfailedQualityFiltering\tused";
	extractionStatsFile << "\tcontamination";
	extractionStatsFile << std::endl;
	stats.outTotalStats(extractionStatsFile, "\t");
	std::ofstream failedForwadFile;
	openTextFile(failedForwadFile, setUp.pars_.directoryName_ + "failedForwad.tab.txt",
			".txt", false, false);
	failedForwadFile << "MidName\ttotalFailed\tfailedInFor\tfailedInRev"
			<< std::endl;
	stats.outFailedForwardStats(failedForwadFile, "\t");

	if (!setUp.pars_.debug_) {
		bib::files::rmDirForce(unfilteredReadsDir);
	}
	if (setUp.pars_.writingOutAlnInfo_) {
		setUp.rLog_ << "Number of alignments done" << "\n";
		alignObj.alnHolder_.write(setUp.pars_.outAlnInfoDirName_, setUp.pars_.verbose_);
	}
	setUp.logRunTime(std::cout);
	return 0;
}

int SeekDeepRunner::qluster(const bib::progutils::CmdArgs & inputCommands) {
	SeekDeepSetUp setUp(inputCommands);
	// parameters
	clusterDownPars pars;

	setUp.setUpClusterDown(pars);
	Json::Value metaData;
	auto analysisDirPath = bib::files::bfs::canonical(setUp.pars_.directoryName_);
	metaData["analysisDirPath"] = analysisDirPath.string();
	// print out the parameters read in
	if (pars.useNucComp && setUp.pars_.verbose_) {
		std::cout << "Nucleotide Composition Binning cut offs" << std::endl;
		printVector(pars.diffCutOffVec, ", ", std::cout);
	}
	// read in the sequences
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto reads = reader.readAllReads<readObject>();
	auto splitOnSize = readVecSplitter::splitVectorBellowLength(reads,
			pars.smallReadSize);
	reads = splitOnSize.first;
	if (!splitOnSize.second.empty()) {
		SeqOutput::write(splitOnSize.second,SeqIOOptions(setUp.pars_.directoryName_ + "smallReads",
				setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
	}
	if (pars.removeLowQualBases) {
		readVec::allRemoveLowQualityBases(reads, pars.lowQualityCutOff);
	}
	if (setUp.pars_.adjustHomopolyerRuns_) {
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
	processRunCutoff(setUp.pars_.runCutoff_, setUp.pars_.runCutOffString_,
			readVec::getTotalReadCount(reads));
	if (setUp.pars_.verbose_ && !pars.onPerId) {
		std::cout << "Kmer Low Frequency Error Cut off Is: " << setUp.pars_.runCutoff_
				<< std::endl;
	}
	// make the runLog, this is what is seen on the terminal screen at run time
	setUp.startARunLog(setUp.pars_.directoryName_);
	// parameter file
	setUp.writeParametersFile(setUp.pars_.directoryName_ + "parametersUsed.txt", false,
			true);
	// read in the paramteres from the parameters file
	setUp.rLog_ << "Parameters used" << "\n";
	table parTab(pars.parameters, ":");
	parTab.outPutContents(setUp.rLog_.runLogFile_, ":");
	if (setUp.pars_.verbose_) {
		std::cout << "Reading clusters from " << setUp.pars_.ioOptions_.firstName_ << " "
				<< setUp.pars_.ioOptions_.secondName_ << std::endl;
		std::cout << "Read in " << counter << " reads" << std::endl;
	}
	setUp.rLog_ << "Reading clusters from " << setUp.pars_.ioOptions_.firstName_ << " "
			<< setUp.pars_.ioOptions_.secondName_ << "\n";
	setUp.rLog_ << "Read in " << counter << " reads" << "\n";
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

	KmerMaps kMaps = indexKmers(clusters, setUp.pars_.kLength_, setUp.pars_.runCutoff_,
			setUp.pars_.kmersByPosition_, setUp.pars_.expandKmerPos_, setUp.pars_.expandKmerSize_);
	// create aligner class object
	aligner alignerObj = aligner(maxSize, setUp.pars_.gapInfo_, setUp.pars_.scoring_, kMaps,
			setUp.pars_.qScorePars_,
			setUp.pars_.countEndGaps_);
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
		;
	}

	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);

	collapser collapserObj = collapser(!setUp.pars_.largestFirst_,
			setUp.pars_.bestMatchCheck_, setUp.pars_.local_, true, setUp.pars_.kmersByPosition_,
			setUp.pars_.runCutoff_, setUp.pars_.kLength_, setUp.pars_.verbose_, !setUp.pars_.largestFirst_,
			true, setUp.pars_.weightHomopolymers_, setUp.pars_.skipOnLetterCounterDifference_,
			setUp.pars_.fractionDifferenceCutOff_, setUp.pars_.adjustHomopolyerRuns_);
	collapserObj.opts_.lowQualityBaseTrim_ = pars.lowQualityCutOff;
	collapserObj.opts_.removeLowQualityBases_ = pars.removeLowQualBases;
	collapserObj.opts_.debug_ = setUp.pars_.debug_;
	collapserObj.opts_.noAlign_ = pars.noAlign_;

	uint32_t singletonNum = 0;
	std::vector<cluster> singletons;
	if (!pars.startWithSingles) {
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
	}

	//run clustering
	collapserObj.runFullClustering(clusters, pars.onPerId, pars.iteratorMap,
			pars.binIteratorMap, pars.useNucComp, pars.useMinLenNucComp,
			pars.findBest, pars.diffCutOffVec, pars.useKmerBinning, pars.kCompareLen,
			pars.kmerCutOff, alignerObj, setUp.pars_, pars.snapShots, "firstSnaps");
	//run again with singlets if needed
	if (!pars.startWithSingles && !pars.leaveOutSinglets) {
		addOtherVec(clusters, singletons);
		collapserObj.runFullClustering(clusters, pars.onPerId, pars.iteratorMap,
				pars.binIteratorMap, pars.useNucComp, pars.useMinLenNucComp,
				pars.findBest, pars.diffCutOffVec, pars.useKmerBinning,
				pars.kCompareLen, pars.kmerCutOff, alignerObj, setUp.pars_, pars.snapShots,
				"secondSnaps");
	}

	//remove reads if they are made up of reads only in one direction
	if (containsCompReads && pars.useCompPerCutOff) {
		for (auto& clus : clusters) {
			if (containsCompReads) {
				uint32_t currentCompAmount = 0;
				for(const auto & read : clus.reads_){
					if(containsSubString(read->seqBase_.name_, "_Comp")){
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
	std::string seqName = bib::files::getFileName(setUp.pars_.ioOptions_.firstName_);
	readVecSorter::sort(clusters);
	renameReadNames(clusters, seqName, true, false, false);


	if (pars.markChimeras) {
		std::ofstream chimerasInfoFile;
		openTextFile(chimerasInfoFile,
				setUp.pars_.directoryName_ + "chimeraNumberInfo.txt", ".txt", false, false);
		chimerasInfoFile << "#chimericClusters\t#chimericReads" << std::endl;

		comparison chiOverlap;
		chiOverlap.oneBaseIndel_ = 2;
		chiOverlap.twoBaseIndel_ = 1;
		chiOverlap.largeBaseIndel_ = .99;
		chiOverlap.lowKmerMismatches_ = 1;
		uint32_t overLapSizeCutoff = 5;
		uint32_t allowableError = 0;
		uint32_t chiCount = 0;
		uint32_t chiRunCutOff = 1;

		collapserObj.markChimerasAdvanced(clusters, alignerObj, pars.parFreqs,
				chiRunCutOff, chiOverlap, overLapSizeCutoff, chiCount, allowableError);
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
		SeqOutput::write(clusters,
				SeqIOOptions(
						setUp.pars_.directoryName_
								+ setUp.pars_.ioOptions_.out_.outFilename_ + "_befroeTanCol",
								setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
		std::cout << "Collapsing on tandem repeat gaps" << std::endl;
		std::cout << "Starting with " << clusters.size() << " clusters"
				<< std::endl;
		clusterCollapser::collapseTandems(clusters, alignerObj, setUp.pars_.runCutoff_,
				setUp.pars_.kLength_, setUp.pars_.kmersByPosition_, pars.parFreqs, setUp.pars_.local_,
				true);
		clusters = readVecSplitter::splitVectorOnRemove(clusters).first;
		std::cout << "Collapsed down to " << clusters.size() << " clusters"
				<< std::endl;
	}
	if (setUp.pars_.refIoOptions_.firstName_ == "") {
		profiler::getFractionInfoCluster(clusters, setUp.pars_.directoryName_,
				"outputInfo");
	} else {
		profiler::getFractionInfoCluster(clusters, setUp.pars_.directoryName_,
				"outputInfo", setUp.pars_.refIoOptions_.firstName_, alignerObj, setUp.pars_.local_,
				true);
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
				pars.additionalOutLocationFile, setUp.pars_.ioOptions_.firstName_);
		if (additionalOutDir == "") {
			std::cerr << bib::bashCT::red << bib::bashCT::bold;
			std::cerr << "No additional out directory found for: "
					<< setUp.pars_.ioOptions_.firstName_ << std::endl;
			std::cerr << bib::bashCT::reset;
		} else {
			SeqOutput::write(clusters, SeqIOOptions(additionalOutDir + setUp.pars_.ioOptions_.out_.outFilename_,
					setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
			std::ofstream metaDataFile;
			openTextFile(metaDataFile, additionalOutDir + "/" + "metaData", ".json",
					setUp.pars_.ioOptions_.out_);
			metaDataFile << metaData;
		}

	}

	SeqOutput::write(clusters,
			SeqIOOptions(
					setUp.pars_.directoryName_ + setUp.pars_.ioOptions_.out_.outFilename_,
					setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
	{
		std::string snpDir = bib::files::makeDir(setUp.pars_.directoryName_,
				"internalSnpInfo", false);
		for (const auto & readPos : iter::range(clusters.size())) {
			std::unordered_map<uint32_t,
					std::unordered_map<char, std::vector<baseReadObject>>>mismatches;

			for (const auto & subReadPos : iter::range(
							clusters[readPos].reads_.size())) {
				const auto & subRead = clusters[readPos].reads_[subReadPos];
				alignerObj.alignCacheGlobal(clusters[readPos], subRead);
				//count gaps and mismatches and get identity
				alignerObj.profilePrimerAlignment(clusters[readPos], subRead,
						setUp.pars_.weightHomopolymers_);
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
		std::string minTreeDirname = bib::files::makeDir(setUp.pars_.directoryName_,
				"minTree", false);
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
				numThreads, setUp.pars_.weightHomopolymers_);
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

	std::string clusterDirectoryName = bib::files::makeDir(setUp.pars_.directoryName_,
			"clusters", false);
	clusterVec::allWriteClustersInDir(clusters, clusterDirectoryName,
			setUp.pars_.ioOptions_);

	if (!setUp.pars_.ioOptions_.processed_) {
		std::ofstream compStats;
		if (containsCompReads) {
			openTextFile(compStats, setUp.pars_.directoryName_ + "compStats.tab.txt",
					".txt", false, false);
			compStats << "cluster\tcompAmount" << std::endl;
		}
		std::string allInputReadsDir = bib::files::makeDir(setUp.pars_.directoryName_,
				"allInputReadsForEachCluster", false);
		for (const auto& clus : clusters) {
			SeqOutput writer(
					SeqIOOptions(allInputReadsDir + clus.seqBase_.name_,
							setUp.pars_.ioOptions_.outFormat_,setUp.pars_.ioOptions_.out_));
			writer.openOut();
			uint32_t currentCompAmount = 0;
			for (const auto & read : clus.reads_) {
				auto input = readVec::getReadByName(identicalClusters,read->seqBase_.name_);
				for(const auto & inputRead : input.reads_){
					writer.write(inputRead);
				}
				if (containsSubString(read->seqBase_.name_, "_Comp")) {
					currentCompAmount += read->seqBase_.cnt_;
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

int SeekDeepRunner::processClusters(const bib::progutils::CmdArgs & inputCommands) {
	// parameters
	SeekDeepSetUp setUp(inputCommands);
	processClustersPars pars;
	setUp.setUpMultipleSampleCluster(pars);

	if (containsSubString(pars.experimentName, ".")) {
		throw std::runtime_error {
				"Error in populationCollapse::populationCollapse, populationName can't contain '.', "
						+ pars.experimentName };
	}

	table customCutOffsTab;
	std::map<std::string, double> customCutOffsMap;
	if (pars.customCutOffs != "") {
		customCutOffsTab = table(pars.customCutOffs, "whitespace", true);
		if (!bib::in(std::string("sample"), customCutOffsTab.columnNames_)
				|| !bib::in(std::string("cutOff"), customCutOffsTab.columnNames_)) {
			throw std::runtime_error { "Error in loading custom cut off file, "
					+ pars.customCutOffs + ", missing sample or cutOff\n"
					+ vectorToString(customCutOffsTab.columnNames_, ",") };
		}
		for (const auto & rowPos : iter::range(customCutOffsTab.content_.size())) {
			customCutOffsMap[customCutOffsTab.content_[rowPos][customCutOffsTab.getColPos(
					"sample")]] =
					bib::lexical_cast<double>(
							customCutOffsTab.content_[rowPos][customCutOffsTab.getColPos(
									"cutOff")]);
		}
	}

	//read in the files in the corresponding sample directories
	std::map<std::string, std::pair<std::string, bool>> files = listFilesInDir(
			".", true);
	VecStr specificFiles;
	for (const auto& fileIter : files) {
		if (fileIter.first.find(setUp.pars_.ioOptions_.firstName_) != std::string::npos
				&& fileIter.first.find(".qual") == std::string::npos) {
			specificFiles.push_back(fileIter.first);
		}
	}
	std::cout << "Reading from" << std::endl;
	for (const auto& sfIter : specificFiles) {
		std::cout << sfIter << std::endl;
	}

	setUp.startARunLog(setUp.pars_.directoryName_);
	// parameters file
	setUp.writeParametersFile(setUp.pars_.directoryName_ + "parametersUsed.txt", false,
			false);
	// reading in reads

	uint64_t maxSize = 0;
	std::map<std::pair<std::string, std::string>, std::vector<readObject>> outputReads;
	for (const auto& strIter : specificFiles) {
		SeqIOOptions inOpts;
		inOpts.firstName_ = strIter;
		inOpts.secondName_ = strIter + ".qual";
		inOpts.processed_ = true;
		inOpts.inFormat_ = setUp.pars_.ioOptions_.inFormat_;
		inOpts.lowerCaseBases_ = setUp.pars_.ioOptions_.lowerCaseBases_;
		SeqInput reader(inOpts);
		reader.openIn();
		auto reads = reader.readAllReads<readObject>();
		readVec::getMaxLength(reads, maxSize);
		VecStr toks = tokenizeString(strIter, "/");
		outputReads.insert(
				std::make_pair(std::make_pair(toks[1], toks[2]), reads));
	}
	// output info about the read In reads
	for (const auto& readsIter : outputReads) {
		std::cout << readsIter.first.first << ": " << readsIter.first.second
				<< " clustersNum: " << readsIter.second.size() << " readsNum: "
				<< readVec::getTotalReadCount(readsIter.second) << std::endl;
		setUp.rLog_ << readsIter.first.first << ": " << readsIter.first.second
				<< " clustersNum: " << readsIter.second.size() << " readsNum: "
				<< readVec::getTotalReadCount(readsIter.second) << "\n";
	}
	// reading expected sequences to compare to
	bool checkingExpected = setUp.pars_.refIoOptions_.firstName_ != "";
	std::vector<readObject> expectedSeqs;
	if (checkingExpected) {
		expectedSeqs = SeqInput::getReferenceSeq(setUp.pars_.refIoOptions_, maxSize);
	}
	// create aligner class object
	aligner alignerObj(maxSize, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
			KmerMaps(), setUp.pars_.qScorePars_, setUp.pars_.countEndGaps_);
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);

	collapser collapserObj = collapser(!setUp.pars_.firstMatch_, setUp.pars_.bestMatchCheck_,
			setUp.pars_.local_, true, setUp.pars_.kmersByPosition_, setUp.pars_.runCutoff_,
			setUp.pars_.kLength_, setUp.pars_.verbose_, !setUp.pars_.largestFirst_, true,
			setUp.pars_.weightHomopolymers_, setUp.pars_.skipOnLetterCounterDifference_,
			setUp.pars_.fractionDifferenceCutOff_, setUp.pars_.adjustHomopolyerRuns_);

	collapserObj.opts_.lowQualityBaseTrim_ = pars.lowQualityCutOff;
	collapserObj.opts_.removeLowQualityBases_ = pars.removeLowQualBases;
	collapserObj.opts_.debug_ = setUp.pars_.debug_;

	// collect all the reads together for each sample
	std::map<std::string, std::vector<std::vector<cluster>>>clustersBySample;

	for (auto& readsIter : outputReads) {
		std::vector<cluster> clusters = baseCluster::convertVectorToClusterVector<
				cluster>(readsIter.second);
		readVecSorter::sortReadVector(clusters, "totalCount");
		// consider adding the sample name in the name as well
		renameReadNamesNewClusters(clusters, readsIter.first.second, true, true,
				false);
		if (pars.checkChimeras) {
			// readVec::allSetFractionByTotalCount(clusters);
			clusterVec::allSetFractionClusters(clusters);
			comparison chiOverlap;
			chiOverlap.oneBaseIndel_ = 2;
			chiOverlap.twoBaseIndel_ = 1;
			chiOverlap.lowKmerMismatches_ = 1;
			uint32_t overLapSizeCutoff = 5;
			uint32_t allowableError = 0;
			uint32_t chiCount = 0;
			uint32_t chiRunCutOff = 1;
			collapserObj.markChimerasAdvanced(clusters, alignerObj, pars.parFreqs,
					chiRunCutOff, chiOverlap, overLapSizeCutoff, chiCount,
					allowableError);
		}
		clusterVec::allSetFractionClusters(clusters);
		clustersBySample[readsIter.first.first].emplace_back(clusters);
	}
	std::map<std::string, collapse::sampleCollapse> sampleCollapses;

	collapserObj.opts_.checkKmers_ = false;

	for (const auto& samp : clustersBySample) {
		std::cout << "Currently on : " << samp.first << std::endl;
		sampleCollapses[samp.first] = collapse::sampleCollapse(samp.second,
				samp.first, pars.clusterCutOff);
		// std::cout << "made sampleCollapse fine" << std::endl;
		// sampleCollapses[samp.first].updateInitialInfos(false);
		if (pars.onPerId) {
			sampleCollapses[samp.first].clusterOnPerId(collapserObj, pars.iteratorMap,
					pars.sortBy, alignerObj);
		} else {
			sampleCollapses[samp.first].cluster(collapserObj, pars.iteratorMap,
					pars.sortBy, alignerObj);
		}

		// std::cout << "clustered fine" << std::endl;
		if (checkingExpected) {
			sampleCollapses[samp.first].collapsed_.checkAgainstExpected(expectedSeqs,
					alignerObj, setUp.pars_.local_, setUp.pars_.weightHomopolymers_);
		}

		sampleCollapses[samp.first].updateCollapsedInfos();
		// std::cout << "updated collapse infos fine" << std::endl;
		sampleCollapses[samp.first].updateExclusionInfos();
		// std::cout << "updated excusesion fine" << std::endl;
		sampleCollapses[samp.first].renameClusters(pars.sortBy);
		// std::cout << "ranmed clusters fine" << std::endl;
	}
	collapse::populationCollapse popCollapse;
	infoPrinter::printSampleCollapseInfo(sampleCollapses, checkingExpected,
			setUp.pars_.directoryName_ + "allClustersInfo.tab.txt", popCollapse, false);

	std::vector<sampleCluster> allSamples;
	std::string originalsDir = bib::files::makeDir(setUp.pars_.directoryName_,
			"originals", false);
	std::string initialDir = bib::files::makeDir(setUp.pars_.directoryName_, "initial", false);
	std::string excludedDir = bib::files::makeDir(setUp.pars_.directoryName_,
			"excluded", false);
	std::string excludedInitialDir = "";
	if (pars.writeExcludedOriginals) {
		excludedInitialDir = bib::files::makeDir(setUp.pars_.directoryName_,
				"excludedInitial", false);
	}

	comparison onlyHpErrors;
	onlyHpErrors.oneBaseIndel_ = 1;
	onlyHpErrors.twoBaseIndel_ = 1;
	onlyHpErrors.largeBaseIndel_ = .99;

	std::string finalDir = bib::files::makeDir(setUp.pars_.directoryName_, "final", false);
	table chiInfoTab(VecStr { "sample", "numClustersSaved",
			"totalClustersChecked", "clusterSavedNames" });
	for (auto& sampCollapse : sampleCollapses) {
		// exclude
		if (0 != pars.runsRequired) {
			sampCollapse.second.excludeBySampNum(pars.runsRequired, false);
		} else {
			sampCollapse.second.excludeBySampNum(
					sampCollapse.second.input_.info_.infos_.size(), false);
		}
		VecStr clustersSavedFromChi;
		uint32_t clustersNotSaved = 0;
		if (!pars.keepChimeras) {
			if (pars.investigateChimeras) {
				//first mark the suspicous clusters as being chimeric
				for (auto &clus : sampCollapse.second.collapsed_.clusters_) {
					if (clus.isClusterAtLeastChimericCutOff(pars.chiCutOff)) {
						clus.seqBase_.markAsChimeric();
					}
				}
				//now check to see if it is ever the top two variants of a sample and if it is unmark it
				for (auto & clus : sampCollapse.second.collapsed_.clusters_) {
					if (clus.seqBase_.frac_ < pars.fracCutoff) {
						continue;
					}
					if (clus.seqBase_.name_.find("CHI_") != std::string::npos) {
						bool saved = false;
						for (const auto & otherCollapse : sampleCollapses) {

							if (otherCollapse.first != sampCollapse.first
									&& !containsSubString(
											stringToLowerReturn(otherCollapse.first), "control")
									&& !containsSubString(stringToLowerReturn(sampCollapse.first),
											"control")) {
								for (const auto & pos : iter::range(
										std::min<size_t>(2,
												otherCollapse.second.collapsed_.clusters_.size()))) {
									alignerObj.alignCacheGlobal(
											otherCollapse.second.collapsed_.clusters_[pos], clus);
									alignerObj.profilePrimerAlignment(
											otherCollapse.second.collapsed_.clusters_[pos], clus,
											setUp.pars_.weightHomopolymers_);
									if (onlyHpErrors.passErrorProfile(alignerObj.comp_)) {

										clus.seqBase_.unmarkAsChimeric();
										for (auto & subRead : clus.reads_) {
											subRead->seqBase_.unmarkAsChimeric();
										}
										clus.resetInfos();
										clustersSavedFromChi.emplace_back(clus.seqBase_.name_);
										saved = true;
										break;
									}
								}
							}
						}
						if (!saved) {
							++clustersNotSaved;
						}
					}
					sampCollapse.second.collapsed_.setSetInfo();
				}
			}
			//now exclude all marked chimeras, currently this will also remark chimeras unneccessarily
			sampCollapse.second.excludeChimeras(false, pars.chiCutOff);
		}
		chiInfoTab.content_.emplace_back(
				toVecStr(sampCollapse.first,
						getPercentageString(clustersSavedFromChi.size(),
								clustersSavedFromChi.size() + clustersNotSaved),
						clustersSavedFromChi.size() + clustersNotSaved,
						vectorToString(clustersSavedFromChi, ",")));
		if (bib::in(sampCollapse.first, customCutOffsMap)) {
			if (setUp.pars_.debug_) {
				std::cout << "Custom Cut off for " << sampCollapse.first << " : "
						<< customCutOffsMap[sampCollapse.first] << std::endl;
			}
			sampCollapse.second.excludeFraction(customCutOffsMap[sampCollapse.first],
					true);
		} else {
			sampCollapse.second.excludeFraction(pars.fracCutoff, true);
		}

		std::string sortBy = "fraction";
		sampCollapse.second.renameClusters(sortBy);
		// write
		sampCollapse.second.writeInitial(originalsDir, setUp.pars_.ioOptions_);
		sampCollapse.second.writeExcluded(excludedDir, setUp.pars_.ioOptions_);
		if (pars.writeExcludedOriginals) {
			sampCollapse.second.writeExcludedOriginalClusters(excludedInitialDir,
					setUp.pars_.ioOptions_);
		}
		sampCollapse.second.writeFinal(finalDir, setUp.pars_.ioOptions_);
		sampCollapse.second.writeFinalOrignalClusters(initialDir, setUp.pars_.ioOptions_);
		// add to all sample cluster

		addOtherVec(allSamples, sampCollapse.second.createOutput(false, sortBy));
		// sampCollapse.second.updateCollapsedInfos(true);
		// sampCollapse.second.updateExclusionInfos(true);
	}
	TableIOOpts chiOutOptions(OutOptions(setUp.pars_.directoryName_ + "chiInfo", ".tab.txt"),
			"\t", chiInfoTab.hasHeader_);
	chiInfoTab.outPutContents(chiOutOptions);
	std::unordered_map<std::string, std::unique_ptr<aligner>> aligners;
	uint32_t numThreads = 2;
	std::mutex alignerLock;
	if (!pars.noPopulation) {
		std::cout << bib::bashCT::boldGreen("Pop Clustering") << std::endl;
		popCollapse = collapse::populationCollapse(allSamples, pars.experimentName);

		std::string popDir = bib::files::makeDir(setUp.pars_.directoryName_,
				"population", false);
		if (pars.onPerId) {
			popCollapse.popClusterOnId(collapserObj, pars.popIteratorMap, "fraction",
					alignerObj);
		} else {
			popCollapse.popCluster(collapserObj, pars.popIteratorMap, "fraction",
					alignerObj);
		}
		popCollapse.updateInfoWithSampCollapses(sampleCollapses);

		if (pars.previousPopFilename != "") {
			SeqIOOptions popOpts;
			popOpts.firstName_ = pars.previousPopFilename;
			popOpts.inFormat_ = SeqIOOptions::getInFormat(bib::files::getExtension(pars.previousPopFilename));
			SeqInput reader(popOpts);
			reader.openIn();
			popCollapse.renameToOtherPopNames(reader.readAllReads<readObject>(), alignerObj);
		}

		if (checkingExpected) {
			popCollapse.collapsed_.checkAgainstExpected(expectedSeqs, alignerObj,
					setUp.pars_.local_, setUp.pars_.weightHomopolymers_);
		}

		infoPrinter::printPopulationCollapseInfo(popCollapse,
				popDir + "populationCluster.tab.txt", checkingExpected);
		std::string popInitialDir = bib::files::makeDir(popDir, "initial", false);
		popCollapse.writeFinalInitial(popInitialDir, setUp.pars_.ioOptions_);
		popCollapse.writeFinal(popDir, setUp.pars_.ioOptions_);

		std::unordered_map<std::string, bib::color> colorsForGraph;
		colorsForGraph = getColorsForNames(popCollapse.collapsed_.clusters_,
				pars.sat, pars.lum);

		std::string dotDir = bib::files::makeDir(setUp.pars_.directoryName_, "dotFiles", false);
		for (const auto& samp : sampleCollapses) {
			std::vector<readObject> readsForGraph;
			if (samp.second.collapsed_.clusters_.empty()) {
				continue;
			}
			for (const auto& clus : samp.second.collapsed_.clusters_) {
				std::string popName =
						popCollapse.collapsed_.clusters_[popCollapse.collapsed_.subClustersPositions_.at(
								clus.getStubName(true))].seqBase_.name_;
				readsForGraph.emplace_back(
						readObject(
								seqInfo(popName, clus.seqBase_.seq_, clus.seqBase_.qual_,
										clus.seqBase_.cnt_, clus.seqBase_.frac_)));
				readsForGraph.back().seqBase_.cnt_ = clus.seqBase_.cnt_;
				readsForGraph.back().seqBase_.frac_ = clus.seqBase_.frac_;
			}
			std::set<std::string> alreadyAdded;
			for (auto & read : readsForGraph) {
				std::string nName = read.getReadId();
				uint32_t num = 1;
				while (bib::in(nName, alreadyAdded)) {
					nName = read.getReadId() + "_" + estd::to_string(num);
					++num;
				}
				alreadyAdded.emplace(nName);
				read.seqBase_.name_ = nName;
			}

			auto readCompGraph = genReadComparisonGraph(readsForGraph, alignerObj,
					aligners, alignerLock, numThreads, setUp.pars_.weightHomopolymers_);
			readCompGraph.setMinimumEventConnections();
			std::vector<std::string> popNames;
			for (const auto & n : readCompGraph.nodes_) {
				if (n->on_) {
					popNames.emplace_back(n->name_);
				}
			}
			auto nameColors = getColorsForNames(popNames);
			Json::Value minTreeData = readCompGraph.toD3Json(bib::color("#000000"),
					nameColors);
			std::ofstream outJson(dotDir + samp.first + ".json");
			std::ofstream outHtml(dotDir + samp.first + ".html");
			std::ofstream outDot(dotDir + samp.first + ".dot");
			outJson << minTreeData;
			outHtml
					<< genHtmlStrForPsuedoMintree(samp.first + ".json",
							"http://bib8.umassmed.edu/~hathawan/js/psuedoMinTreeWithIndels.js");
			jsonTreeToDot(minTreeData, outDot);
		}
		for (const auto & alnObj : aligners) {
			if (setUp.pars_.debug_) {
				std::cout << alnObj.first << ": "
						<< alnObj.second->numberOfAlingmentsDone_ << std::endl;
			}
			if (alnObj.second->numberOfAlingmentsDone_ > 0) {
				alignerObj.alnHolder_.mergeOtherHolder(alnObj.second->alnHolder_);
			}
		}
		if (pars.groupingsFile != "") {
			table groupsTab(pars.groupingsFile, "\t", true);
			/**@todo add in safety checks */
			VecStr sampleNames = groupsTab.getColumn(0);
			bib::for_each(sampleNames, [](std::string & str) {stringToUpper(str);});
			std::unordered_map<std::string, std::unordered_map<std::string, VecStr>> groups;
			std::set<std::string> clustered;
			std::set<std::string> notClustered;
			for (const auto & colPos : iter::range<uint32_t>(1,
					groupsTab.columnNames_.size())) {
				VecStr col = groupsTab.getColumn(colPos);
				for (const auto & pos : iter::range(col.size())) {
					bool found = false;
					for (const auto & samp : sampleCollapses) {
						if (stringToUpperReturn(samp.first) == sampleNames[pos]) {
							groups[groupsTab.columnNames_[colPos]][col[pos]].emplace_back(
									samp.first);
							clustered.emplace(samp.first);
							found = true;
							break;
						}
					}
					if (!found) {
						notClustered.emplace(sampleNames[pos]);
					}
				}
			}
			for (const auto & samp : sampleCollapses) {
				if (!bib::in(samp.first, clustered)) {
					std::cout << "Sample: " << samp.first << " No data supplied for "
							<< std::endl;
				}
			}
			for (const auto nc : notClustered) {
				std::cout << "Sample: " << nc << " was not clustered" << std::endl;
			}
			std::string groupDir = bib::files::makeDir(setUp.pars_.directoryName_,
					"groups", false);
			for (const auto & group : groups) {
				std::string currentGroupDir = bib::files::makeDir(groupDir,
						group.first, false);
				for (const auto & subGroup : group.second) {
					std::string currentSubGroupDir = bib::files::makeDir(currentGroupDir,
							subGroup.first, false);
					std::ofstream popFile;
					std::ofstream sampFile;
					openTextFile(popFile, currentSubGroupDir + "popFile.tab.txt", ".txt",
							false, false);

					openTextFile(sampFile, currentSubGroupDir + "sampFile.tab.txt",
							".txt", false, false);
					std::set<uint32_t> otherPopPositions;
					for (const auto & otherGroup : group.second) {
						if (otherGroup.first != subGroup.first) {
							for (const auto & samp : otherGroup.second) {
								/**@todo should check to see if the sample exist*/
								for (const auto & subClus : sampleCollapses.at(samp).collapsed_.clusters_) {
									otherPopPositions.emplace(
											popCollapse.collapsed_.subClustersPositions_.at(
													subClus.getStubName(true)));
								}
							}
						}
					}
					infoPrinter::printInfoForSamps(sampleCollapses, checkingExpected,
							sampFile, popFile, popCollapse, !pars.noPopulation,
							subGroup.second, group.first + ":" + subGroup.first,
							otherPopPositions);
				}
			}
		}
	}

	infoPrinter::printSampleCollapseInfo(sampleCollapses, checkingExpected,
			setUp.pars_.directoryName_ + "selectedClustersInfo.tab.txt", popCollapse,
			!pars.noPopulation);

	std::string sampRepInfoDir = bib::files::makeDir(setUp.pars_.directoryName_,
			"sampRepAgreementInfo", false);
	table rmseTab(VecStr { "sampleName", "RMSE" });
	table repInfoTab(
			VecStr { "sampleName", "clusterName", "repName", "fraction" });
	for (const auto & samp : sampleCollapses) {
		auto repInfoforSample = samp.second.collapsed_.getReplicateInfo();
		for (const auto & row : repInfoforSample.content_) {
			VecStr addingRow { samp.first };
			addOtherVec(addingRow, row);
			repInfoTab.content_.emplace_back(addingRow);
		}

		//std::cout << samp.second.collapsed_.clusters_.size() << std::endl;
		if (!samp.second.collapsed_.clusters_.empty()) {
			rmseTab.content_.emplace_back(
					toVecStr(samp.first, samp.second.collapsed_.getRMSE()));
		}
	}
	repInfoTab.outPutContents(
			TableIOOpts(OutOptions(sampRepInfoDir + "replicatesFractions", ".tab.txt"), "\t",
					repInfoTab.hasHeader_));
	rmseTab.sortTable("RMSE", true);
	rmseTab.outPutContents(
			TableIOOpts(OutOptions(sampRepInfoDir + "RMSE.tab.txt", ".tab.txt"), "\t", rmseTab.hasHeader_));

	if (setUp.pars_.writingOutAlnInfo_) {
		alignerObj.alnHolder_.write(setUp.pars_.outAlnInfoDirName_);
	}
	std::cout << alignerObj.numberOfAlingmentsDone_ << std::endl;
	setUp.rLog_ << alignerObj.numberOfAlingmentsDone_ << "\n";

	setUp.logRunTime(std::cout);
	return 0;
}

int SeekDeepRunner::makeSampleDirectories(const bib::progutils::CmdArgs & inputCommands) {
	SeekDeepSetUp setUp(inputCommands);
	std::string sampleNameFilename = "";
	bool separatedDirs = false;
	setUp.setOption(separatedDirs, "--separatedDirs",
			"Create a separate directory for each index");
	setUp.setUpMakeSampleDirectories(sampleNameFilename);
	setUpSampleDirs(sampleNameFilename, setUp.pars_.directoryName_, separatedDirs);
	setUp.startARunLog(setUp.pars_.directoryName_);
	setUp.logRunTime(std::cout);
	return 0;
}

}  // namespace bib
