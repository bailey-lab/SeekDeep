//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012, 2015 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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

#include "SeekDeepRunner.hpp"
#include <bibcpp.h>
#include <njhRInside.h>
#include "SeekDeep/programs/SeekDeepUtils.h"


namespace bibseq {


SeekDeepRunner::SeekDeepRunner()
    : bib::progutils::oneRing({addRing<SeekDeepUtilsRunner>()},
    								{addFunc("extractor", extractor, false),
										 addFunc("sffExtractor", sffExtractor, false),
                     addFunc("processClusters", processClusters, false),
                     addFunc("qluster", qluster, false),
										 addFunc("clusterDown", qluster, true),
                     addFunc("makeSampleDirectories", makeSampleDirectories, false),
										 },
                    "SeekDeep") {}

template <typename T>
std::vector<T> splitVectorOnKDistSimpleAdd(
    const std::vector<T>& reads, const T& compareObject,uint32_t kLength, double kmerCutoff,
    std::vector<T>& badReads, uint32_t& splitCount) {
	std::vector<T> ret;
	kmerInfo compareInfo(compareObject.seqBase_.seq_, kLength, false);
  for (const auto& read : reads) {
  	//readChecker::checkReadOnKmerComp(read.seqBase_, compareInfo,pars.contaminationKLen, pars.kmerCutOff, true);
  	kmerInfo currentInfo(read.seqBase_.seq_, kLength, false);
  	auto dist = compareInfo.compareKmers(currentInfo);
  	if(dist.second < kmerCutoff){
      ++splitCount;
      badReads.emplace_back(read);
  	}else{
  		 ret.emplace_back(read);
  	}
  }
  return ret;
};

int SeekDeepRunner::sffExtractor(MapStrStr inputCommands) {

	SeekDeepSetUp setUp(inputCommands);
	std::string idFilename = "", idFileDelim = "whitespace";
  int numberOfNs = 1;
  bool findReversePrimer = true;
  bool screenForPossibleContamination = false;
  std::string compareSeq = "";
  bool flowFiltering = false;
  int maxFlowCutoff = 0;
  bool qualWindowTrim = false;
  bool condensed = false;
  bool unknown = false;
  int unknownPrimerSize = 10;
  int unknownMinSize = 100;
  bool multiplex = false;
  int qualityWindowLength = 0, qualityWindowStep = 0, qualityWindowThres = 0;
  int minLen = 200, within = 20;
  int maxLength = 300;
  bool mothurExtract = false;
  bool pyroExtract = false;
  bool reversePrimerToLowerCase = true;
  bool checkComplement = false;
  int smallFragmentCutoff = 50;
  double queryCoverageCutoff = 0.5;
  double percentIdentityCutoff = 0.8;
  bool HMP = false;
  uint32_t primerLen = 30;
  bool trimTcag = false;
  bool noForwardPrimer = false;
  bool trimForwardPrimer = false;
  uint32_t barcodeErrors = 0;
  bool rename = false;
  double forwardCoverage = 1;
  double forwardPercentageGaps = 0.01;
  uint32_t forwardNumOfMismatches = 0;
  double kmerCutOff = .20;
  uint32_t contaminationKLen = 7;
  setUp.setOption(contaminationKLen, "-contaminationKLen", "Contamination Kmer Length");
  setUp.setOption(kmerCutOff, "-kmerCutOff", "Kmer cut off for contamination check");
  setUp.setOption(mothurExtract, "-mothurExtract", "Extract Data Files for mothur analysis");
  setUp.setOption(forwardCoverage, "-forwardCoverage", "Amount of Foward Primer to find");
  setUp.setOption(forwardPercentageGaps, "-forwardPercentageGaps", "Percentage of gaps to allow in the Forward Primer");
  setUp.setOption(forwardNumOfMismatches, "-forwardNumOfMismatches", "Number of Mismatches to allow in Foward Primer");

  setUp.setOption(rename, "-rename", "Rename With Barcode Names");
  setUp.setOption(barcodeErrors, "-barcodeErrors", "Errors Allowed in Barcode");
  setUp.setOption(trimForwardPrimer, "-trimForwardPrimer", "Trim Forward Primer");
  setUp.setOption(noForwardPrimer, "-noForwardPrimer", "No Forward Primer Required");
  setUp.setOption(trimTcag, "-trimTcag", "trim TCAG(454 tag) if it is present ");
  if(setUp.setOption(HMP, "-HMP", "HMP")){
  	setUp.setOption(primerLen, "-len,-primerLen", "PrimerLen");
  }

  bool barcodesBothEnds = false;
  setUp.setOption(barcodesBothEnds, "-barcodeBothEnds", "Look for Barcodes in Both Primers");
  setUp.processDebug();
  bool variableStart = false;
  uint32_t variableStop = 50;
  variableStart = setUp.setOption(variableStop, "-variableStart", "variableStart");
  uint32_t qualCheck = 25;
  bool checkingQCheck = false;
  if(setUp.setOption(qualCheck, "-qualCheck", "Qual Check Level")){
  	checkingQCheck = true;
  }
  double qualCheckCutOff = 0.75;
  if(setUp.setOption(qualCheckCutOff, "-qualCheckCutOff",
  		"Cut Off for fraction of bases above qual check of " + estd::to_string(qualCheck))){
  	checkingQCheck = true;
  }
  setUp.setUpExtractor(
      idFilename, multiplex, condensed, minLen, maxLength, within,
      qualityWindowLength, qualityWindowStep, qualityWindowThres, unknown,
      unknownPrimerSize, unknownMinSize, findReversePrimer, queryCoverageCutoff,
      percentIdentityCutoff, numberOfNs, pyroExtract, reversePrimerToLowerCase,
      flowFiltering, maxFlowCutoff, checkComplement,
      screenForPossibleContamination, compareSeq, idFileDelim,
      smallFragmentCutoff, qualWindowTrim);

  int readsExtracted = 0;
  int readsNotMatchedToBarcode = 0;
  int readsFailedForwardPrimer = 0;
  int contaminationAmount = 0;
  int readsFinal = 0;
  int failedQualityFiltering = 0;

  // int readsExtractedComp=0; cumulativeFraction
  // int readsNotMatchedToBarcodeComp=0;
  int readsFailedForwardPrimerComp = 0;
  int contaminationAmountComp = 0;
  int readsFinalComp = 0;
  int failedQualityFilteringComp = 0;
  std::cout << "wl:" << qualityWindowLength << " ws:" << qualityWindowStep
            << " wt:" << qualityWindowThres << std::endl;
  // run log
  setUp.startARunLog(setUp.directoryName_);
  // parameter file
  setUp.writeParametersFile(setUp.directoryName_ + "parametersUsed.txt", false,
                            false);
  readObject compareObject;
  if (screenForPossibleContamination) {
    compareObject = readObject(seqInfo("Compare", compareSeq));
  }
  // read in reads and remove lower case bases indicating tech low quality like
  // tags and such
  std::cout << "Reading in reads:" << std::endl;
  readObjectIO reader;
  reader.read(setUp.ioOptions_);
  if(HMP){
  	bib::for_each(reader.reads, [&](readObject & read){
  		subStrToUpper(read.seqBase_.seq_, 4, primerLen);});
  }
  std::string tcagLower = "tcag";
  std::string tcagUpper = "TCAG";
  if(trimTcag){
  	for(auto & read : reader.reads){
  		if(std::equal(tcagLower.begin(), tcagLower.end(), read.seqBase_.seq_.begin()) ||
  				std::equal(tcagUpper.begin(), tcagUpper.end(), read.seqBase_.seq_.begin())){
  			read.trimFront(4);
  		}
  	}
  }
  //letterCounter mainCounter = readVec::getAverageLetterFraction(reader.reads);
  std::cout << "Read in " << reader.reads.size() << " reads from "
            << setUp.ioOptions_.firstName_ << std::endl;
  setUp.rLog_ << "Read in " << reader.reads.size() << " reads from "
         << setUp.ioOptions_.firstName_ << "\n";
  readsExtracted = reader.reads.size();
  readVec::removeLowerCaseBases(reader.reads);
  uint64_t maxReadSize = 0;
  readVec::getMaxLength(reader.reads, maxReadSize);
  std::cout << "Extracting" << std::endl;
  std::map<std::string, std::vector<readObject>> reads;

  // create aligner for primer identification
  auto scoreMatrixMap = substituteMatrix::createDegenScoreMatrix(1,-1);
  gapScoringParameters gapPars(setUp.gapInfo_);
  substituteMatrix scoreMatrix(scoreMatrixMap);
  kmerMaps emptyMaps;
  int primaryQual = 20, secondayQual = 15;
  bool countEndGaps = true;
	aligner alignObj = aligner(maxReadSize, gapPars, scoreMatrix, emptyMaps,
			primaryQual, setUp.qualThresWindow_, secondayQual, countEndGaps);
  // Read in the primers and the mids if present, current just assuming they are
  // will fix in future nick 10.24.2013
  table primerTable;
  if (!unknown) {
    primerTable = seqUtil::readPrimers(idFilename, idFileDelim, false);
  }

  //remove small fragment reads
  std::vector<readObject> smallFragments;
  uint32_t smallFragmentCount = 0;
  std::vector<readObject> largeReads;
  largeReads = readVecSplitter::splitVectorBellowLengthAdd(
      reader.reads, smallFragmentCutoff, smallFragments,
      smallFragmentCount);

  std::vector<readObject> unrecognizedBarcode;
  if (unknown) {
    sortReadsBySeqFront(largeReads, reads, unknownPrimerSize,
                        unrecognizedBarcode, unknownMinSize);
    readsNotMatchedToBarcode = unrecognizedBarcode.size();
  } else {
    if (!multiplex) {
      reads.insert({"", largeReads});
    } else {
      unrecognizedBarcode = readVecExtractor::extractReadsByMIDFull(
          largeReads, reads, idFilename, idFileDelim, 0, checkComplement,
          false, true,variableStart, variableStop, barcodeErrors, barcodesBothEnds);
      readsNotMatchedToBarcode = unrecognizedBarcode.size();
    }
  }
  std::map<std::string, std::pair<std::string, std::string>> primers;
  int maxPrimerSize = 0;
  for (const auto& ptIter : primerTable.content_) {
    primers.insert(
        std::make_pair(ptIter[0], std::make_pair(ptIter[1], ptIter[2])));
    if ((int)ptIter[1].size() > maxPrimerSize) {
      maxPrimerSize = (int)ptIter[1].size();
    }
  }

  if(setUp.debug_){
  	std::cout << "primerName\tforwardPrimer\trevPrimer"<< std::endl;
  	for(const auto & p : primers){
  		std::cout << p.first << "\t" << p.second.first << "\t" << p.second.second << std::endl;
  	}
  	int32_t holdMidSize = 0;
  	table mids = seqUtil::readBarcodes(idFilename, idFileDelim, holdMidSize);
  	mids.outPutContentOrganized(std::cout);
  }
  // make some directories for outputs
  std::string unfilteredReadsDir =
      bib::files::makeDir(setUp.directoryName_, "unfilteredReads");
  std::string unfilteredByBarcodesDir =
      bib::files::makeDir(unfilteredReadsDir, "byBarcodes");
  std::string unfilteredByPrimersDir =
      bib::files::makeDir(unfilteredReadsDir, "byPrimers");
  std::string filteredOffDir =
      bib::files::makeDir(setUp.directoryName_, "filteredOff");
  std::string badDir = bib::files::makeDir(filteredOffDir, "bad");
  std::string unrecognizedPrimerDir =
      bib::files::makeDir(filteredOffDir, "unrecognizedPrimer");
  std::string contaminationDir = "";
  if (screenForPossibleContamination) {
    contaminationDir = bib::files::makeDir(filteredOffDir, "contamination");
  }
  if(setUp.debug_){
    for (const auto& unfilReadIter : reads) {
      reader.write(unfilReadIter.second,
                   unfilteredByBarcodesDir + unfilReadIter.first,
                   setUp.ioOptions_);
    }
  }


  std::ofstream profileLog;
  openTextFile(profileLog, setUp.directoryName_ + "extractionProfile.tab.txt",
               ".txt", false, false);

  profileLog << "name\ttotalReadsExtracted\tgoodReadsExtracted\tforGood\trevGood\ttotalBadReads\t"
                "badReverse\tcontainsNs\tlen<" << minLen << "\tlen>"
             << maxLength;
  if (checkingQCheck){
  	profileLog << "\tq" + estd::to_string(qualCheck) + "<" << qualCheckCutOff;
  }else{
  	profileLog << "\tbadQualityWindow";
  }
  if (flowFiltering) {
    profileLog << "\tbadFlowNoise";
  }
  if (screenForPossibleContamination) {
    profileLog << "\tcontamination";
  }
  profileLog << std::endl;

  int count = 0;
  std::ofstream renameKeyFile;
  if(rename){
  	openTextFile(renameKeyFile, setUp.directoryName_ + "renameKey.tab.txt", ".tab.txt", false, false);
  	renameKeyFile << "originalName\tnewName\n";
  }
  for (auto& rIter : reads) {
    ++count;
    std::cout << "Currently on " << rIter.first << " , number " << count
              << " of " << reads.size() << std::endl;

    // sort reads by forward primer
    std::map<std::string, std::vector<readObject>> readsByPrimer;
    std::vector<readObject> unrecognizedPrimer;
    if (!unknown) {
      bool anyPrimersWithDegenerativeBases = false;
      for (const auto& pIter : primers) {
        if (seqUtil::doesSequenceContainDegenerativeBase(pIter.second.first)) {
          anyPrimersWithDegenerativeBases = true;
          break;
        }
      }
      if (noForwardPrimer){
      	readsByPrimer[primers.begin()->first] = rIter.second;
      	if (trimForwardPrimer){
      		readObject forSeq(seqInfo("forwardPrimer", primers.begin()->second.first));
      	  runningParameters runParams;
      	  runParams.errors_.oneBaseIndel_ = 1;
      		readVecTrimmer::trimBeforeSequence(readsByPrimer[primers.begin()->first], forSeq, runParams,
      				alignObj, false, false, setUp.weightHomopolymers_, false);
      		for(auto & read : readsByPrimer[primers.begin()->first] ){
      			read.remove = false;
      		}
      	}
      } else {
        if (anyPrimersWithDegenerativeBases) {
          findForwardPrimerAndSortDegenerative(rIter.second, primers, readsByPrimer,
                                               unrecognizedPrimer, variableStop,
                                               forwardCoverage, forwardNumOfMismatches,
                                               forwardPercentageGaps);
        } else if (!multiplex && variableStart){
          findForwardPrimerAndSortDegenerative(rIter.second, primers, readsByPrimer,
                                               unrecognizedPrimer, variableStop,
                                               forwardCoverage, forwardNumOfMismatches,
                                               forwardPercentageGaps);
          if(checkComplement){
          	std::vector<readObject> unrecogComp;
          	for(auto & read : unrecognizedPrimer){
          		read.seqBase_.reverseComplementRead(true, false);
          	}
          	findForwardPrimerAndSortDegenerative(unrecognizedPrimer, primers, readsByPrimer,
          			unrecogComp, variableStop,
          	                                               forwardCoverage, forwardNumOfMismatches,
          	                                               forwardPercentageGaps);
          	unrecognizedPrimer = unrecogComp;
          }
        } else {
          findForwardPrimerAndSort(rIter.second, primers, readsByPrimer,
                                   unrecognizedPrimer, condensed);
        }
      }
    } else {
      readsByPrimer.insert({"", rIter.second});
    }
    reader.write(unrecognizedPrimer,
                 unrecognizedPrimerDir + rIter.first + "_unrecognizedPrimer",
                 setUp.ioOptions_);
    readsFailedForwardPrimer += (int)unrecognizedPrimer.size();
    readVec::getCountOfReadNameContaining(unrecognizedPrimer, "Comp",
                                          readsFailedForwardPrimerComp);
    for (auto& readsIter : readsByPrimer) {
    	if(setUp.debug_){
        reader.write(largeReads,
                     unfilteredByPrimersDir + readsIter.first + rIter.first,
                     setUp.ioOptions_);
    	}
      bool useClipped = false;
      uint32_t badReverse = 0;
      uint32_t badReadsMinLength = 0;
      uint32_t badReadsMaxLength = 0;
      uint32_t containsNs = 0;
      uint32_t badQualityWindow = 0;
      uint32_t flowNoise = 0;
      // std::cout<<"mark 1"<<std::endl;
      auto filteredReads =
              readVecSplitter::splitVectorOnRemove(readsIter.second);
      // std::cout<<"firstSize: "<<filteredReads.first.size()<<" secondSize:
      // "<<filteredReads.second.size()<<std::endl;
      // std::cout<<"mark 1"<<std::endl;
      // filteredReads.first.front().outPutSeq(std::cout);
      /*for(const auto & testRead : filteredReads.first){
        testRead.outPutSeq(std::cout);
      }*/
      std::vector<readObject> contamination;
      uint32_t contaminationNumber = 0;
      // std::cout<<"mark 1.5"<<std::endl;
      if (screenForPossibleContamination) {
        // filteredReads.first=readVectorSplitter::splitVectorOnNucleotideCompositionAdd(filteredReads.first,
        // contamination,mainCounter, contaminationNumber);
       /* filteredReads.first = readVecSplitter::splitVectorOnSeqSimilarityAdd(
            filteredReads.first, compareObject, contamination,
            contaminationNumber);*/
      	filteredReads.first = splitVectorOnKDistSimpleAdd(
      	            filteredReads.first, compareObject,contaminationKLen, kmerCutOff, contamination,
      	            contaminationNumber);
        reader.write(contamination, contaminationDir + readsIter.first +
                                        rIter.first + "_contamination",
                     setUp.ioOptions_);
        contaminationAmount += contamination.size();
        readVec::getCountOfReadNameContaining(contamination, "Comp",
                                              contaminationAmountComp);
      }
      // std::cout<<"mark 2"<<std::endl;
      // std::cout<<"firstSize: "<<filteredReads.first.size()<<" secondSize:
      // "<<filteredReads.second.size()<<std::endl;
      // filteredReads.first.front().outPutSeq(std::cout);
      filteredReads.first = readVecSplitter::splitVectorBellowLengthAdd(
          filteredReads.first, minLen, filteredReads.second, badReadsMinLength);
      // std::cout<<"mark 3"<<std::endl;
      // std::cout<<"firstSize: "<<filteredReads.first.size()<<" secondSize:
      // "<<filteredReads.second.size()<<std::endl;
      // filteredReads.first.front().outPutSeq(std::cout);
      if (findReversePrimer && !unknown) {
        readObject rPrimer(
            seqInfo("rPrimer", seqUtil::reverseComplement(
                                   primers[readsIter.first].second, "DNA")));

        runningParameters runParams;
        runParams.errors_.largeBaseIndel_ = .99;
        runParams.errors_.oneBaseIndel_ = 2;
        runParams.errors_.twoBaseIndel_ = 1;
        checkForReversePrimer(filteredReads.first, rPrimer, runParams, alignObj,
                              reversePrimerToLowerCase, true,
                              queryCoverageCutoff, percentIdentityCutoff);
      }
      // std::cout<<"mark 4"<<std::endl;
      // std::cout<<"firstSize: "<<filteredReads.first.size()<<" secondSize:
      // "<<filteredReads.second.size()<<std::endl;
      // filteredReads.first.front().outPutSeq(std::cout);
      filteredReads.first = readVecSplitter::splitVectorOnRemoveAdd(
          filteredReads.first, filteredReads.second, badReverse, "_badReverse",
          true);
      // std::cout<<"mark 5"<<std::endl;
      // std::cout<<"firstSize: "<<filteredReads.first.size()<<" secondSize:
      // "<<filteredReads.second.size()<<std::endl;
      // filteredReads.first.front().outPutSeq(std::cout);

      filteredReads.first = readVecSplitter::splitVectorBellowLengthAdd(
          filteredReads.first, minLen, filteredReads.second, badReadsMinLength);
      filteredReads.first = readVecSplitter::splitVectorAboveLengthAdd(
          filteredReads.first, maxLength, filteredReads.second,
          badReadsMaxLength);
      /*filteredReads.first = readVecSplitter::splitVectorSeqContainingAdd(
          filteredReads.first, "N", numberOfNs, filteredReads.second,
          containsNs);*/
      filteredReads.first = readVecSplitter::splitVectorOnNsPlusLengthAdd(
          filteredReads.first, minLen, filteredReads.second,
          containsNs);
      if (checkingQCheck) {
      	readVec::allSetQualCheck(filteredReads.first, qualCheck);
        filteredReads.first =
            readVecSplitter::splitVectorOnQualCheckAdd(
                filteredReads.first, qualCheck, qualCheckCutOff,
                filteredReads.second,
                badQualityWindow);
      } else if (qualWindowTrim) {
        filteredReads.first =
            readVecSplitter::splitVectorOnQualityWindowPlusLengthAdd(
                filteredReads.first, qualityWindowLength, qualityWindowStep,
                qualityWindowThres, useClipped, minLen, filteredReads.second,
                badQualityWindow);
      } else {
        filteredReads.first = readVecSplitter::splitVectorOnQualityWindowAdd(
            filteredReads.first, qualityWindowLength, qualityWindowStep,
            qualityWindowThres, useClipped, filteredReads.second,
            badQualityWindow);
      }

      // std::cout<<"mark 9"<<std::endl;
      if (flowFiltering) {
        filteredReads.first = readVecSplitter::splitVectorOnFlowNoiseAdd(
            filteredReads.first, maxFlowCutoff, filteredReads.second,
            flowNoise);
      }
      // std::cout<<"mark 10"<<std::endl;
      filteredReads.first = readVecSplitter::splitVectorBellowLengthAdd(
          filteredReads.first, minLen, filteredReads.second, badReadsMinLength);
      // write out the profile
      profileLog << readsIter.first << rIter.first << "\t"
                 << readsIter.second.size() << "\t"
                 << filteredReads.first.size() << "("
                 << 100 * filteredReads.first.size() / readsIter.second.size()
                 << "%)\t";
      int32_t currentCompCount = 0;
      readVec::getCountOfReadNameContaining(filteredReads.first,"_Comp", currentCompCount);
      profileLog << getPercentageString(filteredReads.first.size() - currentCompCount, filteredReads.first.size())
      					 << "\t" << getPercentageString(currentCompCount, filteredReads.first.size())
      					 << "\t";
      uint32_t badReadNumber = badReverse + badQualityWindow + badReadsMinLength +
                             badReadsMaxLength + containsNs;
      badReadNumber += flowNoise;
      profileLog << getPercentageString(badReadNumber, readsIter.second.size())
                 << "\t";
      // std::cout<<"mark 12"<<std::endl;
      if (badReadNumber == 0.00) {
        profileLog << badReverse << "\t";
        profileLog << containsNs << "\t";
        profileLog << badReadsMinLength << "\t";
        profileLog << badReadsMaxLength << "\t";
        profileLog << badQualityWindow << "\t";
        if (flowFiltering) {
          profileLog << flowNoise;
        }
        if (screenForPossibleContamination) {
          profileLog << getPercentageString(contaminationNumber,
                                            readsIter.second.size())
                     << std::endl;
        } else {
          profileLog << std::endl;
        }
      } else {
        profileLog << getPercentageString(badReverse, badReadNumber) << "\t";
        profileLog << getPercentageString(containsNs, badReadNumber) << "\t";
        profileLog << getPercentageString(badReadsMinLength, badReadNumber)
                   << "\t";
        profileLog << getPercentageString(badReadsMaxLength, badReadNumber)
                   << "\t";
        profileLog << getPercentageString(badQualityWindow, badReadNumber)
                   << "\t";
        if (flowFiltering) {
          profileLog << getPercentageString(flowNoise, badReadNumber);
        }
        if (screenForPossibleContamination) {
          profileLog << getPercentageString(contaminationNumber,
                                            readsIter.second.size())
                     << std::endl;
        } else {
          profileLog << std::endl;
        }
      }
      // std::cout<<"mark 13"<<std::endl;

      // write out the reads
      std::string currentMidGeneName = readsIter.first + rIter.first;
      std::string tempNameGood = setUp.directoryName_ + currentMidGeneName;
      if(rename){
				auto originalNames = readVec::getNames(filteredReads.first);
				renameReadNames(filteredReads.first, currentMidGeneName, false, true,
						true, "none");
				for (const auto & renamedReadsPos : iter::range(
						filteredReads.first.size())) {
					renameKeyFile << originalNames[renamedReadsPos] << "\t"
							<< filteredReads.first[renamedReadsPos].seqBase_.name_ << "\n";
				}
      }
      reader.write(filteredReads.first, tempNameGood,
                   setUp.ioOptions_);
      readsFinal += (int)filteredReads.first.size();
      readVec::getCountOfReadNameContaining(filteredReads.first, "Comp",
                                            readsFinalComp);
      if (pyroExtract) {
      	if(mothurExtract){
      		uint32_t mothurExtractFlowNum = 800;
        	if(maxFlowCutoff == 720){
        		mothurExtractFlowNum = 800;
        	}else if (maxFlowCutoff <= 400){
        		mothurExtractFlowNum = 400;
        	}
        	readObjectIOOptions options = setUp.ioOptions_;
        	options.outFilename_ = tempNameGood;
        	options.outFormat_ = "mothurData";
        	options.extra_ = mothurExtractFlowNum;
          reader.write(filteredReads.first, options);
      	}
      	readObjectIOOptions options = setUp.ioOptions_;
      	options.outFilename_ = tempNameGood;
      	options.outFormat_ = "pyroData";
      	options.extra_ = maxFlowCutoff;
        reader.write(filteredReads.first, options);
      }

      reader.write(filteredReads.second, badDir + currentMidGeneName + "_bad",
                   setUp.ioOptions_);
      failedQualityFiltering += (int)filteredReads.second.size();
      readVec::getCountOfReadNameContaining(filteredReads.second, "Comp",
                                            failedQualityFilteringComp);
      std::cout << "\t";
      setUp.logRunTime(std::cout);
    }
  }

  reader.write(unrecognizedBarcode, badDir + "unrecognizedBarcode",
               setUp.ioOptions_);
  reader.write(smallFragments,badDir +"smallFragments",setUp.ioOptions_);
  std::ofstream extractionStatsFile;
  openTextFile(extractionStatsFile,
               setUp.directoryName_ + "extractionStats.tab.txt", ".txt",
               setUp.ioOptions_.overWriteFile_,
               setUp.ioOptions_.exitOnFailureToWrite_);
  extractionStatsFile
      << "TotalReads\tReadsNotMatchedBarcodes\tSmallFragments(len<"
      << smallFragmentCutoff
      << ")\tfailedForwardPrimer\tfailedQualityFiltering\tused";
  if (screenForPossibleContamination) {
    extractionStatsFile << "\tcontamination";
  }
  extractionStatsFile << std::endl;
  extractionStatsFile
      << readsExtracted << "\t"
      << getPercentageString(readsNotMatchedToBarcode, readsExtracted) << "\t"
      << getPercentageString(smallFragmentCount, readsExtracted) << "\t"
      << getPercentageString(readsFailedForwardPrimer, readsExtracted) << "\t"
      << getPercentageString(failedQualityFiltering, readsExtracted) << "\t"
      << getPercentageString(readsFinal, readsExtracted);
  if (screenForPossibleContamination) {
    extractionStatsFile << "\t" << getPercentageString(contaminationAmount,
                                                       readsExtracted);
  }
  extractionStatsFile << std::endl;
  // get complement stats
  if (checkComplement) {
    std::string conplementStatsDir =
        bib::files::makeDir(setUp.directoryName_, "complementStats");
    // complement stats
    std::ofstream extractionStatsFileComp;
    openTextFile(extractionStatsFileComp,
                 conplementStatsDir + "complementExtractionStats.tab.txt",
                 ".txt", setUp.ioOptions_.overWriteFile_,
               setUp.ioOptions_.exitOnFailureToWrite_);
    int totalComplementCount =
         failedQualityFilteringComp +
        readsFailedForwardPrimerComp + readsFinalComp + contaminationAmountComp;

    extractionStatsFileComp
        << "TotalReads\tfailedForwardPrimer\tfailedQualityFiltering\tused";
    if (screenForPossibleContamination) {
      extractionStatsFileComp << "\tcontamination";
    }
    extractionStatsFileComp << std::endl;
    extractionStatsFileComp
        << totalComplementCount
        << "\t" << getPercentageString(readsFailedForwardPrimerComp,
                                       totalComplementCount) << "\t"
        << getPercentageString(failedQualityFilteringComp, totalComplementCount)
        << "\t" << getPercentageString(readsFinalComp, totalComplementCount);
    if (screenForPossibleContamination) {
      extractionStatsFileComp << "\t"
                              << getPercentageString(contaminationAmountComp,
                                                     totalComplementCount);
    }
    extractionStatsFileComp << std::endl;
    // Regular
    std::ofstream extractionStatsFileRegular;
    openTextFile(extractionStatsFileRegular,
                 conplementStatsDir + "regularExtractionStats.tab.txt", ".txt",
                 setUp.ioOptions_.overWriteFile_,
               setUp.ioOptions_.exitOnFailureToWrite_);

    int regularCount = smallFragmentCount + readsFailedForwardPrimer +
                       failedQualityFiltering + readsFinal +
                       contaminationAmount - totalComplementCount;
    extractionStatsFileRegular
        << "TotalReads\tSmallFragments(len<" << smallFragmentCutoff
        << ")\tfailedForwardPrimer\tfailedQualityFiltering\tused";
    if (screenForPossibleContamination) {
      extractionStatsFileRegular << "\tcontamination";
    }
    extractionStatsFileRegular << std::endl;
    extractionStatsFileRegular
        << regularCount << "\t"
        << getPercentageString(smallFragmentCount,
                               regularCount) << "\t"
        << getPercentageString(
               readsFailedForwardPrimer - readsFailedForwardPrimerComp,
               regularCount) << "\t"
        << getPercentageString(
               failedQualityFiltering - failedQualityFilteringComp,
               regularCount) << "\t"
        << getPercentageString(readsFinal - readsFinalComp, regularCount);
    if (screenForPossibleContamination) {
      extractionStatsFileRegular << "\t" << getPercentageString(
                                                contaminationAmount -
                                                    contaminationAmountComp,
                                                regularCount);
    }
    extractionStatsFileRegular << std::endl;
  }
  setUp.logRunTime(std::cout);
  return 0;
}



int SeekDeepRunner::extractor(MapStrStr inputCommands) {
	SeekDeepSetUp setUp(inputCommands);
	extractorPars pars;
  setUp.setUpExtractor(pars);

  uint32_t readsNotMatchedToBarcode = 0;
  if(pars.checkingQCheck){
    std::cout << "Quality Check: " << pars.qualCheck << "\n";
    std::cout << "Quality Check Cut Off: " << pars.qualCheckCutOff << "\n";
    std::cout << "Q" << pars.qualCheck << ">" << pars.qualCheckCutOff << "\n";
  }else{
    std::cout << "Quality Window Length: " << pars.qualityWindowLength << "\n";
    std::cout << "Quality Window Step: " << pars.qualityWindowStep << "\n";
    std::cout << "Quality Window Threshold: " << pars.qualityWindowThres << "\n";
  }
  // run log
  setUp.startARunLog(setUp.directoryName_);
  // parameter file
  setUp.writeParametersFile(setUp.directoryName_ + "parametersUsed.txt", false,
                            false);
  readObject compareObject =readObject(seqInfo("Compare", pars.compareSeq));
  table primerTable = seqUtil::readPrimers(pars.idFilename, pars.idFileDelim, false);
  int midSize = 10;
  table mids = seqUtil::readBarcodes(pars.idFilename, pars.idFileDelim, midSize);
  std::unique_ptr<MidDeterminator> determinator;
  if(pars.multiplex){
  	determinator = std::make_unique<MidDeterminator>(mids);
  }
  PrimerDeterminator pDetermine(primerTable);


  std::string tcagLower = "tcag";
  std::string tcagUpper = "TCAG";
  auto checkForTcag = [&tcagLower, &tcagUpper](const std::string & seq)->bool{
  	return std::equal(tcagLower.begin(), tcagLower.end(), seq.begin()) ||
  	  				std::equal(tcagUpper.begin(), tcagUpper.end(), seq.begin());
  };
  // make some directories for outputs
  std::string unfilteredReadsDir =
      bib::files::makeDir(setUp.directoryName_, "unfilteredReads");
  std::string unfilteredByBarcodesDir =
      bib::files::makeDir(unfilteredReadsDir, "byBarcodes");
  std::string unfilteredByPrimersDir =
      bib::files::makeDir(unfilteredReadsDir, "byPrimers");
  std::string filteredOffDir =
      bib::files::makeDir(setUp.directoryName_, "filteredOff");
  std::string badDir = bib::files::makeDir(filteredOffDir, "bad");
  std::string unrecognizedPrimerDir =
      bib::files::makeDir(filteredOffDir, "unrecognizedPrimer");
  std::string contaminationDir = "";
  if (pars.screenForPossibleContamination) {
    contaminationDir = bib::files::makeDir(filteredOffDir, "contamination");
  }
  readObject read;
  // read in reads and remove lower case bases indicating tech low quality like
  // tags and such
  std::cout << "Reading in reads:" << std::endl;
  readObjectIOOpt reader(setUp.ioOptions_);
  reader.openIn();
  auto smallOpts = setUp.ioOptions_;
  smallOpts.outFilename_ = badDir +"smallFragments";
  readObjectIOOpt smallFragMentOut(smallOpts);
  smallFragMentOut.openOut();
  uint32_t smallFragmentCount = 0;
  uint64_t maxReadSize = 0;
  uint32_t count = 0;

  std::map<std::string, uint32_t> midOutsPostions;
  std::vector<std::unique_ptr<readObjectIOOpt>> midOuts;
  std::map<std::string, std::pair<uint32_t, uint32_t>> counts;
	if(pars.multiplex && pars.barcodeErrors > 0){
		if(setUp.debug_){
			std::cout << "Allowing " << pars.barcodeErrors << " errors in barcode" << std::endl;
		}
	}
  if(pars.multiplex){
    for(const auto & mid : determinator->mids_){
      auto midOpts = setUp.ioOptions_;
      midOpts.outFilename_ = unfilteredByBarcodesDir + mid.first;
      if(setUp.debug_) std::cout << "Inserting: " << mid.first << std::endl;
      midOutsPostions[mid.first] = midOuts.size();
      midOuts.emplace_back(std::make_unique<readObjectIOOpt>(midOpts));
    }
  }else{
    auto midOpts = setUp.ioOptions_;
    midOpts.outFilename_ = unfilteredByBarcodesDir + "all";
    if(setUp.debug_) std::cout << "Inserting: " << "all" << std::endl;
    midOutsPostions["all"] = midOuts.size();
    midOuts.emplace_back(std::make_unique<readObjectIOOpt>(midOpts));
  }

  if(pars.multiplex){
    auto midOpts = setUp.ioOptions_;
    midOpts.outFilename_ = badDir + "unrecognizedBarcode";
    if(setUp.debug_) std::cout << "Inserting: " << "unrecognized" << std::endl;
    midOutsPostions["unrecognized"] = midOuts.size();
    midOuts.emplace_back(std::make_unique<readObjectIOOpt>(midOpts));
  }
  for(const auto & midOut : midOuts){
  	midOut->openOut();
  }
  std::cout << bib::bashCT::boldGreen("Extracting on MIDs") << std::endl;
  while(reader.readNextRead(read)){
  	++count;
  	if(count % 50 == 0){
    	std::cout << count << "\r";
    	std::cout.flush();
  	}

    if(pars.HMP){
    	subStrToUpper(read.seqBase_.seq_, 4, pars.primerLen);
    }
    if(pars.trimTcag){
  		if(checkForTcag(read.seqBase_.seq_)){
  			read.trimFront(4);
  		}
    }
    readVec::handelLowerCaseBases(read,setUp.ioOptions_.lowerCaseBases_);

    if(len(read) < pars.smallFragmentCutoff){
    	smallFragMentOut.write(read);
    	++smallFragmentCount;
    	continue;
    }
		readVec::getMaxLength(read, maxReadSize);

		midPos currentMid;
		if(pars.multiplex){
			if(pars.barcodeErrors > 0){
				currentMid = determinator->fullDetermine(read,
									pars.variableStart, pars.variableStop, pars.checkComplement, pars.barcodesBothEnds,
									pars.barcodeErrors);
			}else{
				currentMid = determinator->fullDetermine(read,
									pars.variableStart, pars.variableStop, pars.checkComplement, pars.barcodesBothEnds);
			}
		}else{
			currentMid = midPos("all", 0,0);
		}

    if(!currentMid){
    	++readsNotMatchedToBarcode;
    }

    if(read.seqBase_.name_.find("_Comp") != std::string::npos){
    	++counts[currentMid.midName_].second;
    }else{
    	++counts[currentMid.midName_].first;
    }
    /**@todo need to reorient the reads here before outputing if that's needed*/
    midOuts[midOutsPostions[currentMid.midName_]]->write(read);
  }
  std::cout << std::endl;
  //close mid outs;
  for(const auto & midOut : midOuts){
  	midOut->closeOut();
  }

  if(setUp.debug_){
    table midCounts{VecStr{"MidName", "For", "Rev"}};
    for(const auto & mCount : counts){
    	midCounts.content_.emplace_back(toVecStr(mCount.first, mCount.second.first, mCount.second.second));
    }
    midCounts.outPutContentOrganized(std::cout);
  }


  std::ofstream renameKeyFile;
  if(pars.rename){
  	openTextFile(renameKeyFile, setUp.directoryName_ + "renameKey.tab.txt", ".tab.txt", false, false);
  	renameKeyFile << "originalName\tnewName\n";
  }

  // create aligner for primer identification
  auto scoreMatrix = substituteMatrix::createDegenScoreMatrixNoNInRef();
  gapScoringParameters gapPars(setUp.gapInfo_);
  kmerMaps emptyMaps;
  int primaryQual = 20, secondayQual = 15;
  bool countEndGaps = true;
	aligner alignObj = aligner(maxReadSize, gapPars, scoreMatrix, emptyMaps,
			primaryQual, setUp.qualThresWindow_, secondayQual, countEndGaps);
	processAlnInfoInput(alignObj, setUp.alnInfoDirName_);
	if(pars.noForwardPrimer){
		if(primerTable.content_.size() > 1){
			std::cerr << "Error, if noForwardPrimer is turned on can only supply one gene name, curently have: " << primerTable.content_.size() << std::endl;
			std::cerr << bib::conToStr(primerTable.getColumn("geneName"), ",") << std::endl;
			exit(1);
		}
	}
  auto barcodeFiles = bib::files::listAllFiles(unfilteredByBarcodesDir,false, VecStr{});
	kmerInfo compareInfo;
	kmerInfo compareInfoRev;
	std::map<std::string, kmerInfo> compareInfos;
	std::map<std::string, kmerInfo> compareInfosRev;
  if (pars.screenForPossibleContamination) {
  	compareInfo = kmerInfo(compareObject.seqBase_.seq_, pars.contaminationKLen, false);
  	auto rev = compareObject;
  	rev.seqBase_.reverseComplementRead(true, true);
  	compareInfoRev = kmerInfo(rev.seqBase_.seq_, pars.contaminationKLen, false);
  }
  if(pars.multipleTargets && pars.screenForPossibleContamination){
  	readObjectIO readerCon;
  	readerCon.read(bib::files::getExtension(pars.compareSeqFilename), pars.compareSeqFilename);
  	for(const auto & read : readerCon.reads){
  		compareInfos[read.seqBase_.name_] = kmerInfo(read.seqBase_.seq_, pars.contaminationKLen, false);
    	auto rev = read;
    	rev.seqBase_.reverseComplementRead(true, true);
    	compareInfosRev[read.seqBase_.name_] = kmerInfo(rev.seqBase_.seq_, pars.contaminationKLen, false);
  	}
  }



  ExtractionStator stats = ExtractionStator(count, readsNotMatchedToBarcode, smallFragmentCount);
  std::map<std::string, uint32_t> goodCounts;
  for(const auto & f : barcodeFiles){

  	auto barcodeName = bib::files::getFileName(f.first.string());
  	if((counts[barcodeName].first + counts[barcodeName].second) == 0 && pars.multiplex){
  		//no reads extracted for barcode so skip filtering step
  		continue;
  	}
  	if(pars.multiplex){
  		std::cout << bib::bashCT::boldGreen("Filtering on barcode: "+ barcodeName) << std::endl;
  	}else{
  		std::cout << bib::bashCT::boldGreen("Filtering") << std::endl;
  	}

  	auto format = bib::files::getExtension(f.first.string());
  	auto barcodeOpts = setUp.ioOptions_;
  	barcodeOpts.firstName_ = f.first.string();
  	barcodeOpts.inFormat_ = format;
  	readObjectIOOpt barcodeIn(barcodeOpts);
  	barcodeIn.openIn();

  	std::vector<std::unique_ptr<readObjectIOOpt>> readOuts;
  	std::map<std::string, uint32_t> outPos;

  	{
  		auto unrecogPrimerOutOpts = setUp.ioOptions_;
  		unrecogPrimerOutOpts.outFilename_ = unrecognizedPrimerDir + barcodeName;
  		outPos["unrecognized"] = readOuts.size();
  		readOuts.emplace_back(std::make_unique<readObjectIOOpt>(unrecogPrimerOutOpts));
  	}
  	uint32_t barcodeCount = 1;
  	for(const auto & out : readOuts){
  		out->openOut();
  	}
  	bib::ProgressBar pbar(counts[barcodeName].first + counts[barcodeName].second);
  	pbar.progColors_ = pbar.RdYlGn_;
  	while(barcodeIn.readNextRead(read)){
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
				if(pars.multiplex){
					primerName = pDetermine.determineForwardPrimer(read, 0, alignObj, pars.fPrimerErrors,!pars.forwardPrimerToUpperCase,setUp.weightHomopolymers_);
				}else{
					uint32_t start = 0;
					if(pars.variableStart){
						start = pars.variableStop;
					}
					primerName = pDetermine.determineForwardPrimer(read, start, alignObj, pars.fPrimerErrors,!pars.forwardPrimerToUpperCase,setUp.weightHomopolymers_);
					if(primerName == "unrecognized" && pars.checkComplement){
						primerName = pDetermine.determineWithReversePrimer(read, start, alignObj, pars.fPrimerErrors,!pars.forwardPrimerToUpperCase,setUp.weightHomopolymers_);
						if(read.seqBase_.on_){
							foundInReverse = true;
						}
					}
				}

	    	if(!read.seqBase_.on_){
	    		stats.increaseFailedForward(barcodeName, read.seqBase_.name_);
	    		readOuts[outPos["unrecognized"]]->write(read);
	    		continue;
	    	}
			}
			std::string fullname = primerName;
			if(pars.multiplex){
				fullname += barcodeName;
			}else if (pars.sampleName != ""){
				fullname += pars.sampleName;
			}
	  	{
	  		auto badDirOutOpts = setUp.ioOptions_;
	  		badDirOutOpts.outFilename_ = badDir + fullname;
	  		if(!bib::in(fullname + "bad", outPos)){
		  		outPos[fullname + "bad"] = readOuts.size();
		  		readOuts.emplace_back(std::make_unique<readObjectIOOpt>(badDirOutOpts));
		  		readOuts.back()->openOut();
	  		}
	  	}
	  	{
	  		auto goodDirOutOpts = setUp.ioOptions_;
	  		goodDirOutOpts.outFilename_ = setUp.directoryName_ + fullname;
	  		if(!bib::in(fullname + "good", outPos)){
		  		outPos[fullname + "good"] = readOuts.size();
		  		readOuts.emplace_back(std::make_unique<readObjectIOOpt>(goodDirOutOpts));
		  		readOuts.back()->openOut();
	  		}
	  	}
			if (pars.screenForPossibleContamination) {
				auto contamOutOpts = setUp.ioOptions_;
				contamOutOpts.outFilename_ = contaminationDir + fullname;
	  		if(!bib::in(fullname + "contamination", outPos)){
		  		outPos[fullname + "contamination"] = readOuts.size();
		  		readOuts.emplace_back(std::make_unique<readObjectIOOpt>(contamOutOpts));
		  		readOuts.back()->openOut();
	  		}
			}
  		//look for possible contamination
      if (pars.screenForPossibleContamination) {
      	kmerInfo compareInfo(compareObject.seqBase_.seq_, pars.contaminationKLen, false);
      	if(pars.multipleTargets){
      		if(compareInfos.find(primerName) != compareInfos.end()){
      			if(foundInReverse){
      				readChecker::checkReadOnKmerComp(read.seqBase_, compareInfosRev[primerName],pars.contaminationKLen, pars.kmerCutOff, true);
      			}else{
      				readChecker::checkReadOnKmerComp(read.seqBase_, compareInfos[primerName],pars.contaminationKLen, pars.kmerCutOff, true);
      			}

          	if(!read.seqBase_.on_){
          		stats.increaseCounts(fullname, read.seqBase_.name_, ExtractionStator::extractCase::CONTAMINATION);
          		readOuts[outPos[fullname + "contamination"]]->write(read);
          		continue;
          	}
      		}else{
      			std::cerr << "Error in screening for contamination, multiple targets turned on but no contamination found for " << primerName << std::endl;
      			std::cerr << "Options are: " << vectorToString(getVectorOfMapKeys(compareInfos), ",") << std::endl;
      		}
      	}else{
      		if(foundInReverse){
      			readChecker::checkReadOnKmerComp(read.seqBase_, compareInfoRev,pars.contaminationKLen, pars.kmerCutOff, true);
      		}else{
      			readChecker::checkReadOnKmerComp(read.seqBase_, compareInfo,pars.contaminationKLen, pars.kmerCutOff, true);
      		}
        	if(!read.seqBase_.on_){
        		stats.increaseCounts(fullname, read.seqBase_.name_, ExtractionStator::extractCase::CONTAMINATION);
        		readOuts[outPos[fullname + "contamination"]]->write(read);
        		continue;
        	}
      	}

      }

      //min len
    	readChecker::checkReadLenAbove(read.seqBase_,pars.minLen, true);
    	if(!read.seqBase_.on_){
    		stats.increaseCounts(fullname, read.seqBase_.name_, ExtractionStator::extractCase::MINLENBAD);
    		readOuts[outPos[fullname + "bad"]]->write(read);
    		continue;
    	}


  		//reverse
      if (!pars.noReversePrimer) {
      	if(foundInReverse){
      		pDetermine.checkForForwardPrimerInRev(read, primerName, alignObj, pars.rPrimerErrors,!pars.reversePrimerToUpperCase, setUp.weightHomopolymers_,pars.variableStop,false);
      	}else{
      		pDetermine.checkForReversePrimer(read, primerName, alignObj, pars.rPrimerErrors,!pars.reversePrimerToUpperCase, setUp.weightHomopolymers_, pars.variableStop,false);
      	}

      	if(!read.seqBase_.on_){
      		stats.increaseCounts(fullname, read.seqBase_.name_, ExtractionStator::extractCase::BADREVERSE);
      		read.seqBase_.name_.append("_badReverse");
      		readOuts[outPos[fullname + "bad"]]->write(read);
      		continue;
      	}
      }
      //if found in the reverse direction need to re-orient now
      if(foundInReverse){
      	read.seqBase_.reverseComplementRead(true,true);
      }

  		//contains n
      readChecker::checkReadOnSeqContaining(read.seqBase_, "N", pars.numberOfNs, true);
    	if(!read.seqBase_.on_){
    		stats.increaseCounts(fullname, read.seqBase_.name_, ExtractionStator::extractCase::CONTAINSNS);
    		readOuts[outPos[fullname + "bad"]]->write(read);
    		continue;
    	}

      //max len
    	readChecker::checkReadLenBellow(read.seqBase_, pars.maxLength, true);
    	if(!read.seqBase_.on_){
    		stats.increaseCounts(fullname, read.seqBase_.name_, ExtractionStator::extractCase::MAXLENBAD);
    		readOuts[outPos[fullname + "bad"]]->write(read);
    		continue;
    	}
      //quality
    	if(pars.checkingQCheck){
      	readChecker::checkReadQualCheck(read.seqBase_, pars.qualCheck, pars.qualCheckCutOff, true);
    	}else{
    		if(pars.qualWindowTrim){
    			readChecker::checkReadOnQualityWindowTrim(read.seqBase_, pars.qualityWindowLength, pars.qualityWindowStep, pars.qualityWindowThres, pars.minLen, true);
    		}else{
    			readChecker::checkReadOnQualityWindow(read.seqBase_, pars.qualityWindowLength, pars.qualityWindowStep, pars.qualityWindowThres, true);
    		}
    	}
    	if(!read.seqBase_.on_){
    		stats.increaseCounts(fullname, read.seqBase_.name_, ExtractionStator::extractCase::QUALITYFAILED);
    		readOuts[outPos[fullname + "bad"]]->write(read);
    		continue;
    	}
    	if(read.seqBase_.on_){
    		stats.increaseCounts(fullname, read.seqBase_.name_, ExtractionStator::extractCase::GOOD);
        if(pars.rename){
        	std::string oldName = replaceString(read.seqBase_.name_, "_Comp", "");
        	read.seqBase_.name_ = fullname + "." + leftPadNumStr(goodCounts[fullname],counts[barcodeName].first + counts[barcodeName].second);
        	if(containsSubString(oldName, "_Comp")){
        		read.seqBase_.name_.append("_Comp");
        	}
        	renameKeyFile << oldName << "\t"
        	  							<< read.seqBase_.name_ << "\n";

        }
    		readOuts[outPos[fullname + "good"]]->write(read);

    		++goodCounts[fullname];
    	}
  	}
  	std::cout << std::endl;
  }
  std::ofstream profileLog;
  openTextFile(profileLog, setUp.directoryName_ + "extractionProfile.tab.txt",
               ".txt", false, false);

  profileLog << "name\ttotalReadsExtracted\tgoodReadsExtracted\tforGood\trevGood\ttotalBadReads\t"
                "badReverse\tcontainsNs\tlen<" << pars.minLen << "\tlen>"
             << pars.maxLength;
  if (pars.checkingQCheck){
  	profileLog << "\tq" + estd::to_string(pars.qualCheck) + "<" << pars.qualCheckCutOff;
  }else{
  	profileLog << "\tbadQualityWindow";
  }
  profileLog << "\tcontamination\n";
  stats.outStatsPerName(profileLog, "\t");
  std::ofstream extractionStatsFile;
  openTextFile(extractionStatsFile,
               setUp.directoryName_ + "extractionStats.tab.txt", ".txt",
               setUp.ioOptions_.overWriteFile_,
               setUp.ioOptions_.exitOnFailureToWrite_);
  extractionStatsFile
      << "TotalReads\tReadsNotMatchedBarcodes\tSmallFragments(len<"
      << pars.smallFragmentCutoff
      << ")\tfailedForwardPrimer\tfailedQualityFiltering\tused";
  extractionStatsFile << "\tcontamination";
  extractionStatsFile << std::endl;
  stats.outTotalStats(extractionStatsFile, "\t");
  std::ofstream failedForwadFile;
  openTextFile(failedForwadFile, setUp.directoryName_ + "failedForwad.tab.txt",
               ".txt", false, false);
  failedForwadFile << "MidName\ttotalFailed\tfailedInFor\tfailedInRev" << std::endl;
  stats.outFailedForwardStats(failedForwadFile, "\t");

  if(!setUp.debug_){
  	bib::files::rmDirForce(unfilteredReadsDir);
  }
  if(setUp.writingOutAlnInfo_){
  	setUp.rLog_ << "Number of alignments done" << "\n";
  	alignObj.alnHolder_.write(setUp.outAlnInfoDirName_, setUp.verbose_);
  }
  setUp.logRunTime(std::cout);
  return 0;
}



int SeekDeepRunner::qluster(MapStrStr inputCommands) {
  SeekDeepSetUp setUp(inputCommands);
  // parameters
  clusterDownPars pars;

  setUp.setUpClusterDown(pars);
  Json::Value metaData;
  auto analysisDirPath = bib::files::bfs::canonical(setUp.directoryName_);
  metaData["analysisDirPath"] = analysisDirPath.string();
  // print out the parameters read in
  if(pars.useNucComp && setUp.verbose_){
  	std::cout << "Nucleotide Composition Binning cut offs" << std::endl;
  	printVector(pars.diffCutOffVec, ", ", std::cout);
  }
  // read in the sequences
  readObjectIO reader;
  reader.read(setUp.ioOptions_);
  auto splitOnSize = readVecSplitter::splitVectorBellowLength(reader.reads, pars.smallReadSize);
  reader.reads = splitOnSize.first;
  if(!splitOnSize.second.empty()){
  	readObjectIO::write(splitOnSize.second, setUp.directoryName_ + "smallReads", setUp.ioOptions_);
  }
  if (pars.removeLowQualBases) {
    readVec::allRemoveLowQualityBases(reader.reads, pars.lowQualityCutOff);
  }
  if (setUp.adjustHomopolyerRuns_) {
    readVec::allAdjustHomopolymerRunsQualities(reader.reads);
  }
  bool containsCompReads = false;
  int compCount = 0;
  readVec::getCountOfReadNameContaining(reader.reads, "_Comp", compCount);
  if (compCount > 0) {
    containsCompReads = true;
  }

  readVecSorter::sortReadVector(reader.reads, pars.sortBy);
  // get the count of reads read in and the max length so far
  int counter = readVec::getTotalReadCount(reader.reads);
  uint64_t maxSize = 0;
  readVec::getMaxLength(reader.reads, maxSize);
  std::vector<readObject> refSequences;
  if (setUp.refFilename_ != "") {
    refSequences = readObjectIO::getReferenceSeq(
        setUp.refFilename_, setUp.refFormat_, setUp.refProcessed_, maxSize);
  }
  maxSize = maxSize * 2;
  // calculate the runCutoff if necessary
  processRunCutoff(setUp.runCutoff_, setUp.runCutOffString_,
                   readVec::getTotalReadCount(reader.reads));
  if(setUp.verbose_ && !pars.onPerId ){
  	std::cout << "Kmer Low Frequency Error Cut off Is: " << setUp.runCutoff_ << std::endl;
  }
  // make the runLog, this is what is seen on the terminal screen at run time
  setUp.startARunLog(setUp.directoryName_);
  // parameter file
  setUp.writeParametersFile(setUp.directoryName_ + "parametersUsed.txt", false,true);
  // read in the paramteres from the parameters file
  setUp.rLog_ << "Parameters used" << "\n";
  table parTab(pars.parameters, ":");
  parTab.outPutContents(setUp.rLog_.runLogFile_, ":");
  if(setUp.verbose_){
    std::cout << "Reading clusters from " << setUp.ioOptions_.firstName_ << " "
              << setUp.ioOptions_.secondName_ << std::endl;
    std::cout << "Read in " << counter << " reads" << std::endl;
  }
  setUp.rLog_ << "Reading clusters from " << setUp.ioOptions_.firstName_
                  << " " << setUp.ioOptions_.secondName_ << "\n";
  setUp.rLog_ << "Read in " << counter << " reads" << "\n";
  // create cluster vector
  std::vector<identicalCluster> identicalClusters;
  std::vector<cluster> clusters;
  if (setUp.ioOptions_.processed_) {
    clusters = baseCluster::convertVectorToClusterVector<cluster>(reader.reads);
  } else {
  	identicalClusters = clusterCollapser::collapseIdenticalReads(
        reader.reads, pars.qualRep, setUp.ioOptions_.lowerCaseBases_);
    for(const auto & read : identicalClusters){
    	clusters.push_back(cluster(readObject(read.seqBase_)));
    }
  	//clusters = baseCluster::convertVectorToClusterVector<cluster>(identicalClusters);
  }

  if(setUp.verbose_){
  	std::cout << "Unique clusters numbers: " << clusters.size() << std::endl;
  }
  setUp.rLog_ << "Unique clusters numbers: " << clusters.size() << "\n";
  std::sort(clusters.begin(), clusters.end());
  //readVecSorter::sortReadVector(clusters, sortBy);

  kmerMaps kMaps = kmerCalculator::indexKmerMpas(
      clusters, setUp.kLength_, setUp.runCutoff_,
      setUp.expandKmerPos_, setUp.expandKmerSize_);
  // create aligner class object
  aligner alignerObj = aligner(maxSize, setUp.gapInfo_, setUp.scoring_, kMaps, setUp.primaryQual_,
                               setUp.secondaryQual_, setUp.qualThresWindow_,
                               setUp.countEndGaps_);
  if(setUp.verbose_ && !pars.onPerId){
  	std::cout << bib::bashCT::bold << "Primary Qual: " << alignerObj.primaryQual_ << std::endl;
  	std::cout << "Secondary Qual: " << alignerObj.secondaryQual_ << bib::bashCT::reset << std::endl;
  }
  if(setUp.debug_){
  	std::ofstream scoreArrayFile;
  	openTextFile(scoreArrayFile, "scoreArrayTest.tab.txt", ".txt", true, false);
  	std::ofstream scoreMapFile;
  	openTextFile(scoreMapFile, "scoreMapFile.tab.txt", ".txt", true, false);
  	for(const auto & row : iter::range(len(alignerObj.parts_.scoring_.mat_))){
  		std::vector<int32_t> currentRow;
			for(const auto & col : iter::range(len(alignerObj.parts_.scoring_.mat_[row]))){
				currentRow.emplace_back(alignerObj.parts_.scoring_.mat_[row][col]);
			}
			printVector(currentRow, "\t", scoreArrayFile);
  	}
  	alignerObj.parts_.gapScores_.printDescription(std::cout, true);
  }

  processAlnInfoInput(alignerObj, setUp.alnInfoDirName_);


  collapser collapserObj = collapser(
  		!setUp.largestFirst_, setUp.bestMatchCheck_, setUp.local_, true,
      setUp.kmersByPosition_, setUp.runCutoff_,
      setUp.kLength_, setUp.verbose_, !setUp.largestFirst_, true,
      setUp.weightHomopolymers_, setUp.skipOnLetterCounterDifference_,
			setUp.fractionDifferenceCutOff_, setUp.adjustHomopolyerRuns_);
  collapserObj.opts_.lowQualityBaseTrim_ = pars.lowQualityCutOff;
  collapserObj.opts_.removeLowQualityBases_ = pars.removeLowQualBases;
  collapserObj.opts_.debug_ = setUp.debug_;
  collapserObj.opts_.noAlign_ = pars.noAlign_;

  uint32_t singletonNum = 0;
  std::vector<cluster> singletons;
  if (!pars.startWithSingles) {
    for (auto& clus : clusters) {
      if (clus.seqBase_.cnt_ <= 1.01) {
        clus.remove = true;
      }
    }
    clusters = readVecSplitter::splitVectorOnRemoveAdd(
        clusters, singletons, singletonNum, "none", false);
    for (auto& clus : singletons) {
      clus.remove = false;
    }
  }

  //run clustering
	collapserObj.runFullClustering(clusters, pars.onPerId, pars.iteratorMap,
			pars.useNucComp, pars.useMinLenNucComp, pars.findBest, pars.diffCutOffVec,
			pars.useKmerBinning, pars.kCompareLen, pars.kmerCutOff, alignerObj, setUp,
			pars.snapShots, "firstSnaps");
	//run again with singlets if needed
  if (!pars.startWithSingles && !pars.leaveOutSinglets) {
    addOtherVec(clusters, singletons);
  	collapserObj.runFullClustering(clusters, pars.onPerId, pars.iteratorMap,
  			pars.useNucComp, pars.useMinLenNucComp, pars.findBest, pars.diffCutOffVec,
  			pars.useKmerBinning, pars.kCompareLen, pars.kmerCutOff, alignerObj, setUp,
  			pars.snapShots, "secondSnaps");
  }

  //remove reads if they are made up of reads only in one direction
  if(containsCompReads && pars.useCompPerCutOff){
    for ( auto& clus : clusters) {
      std::vector<readObject> collection =
          clus.getAllBeginingClusters(identicalClusters);
      readVec::allSetFractionByTotalCount(collection);
      if (containsCompReads) {
        int currentCompAmount = 0;
        readVec::getCountOfReadNameContaining(collection, "_Comp",
                                              currentCompAmount);
        double currentCompPer = currentCompAmount/clus.seqBase_.cnt_;
        if(currentCompPer > pars.compPerCutOff && clus.seqBase_.cnt_ != 1){
        	clus.remove = true;
        }
      }
    }
    auto splitClus = readVecSplitter::splitVectorOnRemove(clusters);
    if(!splitClus.second.empty()){
    	clusters = splitClus.first;
    	clusterVec::allSetFractionClusters(clusters);
    	reader.write(splitClus.second, setUp.directoryName_ + "allOneDirection",setUp.ioOptions_);
    }
  }

  readVecSorter::sortReadVector(clusters, "totalCount");
  std::string seqName = bib::files::getFileName(setUp.ioOptions_.firstName_);
  readVecSorter::sort(clusters);
  renameReadNames(clusters, seqName, true, false, false);
  //std::vector<cluster> rejectedClusters;
  std::unordered_map<double, double> bestLikelihood;
  bib::randomGenerator gen;
  simulation::mismatchProfile eProfile({'A', 'C', 'G', 'T'});
  if (pars.sim) {

    // create R session
    njhRInside::OwnRInside ownSession;
    // create error profile to get the error profile of the current clusters
    for (const auto& clus : clusters) {
      clus.updateErrorProfile(eProfile, alignerObj, setUp.local_);
    }

    eProfile.setFractions();
    if (setUp.verbose_) {
      // Print out info about the errors seen
      std::cout << "eProfile" << std::endl;
      std::cout << "counts" << std::endl;
      eProfile.quickPrintCounts(std::cout);
      std::cout << "profile" << std::endl;
      eProfile.quickPrintProfile(std::cout);
      std::cout << "probs" << std::endl;
      eProfile.quickPrintProbs(std::cout);
    }

    // various directories to store the simulation results
    std::string simulationDir =
        bib::files::makeDir(setUp.directoryName_, "simulation");
    std::string simCorrectedDir =
        bib::files::makeDir(setUp.directoryName_, "simClusters");
    std::string qualErrorDir =
        bib::files::makeDir(setUp.directoryName_, "qualErrors");
    // not sure this is needed anymore
    readVec::allSetFractionByTotalCount(clusters);

    // counts of all qualities simply per base
    std::map<double, uint32_t> baseCounts;
    std::map<double, uint32_t> misMatchBaseCounts;
    bib::for_each(reader.reads, [&](const readObject& read) {
      read.updateBaseQualCounts(baseCounts);
    });
    for (const auto& clus : clusters) {
      // collect all the original sequences
      std::vector<readObject> collection =
          clus.getAllBeginingClusters(identicalClusters);
      for (const auto& subRead : collection) {
        alignerObj.alignVec(clus, subRead, setUp.local_);
        alignerObj.profileAlignment(clus, subRead, setUp.kLength_,
                                    setUp.kmersByPosition_, false, false, false,
                                    setUp.weightHomopolymers_);
        for (const auto& mis : alignerObj.mismatches_) {
          subRead.updateBaseQualCounts(misMatchBaseCounts,
                                       mis.second.seqBasePos);
        }
      }
    }
    // calc the total error observed and what the overall error rate per base is
    uint32_t totalBaseCount = 0;
    uint32_t totalErrorCount = 0;
    double basalErrorRate = 0;
    for (const auto& qCount : baseCounts) {
      totalBaseCount += qCount.second;
    }
    for (const auto& mCount : misMatchBaseCounts) {
      totalErrorCount += mCount.second;
    }
    basalErrorRate = (double)totalErrorCount / totalBaseCount;
    std::cout << "totalCount: " << totalBaseCount << std::endl;
    std::cout << "totalErrorCount: " << totalErrorCount << std::endl;
    std::cout << "basalError: " << basalErrorRate << std::endl;
    // best likelihood predictions

    std::unordered_map<double, double> bestErrorRate;

    // std::unordered_map<std::string, std::unordered_map<std::string,
    // std::vector<double>>> predictions;
    // get the weights and then use R to predict error from the quals
    std::unordered_map<std::string, std::vector<double>> countsForModelBase =
        seqUtil::getCountsForSpecificModel(baseCounts, misMatchBaseCounts);
    std::unordered_map<std::string, std::vector<double>> basePredictions =
        njhRInside::getModelQualError(countsForModelBase, 0, 50, 1, "logit");

    // store the best area under the curve
    double bestArea = 0;
    std::string bestAreaPar = "";
    std::string bestPar = "";
    std::unordered_map<double, double> parLikelihood =
        njhRInside::createLikelihoodMap(basePredictions);
    auto errorRatesBase =
        seqUtil::getTrueErrorRateSpecific(baseCounts, misMatchBaseCounts);
    double baseSum = 0;
    for (const auto& subCount : baseCounts) {
      baseSum += subCount.second *
                 std::abs(basalErrorRate - parLikelihood[subCount.first]);
    }
    double difference = parLikelihood[10.00] - parLikelihood[40.00];
    if (baseSum > bestArea && difference > 0) {
      bestArea = baseSum;
      bestLikelihood = parLikelihood;
      bestErrorRate = errorRatesBase;
      bestAreaPar = "base";
      bestPar = "base";
    }
    std::cout << "baseArea: " << baseSum << std::endl;
    // counts for each window size
    uint32_t windowStart = 2;
    uint32_t windowStop = 9;
    std::vector<uint32_t> qualWindows(windowStop - windowStart + 1);
    std::iota(qualWindows.begin(), qualWindows.end(), windowStart);

    std::unordered_map<
        uint32_t, std::map<std::string, std::map<double, uint32_t>>> countsStep;
    std::unordered_map<uint32_t,
                       std::map<std::string, std::map<double, uint32_t>>>
        misMatchCountsStep;
    // get the total counts for each window
    {
      bib::scopedStopWatch qualCounting("qualCounting");
      for (auto i : qualWindows) {
        bib::for_each(reader.reads, [&countsStep,&i](const readObject& read) {
          read.updateQualWindowCounts(countsStep[i], i);
        });
      }
      // get the mismatch counts for each window
      for (const auto& clus : clusters) {
        std::vector<readObject> collection =
            clus.getAllBeginingClusters(identicalClusters);
        for (const auto& subRead : collection) {
          alignerObj.alignVec(clus, subRead, setUp.local_);
          alignerObj.profileAlignment(clus, subRead, setUp.kLength_,
                                      setUp.kmersByPosition_, false, false,
                                      false, setUp.weightHomopolymers_);
          for (const auto& mis : alignerObj.mismatches_) {
            for (auto i : qualWindows) {
              subRead.updateQualWindowCounts(mis.second.seqBasePos,
                                             misMatchCountsStep[i], i);
            }
          }
        }
      }
    }

    for (auto i : qualWindows) {
      std::unordered_map<std::string,
                         std::unordered_map<std::string, std::vector<double>>>
          predictionsStep;
      std::map<std::string,
               std::unordered_map<std::string, std::vector<double>>>
          countsForModelStep =
              seqUtil::getCountsForModel(countsStep[i], misMatchCountsStep[i]);
      auto errorRates =
          seqUtil::getTrueErrorRate(countsStep[i], misMatchCountsStep[i]);
      for (const auto& cm : countsForModelStep) {
        if (cm.first == "base" || cm.first == "min") {
          predictionsStep[cm.first] =
              njhRInside::getModelQualError(cm.second, 0, 50, 1, "logit");
        } else if (cm.first == "mean" || cm.first == "median") {
          predictionsStep[cm.first] =
              njhRInside::getModelQualError(cm.second, 0, 50, 0.01, "logit");
        }
      }
      // std::cout << "made it here" << std::endl;
      std::cout << "window: " << i << std::endl;
      std::unordered_map<std::string, double> areaInCurveStep;
      for (const auto& qCount : countsStep[i]) {
        // std::cout << qCount.first << std::endl;
        if (qCount.first != "base" && qCount.first != "min" &&
            qCount.first != "mean" && qCount.first != "median") {
          continue;
        }
        std::unordered_map<double, double> parLikelihood =
            njhRInside::createLikelihoodMap(predictionsStep[qCount.first]);
        double sum = 0;
        for (const auto& subCount : qCount.second) {
          sum += subCount.second *
                 std::abs(basalErrorRate - parLikelihood[subCount.first]);
        }
        double difference = parLikelihood[10.00] - parLikelihood[40.00];
        if (sum > bestArea && difference > 0) {
          bestArea = sum;
          bestLikelihood = parLikelihood;
          bestErrorRate = errorRates[qCount.first];
          bestAreaPar = qCount.first + "_qw_" + estd::to_string(i);
          bestPar = qCount.first;
        }
        areaInCurveStep[qCount.first] = sum;
      }
      for (const auto& aic : areaInCurveStep) {
        std::cout << aic.first << ": " << aic.second << std::endl;
      }
    }
    // plot the best likelihood vs actual rates
    ownSession.installLib("ggplot2");
    // put in the model error rate
    Rcpp::DataFrame bestLikeFrame = Rcpp::DataFrame::create(
        Rcpp::Named("qual") = getVectorOfMapKeys(bestLikelihood),
        Rcpp::Named("errorRate") = getVectorOfMapValues(bestLikelihood),
        Rcpp::Named("par") = repeatVector(
            VecStr{"model"}, std::vector<uint64_t>{bestLikelihood.size()}));
    auto& r = ownSession.get();
    r["bestPred"] = bestLikeFrame;
    // put in the observed error rate;
    Rcpp::DataFrame bestErrorFrame = Rcpp::DataFrame::create(
        Rcpp::Named("qual") = getVectorOfMapKeys(bestErrorRate),
        Rcpp::Named("errorRate") = getVectorOfMapValues(bestErrorRate),
        Rcpp::Named("par") = repeatVector(
            VecStr{"observed"}, std::vector<uint64_t>{bestErrorRate.size()}));
    r["bestError"] = bestErrorFrame;
    // put in the expected for the quality
    std::unordered_map<double, double> perQualError;
    for (auto q : iter::range(10, 40)) {
      perQualError[q] = std::pow(10, ((double)-q / 10.0));
    }
    // table perQualTable = table(perQualError, {"qual", "ErrorRate"});
    // perQualTable.outPutContentOrganized(std::cout);
    Rcpp::DataFrame expectedFrame = Rcpp::DataFrame::create(
        Rcpp::Named("qual") = getVectorOfMapKeys(perQualError),
        Rcpp::Named("errorRate") = getVectorOfMapValues(perQualError),
        Rcpp::Named("par") = repeatVector(
            VecStr{"expected"}, std::vector<uint64_t>{perQualError.size()}));
    Rcpp::DataFrame qualCountsFrame = Rcpp::DataFrame::create(
        Rcpp::Named("qual") = getVectorOfMapKeys(baseCounts),
        Rcpp::Named("freq") = getVectorOfMapValues(baseCounts));
    r["qualCounts"] = qualCountsFrame;
    r["expected"] = expectedFrame;
    // put in basal error rate
    r["basalError"] = basalErrorRate;
    std::string combineCmd = "all = rbind(bestPred, bestError, expected);";
    std::string combinePredErrorCmd = "two = rbind(bestPred, bestError);";
    r.parseEvalQ(combineCmd);
    r.parseEvalQ(combinePredErrorCmd);
    // R.parseEvalQ("print(all);");
    std::string plotCmd =
        "pdf(\"" + setUp.directoryName_ + seqName +
        "_qualErrorRate.pdf\"); print(ggplot(all, aes(x = qual, y = errorRate, "
        "group = par, color = par)) + geom_line() + ggtitle(\"" +
        seqName + "_" + bestAreaPar +
        "\") +geom_abline(intercept = basalError , slope = 0));";
    std::string secondPlotCmd =
        "print(ggplot(two, aes(x = qual, y = errorRate, group = par, color = "
        "par)) + geom_line() + ggtitle(\"" +
        seqName + "_" + bestAreaPar +
        "\") +geom_abline(intercept = basalError , slope = 0));";
    std::string thirdPlotCmd =
        "print(ggplot(qualCounts, aes(x = qual, y = freq) ) + geom_bar(stat = "
        "\"identity\") );";

    r.parseEvalQ(plotCmd);
    r.parseEvalQ(secondPlotCmd);
    r.parseEvalQ(thirdPlotCmd);
    std::cout << "bestAreaPar: " << bestAreaPar << std::endl;
    std::cout << "bestArea: " << bestArea << std::endl;
    // simulate and remove on FDR cut
    uint32_t clusCount = 1;
    for (auto& clus : clusters) {
    	std::cout << "on " << clusCount << " of " << len(clusters) << std::endl;
    	++clusCount;
      // if cluster is only one read no need to simulate
      if (clus.seqBase_.cnt_ < 2) {
        continue;
      }
      bib::scopedStopWatch t(clus.seqBase_.name_ + "_sim", true);
      // simulate the cluster using the qualities of the original sequences
      clus.simOnQual(identicalClusters, alignerObj, pars.runTimes, bestLikelihood,
                     eProfile, gen);
      // remove clusters on fdr cut off
      //clus.removeClustersOnFDR(alignerObj, pars.fdrCutOff, rejectedClusters);
    }

  }

  if (pars.sim) {
    std::ofstream pValueFile;
    openTextFile(pValueFile, setUp.directoryName_ + "pValues.tab.txt",
                 ".txt", true, false);
    pValueFile << "read1\tread2\tpValue" << std::endl;
    bib::scopedStopWatch pTimer("comparing p value");
    for(const auto & pos : iter::range(len(clusters))){
    	double bestPValue = 0;
    	double bestScore = 0;
    	std::string bestRead = "none";
    	for(const auto & subPos : iter::range<uint32_t>(0, pos)){
    		if(clusters[subPos].seqBase_.cnt_ > clusters[pos].seqBase_.cnt_){
    			alignerObj.alignVec(clusters[subPos], clusters[pos], false);
    			alignerObj.profilePrimerAlignment(clusters[subPos], clusters[pos],setUp.weightHomopolymers_);
    			if(alignerObj.comp_.distances_.eventBasedIdentity_ > bestScore){
    				bestScore = alignerObj.comp_.distances_.eventBasedIdentity_;
    				bestPValue = clusters[subPos].getPValue(clusters[pos].seqBase_, alignerObj, false);
    				bestRead = clusters[subPos].seqBase_.name_;
    			}
    		}
    	}
    	pValueFile << clusters[pos].seqBase_.name_
    			<< "\t" << bestRead
					<< "\t" << std::setprecision(10) << bestPValue << std::endl;
    }
  }



  if (pars.markChimeras) {
    std::ofstream chimerasInfoFile;
    openTextFile(chimerasInfoFile,
                 setUp.directoryName_ + "chimeraNumberInfo.txt", ".txt", false,
                 false);
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
		 clusters, alignerObj,pars.parFreqs, 1, setUp.local_, chiOverlap,
		 overLapSizeCutoff, setUp.weightHomopolymers_, chiCount, allowableError);*/

    int clusterCount = 0;
    int readCount = 0;
    readVec::getCountOfReadNameContaining(clusters, "CHI", clusterCount);
    readVec::getReadCountOfReadNameContaining(clusters, "CHI", readCount);
    chimerasInfoFile << getPercentageString(clusterCount, clusters.size())
                     << "\t"
                     << getPercentageString(
                            readCount, readVec::getTotalReadCount(clusters))
                     << std::endl;
    if(setUp.verbose_){
    	std::cout << "Marked " << clusterCount << " as chimeric" << std::endl;
    }
  }


  if (pars.collapsingTandems) {
    reader.write(clusters, setUp.directoryName_ + setUp.ioOptions_.outFilename_ +
                               "_befroeTanCol",
                 setUp.ioOptions_);
    std::cout << "Collapsing on tandem repeat gaps" << std::endl;
    std::cout << "Starting with " << clusters.size() << " clusters"
              << std::endl;
    clusterCollapser::collapseTandems(clusters, alignerObj, setUp.runCutoff_,
                                      setUp.kLength_, setUp.kmersByPosition_,
																			pars.parFreqs, setUp.local_, true);
    clusters = readVecSplitter::splitVectorOnRemove(clusters).first;
    std::cout << "Collapsed down to " << clusters.size() << " clusters"
              << std::endl;
  }
  if (setUp.refFilename_ == "") {
    profiler::getFractionInfoCluster(clusters, setUp.directoryName_, "outputInfo");
  } else {
    profiler::getFractionInfoCluster(clusters, setUp.directoryName_, "outputInfo",
                              setUp.refFilename_, alignerObj, setUp.local_,
                              true);
  }
  std::ofstream startingInfo;
  openTextFile(startingInfo, setUp.directoryName_ + "startingInfo.txt", ".txt",
               false, false);
  startingInfo << "cluster\tstartingClusterName\tstartingSize" << std::endl;
	for (const auto& clus : clusters) {
		VecStr toks = tokenizeString(clus.firstReadName_, "_t");
		startingInfo << clus.seqBase_.name_ << "\t" << clus.firstReadName_ << "\t"
				<< toks.back() << std::endl;
	}
  if (pars.additionalOut) {
    std::string additionalOutDir = findAdditonalOutLocation(
    		pars.additionalOutLocationFile, setUp.ioOptions_.firstName_);
    if(additionalOutDir == ""){
    	std::cerr << bib::bashCT::red << bib::bashCT::bold;
    	std::cerr << "No additional out directory found for: " << setUp.ioOptions_.firstName_ << std::endl;
    	std::cerr << bib::bashCT::reset;
    }else{
      reader.write(clusters, additionalOutDir + setUp.ioOptions_.outFilename_,
                   setUp.ioOptions_);
      std::ofstream metaDataFile;
      openTextFile(metaDataFile, additionalOutDir + "/" + "metaData", ".json", setUp.ioOptions_);
      metaDataFile << metaData;
    }

  }

  reader.write(clusters, setUp.directoryName_ + setUp.ioOptions_.outFilename_,
                   setUp.ioOptions_);
  {
		std::string snpDir = bib::files::makeDir(setUp.directoryName_,
				"internalSnpInfo");
		for (const auto & readPos : iter::range(clusters.size())) {
			std::unordered_map<uint32_t, std::unordered_map<char, std::vector<baseReadObject>>> mismatches;

			for (const auto & subReadPos : iter::range(
					clusters[readPos].reads_.size())) {
				const auto & subRead = clusters[readPos].reads_[subReadPos];
				alignerObj.alignVec(clusters[readPos], subRead, false);
				//count gaps and mismatches and get identity
				alignerObj.profilePrimerAlignment(clusters[readPos], subRead,
						setUp.weightHomopolymers_);
				for (const auto & m : alignerObj.mismatches_) {
					mismatches[m.second.refBasePos][m.second.seqBase].emplace_back(
							subRead.seqBase_);
				}
			}
			table misTab { VecStr { "refPos", "refBase", "seqBase", "freq",
					"fraction", "seqs" } };
			for (const auto & m : mismatches) {
				for (const auto & seqM : m.second) {
					double totalCount = 0;
					for(const auto & read : seqM.second){
						totalCount += read.seqBase_.cnt_;
					}
					misTab.content_.emplace_back(
							toVecStr(m.first, clusters[readPos].seqBase_.seq_[m.first],
									seqM.first, totalCount,
									totalCount / clusters[readPos].seqBase_.cnt_,
									vectorToString(readVec::getNames(seqM.second), ",")));
				}
			}
			misTab.sortTable("seqBase", false);
			misTab.sortTable("refPos", false);
			misTab.outPutContents(
					CsvIOOptions(snpDir + clusters[readPos].seqBase_.name_,
							".tab.txt", "\t"));
		}
  }
  if(pars.createMinTree){
    std::string minTreeDirname = bib::files::makeDir(setUp.directoryName_, "minTree");
  	auto clusSplit = readVecSplitter::splitVectorOnReadFraction(clusters, 0.02);
    std::vector<readObject> tempReads;
    for(const auto & clus : clusSplit.first){
    	tempReads.emplace_back(readObject(clus.seqBase_));
    }
    auto minTreeJson = genMinTreeData(tempReads);
  	std::ofstream outFile(minTreeDirname + "psudoMinTree.json");
  	outFile << minTreeJson << std::endl;
  	std::string htmlOut = genHtmlStrForPsuedoMintree("psudoMinTree.json");
  	std::ofstream outHtmlFile(minTreeDirname + "psudoMinTree.html");
  	outHtmlFile << htmlOut << std::endl;
  }

  std::string clusterDirectoryName =
      bib::files::makeDir(setUp.directoryName_, "clusters");
  clusterVec::allWriteOutClusters(clusters, clusterDirectoryName, setUp.ioOptions_);

  if (!setUp.ioOptions_.processed_) {
    std::ofstream compStats;
    if (containsCompReads) {
      openTextFile(compStats, setUp.directoryName_ + "compStats.tab.txt",
                   ".txt", false, false);
      compStats << "cluster\tcompAmount" << std::endl;
    }
    std::string allInputReadsDir =
        bib::files::makeDir(setUp.directoryName_, "allInputReadsForEachCluster");
    for (const auto& clus : clusters) {
      std::vector<readObject> collection =
          clus.getAllBeginingClusters(identicalClusters);
      readVec::allSetFractionByTotalCount(collection);
      if (containsCompReads) {
        int currentCompAmount = 0;
        readVec::getCountOfReadNameContaining(collection, "_Comp",
                                              currentCompAmount);
        compStats << clus.seqBase_.name_ << "\t"
                  << getPercentageString(currentCompAmount, clus.seqBase_.cnt_)
                  << std::endl;
      }
      reader.write(collection, allInputReadsDir + clus.seqBase_.name_,setUp.ioOptions_);
    }
  }

  if (!pars.startWithSingles && pars.leaveOutSinglets) {
  	readObjectIO::write(singletons, setUp.directoryName_ + "singletons", setUp.ioOptions_);
  	std::ofstream singletonsInfoFile;
  	openTextFile(singletonsInfoFile, setUp.directoryName_ + "singletonsInfo", ".tab.txt", false, true);
  	singletonsInfoFile << "readCnt\treadFrac\n";
  	singletonsInfoFile << singletons.size() << "\t" << singletons.size()/static_cast<double>(reader.reads.size()) << std::endl;
  }

  if (setUp.writingOutAlnInfo_) {
    alignerObj.alnHolder_.write(setUp.outAlnInfoDirName_);
  }
  //log number of alignments done
  setUp.rLog_ << "Number of Alignments Done: " <<  alignerObj.numberOfAlingmentsDone_ << "\n";
  if(setUp.verbose_){
  	std::cout << "Number of Alignments Done: " << alignerObj.numberOfAlingmentsDone_ << std::endl;
  	//log time
  	setUp.logRunTime(std::cout);
  }

  return 0;
}


int SeekDeepRunner::processClusters(MapStrStr inputCommands) {
  // parameters
  SeekDeepSetUp setUp(inputCommands);
  processClustersPars pars;
  setUp.setUpMultipleSampleCluster(pars);

	if(containsSubString(pars.experimentName, ".")){
		throw std::runtime_error{"Error in populationCollapse::populationCollapse, populationName can't contain '.', " + pars.experimentName};
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
  std::map<std::string, std::pair<std::string, bool>> files =
      listFilesInDir(".", true);
  VecStr specificFiles;
  for (const auto& fileIter : files) {
    if (fileIter.first.find(setUp.ioOptions_.firstName_) != std::string::npos &&
        fileIter.first.find(".qual") == std::string::npos) {
      specificFiles.push_back(fileIter.first);
    }
  }
  std::cout << "Reading from" << std::endl;
  for (const auto& sfIter : specificFiles) {
    std::cout << sfIter << std::endl;
  }

  setUp.startARunLog(setUp.directoryName_);
  // parameters file
  setUp.writeParametersFile(setUp.directoryName_ + "parametersUsed.txt", false,
                            false);
  // reading in reads
  readObjectIO reader;
  uint64_t maxSize = 0;
  std::map<std::pair<std::string, std::string>, std::vector<readObject>>
      outputReads;
  for (const auto& strIter : specificFiles) {
    reader.read(setUp.ioOptions_.inFormat_, strIter, strIter + ".qual", true);
    readVec::handelLowerCaseBases(reader.reads, setUp.ioOptions_.lowerCaseBases_);
    readVec::getMaxLength(reader.reads, maxSize);
    VecStr toks = tokenizeString(strIter, "/");
    outputReads.insert(
        std::make_pair(std::make_pair(toks[1], toks[2]), reader.reads));
  }
  // output info about the read In reads
  for (const auto& readsIter : outputReads) {
    std::cout << readsIter.first.first << ": " << readsIter.first.second
              << " clustersNum: " << readsIter.second.size()
              << " readsNum: " << readVec::getTotalReadCount(readsIter.second)
              << std::endl;
    setUp.rLog_ << readsIter.first.first << ": " << readsIter.first.second
           << " clustersNum: " << readsIter.second.size()
           << " readsNum: " << readVec::getTotalReadCount(readsIter.second)
           << "\n";
  }
  // reading expected sequences to compare to
  bool checkingExpected = setUp.refFilename_ != "";
  std::vector<readObject> expectedSeqs;
  if (checkingExpected) {
    expectedSeqs = readObjectIO::getReferenceSeq(
        setUp.refFilename_, setUp.refFormat_, setUp.refProcessed_, maxSize);
  }
  // create aligner class object
  // Hard-coded parameters for alignment parameters follow follow:
  kmerMaps emptyMaps = kmerMaps();
  gapScoringParameters  gapPars(setUp.gapInfo_);
  aligner alignerObj = aligner(maxSize, gapPars, setUp.scoring_,
  		emptyMaps, setUp.primaryQual_,
                               setUp.secondaryQual_, setUp.qualThresWindow_,
                               setUp.countEndGaps_);
  if (setUp.alnInfoDirName_ != "") {
    processAlnInfoInput(alignerObj, setUp.alnInfoDirName_);
  }

  collapser collapserObj = collapser(
      !setUp.firstMatch_, setUp.bestMatchCheck_, setUp.local_, true,
      setUp.kmersByPosition_, setUp.runCutoff_,
      setUp.kLength_, setUp.verbose_, !setUp.largestFirst_, true,
      setUp.weightHomopolymers_, setUp.skipOnLetterCounterDifference_,
      setUp.fractionDifferenceCutOff_, setUp.adjustHomopolyerRuns_);

  collapserObj.opts_.lowQualityBaseTrim_ = pars.lowQualityCutOff;
  collapserObj.opts_.removeLowQualityBases_ = pars.removeLowQualBases;
  collapserObj.opts_.debug_ = setUp.debug_;

  // collect all the reads together for each sample
  std::map<std::string, std::vector<std::vector<cluster>>>
      clustersBySample;

  for (auto& readsIter : outputReads) {
    std::vector<cluster> clusters =
        baseCluster::convertVectorToClusterVector<cluster>(readsIter.second);
    readVecSorter::sortReadVector(clusters, "totalCount");
    // consider adding the sample name in the name as well
    renameReadNamesNewClusters(clusters, readsIter.first.second, true, true, false);
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
  				chiRunCutOff, chiOverlap, overLapSizeCutoff, chiCount, allowableError);
    	/*
    	clusterCollapser::markChimerasAdvanced(
    			clusters, alignerObj, pars.parFreqs, 1, setUp.local_, chiOverlap,
    			overLapSizeCutoff, setUp.weightHomopolymers_, chiCount, allowableError);*/
      //kmerMaps currentKmerMaps = kmerCalculator::indexKmerMpas(
        //  clusters, setUp.kLength_, setUp.runCutoff_, setUp.qualRunCutoff_);
      //alignerObj.setKmerMpas(currentKmerMaps);

      /*
      clusters = clusterCollapser::markChimeras(
          clusters, alignerObj, parFreqs, setUp.kmersByPosition_,
          setUp.kLength_, setUp.runCutoff_, setUp.local_,
          setUp.weightHomopolymers_);*/
    }
    // readVec::allSetFractionByTotalCount(clusters);
    clusterVec::allSetFractionClusters(clusters);
    clustersBySample[readsIter.first.first].emplace_back(clusters);
  }
  std::map<std::string, collapse::sampleCollapse> sampleCollapses;

  collapserObj.opts_.checkKmers_ = false;

  for (const auto& samp : clustersBySample) {
    std::cout << "Currently on : " << samp.first << std::endl;
    sampleCollapses[samp.first] =
        collapse::sampleCollapse(samp.second, samp.first, pars.clusterCutOff);
    // std::cout << "made sampleCollapse fine" << std::endl;
    // sampleCollapses[samp.first].updateInitialInfos(false);
    if(pars.onPerId){
      sampleCollapses[samp.first]
          .clusterOnPerId(collapserObj, pars.iteratorMap, pars.sortBy, alignerObj);
    }else{
      sampleCollapses[samp.first]
          .cluster(collapserObj, pars.iteratorMap, pars.sortBy, alignerObj);
    }

    // std::cout << "clustered fine" << std::endl;
    if (checkingExpected) {
      sampleCollapses[samp.first].collapsed_.checkAgainstExpected(
          expectedSeqs, alignerObj, setUp.local_, setUp.weightHomopolymers_);
    }

    sampleCollapses[samp.first].updateCollapsedInfos();
    // std::cout << "updated collapse infos fine" << std::endl;
    sampleCollapses[samp.first].updateExclusionInfos();
    // std::cout << "updated excusesion fine" << std::endl;
    sampleCollapses[samp.first].renameClusters(pars.sortBy);
    // std::cout << "ranmed clusters fine" << std::endl;
  }
  collapse::populationCollapse popCollapse;
  infoPrinter::printSampleCollapseInfo(
      sampleCollapses, checkingExpected,
      setUp.directoryName_ + "allClustersInfo.tab.txt", popCollapse, false);

  std::vector<sampleCluster> allSamples;
  std::string originalsDir = bib::files::makeDir(setUp.directoryName_, "originals");
  std::string initialDir = bib::files::makeDir(setUp.directoryName_, "initial");
  std::string excludedDir = bib::files::makeDir(setUp.directoryName_, "excluded");
  std::string excludedInitialDir = "";
  if(pars.writeExcludedOriginals){
  	 excludedInitialDir = bib::files::makeDir(setUp.directoryName_, "excludedInitial");
  }

  comparison onlyHpErrors;
  onlyHpErrors.oneBaseIndel_ = 1;
  onlyHpErrors.twoBaseIndel_ = 1;
  onlyHpErrors.largeBaseIndel_ = .99;

  std::string finalDir = bib::files::makeDir(setUp.directoryName_, "final");
  table chiInfoTab(VecStr{"sample", "numClustersSaved", "totalClustersChecked", "clusterSavedNames"});
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
  	if(!pars.keepChimeras){
  		if(pars.investigateChimeras){
    		//first mark the suspicous clusters as being chimeric
    	  for (auto &clus : sampCollapse.second.collapsed_.clusters_) {
    	    if (clusterVec::isClusterAtLeastChimericCutOff(clus.reads_, pars.chiCutOff)) {
    	      clus.seqBase_.markAsChimeric();
    	    }
    	  }
    	  //now check to see if it is ever the top two variants of a sample and if it is unmark it
  			for (auto & clus : sampCollapse.second.collapsed_.clusters_) {
  				if(clus.seqBase_.frac_ < pars.fracCutoff){
  					continue;
  				}
  				if (clus.seqBase_.name_.find("CHI_") != std::string::npos) {
  					bool saved = false;
  					for (const auto & otherCollapse : sampleCollapses) {

  						if (otherCollapse.first != sampCollapse.first
  								&& !containsSubString(stringToLowerReturn(otherCollapse.first), "control")
  								&& !containsSubString(stringToLowerReturn(sampCollapse.first), "control")) {
  							for(const auto & pos : iter::range(std::min(2UL,otherCollapse.second.collapsed_.clusters_.size()))){
  								alignerObj.alignVec(otherCollapse.second.collapsed_.clusters_[pos], clus, false);
  								alignerObj.profilePrimerAlignment(otherCollapse.second.collapsed_.clusters_[pos], clus, setUp.weightHomopolymers_);
  								if(onlyHpErrors.passErrorProfile(alignerObj.comp_)){

  									clus.seqBase_.unmarkAsChimeric();
  									for(auto & subRead : clus.reads_){
  										subRead.seqBase_.unmarkAsChimeric();
  									}
  									clus.resetInfos();
  									clustersSavedFromChi.emplace_back(clus.seqBase_.name_);
  									saved = true;
  									break;
  								}
  							}
  						}

  					}
  					if(!saved){
  						++clustersNotSaved;
  					}
  				}
  				sampCollapse.second.collapsed_.setSetInfo();
  			}
  		}
			//now exclude all marked chimeras, currently this will also remark chimeras unneccessarily
  		sampCollapse.second.excludeChimeras(false,pars.chiCutOff);
  	}
		chiInfoTab.content_.emplace_back(
				toVecStr(sampCollapse.first,
						getPercentageString(clustersSavedFromChi.size(),
								clustersSavedFromChi.size() + clustersNotSaved),
						clustersSavedFromChi.size() + clustersNotSaved,
						vectorToString(clustersSavedFromChi, ",")));
		if (bib::in(sampCollapse.first, customCutOffsMap)) {
			if (setUp.debug_) {
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
    sampCollapse.second.writeInitial(originalsDir, setUp.ioOptions_);
    sampCollapse.second.writeExcluded(excludedDir, setUp.ioOptions_);
    if(pars.writeExcludedOriginals){
    	sampCollapse.second.writeExcludedOriginalClusters(excludedInitialDir, setUp.ioOptions_);
    }
    sampCollapse.second.writeFinal(finalDir, setUp.ioOptions_);
    sampCollapse.second.writeFinalOrignalClusters(initialDir, setUp.ioOptions_);
    // add to all sample cluster

    addOtherVec(allSamples, sampCollapse.second.createOutput(false, sortBy));
    // sampCollapse.second.updateCollapsedInfos(true);
    // sampCollapse.second.updateExclusionInfos(true);
  }
  CsvIOOptions chiOutOptions(setUp.directoryName_ + "chiInfo", ".tab.txt", "\t");
  chiInfoTab.outPutContents(chiOutOptions);

  if (!pars.noPopulation) {
		std::cout << bib::bashCT::boldGreen("Pop Clustering") << std::endl;
		popCollapse = collapse::populationCollapse(allSamples, pars.experimentName);

    std::string popDir = bib::files::makeDir(setUp.directoryName_, "population");
		if (pars.onPerId) {
			popCollapse.popClusterOnId(collapserObj, pars.popIteratorMap, "fraction",
					alignerObj);
		} else {
			popCollapse.popCluster(collapserObj, pars.popIteratorMap, "fraction",
					alignerObj);
		}
		popCollapse.updateInfoWithSampCollapses(sampleCollapses);
    readObjectIO popReader;
    if (pars.previousPopFilename != "") {
      std::string format = bib::files::getExtension(pars.previousPopFilename);
      popReader.read(format, pars.previousPopFilename);
      popCollapse.renameToOtherPopNames(popReader.reads, alignerObj);
    }

    if (checkingExpected) {
      popCollapse.collapsed_.checkAgainstExpected(
          expectedSeqs, alignerObj, setUp.local_, setUp.weightHomopolymers_);
    }

    infoPrinter::printPopulationCollapseInfo(
        popCollapse, popDir + "populationCluster.tab.txt", checkingExpected);
    std::string popInitialDir = bib::files::makeDir(popDir, "initial");
    popCollapse.writeFinalInitial(popInitialDir, setUp.ioOptions_);
    popCollapse.writeFinal(popDir, setUp.ioOptions_.outFormat_,
                           setUp.ioOptions_.overWriteFile_,
                           setUp.ioOptions_.exitOnFailureToWrite_);

    std::unordered_map<std::string, bib::color> colorsForGraph;
    colorsForGraph =
        getColorsForNames(popCollapse.collapsed_.clusters_, pars.sat, pars.lum);

    std::string dotDir = bib::files::makeDir(setUp.directoryName_, "dotFiles");
    for (const auto& samp : sampleCollapses) {
      std::vector<readObject> readsForGraph;
      if (samp.second.collapsed_.clusters_.empty()) {
        continue;
      }
      for (const auto& clus : samp.second.collapsed_.clusters_) {
        std::string popName =
            popCollapse.collapsed_.clusters_
                [popCollapse.collapsed_.subClustersPositions_.at(
                     clus.getStubName(true))].seqBase_.name_;
        readsForGraph.emplace_back(
            readObject(seqInfo(popName, clus.seqBase_.seq_, clus.seqBase_.qual_,
                               clus.seqBase_.cnt_, clus.seqBase_.frac_)));
        readsForGraph.back().seqBase_.cnt_ = clus.seqBase_.cnt_;
        readsForGraph.back().seqBase_.frac_ = clus.seqBase_.frac_;
      }
      bestDistGraph currentGraphDist(readsForGraph, alignerObj, setUp.local_,
                                     samp.first);
      currentGraphDist.createDotBestConnectedFile(dotDir, colorsForGraph,
                                                  false);
      for(auto & read : readsForGraph){
      	read.seqBase_.name_ = read.getReadId();
      }
      //auto treeData = genMinTreeData(readsForGraph, alignerObj);
      //std::ofstream outDot(dotDir + samp.first + "_2.dot");
      //jsonTreeToDot(treeData, outDot);
    }
    if(pars.groupingsFile != ""){
    	table groupsTab(pars.groupingsFile, "\t", true);
    	/**@todo add in safety checks */
    	VecStr sampleNames = groupsTab.getColumn(0);
    	std::unordered_map<std::string, std::unordered_map<std::string, VecStr>> groups;
    	for(const auto & colPos : iter::range<uint32_t>(1,groupsTab.columnNames_.size())){
    		VecStr col = groupsTab.getColumn(colPos);
    		for(const auto & pos : iter::range(col.size())){
    			if(sampleCollapses.find(sampleNames[pos]) == sampleCollapses.end()){
    				std::cout << sampleNames[pos] << " was not clustered" << "\n";
    			}else{
    				groups[groupsTab.columnNames_[colPos]][col[pos]].emplace_back(sampleNames[pos]);
    			}
    		}
    	}
    	std::string groupDir = bib::files::makeDir(setUp.directoryName_, "groups");
			for (const auto & group : groups) {
				std::string currentGroupDir = bib::files::makeDir(groupDir,
						group.first);
				for (const auto & subGroup : group.second) {
					std::string currentSubGroupDir = bib::files::makeDir(currentGroupDir,
							subGroup.first);
					std::ofstream popFile;
					std::ofstream sampFile;
					openTextFile(popFile, currentSubGroupDir + "popFile.tab.txt", ".txt",
							false, false);

					openTextFile(sampFile, currentSubGroupDir + "sampFile.tab.txt",
							".txt", false, false);
					std::set<uint32_t> otherPopPositions;
					for (const auto & otherGroup : group.second) {
						if(otherGroup.first != subGroup.first){
						  for(const auto & samp : otherGroup.second){
						  	/**@todo should check to see if the sample exist*/
						  	for(const auto & subClus : sampleCollapses.at(samp).collapsed_.clusters_){
						  		otherPopPositions.emplace(popCollapse.collapsed_.subClustersPositions_.at(subClus.getStubName(true)));
						  	}
						  }
						}
					}
					infoPrinter::printInfoForSamps(sampleCollapses, checkingExpected,
							sampFile, popFile, popCollapse, !pars.noPopulation, subGroup.second,
							group.first + ":" + subGroup.first,otherPopPositions);
				}
			}
    }
  }

  infoPrinter::printSampleCollapseInfo(
      sampleCollapses, checkingExpected,
      setUp.directoryName_ + "selectedClustersInfo.tab.txt", popCollapse,
			!pars.noPopulation);

  std::string sampRepInfoDir = bib::files::makeDir(setUp.directoryName_, "sampRepAgreementInfo");
	njhRInside::OwnRInside rSes;
	auto & R = rSes.get();
	std::string colFunc= "nickColors<-function(n, h = c(120,420), l = c(.40,.70), s = c(.8,1), alpha = 1){\n"
	"require(colorspace)\n"
	"return (hex(HLS(seq(h[1],h[2],length.out = n), seq(l[1],l[2],length.out = n), seq(s[1],s[2],length.out=n))))\n"
  "}";
	R.parseEvalQ(colFunc);
	rSes.openPdfDevice(sampRepInfoDir + "repComparisonGraphs.pdf", 11.5, 8);
	table rmseTab(VecStr{"sampleName", "RMSE"});
	table repInfoTab(VecStr{"sampleName", "clusterName","repName", "fraction"});
	for (const auto & samp : sampleCollapses) {
		auto repInfoforSample = samp.second.collapsed_.getReplicateInfo();
		for(const auto & row : repInfoforSample.content_){
			VecStr addingRow{samp.first};
			addOtherVec(addingRow, row);
			repInfoTab.content_.emplace_back(addingRow);
		}

		//std::cout << samp.second.collapsed_.clusters_.size() << std::endl;
		if(!samp.second.collapsed_.clusters_.empty()){
			rSes.installLib("ggplot2");
			Rcpp::DataFrame inData = Rcpp::DataFrame::create(
					Rcpp::Named("clusterName") = repInfoforSample.getColumn("clusterName"),
					Rcpp::Named("repName") = repInfoforSample.getColumn("repName"),
					Rcpp::Named("fraction") = vecStrToVecNum<double>(repInfoforSample.getColumn("fraction")));

			R["inData"] = inData;
			std::string plotCmd =
					"print(ggplot(inData, aes(x = clusterName, y = fraction, fill = repName)) + geom_bar(stat = \"identity\", position =\"dodge\")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(\""
							+ samp.first
							+ "\") + scale_fill_manual(values = nickColors(length(levels(inData$repName)) + 1)))";
			R.parseEvalQ(plotCmd);
			rmseTab.content_.emplace_back(toVecStr(samp.first, samp.second.collapsed_.getRMSE()));
		}
	}
	repInfoTab.outPutContents(
					CsvIOOptions(sampRepInfoDir + "replicatesFractions", ".tab.txt", "\t", false));
	rmseTab.sortTable("RMSE", true);
	rmseTab.outPutContents(CsvIOOptions(sampRepInfoDir + "RMSE.tab.txt", ".tab.txt", "\t"));
  rSes.closeCurrentDevice();
  if (setUp.writingOutAlnInfo_) {
    alignerObj.alnHolder_.write(setUp.outAlnInfoDirName_);
  }
  std::cout << alignerObj.numberOfAlingmentsDone_ << std::endl;
  setUp.rLog_ << alignerObj.numberOfAlingmentsDone_ << "\n";

  setUp.logRunTime(std::cout);
  return 0;
}


int SeekDeepRunner::makeSampleDirectories(MapStrStr inputCommands) {
  SeekDeepSetUp setUp(inputCommands);
  std::string sampleNameFilename = "";
  bool separatedDirs = false;
  setUp.setOption(separatedDirs, "--separatedDirs", "Create a separate directory for each index");
  setUp.setUpMakeSampleDirectories(sampleNameFilename);
  setUpSampleDirs(sampleNameFilename,setUp.directoryName_, separatedDirs);
  setUp.startARunLog(setUp.directoryName_);
  setUp.logRunTime(std::cout);
  return 0;
}


}  // namespace bib
