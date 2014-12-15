//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012, 2014 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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


#include "SeekDeepSetUp.hpp"
#include "SeekDeepRunner.hpp"
#include <bibcpp.h>
#include <functional>

namespace bibseq {

SeekDeepRunner::SeekDeepRunner()
    : bib::progutils::programRunner({addFunc("extractor", extractor, false),
                     addFunc("processClusters", processClusters, false),
                     addFunc("qluster", qluster, false),
                     addFunc("makeSampleDirectories", makeSampleDirectories, false),},
                    "SeekDeep") {}



namespace readVec{


}

int SeekDeepRunner::extractor(MapStrStr inputCommands) {

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
  		"Cut Off for fraction of bases above qual check of " + to_string(qualCheck))){
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

  // int readsExtractedComp=0;
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
  	for_each(reader.reads, [&](readObject & read){
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
  readsExtracted = (int)reader.reads.size();
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
  aligner alignObj = aligner(maxReadSize, gapPars, scoreMatrix,
  		emptyMaps, primaryQual,
                             setUp.qualThresWindow_, secondayQual, countEndGaps
                             );
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
    readsNotMatchedToBarcode = (int)unrecognizedBarcode.size();
  } else {
    if (!multiplex) {
      reads.insert({"", largeReads});
    } else {

      unrecognizedBarcode = readVecExtractor::extractReadsByMIDFull(
          largeReads, reads, idFilename, idFileDelim, 0, checkComplement,
          false, true,variableStart, variableStop, barcodeErrors);
      readsNotMatchedToBarcode = (int)unrecognizedBarcode.size();
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
  std::string smallFragmentDir =
      bib::files::makeDir(filteredOffDir, "smallFragments");
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
  	profileLog << "\tq" + to_string(qualCheck) + "<" << qualCheckCutOff;
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
        filteredReads.first = readVecSplitter::splitVectorOnSeqSimilarityAdd(
            filteredReads.first, compareObject, contamination,
            contaminationNumber);
        reader.write(contamination, contaminationDir + readsIter.first +
                                        rIter.first + "_contamination",
                     setUp.ioOptions_);
        contaminationAmount += (int)contamination.size();
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
        runParams.errors_.largeBaseIndel_ = 0;
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
      	renameReadNames(filteredReads.first, currentMidGeneName, false, true, true, "none");
      }
      reader.write(filteredReads.first, tempNameGood,
                   setUp.ioOptions_);
      readsFinal += (int)filteredReads.first.size();
      readVec::getCountOfReadNameContaining(filteredReads.first, "Comp",
                                            readsFinalComp);
      if (pyroExtract) {
        reader.write(filteredReads.first, tempNameGood, "pyroData",
                     setUp.ioOptions_.overWriteFile_,
                     setUp.ioOptions_.exitOnFailureToWrite_, maxFlowCutoff);
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
/*
auto runCollapse = [&](collapser col, std::vector<cluster> & clusVec,
		runningParameters rPars, aligner & currentAlignerObj ){
	 col.collapseWithParametersRegKmer(clusVec, rPars, currentAlignerObj);

	 return;
};*/

template<typename T>
uint32_t getMismatches(const T & read1,
				const T & read2,
				aligner alignerObj, bool weightHomopolymers){
	alignerObj.alignVec(read1.seqBase_,read2.seqBase_, false);
	alignerObj.profilePrimerAlignment(read1.seqBase_, read2.seqBase_, weightHomopolymers);
	return alignerObj.errors_.hqMismatches_;
};

template<typename T>
Json::Value getMinMismatchTreeJson(const std::vector<T> & reads, aligner & alignerObj,
		uint32_t numThreads, bool weightHomopolymers,double hueStart, double hueStop,
    double lumStart, double lumStop,
    double satStart, double satStop, std::string backgroundColor){
	std::function<uint32_t(const T & ,
	  		const T &, aligner, bool)> misFun = getMismatches<T>;
		auto misDistances = getDistanceCopy(reads, numThreads, misFun,
				alignerObj, weightHomopolymers);
	  readDistGraph<uint32_t> graphMis(misDistances, reads);
		std::vector<std::string> names;
	  for(const auto & n : graphMis.nodes_){
	  	names.emplace_back(n->name_);
	  }
	  if(hueStop == 360 && hueStart == 0){
	  	hueStop = 360 - (360.0/names.size());
	  }
		auto nameColors = bib::getColorsForNames(names, hueStart, hueStop,
				lumStart, lumStop, satStop, satStart);
	  Json::Value graphJson = graphMis.toJsonMismatchGraphAll(bib::color(backgroundColor), nameColors);
	  return graphJson;
}



int SeekDeepRunner::qluster(MapStrStr inputCommands) {
  // parameters
  std::string parameters = "";
  bool snapShots = false, extra = false;
  std::string qualRep = "median", sortBy = "totalCount";
  bool smallestFirst = true, bestMatch = true, markChimeras = false;
  int bestMatchCheck = 10;
  int parFreqs = 2;
  bool collapsingTandems = false;
  bool kmerCheckingOnAlignProfile = false;
  bool additionalOut = false;
  std::string additionalOutLocationFile = "";

  //std::map<char, std::map<char, int>> scoringMatrixMap;
  std::map<int, std::vector<double>> iteratorMap;
  bool removeLowQualBases = false;
  int lowQualityCutOff = 3;
  bool adjustHomopolyerRuns = false;
  bool condensedCollapse = true;
  // bool byIndices = false;
  bool skipOnLetterCounterDifference = false;
  double fractionDifferenceCutOff = 0.05;
  SeekDeepSetUp setUp(inputCommands);
  uint32_t runTimes = 100;
  // const auto& errorLookUp = simulation::constants::QualErrorArr;
  bool sim = false;
  bool regKmers = true;
  bool printSimClusters = false;
  double pValueCutOff = 0.01;
  double fdrCutOff = 0.01;
  bool removeSingletons = false;
  bool debug = false;
  bool createMinTree = false;
  std::string diffCutOffStr = "0.1";
  bool findBest = true;
  std::string alphabetStr = "A,C,G,T";
  setUp.setOption(diffCutOffStr, "-diffCutOffs","diffCutOff");
  std::vector<double> diffCutOffVec = vecStrToVecNum<double>(tokenizeString(diffCutOffStr,","));
  bool keepCondensed = false;
  setUp.setOption(keepCondensed, "-keepCondensed", "keepCondensed");
  setUp.setOption(findBest, "-findBest","findBest");
  setUp.setBoolOptionFalse(findBest, "-findNucFirst,-findFirstNuc,-firstNuc","findBest");
  setUp.setOption(alphabetStr, "-alph,-alphabet","alphabetStr");
  std::vector<char> alphabet = processAlphStrVecChar(alphabetStr, ",");
  bool useNucComp = false;
  setUp.setOption(useNucComp, "-useNucComp", "useNucComp");
  bool useMinLenNucComp = false;
  setUp.setOption(useMinLenNucComp, "-useMinLenNucComp", "useMinLenNucComp");
  setUp.setOption(debug, "-debug", "debug");
  setUp.setOption(removeSingletons,
                  "-startWithoutSingles",
                  "Start Clustering Without Singles");
  setUp.setOption(fdrCutOff, "-fdrCutOff", "False_discovery_cutoff");
  setUp.setOption(pValueCutOff, "-pValueCutOff", "pValueCutOff");
  setUp.setOption(printSimClusters, "-printSimClusters,-printSim",
                  "printSimClusters");

  setUp.setBoolOptionFalse(regKmers, "-qualKmers,-qualK,-qualKmer", "doQualKmer");

  setUp.setOption(sim, "-sim", "sim");
  setUp.setOption(runTimes, "-runTimes", "runTimes");
  setUp.setOption(skipOnLetterCounterDifference, "-skip",
                  "skipOnLetterCounterDifference");
  setUp.setOption(fractionDifferenceCutOff, "-skipCutOff",
                  "fractionDifferenceCutOff");
  setUp.processAlnInfoInput();
  setUp.setOption(createMinTree, "--createMinTree",
  		"Create Psudo minimum Spanning Trees For Mismatches for Final Clusters");

	std::string backgroundColor = "#000000";
	double hueStart = 0;
	double hueStop = 360;
	double satStart = 1.0;
	double satStop = 1.0;
	double lumStart = 0.5;
	double lumStop = 0.5;
	setUp.setOption(backgroundColor, "-b,--backgroundColor", "Hex String for the color for the Background");
	setUp.setOption(hueStart, "--hueStart", "Hue Start for the color of Reads");
	setUp.setOption(hueStop, "--hueStop", "Hue Stop for the color of Reads");
	setUp.setOption(satStart, "--satStart", "Sat Start for the color of Reads");
	setUp.setOption(satStop, "--satStop", "Sat Stop for the color of Reads");
	setUp.setOption(lumStart, "--lumStart", "Lum Start for the color of Reads");
	setUp.setOption(lumStop, "--lumStop", "Lum Stop for the color of Reads");
  setUp.setUpClusterDown(
      qualRep, parameters, extra, iteratorMap, smallestFirst, markChimeras,
      parFreqs, bestMatch, bestMatchCheck, snapShots, sortBy, additionalOut,
      additionalOutLocationFile, collapsingTandems, kmerCheckingOnAlignProfile,
      condensedCollapse, removeLowQualBases, lowQualityCutOff,
      adjustHomopolyerRuns);
  // print out the parameters read in
/*
  for (auto& iteratorMapIter : iteratorMap) {
    outputIterationMapIter(iteratorMapIter.first, iteratorMapIter.second,
                           std::cout);
  }*/
  if(useNucComp){
  	std::cout << "nuc comp cutoffs" << std::endl;
  	printVector(diffCutOffVec, ", ", std::cout);
  }
  // read in the sequences
  readObjectIO reader;
  reader.read(setUp.ioOptions_);
  auto splitOnSize = readVecSplitter::splitVectorBellowLength(reader.reads, setUp.kLength_ * 2);
  reader.reads = splitOnSize.first;
  if(!splitOnSize.second.empty()){
  	readObjectIO::write(splitOnSize.second, setUp.directoryName_ + "smallReads", setUp.ioOptions_);
  }
  /*std::cout << "after reading: " << std::endl;
  uint32_t afterReading = seqQualSizeCheck(reader.reads);
  if(afterReading > 0){
    exit(1);
  }*/

  if (removeLowQualBases) {
    readVec::allRemoveLowQualityBases(reader.reads, lowQualityCutOff);
  }

  /*
  std::cout << "after low quality removal: " << std::endl;
  uint32_t afterQRemoval = seqQualSizeCheck(reader.reads);
  if(afterQRemoval > 0){
    exit(1);
  }*/
  if (adjustHomopolyerRuns) {
    readVec::allAdjustHomopolymerRunsQualities(reader.reads);
  }
  /*
  std::cout << "after hp adjust: " << std::endl;
  uint32_t afterHPAdjust = seqQualSizeCheck(reader.reads);
  if(afterHPAdjust > 0){
    exit(1);
  }*/
  /*
    std::map<std::string, std::map<double, uint32_t>> counts;
    std::map<std::string, std::map<double, uint32_t>> misMatchCounts;
    if (sim) {
      readVec::updateQualCountsMultiple(reader.reads, counts,
                                        setUp.qualThresWindow_, errorLookUp);
    }*/

  // std::string seqName = getFileName(setUp.ioOptions_.firstName_);

  bool containsCompReads = false;

  int compCount = 0;
  readVec::getCountOfReadNameContaining(reader.reads, "_Comp", compCount);
  if (compCount > 0) {
    containsCompReads = true;
  }

  readVecSorter::sortReadVector(reader.reads, sortBy);
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
  processRunCutoff(setUp.qualRunCutoff_, setUp.qualRunCutOffString_,
                   readVec::getTotalReadCount(reader.reads));
  if (setUp.qualRunCutoff_ < setUp.runCutoff_) {
    setUp.qualRunCutoff_ = setUp.runCutoff_;
  }
  std::cout << "run cut off is " << setUp.runCutoff_ << std::endl;
  std::cout << "qaul run cut off is " << setUp.qualRunCutoff_ << std::endl;

  // make the runLog, this is what is seen on the terminal screen at run time
  setUp.startARunLog(setUp.directoryName_);
  // parameter file
  setUp.writeParametersFile(setUp.directoryName_ + "parametersUsed.txt", false,
                            false);
  // read in the paramteres from the parameters file
  setUp.rLog_ << "Parameters used" << "\n";
  textFileReader txtReader(":");
  txtReader.readFile(parameters);
  txtReader.fileContent.outPutContents(setUp.rLog_.runLogFile_, ":");

  std::cout << "Reading clusters from " << setUp.ioOptions_.firstName_ << " "
            << setUp.ioOptions_.secondName_ << std::endl;
  setUp.rLog_ << "Reading clusters from " << setUp.ioOptions_.firstName_
                  << " " << setUp.ioOptions_.secondName_ << "\n";
  std::cout << "Read in " << counter << " reads" << std::endl;
  setUp.rLog_ << "Read in " << counter << " reads" << "\n";

  // create cluster vector
  std::vector<identicalCluster> identicalClusters;
  std::vector<cluster> clusters;
  if (setUp.ioOptions_.processed_) {
    clusters = baseCluster::convertVectorToClusterVector<cluster>(reader.reads);
  } else {
    identicalClusters = clusterCollapser::collapseIdenticalReads(
        reader.reads, qualRep, setUp.ioOptions_.lowerCaseBases_);
    clusters =
        baseCluster::convertVectorToClusterVector<cluster>(identicalClusters);
  }

  /*
  std::cout << "after collapse: " << std::endl;
  uint32_t afterCollapse = seqQualSizeCheck(clusters);
  if(afterCollapse > 0){
    exit(1);
  }*/
  // reader.reads.front().printDescription(std::cout);
  // clusters.front().printDescription(std::cout);
  std::cout << "Identical clusters numbers: " << clusters.size() << std::endl;
  setUp.rLog_ << "Identical clusters numbers: " << clusters.size()
                  << "\n";
  readVecSorter::sortReadVector(clusters, sortBy);
  kmerMaps kMaps = kmerCalculator::indexKmerMpas(
      clusters, setUp.kLength_, setUp.runCutoff_, setUp.qualRunCutoff_,
      setUp.expandKmerPos_, setUp.expandKmerSize_);
  // create aligner class object
  gapScoringParameters gapPars(setUp.gapInfo_);
  gapPars.setIdentifer();
  aligner alignerObj = aligner(maxSize, gapPars, setUp.scoring_, kMaps, setUp.primaryQual_,
                               setUp.secondaryQual_, setUp.qualThresWindow_,
                               setUp.countEndGaps_);
  if(debug){
  	for (const auto& i : iter::range(kMaps.qualRunCutOffs_.size())) {
			std::cout << i << "\t" << std::pow(10.0, -(i / 10.0)) << "\t"
								<< kMaps.qualRunCutOffs_[i] << std::endl;
		}
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
  	/*
  	for(const auto & f : setUp.scoring_){
  		for(const auto & s : f.second){
  			scoreMapFile << f.first << ":" << s.first << " = " << s.second << std::endl;;
  		}
  	}*/
  	alignerObj.parts_.gapScores_.printDescription(std::cout, true);
  }
  processAlnInfoInput(alignerObj, setUp.alnInfoDirName_);


  collapser collapserObj = collapser(
      bestMatch, bestMatchCheck, setUp.local_, setUp.checkKmers_,
      setUp.kmersByPosition_, setUp.runCutoff_, setUp.qualRunCutoff_,
      setUp.kLength_, setUp.verbose_, smallestFirst, condensedCollapse,
      setUp.weightHomopolymers_, skipOnLetterCounterDifference,
      fractionDifferenceCutOff, regKmers, setUp.adjustHomopolyerRuns_);

  std::string snapShotsDirectoryName = "";
  if (snapShots) {
    snapShotsDirectoryName = bib::files::makeDir(setUp.directoryName_, "snapShots");
  }
  uint32_t singletonNum = 0;
  std::vector<cluster> singletons;
  if (removeSingletons) {
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

  clusterVec::allSetFractionClusters(clusters);
  readVec::allSetLetterCount(clusters, alphabet);
  /*for (auto& clus : clusters) {
    clus.counter_.setFractions();
  };*/
  //auto testScore = alignerObj.getGapScoring();
  //testScore.printDescription(std::cout, true);
  uint32_t minLen = readVec::getMinLength(clusters);
  if (useNucComp) {
  	bool oldVerbose = collapserObj.verbose_;
  	collapserObj.verbose_ = debug;
  	for(const auto & diffCutOff : diffCutOffVec){
      clusterVec::allSetFractionClusters(clusters);
      readVec::allSetLetterCount(clusters, alphabet);
    	//std::vector<cluster> processedClusters;
    	std::vector<nucCompCluster> comps;
    	if(useMinLenNucComp){
    		comps = clusterOnNucComp(clusters,minLen, alphabet, diffCutOff, findBest, true, setUp.debug_);
    	}else{
    		comps = clusterOnNucComp(clusters, alphabet, diffCutOff, findBest, true, setUp.debug_);
    	}

    	std::cout << "On nuc comp diff of " << diffCutOff << " in " << vectorToString(diffCutOffVec, ", ") << std::endl;
    	uint32_t compCount = 0;
    	for(auto & comp : comps){
    		compCount++;
    		if(compCount % 50 == 0){
    			std::cout << "On " << compCount << " of " << comps.size() << std::endl;
    		}
    		if(comp.readPositions_.size() < 2){
    			continue;
    		}
    		/*std::vector<cluster> currentClusters;
    		for(const auto & readPos : comp.readPositions_){
    			currentClusters.emplace_back(clusters[readPos]);
    		}*/
    		//addOtherVec(processedClusters, collapserObj.runClustering(currentClusters, iteratorMap, alignerObj));
    		//collapserObj.runClustering(clusters, comp.readPositions_, iteratorMap, alignerObj);
    		//auto oldMaps = alignerObj.getKmerMaps();
    		double currentReadCnt = 0;
    		for(const auto & pos : comp.readPositions_){
    			currentReadCnt += clusters[pos].seqBase_.cnt_;
    		}
    		auto oldRunCutoff = alignerObj.kMaps_.runCutOff_ ;
    		auto oldQualRunCutoff = alignerObj.kMaps_.qualRunCutOff_;
    		uint32_t currentRunCutoff = processCutOffStr(setUp.runCutOffString_,currentReadCnt );
    		uint32_t currentQualRunCutoff = processCutOffStr(setUp.qualRunCutOffString_,currentReadCnt );
    		alignerObj.kMaps_.runCutOff_ = currentRunCutoff;
    		alignerObj.kMaps_.qualRunCutOff_ = currentQualRunCutoff;
    		//auto currentKMaps = kmerCalculator::indexKmerMpas(clusters, kClus.readPositions_,setUp.kLength_,
    			//	currentRunCutoff ,currentQualRunCutoff,setUp.expandKmerPos_, setUp.expandKmerSize_);
    		//alignerObj.setKmerMpas(currentKMaps);

    		collapserObj.runClustering(clusters, comp.readPositions_, iteratorMap, alignerObj);
    		alignerObj.kMaps_.runCutOff_ = oldRunCutoff;
    		alignerObj.kMaps_.qualRunCutOff_ = oldQualRunCutoff;
    		//alignerObj.setKmerMpas(oldMaps);
    	}
      //clusters = processedClusters;
  	}
  	collapserObj.verbose_ = oldVerbose;
  	clusters = readVecSplitter::splitVectorOnRemove(clusters).first;
  }


  //alignerObj.kMaps_.outputKmerInfo(alignerObj.kMaps_, std::cout);
  readVecSorter::sortReadVector(clusters, "totalCount");
  for (uint32_t i = 1; i <= iteratorMap.size(); ++i) {

  	std::cout << std::endl;
    uint32_t sizeOfReadVector = readVec::getReadVectorSize(clusters);
    runningParameters runPars(iteratorMap[i],i,  sizeOfReadVector);
    runPars.printIterInfo(std::cout, true);
    if (condensedCollapse && runPars.errors_.hqMismatches_ == 0 &&
        runPars.errors_.lqMismatches_ == 0 &&
        runPars.errors_.largeBaseIndel_ == 0 &&
        runPars.errors_.twoBaseIndel_ == 0) {
      auto byCondensed = readVec::organizeByCondensedSeq(clusters);
    	if (collapserObj.regKmer_) {
    		for (auto & condensedReads : byCondensed){
    			collapserObj.collapseWithParametersRegKmer(condensedReads.second,
    			                                                  runPars, alignerObj);
    		}
      } else {
      	for (auto & condensedReads : byCondensed) {
      		collapserObj.collapseWithParametersQualKmer(condensedReads.second,
                                                    runPars, alignerObj);
      	}
      }

      clusters.clear();
      for (const auto& condensedReads : byCondensed) {
        addOtherVec(clusters, condensedReads.second);
      }
      int sizeOfReadVector = readVec::getReadVectorSize(clusters);
      std::cout << "collapsed down to " << sizeOfReadVector << std::endl;
    } else {
      if (collapserObj.regKmer_) {
        collapserObj.collapseWithParametersRegKmer(clusters, runPars,
                                                   alignerObj);
      } else {
        collapserObj.collapseWithParametersQualKmer(clusters, runPars,
                                                    alignerObj);
      }
    }
    readVec::allUpdateName(clusters);
    readVecSorter::sortReadVector(clusters, sortBy);
    clusterVec::allCalculateConsensus(clusters, alignerObj, true);
    if (adjustHomopolyerRuns) {
      readVec::allAdjustHomopolymerRunsQualities(clusters);
    }

    if (snapShots) {
      std::string iterDir =
          bib::files::makeDir(snapShotsDirectoryName, std::to_string(i));
      std::vector<cluster> currentClusters =
          readVecSplitter::splitVectorOnRemove(clusters).first;
      std::string seqName = bib::files::getFileName(setUp.ioOptions_.firstName_);
      renameReadNames(currentClusters, seqName, true, false);
      reader.write(currentClusters, snapShotsDirectoryName + std::to_string(i),
                   setUp.ioOptions_);
      clusterVec::allWriteOutClusters(currentClusters, iterDir);
      if (setUp.refFilename_ == "") {
        profiler::getFractionInfoCluster(currentClusters, snapShotsDirectoryName,
                                  std::to_string(i) + ".tab.txt");
      } else {
        profiler::getFractionInfoCluster(currentClusters, snapShotsDirectoryName,
                                  std::to_string(i) + ".tab.txt",
                                  setUp.refFilename_, alignerObj, setUp.local_,
                                  true);
      }
    }
    std::cout << "Current duration: ";
    setUp.logRunTime(std::cout);
    setUp.rLog_ << "Current duration: "; setUp.rLog_.logTotalTime(6);
    setUp.rLog_ << "Current clusters size " << sizeOfReadVector
                    << "\n";
    /*std::cout << "after clustering: " << std::endl;
    uint32_t afterClustering = seqQualSizeCheck(clusters);
    if(afterClustering > 0){
      exit(1);
    }*/
  }

  if (removeSingletons) {
    addOtherVec(clusters, singletons);
    if (useNucComp) {
    	bool oldVerbose = collapserObj.verbose_;
    	collapserObj.verbose_ = debug;
    	for(const auto & diffCutOff : diffCutOffVec){
        clusterVec::allSetFractionClusters(clusters);
        readVec::allSetLetterCount(clusters, alphabet);
      	//std::vector<cluster> processedClusters;
      	std::vector<nucCompCluster> comps;
      	if(useMinLenNucComp){
      		comps = clusterOnNucComp(clusters,minLen, alphabet, diffCutOff, findBest, true, setUp.debug_);
      	}else{
      		comps = clusterOnNucComp(clusters, alphabet, diffCutOff, findBest, true, setUp.debug_);
      	}
      	std::cout << "On nuc comp diff of " << diffCutOff << " in " << vectorToString(diffCutOffVec, ", ") << std::endl;
      	uint32_t compCount = 0;
      	for(auto & comp : comps){
      		compCount++;
      		if(compCount % 50 == 0){
      			std::cout << "On " << compCount << " of " << len(comps) << std::endl;
      		}
      		if(comp.readPositions_.size() < 2){
      			continue;
      		}
      		/*std::vector<cluster> currentClusters;
      		for(const auto & readPos : comp.readPositions_){
      			currentClusters.emplace_back(clusters[readPos]);
      		}*/
      		//addOtherVec(processedClusters, collapserObj.runClustering(currentClusters, iteratorMap, alignerObj));
      		//collapserObj.runClustering(clusters, comp.readPositions_, iteratorMap, alignerObj);
      		//auto oldMaps = alignerObj.getKmerMaps();
      		double currentReadCnt = 0;
      		for(const auto & pos : comp.readPositions_){
      			currentReadCnt += clusters[pos].seqBase_.cnt_;
      		}
      		auto oldRunCutoff = alignerObj.kMaps_.runCutOff_ ;
      		auto oldQualRunCutoff = alignerObj.kMaps_.qualRunCutOff_;
      		uint32_t currentRunCutoff = processCutOffStr(setUp.runCutOffString_,currentReadCnt );
      		uint32_t currentQualRunCutoff = processCutOffStr(setUp.qualRunCutOffString_,currentReadCnt );
      		alignerObj.kMaps_.runCutOff_ = currentRunCutoff;
      		alignerObj.kMaps_.qualRunCutOff_ = currentQualRunCutoff;
      		//auto currentKMaps = kmerCalculator::indexKmerMpas(clusters, kClus.readPositions_,setUp.kLength_,
      			//	currentRunCutoff ,currentQualRunCutoff,setUp.expandKmerPos_, setUp.expandKmerSize_);
      		//alignerObj.setKmerMpas(currentKMaps);

      		collapserObj.runClustering(clusters, comp.readPositions_, iteratorMap, alignerObj);
      		alignerObj.kMaps_.runCutOff_ = oldRunCutoff;
      		alignerObj.kMaps_.qualRunCutOff_ = oldQualRunCutoff;
      		//alignerObj.setKmerMpas(oldMaps);
      	}
        //clusters = processedClusters;
    	}
    	collapserObj.verbose_ = oldVerbose;
    	clusters = readVecSplitter::splitVectorOnRemove(clusters).first;
    }
    for (uint32_t i = 1; i <= iteratorMap.size(); ++i) {
      std::cout << std::endl;
      uint32_t sizeOfReadVector = readVec::getReadVectorSize(clusters);
      runningParameters runPars(iteratorMap[i],i, sizeOfReadVector);
      runPars.printIterInfo(std::cout, true);
      if (condensedCollapse && runPars.errors_.hqMismatches_ == 0 &&
          runPars.errors_.lqMismatches_ == 0 &&
          runPars.errors_.largeBaseIndel_ == 0 &&
          runPars.errors_.twoBaseIndel_ == 0) {
        auto byCondensed = readVec::organizeByCondensedSeq(clusters);
      	if (collapserObj.regKmer_) {
      		for (auto & condensedReads : byCondensed){
      			collapserObj.collapseWithParametersRegKmer(condensedReads.second,
      			                                                  runPars, alignerObj);
      		}
        } else {
        	for (auto & condensedReads : byCondensed) {
        		collapserObj.collapseWithParametersQualKmer(condensedReads.second,
                                                      runPars, alignerObj);
        	}
        }
				clusters.clear();
				for (const auto& condensedReads : byCondensed) {
					addOtherVec(clusters, condensedReads.second);
				}
      } else {
        if (collapserObj.regKmer_) {

          collapserObj.collapseWithParametersRegKmer(clusters, runPars,
                                                     alignerObj);
        } else {
          collapserObj.collapseWithParametersQualKmer(clusters, runPars,
                                                      alignerObj);
        }
      }
      readVec::allUpdateName(clusters);
      readVecSorter::sortReadVector(clusters, sortBy);
      clusterVec::allCalculateConsensus(clusters, alignerObj, true);
      if (adjustHomopolyerRuns) {
        readVec::allAdjustHomopolymerRunsQualities(clusters);
      }
      std::cout << "Current duration: ";
      setUp.logRunTime(std::cout);
      setUp.rLog_ << "Current duration: "; setUp.rLog_.logTotalTime(6);
      setUp.rLog_ << "Current clusters size " << sizeOfReadVector
                      << "\n";
    }
  }
  clusters = readVecSplitter::splitVectorOnRemove(clusters).first;
  readVecSorter::sortReadVector(clusters, "totalCount");

  std::string seqName = bib::files::getFileName(setUp.ioOptions_.firstName_);
  std::vector<cluster> rejectedClusters;
  std::unordered_map<double, double> bestLikelihood;
  randomGenerator gen;
  simulation::errorProfile eProfile({'A', 'C', 'G', 'T'});
  if (sim) {
    if (markChimeras) {
      std::ofstream chimerasInfoFile;
      openTextFile(chimerasInfoFile,
                   setUp.directoryName_ + "chimeraNumberInfoBeforeSim.txt", ".txt", false,
                   false);
      chimerasInfoFile << "#chimericClusters\t#chimericReads" << std::endl;
      std::cout << "Marking chimeras" << std::endl;
      errorProfile chiOverlap;
      chiOverlap.oneBaseIndel_ = 2;
      chiOverlap.twoBaseIndel_ = 1;
      chiOverlap.lowKmerMismatches_ = 1;
      uint32_t overLapSizeCutoff = 5;
      uint32_t allowableError = 0;
      uint32_t chiCount = 0;
      std::vector<cluster> tempClusters;
      for (const auto& clus : clusters) {
        tempClusters.emplace_back(cluster(readObject(clus.seqBase_)));
      }
      readVec::allSetFractionByTotalCount(tempClusters);
      clusterCollapser::markChimerasAdvanced(
          tempClusters, alignerObj, parFreqs, 1, setUp.local_, chiOverlap,
          overLapSizeCutoff, setUp.weightHomopolymers_, chiCount, allowableError);

      int clusterCount = 0;
      int readCount = 0;
      readVec::getCountOfReadNameContaining(tempClusters, "CHI", clusterCount);
      readVec::getReadCountOfReadNameContaining(tempClusters, "CHI", readCount);
      chimerasInfoFile << getPercentageString(clusterCount, tempClusters.size())
                       << "\t"
                       << getPercentageString(
                              readCount, readVec::getTotalReadCount(tempClusters))
                       << std::endl;
      std::cout << "Marked " << clusterCount << " as chimeric before sim" << std::endl;
      reader.write(tempClusters, setUp.directoryName_ + "beforSimOutput.fastq",
                   setUp.ioOptions_);
      if (setUp.refFilename_ == "") {
        profiler::getFractionInfoCluster(tempClusters, setUp.directoryName_, "beforSimOutputInfo");
      } else {
        profiler::getFractionInfoCluster(tempClusters, setUp.directoryName_, "beforSimOutputInfo",
                                  setUp.refFilename_, alignerObj, setUp.local_,
                                  true);
      }

    }
    // create random generator

    // create R session
    ownRInside ownSession;
    // create error profile to get the error profile of the current clusters

    for (const auto& clus : clusters) {
      clus.updateErrorProfile(eProfile, alignerObj, setUp.local_);
    }
    letterCounter masterCounter(std::vector<char>{'A', 'C', 'G', 'T'});
    for_each(reader.reads, [&](const readObject& read) {
      masterCounter.increaseCountByString(read.seqBase_.seq_);
    });

    masterCounter.setFractions();
    table letCount(masterCounter.letters, {"let", "count"});
    table letFrac(masterCounter.fractions, {"let", "frac"});
    letCount.outPutContentOrganized(std::cout);
    letFrac.outPutContentOrganized(std::cout);
    // masterCounter.outPutACGTInfo(std::cout);
    // masterCounter.outPutACGTFractionInfo(std::cout);
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
    for_each(reader.reads, [&](const readObject& read) {
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
        getModelQualError(countsForModelBase, 0, 50, 1, "logit");

    // store the best area under the curve
    double bestArea = 0;
    std::string bestAreaPar = "";
    std::string bestPar = "";
    std::unordered_map<double, double> parLikelihood =
        createLikelihoodMap(basePredictions);
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
      simulation::timeTracker qualCounting("qualCounting");
      for (auto i : qualWindows) {
        for_each(reader.reads, [&](const readObject& read) {
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
              getModelQualError(cm.second, 0, 50, 1, "logit");
        } else if (cm.first == "mean" || cm.first == "median") {
          predictionsStep[cm.first] =
              getModelQualError(cm.second, 0, 50, 0.01, "logit");
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
            createLikelihoodMap(predictionsStep[qCount.first]);
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
          bestAreaPar = qCount.first + "_qw_" + to_string(i);
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
        "pdf(\"" + seqName +
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
      simulation::timeTracker t(clus.seqBase_.name_ + "_sim", true);
      // simulate the cluster using the qualities of the original sequences
      clus.simOnQual(identicalClusters, alignerObj, runTimes, bestLikelihood,
                     eProfile, gen);
      // remove clusters on fdr cut off
      //std::cout << "fdr1" << std::endl;
      if(keepCondensed){
      	clus.qualityClip = std::vector<uint32_t> (100);
      }
      clus.removeClustersOnFDR(alignerObj, fdrCutOff, rejectedClusters);
      //std::cout << "fdr2" << std::endl;
    }
    //std::cout << "?" << std::endl;
    //std::cout << "len of rejectedClusters" << len(rejectedClusters) << std::endl;
    /*for(const auto & rej : rejectedClusters){
    	std::cout << rej.seqBase_.name_ << std::endl;
    	std::cout << rej.seqBase_.seq_ << std::endl;
    	std::cout << rej.seqBase_.cnt_ << std::endl;
    }*/
    uint32_t quickCount = 0;
    //printDescriptionVec(rejectedClusters,std::cout, true);
    clusterVec::allCalculateConsensus(rejectedClusters, alignerObj, true);
    readVec::allUpdateName(rejectedClusters);
    for (auto& clus : rejectedClusters) {
    	// simulate the rejected reads so they can do p value comparison latter
    	++quickCount;
      std::cout << quickCount << ":" << len(rejectedClusters) << std::endl;

      if (clus.rejected_) {
        // std::cout <<"rejections" << std::endl;

        // if cluster is only one read no need to simulate
        if (clus.seqBase_.cnt_ < 2) {
          continue;
        }
        simulation::timeTracker t(clus.seqBase_.name_ + "_sim", true);
        // simulate the cluster using the qualities of the original sequences
        clus.simOnQual(identicalClusters, alignerObj, runTimes, bestLikelihood,
                       eProfile, gen);

        // remove clusters on fdr cut off
        //clus.removeClustersOnFDR(alignerObj, fdrCutOff, rejectedClusters);
      }
    }
  }
  //std::cout << "made it here " << std::endl;

  addOtherVec(clusters, rejectedClusters);
  readVec::allUpdateName(clusters);
  reader.write(clusters, setUp.directoryName_ + "clusWithRejected",
               setUp.ioOptions_.outFormat_, true, false);
  profiler::getFractionInfoCluster(clusters, setUp.directoryName_,
                            "clusWithRejectedInfo");

  readVecSorter::sort(clusters);
  renameReadNames(clusters, seqName, true, false);
  if (sim) {
  	if(debug){
      std::ofstream pValueFile;
      openTextFile(pValueFile, setUp.directoryName_ + "pValues.tab.txt", ".txt",
                   true, false);
      pValueFile << "read1\tread2\tpValue" << std::endl;
      for (const auto& readPos : iter::range(len(clusters))) {
        uint32_t realPos = len(clusters) - readPos - 1;
        for (const auto& readSubPos : iter::range<uint32_t>(0, realPos)) {
          if (clusters[readSubPos].seqBase_.cnt_ <=
              clusters[realPos].seqBase_.cnt_) {
            continue;
          }
          double currentPvalue = clusters[readSubPos].getPValue(
              clusters[realPos].seqBase_, alignerObj, setUp.local_);
          if (currentPvalue >= pValueCutOff) {
            pValueFile << clusters[realPos].seqBase_.name_ << "\t"
                       << clusters[readSubPos].seqBase_.name_ << "\t"
                       << currentPvalue << std::endl;
          }
        }
      }

  	}

    std::ofstream pValueBestFile;
    openTextFile(pValueBestFile, setUp.directoryName_ + "pValuesBest.tab.txt",
                 ".txt", true, false);
    pValueBestFile << "read1\tread2\tpValue" << std::endl;
    simulation::timeTracker pTimer("comparing p value");
    uint32_t pTenPer = .1 * len(clusters);
    if(pTenPer == 0){
    	pTenPer = 1;
    }
    std::ofstream allPValuesBeforeCollapse;
    openTextFile(allPValuesBeforeCollapse, setUp.directoryName_ + "allPValuesBeforeCollapse.tab.txt", ".txt", false, false);
    allPValuesBeforeCollapse << "clusters";
    for(const auto & pos : iter::range(len(clusters))){
    	allPValuesBeforeCollapse << "\t" << clusters[pos].seqBase_.name_;
    }
    allPValuesBeforeCollapse << std::endl;
    for(const auto & pos : iter::range(len(clusters))){
    	allPValuesBeforeCollapse << clusters[pos].seqBase_.name_;
    	for(const auto & subPos : iter::range<uint32_t>(0, pos + 1)){
    		allPValuesBeforeCollapse << std::setprecision(10);
    		allPValuesBeforeCollapse << "\t" << clusters[subPos].getPValue(clusters[pos].seqBase_, alignerObj, false);
    	}
    	allPValuesBeforeCollapse << std::endl;
    }
    for (const auto& readPos : iter::range(len(clusters))) {
    	if(readPos % pTenPer == 0){
    		std::cout << "On " << readPos << " of " << len(clusters) << " "; pTimer.print(pTimer.prefix_, std::cout, 0, 2);
    	}
      double bestPvalue = 0;
      uint32_t bestReadPos = std::numeric_limits<uint32_t>::max();
      uint32_t realPos = len(clusters) - readPos - 1;
      for (const auto& readSubPos : iter::range<uint32_t>(0, realPos)) {
        if (clusters[readSubPos].seqBase_.cnt_ <=
            clusters[realPos].seqBase_.cnt_ * 2) {
          continue;
        }
        double currentPvalue = clusters[readSubPos].getPValue(
            clusters[realPos].seqBase_, alignerObj, setUp.local_);
        if (currentPvalue >= pValueCutOff) {
          if (currentPvalue > bestPvalue) {
            bestPvalue = currentPvalue;
            bestReadPos = readSubPos;
          } else if (currentPvalue == bestPvalue) {
            if (clusters[realPos].seqBase_.cnt_ >
                clusters[bestReadPos].seqBase_.cnt_) {
              bestReadPos = readSubPos;
            }
          }
        }
      }
      if (bestReadPos != std::numeric_limits<uint32_t>::max()) {
        pValueBestFile << clusters[realPos].seqBase_.name_ << "\t"
                       << clusters[bestReadPos].seqBase_.name_ << "\t"
                       << bestPvalue << std::endl;
        clusters[bestReadPos].addRead(clusters[realPos]);
        clusters[realPos].remove = true;
      }
    }

    clusters = readVecSplitter::splitVectorOnRemove(clusters).first;
    std::vector<cluster> newRejected;
    for (auto& clus : clusters) {
      if(keepCondensed){
      	clus.qualityClip = std::vector<uint32_t> (100);
      }
      clus.removeClustersOnFDR(alignerObj, fdrCutOff, newRejected);
    }
    clusterVec::allCalculateConsensus(newRejected, alignerObj, true);
		readVec::allUpdateName(newRejected);
		for (auto& clus : newRejected) {
			// simulate the rejected reads so they can do p value comparison latter
			//++quickCount;
			//std::cout << quickCount << ":" << len(rejectedClusters) << std::endl;

			if (clus.rejected_) {
				// std::cout <<"rejections" << std::endl;

				// if cluster is only one read no need to simulate
				if (clus.seqBase_.cnt_ < 2) {
					continue;
				}
				simulation::timeTracker t(clus.seqBase_.name_ + "_sim", true);
				// simulate the cluster using the qualities of the original sequences
				clus.simOnQual(identicalClusters, alignerObj, runTimes, bestLikelihood,
											 eProfile, gen);

				// remove clusters on fdr cut off
				//clus.removeClustersOnFDR(alignerObj, fdrCutOff, rejectedClusters);
			}
		}
    addOtherVec(clusters, newRejected);
    readVecSorter::sort(clusters);
    renameReadNames(clusters, seqName, true, false);
    std::ofstream allPValues;
    openTextFile(allPValues, setUp.directoryName_ + "allPvalues.tab.txt", ".txt", false, false);
    std::ofstream allComparisonFile;
    openTextFile(allComparisonFile, setUp.directoryName_ + "allComparison.tab.txt", ".txt", false, false);
    allPValues << "clusters";

    allComparisonFile << "clusters";
    for(const auto & pos : iter::range(len(clusters))){
    	allPValues << "\t" << clusters[pos].seqBase_.name_;
    	allComparisonFile << "\t" << clusters[pos].seqBase_.name_;
    }
    allPValues << std::endl;
    for(const auto & pos : iter::range(len(clusters))){
    	allPValues << clusters[pos].seqBase_.name_;

    	allComparisonFile << clusters[pos].seqBase_.name_;
    	for(const auto & subPos : iter::range<uint32_t>(0, pos + 1)){
    		allPValues << std::setprecision(10);
    		allPValues << "\t" << clusters[subPos].getPValue(clusters[pos].seqBase_, alignerObj, false);
    		allComparisonFile << "\t" << to_string(alignerObj.errors_.oneBaseIndel_) + "," << to_string(alignerObj.errors_.twoBaseIndel_)
    				<< "," << to_string(alignerObj.errors_.largeBaseIndel_) << "," << to_string(alignerObj.errors_.lqMismatches_)
    				<< "," << to_string(alignerObj.errors_.hqMismatches_);
    	}
    	allPValues << std::endl;
    	allComparisonFile << std::endl;
    }
  }

  if (markChimeras) {
    std::ofstream chimerasInfoFile;
    openTextFile(chimerasInfoFile,
                 setUp.directoryName_ + "chimeraNumberInfo.txt", ".txt", false,
                 false);
    chimerasInfoFile << "#chimericClusters\t#chimericReads" << std::endl;
    std::cout << "Marking chimeras" << std::endl;
    errorProfile chiOverlap;
    chiOverlap.oneBaseIndel_ = 2;
    chiOverlap.twoBaseIndel_ = 1;
    chiOverlap.lowKmerMismatches_ = 1;
    uint32_t overLapSizeCutoff = 5;
    uint32_t allowableError = 0;
    uint32_t chiCount = 0;
    std::vector<cluster> tempClusters;
    for (const auto& clus : clusters) {
      tempClusters.emplace_back(cluster(readObject(clus.seqBase_)));
    }
    readVec::allSetFractionByTotalCount(tempClusters);
    clusterCollapser::markChimerasAdvanced(
        tempClusters, alignerObj, parFreqs, 1, setUp.local_, chiOverlap,
        overLapSizeCutoff, setUp.weightHomopolymers_, chiCount, allowableError);

    int clusterCount = 0;
    int readCount = 0;
    readVec::getCountOfReadNameContaining(tempClusters, "CHI", clusterCount);
    readVec::getReadCountOfReadNameContaining(tempClusters, "CHI", readCount);
    chimerasInfoFile << getPercentageString(clusterCount, tempClusters.size())
                     << "\t"
                     << getPercentageString(
                            readCount, readVec::getTotalReadCount(tempClusters))
                     << std::endl;
    std::cout << "Marked " << clusterCount << " as chimeric" << std::endl;
    for (const auto& clusNum : iter::range(clusters.size())) {
      clusters[clusNum].seqBase_.name_ = tempClusters[clusNum].seqBase_.name_;
    }
  }

  if (collapsingTandems) {
    reader.write(clusters, setUp.directoryName_ + setUp.ioOptions_.outFilename_ +
                               "_befroeTanCol",
                 setUp.ioOptions_);
    std::cout << "Collapsing on tandem repeat gaps" << std::endl;
    std::cout << "Starting with " << clusters.size() << " clusters"
              << std::endl;
    clusterCollapser::collapseTandems(clusters, alignerObj, setUp.runCutoff_,
                                      setUp.kLength_, setUp.kmersByPosition_,
                                      parFreqs, setUp.local_, true);
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
    VecStr toks = tokenizeString(clus.firstReadName, "_t");
    startingInfo << clus.seqBase_.name_ << "\t" << clus.firstReadName << "\t"
                 << toks.back() << std::endl;
  }
  if (additionalOut) {
    std::string additionalOutDir = findAdditonalOutLocation(
        additionalOutLocationFile, setUp.ioOptions_.firstName_);
    reader.write(clusters, additionalOutDir + setUp.ioOptions_.outFilename_,
                 setUp.ioOptions_);
  }

  reader.write(clusters, setUp.directoryName_ + setUp.ioOptions_.outFilename_,
               setUp.ioOptions_);
  if(createMinTree){
    auto clusSplit = readVecSplitter::splitVectorOnReadFraction(clusters, 0.02);
    auto minTreeJson = getMinMismatchTreeJson(clusSplit.first, alignerObj, 1, setUp.weightHomopolymers_,
    		hueStart, hueStop, lumStart, lumStop, satStart, satStop, backgroundColor);
  	std::ofstream outFile;
  	openTextFile(outFile, setUp.directoryName_ + "psudoMinTree.json", "json", false, false);
  	outFile << minTreeJson << std::endl;
  	std::string htmlOut = genHtmlStrForPsuedoMintree("psudoMinTree.json");
  	std::ofstream outHtmlFile;
  	openTextFile(outFile, setUp.directoryName_ + "psudoMinTree.html", "html", false, false);
  	outHtmlFile << htmlOut << std::endl;
  }

  std::string clusterDirectoryName =
      bib::files::makeDir(setUp.directoryName_, "clusters");
  clusterVec::allWriteOutClusters(clusters, clusterDirectoryName);

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
      reader.write(collection, allInputReadsDir + clus.seqBase_.name_,
                   "fastaQual", false, false);
    }
  }
  if (extra) {
    std::string compInfo =
        bib::files::makeDir(setUp.directoryName_, "identicalCompositionInfos");
    profiler::allOutputIdenticalReadComp(clusters, compInfo);

    clusterVec::outPutAlignments(clusters, alignerObj, setUp.directoryName_);
    clusterVec::allWriteOutLongestAlignments(clusters, setUp.directoryName_);
    std::string alnProfileDir =
        bib::files::makeDir(setUp.directoryName_, "alnProfiles");
    for (const auto& iter : clusters) {
      iter.alignmentProfile(alnProfileDir, alignerObj, setUp.local_,
                            kmerCheckingOnAlignProfile, setUp.kLength_,
                            setUp.kmersByPosition_, true);
    }
    std::string inputClustersDir =
        bib::files::makeDir(setUp.directoryName_, "inputClusters");
    clusterVec::allWriteAllInputClusters(clusters, inputClustersDir);
    if (!setUp.ioOptions_.processed_) {
      std::string allAlnProfileDir =
          bib::files::makeDir(setUp.directoryName_, "allAlnProfiles");

      for (const auto& clus : clusters) {
        std::vector<readObject> collection =
            clus.getAllBeginingClusters(identicalClusters);
        /*for (auto& colIter : collection) {
          colIter.seqBase_.name_ = colIter.getStubName();
        }*/
        readObject tempObject(seqInfo(clus.seqBase_.name_, clus.seqBase_.seq_,
                                      clus.seqBase_.qual_));
        readVec::allSetFractionByTotalCount(collection);
        alignmentProfiler::getAlignmentInformationForReferenceRawReads(
            collection, std::vector<readObject>(1, tempObject), alignerObj,
            setUp.local_, allAlnProfileDir + clus.seqBase_.name_ + "_profile",
            setUp.kLength_, true);
      }
    }
    // std::string profileAllDir=makeDirectory(setUp.directoryName_,
    // "profilesAll");
    // allProfileAllClusters(output,profileAllDir, alignerObj, setUp.local_,
    // setUp.simple_, kmerMapPos, runCutOff,kLength);
  }

  if (setUp.writingOutAlnInfo_) {
  	simulation::timeTracker writingTimer("Writing aln infos", true);
    alignerObj.alnHolder_.write(setUp.outAlnInfoDirName_);
  }

  std::cout << alignerObj.numberOfAlingmentsDone_ << std::endl;
  setUp.rLog_ << alignerObj.numberOfAlingmentsDone_ << "\n";
  setUp.logRunTime(std::cout);

  return 0;
}


int SeekDeepRunner::processClusters(MapStrStr inputCommands) {
  // parameters
  SeekDeepSetUp setUp(inputCommands);
  bool population = false;

  bool checkChimeras = false;
  int parFreqs = 2;

  std::string parameters = "";
  // clusterSize cutoff
  int cutOff = 1;

  bool extra = false;
  bool condensedCollapse = false;
  double fracCutoff = 0.005;
  uint32_t runsRequired = 0;

  bool smallestFirst = true;
  bool bestMatch = false;
  int bestMatchCheck = 10;
  // bool fractionsByReads=false;
  bool keepChimeras = false;
  std::string experimentName = "PopUID";
  // bool preCollapse=false;
  std::string parametersPopulation = "";
  bool differentPar = false;
  bool popBoth = false;
  std::string sortBy = "fraction";
  std::map<int, std::vector<double>> popIteratorMap;
  std::map<int, std::vector<double>> iteratorMap;
  // bool containsRecurrence = false;

  bool skipOnLetterCounterDifference = false;
  double fractionDifferenceCutOff = 0.05;
  bool grayScale = false;
  double sat = 0.99;
  double lum = 0.5;
  bool regKmers = true;
  bool eventBasedRef = false;
  setUp.setOption(eventBasedRef, "-eventBasedRef", "Do Event Based Ref Count");
  setUp.setBoolOptionFalse(regKmers, "-qualKmers", "doQualKmer");
  //setUp.setOption(regKmers, "-regKmers,-regK", "regularKmerAnalysis");
  setUp.setOption(grayScale, "-gray", "grayScale");
  setUp.setOption(sat, "-sat", "sat");
  setUp.setOption(lum, "-lum", "lum");
  setUp.setOption(skipOnLetterCounterDifference, "-skip",
                  "skipOnLetterCounterDifference");
  setUp.setOption(fractionDifferenceCutOff, "-skipCutOff",
                  "fractionDifferenceCutOff");
  setUp.setOption(condensedCollapse, "-condensedCollapse", "condensedCollapse");

  std::string previousPopFilename = "";
  setUp.setOption(previousPopFilename, "-previousPop", "previousPopFilename");

  setUp.setUpMultipleSampleCluster(
      parameters, extra, cutOff, iteratorMap, population, fracCutoff,
      smallestFirst, bestMatch, bestMatchCheck, checkChimeras, parFreqs,
      parametersPopulation, differentPar, popIteratorMap, popBoth, keepChimeras,
      experimentName, runsRequired);
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
    readVec::getMaxLength(reader.reads, maxSize);
    VecStr toks = tokenizeString(strIter, "/");
    if (toks[2].find("R") != std::string::npos) {
      // containsRecurrence = true;
    }
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
      bestMatch, bestMatchCheck, setUp.local_, setUp.checkKmers_,
      setUp.kmersByPosition_, setUp.runCutoff_, setUp.qualRunCutoff_,
      setUp.kLength_, setUp.verbose_, smallestFirst, condensedCollapse,
      setUp.weightHomopolymers_, skipOnLetterCounterDifference,
      fractionDifferenceCutOff, regKmers, setUp.adjustHomopolyerRuns_);
  // collect all the reads together for each sample
  std::map<std::string, std::vector<std::vector<sampleCluster>>>
      clustersBySample;
  std::map<std::string, std::vector<sampleCluster>>::iterator clusMapIter;
  std::map<std::string, std::vector<sampleCluster>> samplesExcluded;

  for (auto& readsIter : outputReads) {
    std::vector<cluster> clusters =
        baseCluster::convertVectorToClusterVector<cluster>(readsIter.second);
    readVecSorter::sortReadVector(clusters, "totalCount");
    // consider adding the sample name in the name as well
    renameReadNamesNewClusters(clusters, readsIter.first.second, true, true, false);
    if (checkChimeras) {
      // readVec::allSetFractionByTotalCount(clusters);
      clusterVec::allSetFractionClusters(clusters);
    	errorProfile chiOverlap;
    	chiOverlap.oneBaseIndel_ = 2;
    	chiOverlap.twoBaseIndel_ = 1;
    	chiOverlap.lowKmerMismatches_ = 1;
    	uint32_t overLapSizeCutoff = 5;
    	uint32_t allowableError = 0;
    	uint32_t chiCount = 0;
    	clusterCollapser::markChimerasAdvanced(
    			clusters, alignerObj, parFreqs, 1, setUp.local_, chiOverlap,
    			overLapSizeCutoff, setUp.weightHomopolymers_, chiCount, allowableError);
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
    auto currentSampClusters =
        baseCluster::convertVectorToClusterVector<sampleCluster>(clusters);
    clustersBySample[readsIter.first.first].push_back(currentSampClusters);
  }
  std::map<std::string, collapse::sampleCollapse> sampleCollapses;
  collapserObj.checkKmers_ = false;
  for (const auto& samp : clustersBySample) {
    std::cout << "Currently on : " << samp.first << std::endl;
    sampleCollapses[samp.first] =
        collapse::sampleCollapse(samp.second, samp.first, cutOff);
    // std::cout << "made sampleCollapse fine" << std::endl;
    // sampleCollapses[samp.first].updateInitialInfos(false);
    sampleCollapses[samp.first]
        .cluster(collapserObj, iteratorMap, sortBy, alignerObj);
    // std::cout << "clustered fine" << std::endl;
    if (checkingExpected) {
      sampleCollapses[samp.first].collapsed_.checkAgainstExpected(
          expectedSeqs, alignerObj, setUp.local_, setUp.weightHomopolymers_);
    }

    sampleCollapses[samp.first].updateCollapsedInfos(true);
    // std::cout << "updated collapse infos fine" << std::endl;
    sampleCollapses[samp.first].updateExclusionInfos(true);
    // std::cout << "updated excusesion fine" << std::endl;
    sampleCollapses[samp.first].renameClusters(sortBy);
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
  std::string finalDir = bib::files::makeDir(setUp.directoryName_, "final");
  for (auto& sampCollapse : sampleCollapses) {
    // exclude
    sampCollapse.second.excludeChimeras(false);
    if (0 != runsRequired) {
      sampCollapse.second.excludeBySampNum(runsRequired, false);
    } else {
      sampCollapse.second.excludeBySampNum(
          sampCollapse.second.input_.infos_.size(), false);
    }

    sampCollapse.second.excludeFraction(fracCutoff, true);
    std::string sortBy = "fraction";
    sampCollapse.second.renameClusters(sortBy);
    // write
    sampCollapse.second.writeInitial(originalsDir, setUp.ioOptions_.outFormat_,
                                     setUp.ioOptions_.overWriteFile_,
                                     setUp.ioOptions_.exitOnFailureToWrite_);
    sampCollapse.second.writeExcluded(excludedDir, setUp.ioOptions_.outFormat_,
                                      setUp.ioOptions_.overWriteFile_,
                                      setUp.ioOptions_.exitOnFailureToWrite_);
    sampCollapse.second.writeFinal(finalDir, setUp.ioOptions_.outFormat_,
                                   setUp.ioOptions_.overWriteFile_,
                                   setUp.ioOptions_.exitOnFailureToWrite_);
    sampCollapse.second.writeFinalOrginals(initialDir);
    // add to all sample cluster

    addOtherVec(allSamples, sampCollapse.second.createOutput(false, sortBy));
    // sampCollapse.second.updateCollapsedInfos(true);
    // sampCollapse.second.updateExclusionInfos(true);
  }
  popCollapse = collapse::populationCollapse(allSamples, experimentName);
  if (population) {
    std::string popDir = bib::files::makeDir(setUp.directoryName_, "population");
    if (differentPar) {
      popCollapse.popCluster(collapserObj, popIteratorMap, "fraction",
                             alignerObj);
    } else {
      popCollapse.popCluster(collapserObj, iteratorMap, "fraction", alignerObj);
    }
    readObjectIO popReader;
    if (previousPopFilename != "") {

      std::string format = bib::files::getExtension(previousPopFilename);
      popReader.read(format, previousPopFilename);

      popCollapse.renameToOtherPopNames(popReader.reads, alignerObj);
    }
    if (checkingExpected) {
      popCollapse.collapsed_.checkAgainstExpected(
          expectedSeqs, alignerObj, setUp.local_, setUp.weightHomopolymers_);
    }
    infoPrinter::printPopulationCollapseInfo(
        popCollapse, popDir + "populationCluster.tab.txt", checkingExpected,
        true);
    std::string popInitialDir = bib::files::makeDir(popDir, "initial");
    popCollapse.writeFinalInitial(popInitialDir);
    popCollapse.writeFinal(popDir, setUp.ioOptions_.outFormat_,
                           setUp.ioOptions_.overWriteFile_,
                           setUp.ioOptions_.exitOnFailureToWrite_);
    std::unordered_map<std::string, bib::color> colorsForGraph;
    if (previousPopFilename != "") {
      colorsForGraph = getColorsForNames(popReader.reads, sat, lum);
    } else {
      colorsForGraph =
          getColorsForNames(popCollapse.collapsed_.clusters_, sat, lum);
    }

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
      // currentGraphDist.createDotBestConnectedFile(dotDir, colorsForGraph,
      // true);
      currentGraphDist.createDotBestConnectedFile(dotDir, colorsForGraph,
                                                  false);
    }
  }

  infoPrinter::printSampleCollapseInfo(
      sampleCollapses, checkingExpected,
      setUp.directoryName_ + "selectedClustersInfo.tab.txt", popCollapse,
      population);

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
  setUp.setUpMakeSampleDirectories(sampleNameFilename);
  makeSampleDirectoriesWithSubDirectories(sampleNameFilename,
                                          setUp.directoryName_);
  setUp.startARunLog(setUp.directoryName_);
  setUp.logRunTime(std::cout);
  return 0;
}


}  // namespace bib
