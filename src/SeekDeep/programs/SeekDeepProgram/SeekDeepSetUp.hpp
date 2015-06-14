#pragma once
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
//
//  SeekDeepSetUp.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/24/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include <bibseq.h>
#include <bibseq/programUtils/seqSetUp.hpp>

namespace bibseq {

struct extractorPars{
	extractorPars(){
	  rPrimerErrors.hqMismatches_ = 4;
	  rPrimerErrors.distances_.queryCoverage_ = .50;
	  rPrimerErrors.largeBaseIndel_ = .99;
	  rPrimerErrors.oneBaseIndel_ = 2;
	  rPrimerErrors.twoBaseIndel_ = 1;

	  fPrimerErrors.hqMismatches_ = 2;
	  fPrimerErrors.distances_.queryCoverage_ = 1;
	  fPrimerErrors.largeBaseIndel_ = .99;
	  fPrimerErrors.oneBaseIndel_ = 2;
	  fPrimerErrors.twoBaseIndel_ = 1;
	}
	std::string idFilename = "";
	std::string idFileDelim = "whitespace";
  int numberOfNs = 1;

  bool screenForPossibleContamination = false;
  std::string compareSeq = "";
  bool qualWindowTrim = false;
  bool multiplex = false;
  int qualityWindowLength = 0;
  int qualityWindowStep = 0;
  int qualityWindowThres = 0;
  int minLen = 200;
  int maxLength = 300;

  bool checkComplement = false;
  int smallFragmentCutoff = 50;

  bool HMP = false;
  uint32_t primerLen = 30;
  bool trimTcag = false;

  uint32_t barcodeErrors = 0;
  bool rename = false;

  bool noReversePrimer = false;
  bool reversePrimerToUpperCase = false;
  comparison rPrimerErrors;

  bool noForwardPrimer = false;
  bool forwardPrimerToUpperCase = false;
  comparison fPrimerErrors;


  double kmerCutOff = .20;
  uint32_t contaminationKLen = 7;

  bool barcodesBothEnds = false;
  bool variableStart = false;
  uint32_t variableStop = 50;
  uint32_t qualCheck = 25;
  bool checkingQCheck = false;
  double qualCheckCutOff = 0.75;
};

class SeekDeepSetUp : public seqSetUp {

 public:
  // constructors
  SeekDeepSetUp(int argc, char* argv[]) : seqSetUp(argc, argv) {}
  SeekDeepSetUp(const bib::progutils::commandLineArguments& inputCommands)
      : seqSetUp(inputCommands) {}
  SeekDeepSetUp(const MapStrStr& inputCommands)
      : seqSetUp(inputCommands) {}

  void setUpExtractor(std::string& idFileName, bool& multiplex, bool& condensed,
                      int& minLen, int& maxLength, int& within,
                      int& qualityWindowLength, int& qualityWindowStep,
                      int& qualityWindowThres, bool& unknown,
                      int& unknownPrimerSize, int& unknownMinSize,
                      bool& findReversePrimer, double& queryCoverageCutoff,
                      double& percentIdentityCutoff, int& numberOfNs,
                      bool& pyroExtract, bool& reversePrimerToLowerCase,
                      bool& flowFiltering, int& maxFlows, bool& checkComplement,
                      bool& screenForPossibleContaimination,
                      std::string& compareSeq, std::string& idFileDelim,
                      int& smallFragmentCutOff, bool& qualWindowTrim);
  void setUpExtractor(extractorPars & pars);

  void setUpClusterDown(std::string& qualRep, std::string& parameters,
                        bool& extra,
                        std::map<int, std::vector<double>>& iteratorMap,
                        bool& smallestFirst, bool& markChimeras, double& parFreqs,
                        bool& bestMatch, int& bestMatchCheck, bool& snapShots,
                        std::string& sortBy, bool& additionalOut,
                        std::string& additonalOutLocationFile,
                        bool& collapsingTandems,
                        bool& kmerCheckingOnAlignProfile,
                        bool& condesnedCollpase, bool& removeLowQualBases,
                        int& lowQualityCutOff, bool& adjustHomopolyerRuns,
												bool & onPerId);
  void setUpMultipleSampleCluster(
      std::string& parameters, bool& extra, int& cutOff,
      std::map<int, std::vector<double>>& iteratorMap, bool& population,
      double& fracCutoff, bool& smallestFirst, bool& bestMatch,
      int& bestMatchCheck, bool& checkChimeras, double& parFreqs,
      std::string& parametersPopulation, bool& differentPar,
      std::map<int, std::vector<double>>& popIteratorMap, bool& popBoth,
      bool& keepChimeras, std::string& experimentName, uint32_t& runsRequired, bool & onPerId);
  void setUpCollapseTandems(double& freqCutoff, bool& extra,
                            bool& additionalOut,
                            std::string& additionalOutLocationFile);
  void setUpMarkChimeras(double& parentFreqs);
  void setUpMakeSampleDirectories(std::string& sampleNameFilename);
  void setUpPopulationClustering(
      std::string& directory, VecStr& contains, bool& specific, bool& recursive,
      std::string& parameters, bool& extra, int& cutOff,
      std::map<int, std::vector<double>>& popIteratorMap, double& fracCutoff,
      bool& smallestFirst, bool& bestMatch, int& bestMatchCheck,
      bool& checkChimeras, int& parFreqs);
  void setUpCompareTwoReplicates(
      std::string& filename1, std::string& qualFilename1, std::string& format1,
      std::string& filename2, std::string& qualFilename2, std::string& format2,
      std::string& parameters, bool& extra, int& cutOff,
      std::map<int, std::vector<double>>& iteratorMap, bool& smallestFirst,
      bool& markChimeras, int& parFreq, bool& bestMatch, int& bestMatchCheck);
};
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "SeekDeepSetUp.cpp"
#endif
