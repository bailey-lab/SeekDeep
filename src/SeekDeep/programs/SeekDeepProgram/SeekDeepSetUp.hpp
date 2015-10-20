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
  std::string sampleName = "";
  std::string compareSeqFilename = "";
  bool multipleTargets = false;

};

struct clusterDownPars {
  std::string parameters = "";
  bool snapShots = false;
  std::string qualRep = "median", sortBy = "totalCount";
  bool markChimeras = false;
  double parFreqs = 2;
  bool collapsingTandems = false;
  bool additionalOut = false;
  std::string additionalOutLocationFile = "";
  std::map<int, std::vector<double>> iteratorMap;
  bool removeLowQualBases = false;
  int lowQualityCutOff = 3;
  uint32_t runTimes = 100;
  bool sim = false;
  bool printSimClusters = false;
  double pValueCutOff = 0.01;
  double fdrCutOff = 0.01;
  bool startWithSingles = false;
  bool createMinTree = false;
  std::string diffCutOffStr = "0.1";
  bool findBest = true;
  double compPerCutOff = .98;
  bool useCompPerCutOff = false;
  bool ionTorrent = false;
  bool useNucComp = false;
  bool useKmerBinning = false;
  double kmerCutOff = 0.80;
  uint32_t kCompareLen = 10;
  bool leaveOutSinglets = false;
  bool onPerId = false;
	bool extra = false;
	uint32_t smallReadSize = 20;
  bool useMinLenNucComp = false;
  std::vector<double> diffCutOffVec;
  bool noAlign_ = false;
};

struct processClustersPars {
  bool noPopulation = false;
  std::string previousPopFilename = "";
  std::string parameters = "";
  uint32_t clusterCutOff = 1;
  bool extra = false;
  double fracCutoff = 0.005;
  uint32_t runsRequired = 0;
  bool checkChimeras = false;
  double parFreqs = 2;
  bool keepChimeras = false;
  bool investigateChimeras = false;
  bool recheckChimeras = false;
  double chiCutOff = .40;
  std::string experimentName = "PopUID";
  std::string parametersPopulation = "";
  bool differentPar = false;
  bool popBoth = false;
  std::string sortBy = "fraction";
  std::map<int, std::vector<double>> popIteratorMap;
  std::map<int, std::vector<double>> iteratorMap;
  bool grayScale = false;
  double sat = 0.99;
  double lum = 0.5;
  bool eventBasedRef = false;
  bool writeExcludedOriginals = false;
  bool ionTorrent = false;
  bool removeLowQualBases = false;
  uint32_t lowQualityCutOff = 3;
  std::string customCutOffs = "";
  std::string groupingsFile = "";
  bool onPerId = false;
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
	void setUpClusterDown(clusterDownPars & pars);
	void setUpMultipleSampleCluster(processClustersPars & pars);
	void setUpMakeSampleDirectories(std::string& sampleNameFilename);
};
}  // namespace bib

#ifndef NOT_HEADER_ONLY
#include "SeekDeepSetUp.cpp"
#endif
