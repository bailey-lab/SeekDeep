#pragma once
//
//  SeekDeepSetUp.hpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 10/24/13.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//
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
#include <bibseq.h>
#include <bibseq/programUtils/seqSetUp.hpp>

namespace bibseq {

struct extractorPars{
	extractorPars(){
	  rPrimerErrors.hqMismatches_ = 4;
	  rPrimerErrors.distances_.query_.coverage_ = .50;
	  rPrimerErrors.largeBaseIndel_ = .99;
	  rPrimerErrors.oneBaseIndel_ = 2;
	  rPrimerErrors.twoBaseIndel_ = 1;

	  fPrimerErrors.hqMismatches_ = 2;
	  fPrimerErrors.distances_.query_.coverage_ = 1;
	  fPrimerErrors.largeBaseIndel_ = .99;
	  fPrimerErrors.oneBaseIndel_ = 2;
	  fPrimerErrors.twoBaseIndel_ = 1;
	}
	std::string idFilename = "";
	std::string idFileDelim = "whitespace";
	uint32_t numberOfNs = 1;

  bool screenForPossibleContamination = false;
  std::string compareSeq = "";
  bool contaminationMutlipleCompare = false;
  bool qualWindowTrim = false;
  bool multiplex = false;
  uint32_t qualityWindowLength = 0;
  uint32_t qualityWindowStep = 0;
  uint32_t qualityWindowThres = 0;
  uint32_t minLen = 200;
  uint32_t maxLength = 300;

  bool checkComplement = false;
  uint32_t smallFragmentCutoff = 50;

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
  uint32_t qualCheck = 30;
  bool checkingQCheck = false;
  double qualCheckCutOff = 0.75;
  std::string sampleName = "";
  std::string compareSeqFilename = "";
  bool multipleTargets = false;

  std::string multipleLenCutOffFilename = "";

  uint32_t smallExtractReadCount = 5;
  bool filterOffSmallReadCounts = false;

	bool mothurExtract = false;
	bool pyroExtract = false;
	uint32_t maxFlowCutoff = 0;

	uint32_t trimAtQualCutOff = 2;
	bool trimAtQual = false;

};

struct clusterDownPars {
	std::string parameters = "";
	std::string binParameters = "";
	std::string qualRep = "median";
	std::string sortBy = "totalCount";

	bool collapsingTandems = false;
	bool additionalOut = false;
	std::string additionalOutLocationFile = "";
	CollapseIterations iteratorMap;
	CollapseIterations binIteratorMap;


  bool startWithSingles = false;
  bool leaveOutSinglets = false;
  bool mapBackSinglets = false;
  uint32_t singletCutOff = 1;

  bool createMinTree = false;
  std::string diffCutOffStr = "0.1";

  double compPerCutOff = .98;
  bool useCompPerCutOff = false;
  bool ionTorrent = false;
  bool illumina = false;


  bool onPerId = false;
	bool extra = false;
	uint32_t smallReadSize = 20;

	SnapShotsOpts snapShotsOpts_;
};

struct processClustersPars {
	std::string masterDir = ".";
  bool noPopulation = false;
  std::string previousPopFilename = "";
  std::string parameters = "";
  std::string binParameters = "";

  uint32_t clusterCutOff = 1;
  bool extra = false;
  double fracCutoff = 0.005;
  uint32_t runsRequired = 0;

  bool keepChimeras = false;
  bool investigateChimeras = false;
  bool recheckChimeras = false;
  double chiCutOff = .40;
  std::string experimentName = "PopUID";
  std::string parametersPopulation = "";
  bool differentPar = false;
  bool popBoth = false;
  std::string sortBy = "fraction";
  CollapseIterations popIteratorMap;
  CollapseIterations iteratorMap;
  CollapseIterations binIteratorMap;
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
  bool plotRepAgreement = false;

  bool noTrees = false;
};

class SeekDeepSetUp : public seqSetUp {

 public:
	using seqSetUp::seqSetUp;

	void setUpExtractor(extractorPars & pars);
	void setUpClusterDown(clusterDownPars & pars);
	void setUpMultipleSampleCluster(processClustersPars & pars);
	void setUpMakeSampleDirectories(std::string& sampleNameFilename);
};
}  // namespace bibseq


