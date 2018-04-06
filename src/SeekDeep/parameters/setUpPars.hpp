#pragma once

/*
 * setUpPars.hpp
 *
 *  Created on: Nov 25, 2016
 *      Author: nick
 */
//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include <bibseq.h>
#include "SeekDeep/objects/PairedReadProcessor.hpp"
#include "SeekDeep/objects/PrimersAndMids.hpp"

namespace bibseq {

struct CoreExtractorPars{
	std::string idFileDelim = "whitespace";


  uint32_t smallFragmentCutoff = 50;
  bool rename = false;
  QualFilteringPars qPars_;
  uint32_t numberOfNs = 1;


  MidDeterminator::MidDeterminePars mDetPars;
  PrimerDeterminator::PrimerDeterminatorPars pDetPars;
  bool noPrimers_{false};
  PrimersAndMids::InitPars primIdsPars;

  std::string sampleName = "";

  bool keepUnfilteredReads = false;

  void setCorePars(seqSetUp & setUp);

};

struct extractorPars{

	extractorPars();


	CoreExtractorPars corePars_;

  uint32_t minLen = std::numeric_limits<uint32_t>::max();
  uint32_t maxLength = std::numeric_limits<uint32_t>::max();


  bool illumina = false;


  uint32_t smallExtractReadCount = 5;
  bool filterOffSmallReadCounts = false;

	uint32_t trimAtQualCutOff = 2;
	bool trimAtQual = false;
	bool qualWindowTrim = false;

};



struct ExtractorPairedEndPars{

	ExtractorPairedEndPars();
	CoreExtractorPars corePars_;


  uint32_t r1Trim_ = 1;
  uint32_t r2Trim_ = 1;

  PairedReadProcessor::ProcessParams pairProcessorParams_;

};

struct clusterDownPars {

	std::string qualRep = "median";
	std::string sortBy = "totalCount";

	bool collapsingTandems = false;
	bool additionalOut = false;
	std::string additionalOutLocationFile = "";

	bfs::path initalParsFnp = "";
	std::string parameters = "";
	std::string binParameters = "";
	CollapseIterations intialParameters;
	CollapseIterations iteratorMap;
	CollapseIterations binIteratorMap;


  bool startWithSingles = false;
  bool leaveOutSinglets = false;
  bool mapBackSinglets = false;
  uint32_t singletCutOff = 1;

  bool createMinTree = false;
  std::string diffCutOffStr = "0.1";

  bool writeOutFinalInternalSnps = false;

  double compPerCutOff = .98;
  bool useCompPerCutOff = false;
  bool ionTorrent = false;
  bool illumina = false;
  bool tech454 = false;
  uint32_t hq = 0;

  double otuPerc = .99;
  bool onHqPerId = false;
  bool onPerId = false;
	bool extra = false;
	uint32_t smallReadSize = 20;

	bool writeOutInitalSeqs = false;

	SnapShotsOpts snapShotsOpts_;
};

struct processClustersPars {
	bfs::path masterDir = ".";
  bool noPopulation = false;
  std::string previousPopFilename = "";
  comparison previousPopErrors;

  uint32_t numThreads = 1;

  std::string parameters = "";
  std::string binParameters = "";

  uint32_t clusterCutOff = 1;
  bool extra = false;
  double fracCutoff = 0.005;
  uint32_t runsRequired = 0;
  bool fracExcludeOnlyInFinalAverageFrac = false;

  bool collapseLowFreqOneOffs = false;
  double lowFreqMultiplier = 30;

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
  bool eventBasedRef = false;
  bool writeExcludedOriginals = false;

  bool illumina = false;

  bool ionTorrent = false;
  bool removeLowQualBases = false;
  uint32_t lowQualityCutOff = 3;
  std::string customCutOffs = "";
  std::string groupingsFile = "";
  bool onPerId = false;
  bool plotRepAgreement = false;

  bool noTrees = false;

	bool noErrorsSet = false;
	bool strictErrorsSet = false;
	bool strictErrorsSetHq1 = false;
	uint32_t hqMismatches = 0;
	uint32_t stopAfter = 100;

	uint32_t sampleMinTotalReadCutOff = 0;

};


struct makeSampleDirectoriesPars{
	bfs::path sampleNameFilename;
	bool separatedDirs = false;
};


}  // namespace bibseq




