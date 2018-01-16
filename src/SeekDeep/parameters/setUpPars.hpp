#pragma once

/*
 * setUpPars.hpp
 *
 *  Created on: Nov 25, 2016
 *      Author: nick
 */

#include <bibseq.h>
#include "SeekDeep/objects/PairedReadProcessor.hpp"


namespace bibseq {

struct extractorPars{

	extractorPars();

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
  uint32_t minLen = std::numeric_limits<uint32_t>::max();
  uint32_t maxLength = std::numeric_limits<uint32_t>::max();

  uint32_t r1MinLen = std::numeric_limits<uint32_t>::max();
  uint32_t r2MinLen = std::numeric_limits<uint32_t>::max();
  uint32_t r1MaxLen = std::numeric_limits<uint32_t>::max();
  uint32_t r2MaxLen = std::numeric_limits<uint32_t>::max();

  uint32_t smallFragmentCutoff = 50;

  bool HMP = false;
  uint32_t primerLen = 30;
  bool trimTcag = false;

  uint32_t barcodeErrors = 0;
  bool midEndsRevComp = false;
  bool rename = false;

  bool noReversePrimer = false;
  bool reversePrimerToUpperCase = false;
  comparison rPrimerErrors;

  bool noForwardPrimer = false;
  bool forwardPrimerToUpperCase = false;
  comparison fPrimerErrors;


  double kmerCutOff = .20;
  uint32_t contaminationKLen = 7;

  MidDeterminator::MidDeterminePars mDetPars;
  uint32_t qualCheck = 30;
  bool checkingQCheck = false;
  double qualCheckCutOff = 0.75;
  bool illumina = false;
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

struct ExtractorPairedEndPars{

	ExtractorPairedEndPars();

	bfs::path idFilename = "";
	std::string idFileDelim = "whitespace";


  uint32_t smallFragmentCutoff = 50;


  uint32_t barcodeErrors = 0;
  bool midEndsRevComp = false;
  bool rename = false;

  bool primerToUpperCase = false;
  comparison primerErrors;


  MidDeterminator::MidDeterminePars mDetPars;

  uint32_t primerWithinStart_ = 0;

  std::string sampleName = "";

  uint32_t r1Trim_ = 1;
  uint32_t r2Trim_ = 1;
  QualFilteringPars qPars_;
  bfs::path lenCutOffFilename_ = "";

  bfs::path comparisonSeqFnp_ = "";
  uint32_t compKmerLen = 5;
  double compKmerSimCutOff = 0.50;

  bfs::path overlapStatusFnp_ = "";
  bool noOverlapProcessForNoOverlapStatusTargets_ = false;
  uint32_t numberOfNs = 0;

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




