#pragma once
/*
 * extractTargetsFromGenomes.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: nick
 */

#include <njhseq/common.h>
#include <njhseq/GenomeUtils/GenomeMapping/MultiGenomeMapper.hpp>
#include <njhseq/programUtils/seqSetUp.hpp>

#include "SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp"


namespace njhseq {

struct extractBetweenSeqsPars{
	MultiGenomeMapper::inputParameters pars;
	std::string gffExtraAttributesStr = "description";
	bfs::path primersFile = "";
	std::string forwardPrimer;
	std::string reversePrimer;
	std::string targetName;
	uint32_t errors = 0;
	uint32_t sizeLimit = 1000;
	uint32_t maxLenCutOffSizeExpand = 40;
  uint32_t minLenCutOffSizeExpand = 40;
	uint32_t pairedEndLength = std::numeric_limits<uint32_t>::max();
	uint32_t barcodeSize = 0;

  bool useBlast = false;
  uint32_t blastExpandSize = 10;

	bool shortNames = false;
	std::string selectedGenomesStr;
	bool writeOutAllSeqsFile = false;
	bool removeRefAlignments = false;
	njh::files::MkdirPar outputDirPars{"extractedRegions_TODAY"};
	bool verbose_ = false;
	bool debug_ = false;
	bool longRangeAmplicon = false;
	void setUpCoreOptions(seqSetUp & setUp, bool needReadLength = false);

};

void extractBetweenSeqs(const PrimersAndMids & ids,
		const extractBetweenSeqsPars & extractPars);

}  // namespace njhseq


