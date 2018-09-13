#pragma once
/*
 * extractTargetsFromGenomes.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: nick
 */

#include <bibseq.h>
#include "SeekDeep/objects/PrimersAndMids.hpp"


namespace bibseq {

struct extractBetweenSeqsPars{
	MultiGenomeMapper::inputParameters pars;
	std::string gffExtraAttributesStr = "description";
	bfs::path primersFile = "";
	std::string forwardPrimer = "";
	std::string reversePrimer = "";
	std::string targetName = "";
	uint32_t errors = 0;
	uint32_t sizeLimit = 1000;
	uint32_t lenCutOffSizeExpand = 20;
	uint32_t pairedEndLength = std::numeric_limits<uint32_t>::max();
	uint32_t barcodeSize = 0;
	bool shortNames = false;
	std::string selectedGenomesStr = "";
	bool removeRefAlignments = false;
	bib::files::MkdirPar outputDirPars{"extractedRegions_TODAY"};

	void setUpCoreOptions(seqSetUp & setUp, bool needReadLength = false);

};

void extractBetweenSeqs(const PrimersAndMids & ids,
		const extractBetweenSeqsPars & extractPars);

}  // namespace bibseq


