#pragma once

/*
 * SampleFileNameGenerator.hpp
 *
 *  Created on: Jan 16, 2020
 *      Author: nicholashathaway
 */



#include <njhseq/common.h>

#include "SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp"


namespace njhseq {



class SampleFileNameGenerator{
public:

	SampleFileNameGenerator(const bfs::path & idFnp, const bfs::path sampleFnp, bool noBarcodes = false);

	bfs::path idFnp_;
	bfs::path sampleFnp_;
	std::shared_ptr<PrimersAndMids> ids_;

	bool noBarcodes_;
	bool hasRBarcode_{false};
	std::map<std::string, std::map<std::string, VecStr>> barcodesByNamePerLibraryPerSample_;
	std::map<std::string, std::pair<std::string, std::string>> namesToBarcodes_;
	std::map<std::string, VecStr> libraryFilesForSample_;


	void writeSampleNameFile(const OutOptions & sampleNamesOutOpts);
	void writeBarcodePrimerFile(const OutOptions & idFileOutOpts);

};




}  // namespace njhseq



