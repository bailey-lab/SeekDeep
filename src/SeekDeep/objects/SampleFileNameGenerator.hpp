#pragma once

/*
 * SampleFileNameGenerator.hpp
 *
 *  Created on: Jan 16, 2020
 *      Author: nicholashathaway
 */



#include <njhseq.h>

#include "SeekDeep/objects/PrimersAndMids.hpp"


namespace njhseq {



class SampleFileNameGenerator{
public:

	SampleFileNameGenerator(const bfs::path & idFnp, const bfs::path sampleFnp);

	bfs::path idFnp_;
	bfs::path sampleFnp_;
	std::shared_ptr<PrimersAndMids> ids_;

	bool hasRBarcode_{false};
	std::map<std::string, std::map<std::string, VecStr>> barcodesByNamePerLibraryPerSample_;
	std::map<std::string, std::pair<std::string, std::string>> namesToBarcodes_;


	void writeSampleNameFile(const OutOptions & sampleNamesOutOpts);
	void writeBarcodePrimerFile(const OutOptions & idFileOutOpts);

};




}  // namespace njhseq



