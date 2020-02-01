/*
 * SampleFileNameGenerator.cpp
 *
 *  Created on: Jan 16, 2020
 *      Author: nicholashathaway
 */


#include "SampleFileNameGenerator.hpp"

namespace njhseq {

SampleFileNameGenerator::SampleFileNameGenerator(const bfs::path & idFnp,
		const bfs::path sampleFnp) :
		idFnp_(idFnp), sampleFnp_(sampleFnp) {
	ids_ = std::make_shared<PrimersAndMids>(idFnp_);
	ids_->initPrimerDeterminator();



	table sampleTab(sampleFnp, "\t", true);
	njh::for_each(sampleTab.columnNames_, [](std::string & col){
		njh::strToLower(col);
	});
	sampleTab.setColNamePositions();
	sampleTab.checkForColumnsThrow(VecStr{"library", "sample", "fbarcode"}, __PRETTY_FUNCTION__);
	hasRBarcode_ = sampleTab.containsColumn("rbarcode");


	std::map<std::string, std::map<std::string, uint32_t>> barcodes;

	//can't have the sample barcode within library file
	std::unordered_map<std::string, std::set<std::string>> barcodesForLibrary;

	VecStr warnings;

	for(const auto & row : sampleTab){
		std::string library = row[sampleTab.getColPos("library")];
		std::string sample = row[sampleTab.getColPos("sample")];
		std::string fbarcode = row[sampleTab.getColPos("fbarcode")];
		std::string rbarcode = "";
		if(hasRBarcode_){
			rbarcode = row[sampleTab.getColPos("rbarcode")];
		}
		++barcodes[fbarcode][rbarcode];
		std::string barId = njh::pasteAsStr(fbarcode, "-", rbarcode);
		if(njh::in(barId, barcodesForLibrary[library])){
			warnings.emplace_back(njh::pasteAsStr("Already have barcodes ", fbarcode, " and ", rbarcode, " for library ", library));
		}
	}
	if(!warnings.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", errors "<< "\n";
		ss << njh::conToStr(warnings, "\n") << "\n";
		throw std::runtime_error{ss.str()};
	}

	std::map<std::string, std::map<std::string, std::string>> barcodesToNames;

	uint32_t midTotal = 0;
	for(const auto & fBar : barcodes){
		midTotal += fBar.second.size();
	}

	uint32_t midCount = 1;
	for(const auto & fBar : barcodes){
		for(const auto & rBar : fBar.second){
			std::string name = njh::pasteAsStr("MID", njh::leftPadNumStr(midCount, midTotal));
			barcodesToNames[fBar.first][rBar.first] = name;
			namesToBarcodes_[name] = std::make_pair(fBar.first, rBar.first);
			++midCount;
		}
	}

	for(const auto & row : sampleTab){
		std::string library = row[sampleTab.getColPos("library")];
		std::string sample = row[sampleTab.getColPos("sample")];
		std::string fbarcode = row[sampleTab.getColPos("fbarcode")];
		std::string rbarcode = "";
		if(hasRBarcode_){
			rbarcode = row[sampleTab.getColPos("rbarcode")];
		}
		barcodesByNamePerLibraryPerSample_[library][sample].emplace_back(barcodesToNames[fbarcode][rbarcode]);
	}
}

void SampleFileNameGenerator::writeSampleNameFile(const OutOptions & sampleNamesOutOpts){
	uint32_t maxRunCount = 1;
	for(const auto & library : barcodesByNamePerLibraryPerSample_){
		for(const auto & sample : library.second){
			if(sample.second.size() > maxRunCount){
				maxRunCount = sample.second.size();
			}
		}
	}

	VecStr headerForSampNamesFile{"#input", "sample"};
	for(uint32_t run = 1; run <= maxRunCount; ++run){
		headerForSampNamesFile.emplace_back(njh::pasteAsStr("run", run));
	}

	//sample name file
	OutputStream sampleNamesOut(sampleNamesOutOpts);
	sampleNamesOut << njh::conToStr(headerForSampNamesFile, "\t") << std::endl;
	for(const auto & library : barcodesByNamePerLibraryPerSample_){
		for(const auto & sample : library.second){
			sampleNamesOut << library.first << '\t' << sample.first;
			for(const auto & run : sample.second){
				sampleNamesOut << "\t" << run;
			}
			sampleNamesOut << std::endl;
		}
	}
}

void SampleFileNameGenerator::writeBarcodePrimerFile(const OutOptions & idFileOutOpts){
	//barcode/primer file
	OutputStream idFileOut(idFileOutOpts);
	auto tarkeys = getVectorOfMapKeys(ids_->targets_);
	njh::sort(tarkeys);
	idFileOut << "target\tfroward\treverse" << std::endl;

	for(const auto & tar : tarkeys){
		idFileOut << tar
				<< "\t" << ids_->targets_.at(tar).info_.forwardPrimerRaw_
				<< "\t" << ids_->targets_.at(tar).info_.reversePrimerRaw_ << std::endl;
	}
	idFileOut << "id\tbarcode" << std::endl;
	if(hasRBarcode_){
		idFileOut << "\t" << "revBarcode";
	}
	idFileOut << std::endl;

	for(const auto & bar : namesToBarcodes_){
		idFileOut << bar.first
				<< "\t" << bar.second.first;
		if(bar.second.second != ""){
			idFileOut << "\t" << bar.second.second;
		}
		idFileOut << std::endl;
	}
}

}  // namespace njhseq
