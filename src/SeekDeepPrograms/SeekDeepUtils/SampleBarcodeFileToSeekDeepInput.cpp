/*
 * SampleBarcodeFileToSeekDeepInput.cpp
 *
 *  Created on: Jan 14, 2020
 *      Author: nicholashathaway
 */


#include "SeekDeepUtilsRunner.hpp"

#include "SeekDeep/objects.h"
#include "SeekDeep/parameters.h"



namespace njhseq {




int SeekDeepUtilsRunner::SampleBarcodeFileToSeekDeepInput(
		const njh::progutils::CmdArgs & inputCommands) {
	bfs::path idFnp = "";
	bfs::path sampleFnp = "";
	OutOptions outStub(bfs::path("out"));
	seqSetUp setUp(inputCommands);
	setUp.setOption(idFnp, "--id", "SeekDeep primers file, 3 columns 1)target,2)forward,3)reverse", true);
	setUp.setOption(sampleFnp, "--sampleFnp", "Sample file 3 or 4 required columns 1)library,2)sample,3)fbarcode,4(optional))rbarcode."
			"\n\t\t\t1) name of input file without extension/illumina info e.g. Sample1 for Sample1_S2_R1_001.fastq.gz"
			"\n\t\t\t2) sample name to be given to this barcode in this sample"
			"\n\t\t\t3) barcode sequence associated with forward primer"
			"\n\t\t\t4) if sample is dual barcoded, barcode associated with reverse primer", true);
	setUp.processWritingOptions(outStub);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	PrimersAndMids ids(idFnp);
	ids.initPrimerDeterminator();


	table sampleTab(sampleFnp, "\t", true);
	njh::for_each(sampleTab.columnNames_, [](std::string & col){
		njh::strToLowerRet(col);
	});
	sampleTab.setColNamePositions();
	sampleTab.checkForColumnsThrow(VecStr{"library", "sample", "fbarcode"}, __PRETTY_FUNCTION__);
	bool hasRBarcode = sampleTab.containsColumn("rbarcode");


	std::map<std::string, std::map<std::string, uint32_t>> barcodes;

	//can't have the sample barcode within library file

	std::unordered_map<std::string, std::set<std::string>> barcodesForLibrary;

	VecStr warnings;

	for(const auto & row : sampleTab){
		std::string library = row[sampleTab.getColPos("library")];
		std::string sample = row[sampleTab.getColPos("sample")];
		std::string fbarcode = row[sampleTab.getColPos("fbarcode")];
		std::string rbarcode = "";
		if(hasRBarcode){
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
	std::map<std::string, std::pair<std::string, std::string>> namesToBarcodes;

	uint32_t midTotal = 0;
	for(const auto & fBar : barcodes){
		midTotal += fBar.second.size();
	}

	uint32_t midCount = 1;
	for(const auto & fBar : barcodes){
		for(const auto & rBar : fBar.second){
			std::string name = njh::pasteAsStr("MID", njh::leftPadNumStr(midCount, midTotal));
			barcodesToNames[fBar.first][rBar.first] = name;
			namesToBarcodes[name] = std::make_pair(fBar.first, rBar.first);
			++midCount;
		}
	}

	std::map<std::string, std::map<std::string, VecStr>> barcodesByNamePerLibraryPerSample;
	for(const auto & row : sampleTab){
		std::string library = row[sampleTab.getColPos("library")];
		std::string sample = row[sampleTab.getColPos("sample")];
		std::string fbarcode = row[sampleTab.getColPos("fbarcode")];
		std::string rbarcode = "";
		if(hasRBarcode){
			rbarcode = row[sampleTab.getColPos("rbarcode")];
		}
		barcodesByNamePerLibraryPerSample[library][sample].emplace_back(barcodesToNames[fbarcode][rbarcode]);
	}
	uint32_t maxRunCount = 1;
	for(const auto & library : barcodesByNamePerLibraryPerSample){
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
	OutOptions sampleNamesOutOpts(outStub.outFilename_.string(), "_sampleNames.tab.txt");
	OutputStream sampleNamesOut(sampleNamesOutOpts);
	sampleNamesOut << njh::conToStr(headerForSampNamesFile, "\t") << std::endl;
	for(const auto & library : barcodesByNamePerLibraryPerSample){
		for(const auto & sample : library.second){
			sampleNamesOut << library.first << '\t' << sample.first;
			for(const auto & run : sample.second){
				sampleNamesOut << "\t" << run;
			}
			sampleNamesOut << std::endl;
		}
	}

	//barcode/primer file
	OutOptions idFileOutOpts(outStub.outFilename_.string(), "_ids.tab.txt");
	OutputStream idFileOut(idFileOutOpts);

	auto tarkeys = getVectorOfMapKeys(ids.targets_);
	njh::sort(tarkeys);
	idFileOut << "target\tfroward\treverse" << std::endl;

	for(const auto & tar : tarkeys){
		idFileOut << tar
				<< "\t" << ids.targets_.at(tar).info_.forwardPrimerRaw_
				<< "\t" << ids.targets_.at(tar).info_.reversePrimerRaw_ << std::endl;
	}
	idFileOut << "id\tbarcode" << std::endl;
	if(hasRBarcode){
		idFileOut << "\t" << "revBarcode";
	}
	idFileOut << std::endl;

	for(const auto & bar : namesToBarcodes){
		idFileOut << bar.first
				<< "\t" << bar.second.first;
		if(bar.second.second != ""){
			idFileOut << "\t" << bar.second.second;
		}
		idFileOut << std::endl;
	}


	return 0;
}

}  // namespace njhseq

