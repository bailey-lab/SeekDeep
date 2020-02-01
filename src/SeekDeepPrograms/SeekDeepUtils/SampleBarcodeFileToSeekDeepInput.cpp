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

	//
	SampleFileNameGenerator fileGen(idFnp, sampleFnp);

	//sample name file
	OutOptions sampleNamesOutOpts(outStub.outFilename_.string(), "_sampleNames.tab.txt");
	sampleNamesOutOpts.transferOverwriteOpts(outStub);
	fileGen.writeSampleNameFile(sampleNamesOutOpts);

	//barcode/primer file
	OutOptions idFileOutOpts(outStub.outFilename_.string(), "_ids.tab.txt");
	idFileOutOpts.transferOverwriteOpts(outStub);
	fileGen.writeBarcodePrimerFile(idFileOutOpts);


	return 0;
}

}  // namespace njhseq

