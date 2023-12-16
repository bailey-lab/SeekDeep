//
// Created by Nicholas Hathaway on 12/16/23.
//
#include "SeekDeepUtilsRunner.hpp"

#include <njhseq/objects/dataContainers/tables/TableReader.hpp>

namespace njhseq {


int SeekDeepUtilsRunner::variantCallOnSeqAndProtein(
				const njh::progutils::CmdArgs &inputCommands) {


	bfs::path resultsFnp;
	std::string sampleColName = "s_Sample";
	std::string withinSampleReadCntColName = "c_ReadCnt";
	std::string popHapIdColName = "h_popUID";
	std::string popHapSeqColName = "h_Consensus";

	std::string targetNameColName = "p_name";

	bfs::path popSeqsDirFnp = "";

	bfs::path metaFnp = "";

	std::string popSeqsRegexPatRemoval = R"(_([tf])?\d+(\.\d+)?$)";
	uint32_t numThreads = 1;

	CollapseAndCallVariantsPars collapseVarCallPars;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();

	setUp.setOption(collapseVarCallPars.variantCallerRunPars.occurrenceCutOff, "--occurrenceCutOff", "Occurrence Cut Off, don't report variants below this count");
	setUp.setOption(collapseVarCallPars.variantCallerRunPars.lowVariantCutOff, "--lowVariantCutOff", "Low Variant Cut Off, don't report variants below this fraction");
	collapseVarCallPars.calcPopMeasuresPars.lowVarFreq = collapseVarCallPars.variantCallerRunPars.lowVariantCutOff;
	collapseVarCallPars.transPars.setOptions(setUp, true);
	setUp.setOption(collapseVarCallPars.calcPopMeasuresPars.getPairwiseComps, "--getPairwiseComps", "get Pairwise comparison metrics");
	setUp.setOption(collapseVarCallPars.noDiagAlnPairwiseComps, "--noDiagAlnPairwiseComps", "Use diagonal Alignment for Pairwise Comparisons");
	collapseVarCallPars.calcPopMeasuresPars.diagAlnPairwiseComps = !collapseVarCallPars.noDiagAlnPairwiseComps;
	//setOption(collapseVarCallPars.ignoreSubFields, "--ignoreSubFields", "Meta Sub Field values to ignore when calculating variants, e.g. --ignoreSubFields \"isFieldSample:TRUE,PreferredSample:FALSE\"");
	setUp.setOption(collapseVarCallPars.calcPopMeasuresPars.numThreads, "--mappingNumThreads", "Number of threads to use for the alignment portion");

	setUp.setOption(numThreads, "--numThreads", "Number of threads to use");

	
	setUp.setOption(resultsFnp, "--resultsFnp",
									"results tab delimited file, each row is a haplotype, should have at least 5 columns, 1) sample (--sampleColName), 2)within sample freq (--withinSampleFreqColName), 3)within sample read count (--withinSampleReadCntColName), 4)haplotype pop ID (--popHapIdColName), 5)target name column (--targetNameColName), optionally 4th col with hap sequence (--popHapSeqColName) or read in from --popSeqsDirFnp",
									true);

	setUp.setOption(sampleColName, "--sampleColName", "sample Column Name", false, "Results Column Names");
	setUp.setOption(withinSampleReadCntColName, "--withinSampleReadCntColName", "within Sample Read Cnt Col Column Name", false, "Results Column Names");
	setUp.setOption(popHapIdColName, "--popHapIdColName", "popHapIdColName", false, "Results Column Names");
	setUp.setOption(popHapSeqColName, "--popHapSeqColName",
									"population Haplotype Sequence Column Name, the seq to compare to expected", false, "Results Column Names");
	setUp.setOption(targetNameColName, "--targetNameColName",
									"target Name Column Name, the column name in the table which indicates the different targets", false, "Results Column Names");

	setUp.setOption(popSeqsDirFnp, "--popSeqsDirFnp",
									"Population Sequences, in this directory should be a fasta file with the name of each target");

	setUp.setOption(metaFnp, "--metaFnp", "meta data for the control samples");
	setUp.setOption(popSeqsRegexPatRemoval, "--popSeqsRegexPatRemoval",
									"optional regex pattern to process the input pop sequences from --popSeqsFnp to make it match up with the input results folder");

	setUp.processDirectoryOutputName("variantCalling_TODAY", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);


	njh::files::checkExistenceThrow({resultsFnp},
																	__PRETTY_FUNCTION__);


	std::unique_ptr<MultipleGroupMetaData> metaGroupData;
	if (exists(metaFnp)) {
		metaGroupData = std::make_unique<MultipleGroupMetaData>(metaFnp);
	}
	std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, double>>> readCountsPerHapPerSample;
	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> cNameToPopUID;
	std::unordered_map<std::string, std::unordered_map<std::string, std::set<std::string>>> hPopUIDPopSamps;
	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> hPopUID_to_hConsensus;

	std::set<std::string> targetNamesSet;

	VecStr requiredColumns{sampleColName, withinSampleReadCntColName,
												  withinSampleReadCntColName,
												 popHapIdColName,
												 targetNameColName};
	if (!exists(popSeqsDirFnp)) {
		requiredColumns.emplace_back(popHapSeqColName);
	}
	uint64_t maxLen = 0;
	//population seqs;
	std::unordered_map<std::string, std::vector<seqInfo>> popSeqs;
	auto sampInfoFnp = resultsFnp;
	TableReader sampInfoReader(TableIOOpts::genTabFileIn(sampInfoFnp, true));
	sampInfoReader.header_.checkForColumnsThrow(requiredColumns, __PRETTY_FUNCTION__);
	VecStr row;
	std::unordered_map<std::string, std::unordered_set<std::string>> allSamplesInOutput;
	//key1 == target, key2 == sample
	std::unordered_map<std::string, std::unordered_map<std::string, std::vector<seqInfo>>> allResultSeqs;
	while (sampInfoReader.getNextRow(row)) {
		auto target = row[sampInfoReader.header_.getColPos(targetNameColName)];
		targetNamesSet.emplace(target);
		auto sample = row[sampInfoReader.header_.getColPos(sampleColName)];

		auto c_ReadCnt = njh::StrToNumConverter::stoToNum<double>(
						row[sampInfoReader.header_.getColPos(withinSampleReadCntColName)]);
		auto readCnt = njh::StrToNumConverter::stoToNum<double>(
						row[sampInfoReader.header_.getColPos(withinSampleReadCntColName)]);

		auto h_popUID = row[sampInfoReader.header_.getColPos(popHapIdColName)];
		auto hapName = njh::pasteAsStr(sample, "__", h_popUID);
		std::string hapSeq;
		if (popSeqsDirFnp.empty()) {
			hPopUID_to_hConsensus[target][h_popUID] = row[sampInfoReader.header_.getColPos(popHapSeqColName)];
			hapSeq = row[sampInfoReader.header_.getColPos(popHapSeqColName)];
		} else {
			hapSeq = njh::mapAt(njh::mapAt(hPopUID_to_hConsensus, target), h_popUID);
		}
		seqInfo clus(hapName, hapSeq);
		clus.cnt_ = c_ReadCnt;
		bool add = true;
		for (auto &seq: allResultSeqs[target][sample]) {
			if (seq.seq_ == clus.seq_) {
				add = false;
				seq.cnt_ += c_ReadCnt;
				readCountsPerHapPerSample[target][sample][hapName] += readCnt;
				break;
			}
		}
		if (add) {
			readCountsPerHapPerSample[target][sample][hapName] = readCnt;
			allResultSeqs[target][sample].emplace_back(clus);
		}
		readVec::getMaxLength(clus, maxLen);
		cNameToPopUID[target][hapName] = h_popUID;
		hPopUIDPopSamps[target][h_popUID].emplace(sample);
		allSamplesInOutput[target].emplace(sample);
	}


	if (!exists(popSeqsDirFnp)) {
		for (const auto &tarPopHaps: hPopUID_to_hConsensus) {
			for (const auto &popHap: tarPopHaps.second) {
				popSeqs[tarPopHaps.first].emplace_back(popHap.first, popHap.second);
			}
		}
	}

	//write out seqs
	for(const auto & tar : allResultSeqs) {
		auto tarDir = njh::files::make_path(setUp.pars_.directoryName_, tar.first);
		njh::files::makeDir(njh::files::MkdirPar{tarDir});
		auto outSeqOpts = SeqIOOptions::genFastaOutGz(njh::files::make_path(setUp.pars_.directoryName_, tar.first, "inputSeqs.fasta.gz"));
		SeqOutput writer(outSeqOpts);
		writer.openOut();
		for(const auto & sampSeqs : tar.second) {
			for( auto  seq : sampSeqs.second) {
				MetaDataInName meta;
				if(nullptr != metaGroupData) {
					meta = metaGroupData->getMetaForSample(sampSeqs.first);
				}
				meta.addMeta("sample", sampSeqs.first, true);
				meta.addMeta("readCount", std::max(1.0, std::round(seq.cnt_)), true);
				meta.resetMetaInName(seq.name_);
				writer.write(seq);
			}
		}
	}


	njh::concurrent::LockableQueue<std::string> targetNamesQueue(targetNamesSet);

	std::function<void()> callVariants = [&collapseVarCallPars,&targetNamesQueue, &setUp]() {
		std::string target;
		while(targetNamesQueue.getVal(target)) {
			auto inputSeqsOpts = SeqIOOptions::genFastaInGz(njh::files::make_path(setUp.pars_.directoryName_, target, "inputSeqs.fasta.gz"));
			auto inputSeqs = SeqInput::getSeqVec<seqInfo>(inputSeqsOpts);
			const auto varCallDirPath = njh::files::make_path(setUp.pars_.directoryName_, target,  "variantCalling");
			auto collapseVarCallParsForTar = collapseVarCallPars;
			collapseVarCallParsForTar.identifier = target;
			collapseVarCallParsForTar.outputDirectory = varCallDirPath;
			collapseAndCallVariants(collapseVarCallParsForTar, inputSeqs);
		}

	};

	njh::concurrent::runVoidFunctionThreaded(callVariants, numThreads	);


	return 0;

}

}  // namespace njhseq