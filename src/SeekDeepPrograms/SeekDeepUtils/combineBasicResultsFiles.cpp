//
// Created by Nicholas Hathaway on 2/27/24.
//
#include "SeekDeepUtilsRunner.hpp"

#include <njhseq/objects/dataContainers/tables/TableReader.hpp>

namespace njhseq {

int SeekDeepUtilsRunner::combineBasicResultsFiles(
				const njh::progutils::CmdArgs &inputCommands) {


	std::vector<bfs::path> resultFnps;
	std::string sampleColName = "s_Sample";
	std::string withinSampleReadCntColName = "c_ReadCnt";
	std::string popHapIdColName = "h_popUID";
	std::string popHapSeqColName = "h_Consensus";
	std::string targetNameColName = "p_name";

	std::set<std::string> selectTargets;
	std::set<std::string> selectSamples;

	OutOptions out_options("", ".tsv.gz");

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();

	setUp.setOption(selectTargets, "--selectTargets", "Only analzye these select targets");
	setUp.setOption(selectSamples, "--selectSamples", "Only analzye these select samples");


	setUp.setOption(resultFnps, "--resultFnps",
									"results tab delimited files, each row is a haplotype, should have at least 5 columns, 1) sample (--sampleColName), 2)within sample freq (--withinSampleFreqColName), 3)within sample read count (--withinSampleReadCntColName), 4)haplotype pop ID (--popHapIdColName), 5)target name column (--targetNameColName), optionally 4th col with hap sequence (--popHapSeqColName) or read in from --popSeqsDirFnp",
									true);

	setUp.setOption(sampleColName, "--sampleColName", "sample Column Name", false, "Results Column Names");
	setUp.setOption(withinSampleReadCntColName, "--withinSampleReadCntColName", "within Sample Read Cnt Col Column Name", false, "Results Column Names");
	setUp.setOption(popHapIdColName, "--popHapIdColName", "popHapIdColName", false, "Results Column Names");
	setUp.setOption(popHapSeqColName, "--popHapSeqColName",
									"population Haplotype Sequence Column Name, the seq to call variants on", false, "Results Column Names");
	setUp.setOption(targetNameColName, "--targetNameColName",
									"target Name Column Name, the column name in the table which indicates the different targets", false, "Results Column Names");

	setUp.processWritingOptions(out_options);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);


	njh::files::checkExistenceThrow(resultFnps,
																	__PRETTY_FUNCTION__);
	OutputStream out(out_options);

	//key1 == target, key2 == hap seq, value = count
	//used for renaming
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> hapSeqsForTargets;

	//read through to get seqs to count
	for(const auto & fnp : resultFnps) {
		TableReader tabIn(TableIOOpts::genTabFileIn(fnp));
		tabIn.header_.checkForColumnsThrow(VecStr{sampleColName, popHapSeqColName, targetNameColName, withinSampleReadCntColName}, __PRETTY_FUNCTION__);

		VecStr row;
		while(tabIn.getNextRow(row)) {
			auto target = row[tabIn.header_.getColPos(targetNameColName)];
			//filter to just the select targets if filtering for that
			if(!selectTargets.empty() && !njh::in(target, selectTargets)) {
				continue;
			}

			auto sample = row[tabIn.header_.getColPos(sampleColName)];
			//filter to just the select samples if filtering for that
			if(!selectSamples.empty() && njh::notIn(sample, selectSamples)) {
				continue;
			}
			// auto readCnt = njh::StrToNumConverter::stoToNum<double>(row[tabIn.header_.getColPos(withinSampleReadCntColName)]);
			auto seq = row[tabIn.header_.getColPos(popHapSeqColName)];
			++hapSeqsForTargets[target][seq];
		}
	}

	//create names
	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> hapSeqsForTargetsRenamed;
	for(auto & tar : hapSeqsForTargets) {
		auto hapSeqs = getVectorOfMapKeys(tar.second);
		njh::sort(hapSeqs, [&tar](const std::string & seq1, const std::string & seq2) {
			if(tar.second[seq1] == tar.second[seq2]) {
				return seq1 > seq2;
			}
			return tar.second[seq1] > tar.second[seq2];
		});
		for(const auto & hapSeq : iter::enumerate(hapSeqs)) {
			hapSeqsForTargetsRenamed[tar.first][hapSeq.element] = njh::pasteAsStr(tar.first, ".", njh::leftPadNumStr(hapSeq.index, hapSeqs.size()));
		}
	}

	//write out
	out << sampleColName
	<< "\t" << targetNameColName
	<< "\t" << popHapIdColName
	<< "\t" << withinSampleReadCntColName
	<< "\t" << popHapSeqColName << std::endl;

	for(const auto & fnp : resultFnps) {
		TableReader tabIn(TableIOOpts::genTabFileIn(fnp));
		tabIn.header_.checkForColumnsThrow(VecStr{sampleColName, popHapSeqColName, targetNameColName, withinSampleReadCntColName}, __PRETTY_FUNCTION__);

		VecStr row;
		while(tabIn.getNextRow(row)) {
			auto target = row[tabIn.header_.getColPos(targetNameColName)];
			//filter to just the select targets if filtering for that
			if(!selectTargets.empty() && !njh::in(target, selectTargets)) {
				continue;
			}

			auto sample = row[tabIn.header_.getColPos(sampleColName)];
			//filter to just the select samples if filtering for that
			if(!selectSamples.empty() && njh::notIn(sample, selectSamples)) {
				continue;
			}
			auto readCnt = row[tabIn.header_.getColPos(withinSampleReadCntColName)];
			auto seq = row[tabIn.header_.getColPos(popHapSeqColName)];

			out << sample
			<< "\t" << target
			<< "\t" << hapSeqsForTargetsRenamed[target][seq]
			<< "\t" << readCnt
			<< "\t" << seq
			<< std::endl;
		}
	}


	return 0;
}

} //namespace njhseq

