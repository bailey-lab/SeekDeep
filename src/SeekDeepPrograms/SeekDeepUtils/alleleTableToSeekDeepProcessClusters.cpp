//
// Created by Nicholas Hathaway on 10/26/24.
//
#include "SeekDeepUtilsRunner.hpp"

#include <njhseq/objects/dataContainers/tables/TableReader.hpp>

namespace njhseq {

int SeekDeepUtilsRunner::alleleTableToSeekDeepProcessClusters(
				const njh::progutils::CmdArgs &inputCommands) {


	bfs::path alleleTableFnp;
	std::string sampleColName = "sampleID";
	std::string withinSampleReadCntColName = "reads";
	std::string hapSeqColName = "asv";
	std::string targetNameColName = "locus";

	std::set<std::string> selectTargets;
	std::set<std::string> selectSamples;

	std::set<std::string> excludeTargets;
	std::set<std::string> excludeSamples;

	bfs::path targetKeyFnp;
	std::string targetKeyOldTargetNameColName;
	std::string targetKeyNewTargetNameColName;

	bfs::path outputFnp = "output.fasta.gz";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(outputFnp, "--outputFnp", "name of the output directories");
	setUp.setOption(targetKeyFnp, "--targetKeyFnp", "target Key  Fnp, tab delimited assumed first column old name and second column new name, unless columns provided with targetKeyOldTargetNameColName and targetKeyNewTargetNameColName arguments");
	setUp.setOption(targetKeyOldTargetNameColName, "--targetKeyOldTargetNameColName", "if renaming targets, target Key Old Target Name Col Name");
	setUp.setOption(targetKeyNewTargetNameColName, "--targetKeyNewTargetNameColName", "if renaming targets, target Key New Target Name Col Name");

	setUp.setOption(selectTargets, "--selectTargets", "Only analyze these select targets");
	setUp.setOption(selectSamples, "--selectSamples", "Only analyze these select samples");

	setUp.setOption(excludeTargets, "--excludeTargets", "Exclude these select targets from analysis");
	setUp.setOption(excludeSamples, "--excludeSamples", "Exclude these select samples from analysis");

	setUp.setOption(alleleTableFnp, "--alleleTableFnp",
									"results tab delimited file, each row is a haplotype, should have at least 4 columns,  sample (--sampleColName),within sample read count (--withinSampleReadCntColName), haplotype sequnece (--popHapSeqColName), target name column (--targetNameColName)",
									true);
	setUp.setOption(sampleColName, "--sampleColName", "sample Column Name", false, "Results Column Names");
	setUp.setOption(withinSampleReadCntColName, "--withinSampleReadCntColName", "within Sample Read Cnt Col Column Name", false, "Results Column Names");
	setUp.setOption(hapSeqColName, "--hapSeqColName",
									"Haplotype Sequence Column Name", false, "Results Column Names");
	setUp.setOption(targetNameColName, "--targetNameColName",
									"target Name Column Name, the column name in the table which indicates the different targets", false, "Results Column Names");
	setUp.processDirectoryOutputName(njh::pasteAsStr(bfs::basename(alleleTableFnp), "_", setUp.commands_.subProgram_, "_TODAY"), true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);


	std::set<std::string> inputTargets;
	std::set<std::string> inputSamples;

	TableReader alleleTableReader(TableIOOpts::genTabFileIn(alleleTableFnp));
	alleleTableReader.header_.checkForColumnsThrow(VecStr{sampleColName, withinSampleReadCntColName, hapSeqColName, targetNameColName}, __PRETTY_FUNCTION__);
	auto sampleColPos = alleleTableReader.header_.getColPos(sampleColName);
	auto readColPos = alleleTableReader.header_.getColPos(withinSampleReadCntColName);
	auto seqColPos = alleleTableReader.header_.getColPos(hapSeqColName);
	auto targetColPos = alleleTableReader.header_.getColPos(targetNameColName);
	auto samplePasses = [&selectSamples,&excludeSamples](const std::string & sample) {
		return ( selectSamples.empty() && excludeSamples.empty() )  ||
			  (!selectSamples.empty() && njh::in(sample, selectSamples) )||
				(!excludeSamples.empty() && njh::notIn(sample, excludeSamples) )
					;
	};
	auto targetPasses = [&selectTargets,&excludeTargets](const std::string & target) {
		return ( selectTargets.empty() && excludeTargets.empty() )  ||
				(!selectTargets.empty() && njh::in(target, selectTargets) )||
				(!excludeTargets.empty() && njh::notIn(target, excludeTargets))
				;
	};
	{
		//checks
		VecStr warnings;
		for(const auto & t : excludeTargets) {
			if(njh::in(t, selectTargets)) {
				warnings.emplace_back(njh::pasteAsStr("can't have ", t, " in both --selectTargets and in --excludeTargets"));
			}
		}
		for(const auto & s : excludeSamples) {
			if(njh::in(s, selectSamples)) {
				warnings.emplace_back(njh::pasteAsStr("can't have ", s, " in both --selectSamples and in --excludeSamples"));
			}
		}
		if(!warnings.empty()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error "	 << "\n";
			ss << njh::conToStr(warnings, "\n") << "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	//k1 = target, k2 = sample, data = vector of seqs
	std::unordered_map<std::string, std::unordered_map<std::string, std::vector<seqInfo>>> seqsPerTargetPerSample;

	std::unordered_map<std::string, std::string> targetKey;
	if(!targetKeyFnp.empty()) {
		if(!targetKeyNewTargetNameColName.empty() && !targetKeyOldTargetNameColName.empty()) {
			table targetKeyTab(targetKeyFnp, "\t", true);
			targetKeyTab.checkForColumnsThrow(VecStr{targetKeyNewTargetNameColName, targetKeyOldTargetNameColName}, __PRETTY_FUNCTION__);
			auto targetKeyNewTargetNameColPos = targetKeyTab.getColPos(targetKeyNewTargetNameColName);
			auto targetKeyOldTargetNameColPos = targetKeyTab.getColPos(targetKeyOldTargetNameColName);
			for(const auto & row : targetKeyTab) {
				if(njh::in(row[targetKeyOldTargetNameColPos], targetKey)) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "already have target: " << row[targetKeyOldTargetNameColPos] << "\n";
					throw std::runtime_error{ss.str()};
				}
				targetKey[row[targetKeyOldTargetNameColPos]] = row[targetKeyNewTargetNameColPos];
			}
		} else {
			table targetKeyTab(targetKeyFnp, "\t", false);
			if(targetKeyTab.columnNames_.size() != 2) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << targetKeyFnp << " should be two columns, not: " << targetKeyTab.columnNames_.size() << "\n";
				throw std::runtime_error{ss.str()};
			}
			for(const auto & row : targetKeyTab) {
				if(njh::in(row[0], targetKey)) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "already have target: " << row[0] << "\n";
					throw std::runtime_error{ss.str()};
				}
				targetKey[row[0]] = row[1];
			}
		}
	}
	{
		//first pass
		VecStr row;
		while(alleleTableReader.getNextRow(row)) {
			if(targetPasses(row[targetColPos]) && samplePasses(row[sampleColPos])) {
				auto sample = row[sampleColPos];
				auto target = row[targetColPos];
				auto readCnt = njh::StrToNumConverter::stoToNum<double>(row[readColPos]);

				inputTargets.emplace(target);
				inputSamples.emplace(sample);
				MetaDataInName meta;
				meta.addMeta("sample", sample);
				meta.addMeta("target", target);
				meta.addMeta("readCount", readCnt);
				seqsPerTargetPerSample[target][sample].emplace_back(njh::pasteAsStr(meta.createMetaName(), "_t", readCnt), row[seqColPos]);
			}
		}
	}
	{
		//check if targetKey
		if(targetKeyFnp.empty()) {
			for(const auto & tar : inputTargets) {
				targetKey[tar] = tar;
			}
		} else {
			VecStr missingTargets;

			for(const auto & tar : inputTargets) {
				if(njh::notIn(tar, targetKey)) {
					missingTargets.emplace_back(tar);
				}
			}
			if(!missingTargets.empty()) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "missing the following targets from the naming key: " << njh::conToStr(missingTargets, ",") << "\n";
				throw std::runtime_error{ss.str()};
			}
		}
	}
	
	//sort by total count
	for(auto & tar : seqsPerTargetPerSample) {
		for(auto & samp : tar.second) {
			readVecSorter::sortByTotalCount(samp.second, true);
		}
	}
	auto popDir = njh::files::make_path(setUp.pars_.directoryName_, "popClustering");
	njh::files::makeDir(popDir);
	{
		//create directory structure
		for(const auto & t : seqsPerTargetPerSample) {
			auto tarDir = njh::files::makeDir(popDir, bfs::path(targetKey[t.first]));
			for(const auto & s : t.second) {
				auto sampDir = njh::files::make_path(tarDir, s.first, s.first);
				njh::files::makeDirP(sampDir);
				auto outOpts = SeqIOOptions::genFastaOutGz(njh::files::make_path(sampDir, outputFnp));
				SeqOutput::write(s.second, outOpts);
			}
		}
	}


  return 0;
}

}  //namespace njhseq

