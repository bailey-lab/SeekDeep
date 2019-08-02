/*
 * SeekDeepUtilsRunner_benchmarkControlMixtures.cpp
 *
 *  Created on: May 6, 2019
 *      Author: nicholashathaway
 */


#include "SeekDeepUtilsRunner.hpp"
#include "SeekDeep/objects/ControlBenchmarking.h"
namespace njhseq {


int SeekDeepUtilsRunner::benchmarkControlMixtures(
		const njh::progutils::CmdArgs & inputCommands) {
	ControlBencher::ControlBencherPars conBenchPars;
	bfs::path processClustersDir = "";
	bfs::path expectedSeqsFnp = "";
	std::string name = "";
	bfs::path metaFnp = "";
	bool skipMissingSamples = false;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(processClustersDir, "--processClustersDir", "Process Clusters Dir", true);
	setUp.setOption(expectedSeqsFnp, "--expectedSeqsFnp", "Expected Seqs fasta file", true);
	setUp.setOption(name, "--name", "Name to give the current analysis", true);
	setUp.setOption(metaFnp, "--metaFnp", "meta data for the control samples");

	setUp.setOption(conBenchPars.samplesToMixFnp_, "--sampleToMixture", "Sample To Mixture, 2 columns 1)sample, 2)MixName", true);
	setUp.setOption(conBenchPars.mixSetUpFnp_, "--mixtureSetUp", "Mixture Set Up, 3 columns 1)MixName, 2)strain, 3)relative_abundance", true);
	setUp.setOption(skipMissingSamples, "--skipMissingSamples", "Skip Samples if they are missing");
	setUp.processDirectoryOutputName(name + "_benchmarkControlMixtures_TODAY", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	auto coreJsonFnp = njh::files::make_path(processClustersDir,"coreInfo.json");

	njh::files::checkExistenceThrow( { processClustersDir, coreJsonFnp, expectedSeqsFnp,
			conBenchPars.samplesToMixFnp_, conBenchPars.mixSetUpFnp_ },
			__PRETTY_FUNCTION__);

	Json::Value coreJson = njh::json::parseFile(coreJsonFnp.string());

	collapse::SampleCollapseCollection analysisMaster(coreJson);

	ControlBencher bencher(conBenchPars);
	//check for samples in analysis
	auto controlSamples = bencher.getSamples();
	VecStr missingSamples;
	for (const auto & sname : controlSamples) {
		if (!analysisMaster.hasSample(sname)) {
			missingSamples.emplace_back(sname);
		}
	}

	if (!skipMissingSamples && !missingSamples.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error "
				<< "missing the following samples in " << processClustersDir << ": "
				<< njh::conToStrEndSpecial(missingSamples, ", ", " and ") << "\n";
		throw std::runtime_error { ss.str() };
	}else if(skipMissingSamples && !missingSamples.empty()){
		OutOptions outOptsMissing(njh::files::make_path(setUp.pars_.directoryName_, "missingSamples"));
		OutputStream outMissing(outOptsMissing);
		outMissing << njh::conToStr(missingSamples, "\n") << std::endl;
	}


	//read in expected seqs
	SeqInput expSeqsSeqIn(SeqIOOptions::genFastaIn(expectedSeqsFnp));
	auto initialExpSeqs = expSeqsSeqIn.readAllReadsPtrs<seqInfo>();
	//check the needed expected names
	std::set<std::string> expNames;
	std::unordered_map<std::string, uint32_t> initialExpSeqsPositions;
	std::unordered_map<std::string, std::string> expSeqsKey;
	std::unordered_map<std::string, VecStr> matchingExpectedSeqs;

	std::vector<std::shared_ptr<seqInfo>> expSeqs;
	for(const auto & expSeq : initialExpSeqs){
		if(njh::in(expSeq->name_, expNames)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << " already have expected sequence named " << expSeq->name_<< "\n";
			throw std::runtime_error{ss.str()};
		}
		expNames.emplace(expSeq->name_);
		bool found = false;
		for(const auto & otherSeqPos : iter::range(expSeqs.size())){
			const auto & otherSeq = expSeqs[otherSeqPos];
			if(otherSeq->seq_ == expSeq->seq_){
				found = true;
				otherSeq->name_ += "," + expSeq->name_;
				initialExpSeqsPositions[expSeq->name_] = otherSeqPos;
				break;
			}
		}
		if(!found){
			//if not found add
			//add position of the expected seq
			initialExpSeqsPositions[expSeq->name_] = expSeqs.size();
			expSeqs.emplace_back(std::make_shared<seqInfo>(*expSeq));
		}
	}
	std::unordered_map<std::string, uint32_t> finalExpSeqsPositions;

	for(const auto expSeqPos : iter::range(expSeqs.size())){
		expSeqsKey[expSeqs[expSeqPos]->name_] = expSeqs[expSeqPos]->name_;
		finalExpSeqsPositions[expSeqs[expSeqPos]->name_] = expSeqPos;
	}

	bencher.checkForStrainsThrow(expNames, __PRETTY_FUNCTION__);


	std::unordered_map<std::string, std::unordered_map<std::string, double>> readCountsPerHapPerSample;
	auto sampInfoFnp = analysisMaster.getSampInfoPath();
	TableReader sampInfoReader(TableIOOpts::genTabFileIn(sampInfoFnp, true));
	VecStr row;
	while(sampInfoReader.getNextRow(row)){
		auto sample = row[sampInfoReader.header_.getColPos("s_Name")];
		auto hapName = row[sampInfoReader.header_.getColPos("c_name")];
		auto readCnt = njh::StrToNumConverter::stoToNum<double>(row[sampInfoReader.header_.getColPos("c_AveragedFrac")]);
		readCountsPerHapPerSample[sample][hapName] = readCnt;
	}



	OutputStream haplotypesClassified(njh::files::make_path(setUp.pars_.directoryName_, "classifiedHaplotypes.tab.txt"));
	haplotypesClassified << "AnalysisName\tsample\tmix\tseqName\treadCnt\tfrac\tmatchExpcted\texpectedRef\texpectedFrac";
	OutputStream performanceOut(njh::files::make_path(setUp.pars_.directoryName_, "performancePerTarget.tab.txt"));
	performanceOut << "AnalysisName\tsample\tmix\ttotalReads\trecoveredHaps\tfalseHaps\ttotalHaps\ttotalExpectedHaps\thapRecovery\tfalseHapsRate\tRMSE";
	VecStr metalevels;
	if(nullptr != analysisMaster.groupMetaData_){
		metalevels = getVectorOfMapKeys(analysisMaster.groupMetaData_->groupData_);
		njh::sort(metalevels);
		for(const auto & meta : metalevels){
			haplotypesClassified << "\t" << meta;
			performanceOut << "\t" << meta;
		}
	}

	haplotypesClassified << std::endl;
	performanceOut << std::endl;

	for(const auto & sname : controlSamples){
		//skip completely missing
		if(skipMissingSamples && njh::in(sname, missingSamples)){
			continue;
		}
		//read in result sequences
		auto resultsSeqsFnp = analysisMaster.getSampleFinalHapsPath(sname);
		if(!bfs::exists(resultsSeqsFnp) && skipMissingSamples){
			OutOptions outOptsMissing(njh::files::make_path(setUp.pars_.directoryName_, "missingSamples"));
			outOptsMissing.append_ = true;
			OutputStream outMissing(outOptsMissing);
			outMissing << sname << std::endl;
			continue;
		} else if(!bfs::exists(resultsSeqsFnp)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "missing results for the following sample: " << sname << "\n";
			throw std::runtime_error{ss.str()};
		}

		SeqIOOptions resultsSeqsOpts(resultsSeqsFnp, analysisMaster.inputOptions_.inFormat_, true);
		auto resultSeqs = SeqInput::getSeqVec<seqInfo>(resultsSeqsOpts);
		//get current expected seqs
		std::unordered_map<std::string, double> currentExpectedSeqsFrac;
		for(const auto & expSeqFrac : bencher.mixSetups_.at(bencher.samplesToMix_.at(sname)).relativeAbundances_){
			currentExpectedSeqsFrac[expSeqs[initialExpSeqsPositions[expSeqFrac.first]]->name_] += expSeqFrac.second;
		}

		std::vector<std::shared_ptr<seqInfo>> currentExpectedSeqs;
		for(const auto & finalExp : currentExpectedSeqsFrac){
			currentExpectedSeqs.push_back(expSeqs[finalExpSeqsPositions[finalExp.first]]);
		}

		auto res = bencher.benchmark(resultSeqs, currentExpectedSeqs,
				currentExpectedSeqsFrac,
				expSeqsKey);
		double total = 0;
		for(const auto & seq : resultSeqs){
			total += readCountsPerHapPerSample[sname][seq.name_];
			haplotypesClassified << name
					<< "\t" << sname
					<< "\t" << bencher.samplesToMix_[sname]
					<< "\t" << seq.name_
					<< "\t" << readCountsPerHapPerSample[sname][seq.name_]
					<< "\t" << seq.frac_
					<< "\t" << ("" == res.resSeqToExpSeq_[seq.name_] ? "FALSE": "TRUE")
					<< "\t" << res.resSeqToExpSeq_[seq.name_]
					<< "\t" << currentExpectedSeqsFrac[res.resSeqToExpSeq_[seq.name_]];
			if(nullptr != analysisMaster.groupMetaData_){
				for(const auto & meta : metalevels){
					haplotypesClassified << "\t" << analysisMaster.groupMetaData_->groupData_[meta]->getGroupForSample(sname);
				}
			}
			haplotypesClassified << std::endl;
		}
		performanceOut  << name
				<< "\t" << sname
				<< "\t" << bencher.samplesToMix_[sname]
				<< "\t" << total
				<< "\t" << res.recoveredHaps_
				<< "\t" << res.falseHaps_
				<< "\t" << res.totalHaps()
				<< "\t" << res.expectedHapCnt_
				<< "\t" << res.hapRecoveryRate()
				<< "\t" << res.falseHapRate()
				<< "\t" << res.RMSE();
		if(nullptr != analysisMaster.groupMetaData_){
			for(const auto & meta : metalevels){
				performanceOut << "\t" << analysisMaster.groupMetaData_->groupData_[meta]->getGroupForSample(sname);
			}
		}
		performanceOut << std::endl;
	}

	return 0;
}

}  // namespace njhseq

