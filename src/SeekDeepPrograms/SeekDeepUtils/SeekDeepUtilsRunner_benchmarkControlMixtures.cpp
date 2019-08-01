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
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(processClustersDir, "--processClustersDir", "Process Clusters Dir", true);
	setUp.setOption(expectedSeqsFnp, "--expectedSeqsFnp", "Expected Seqs fasta file", true);
	setUp.setOption(name, "--name", "Name to give the current analysis", true);
	setUp.setOption(metaFnp, "--metaFnp", "meta data for the control samples");

	setUp.setOption(conBenchPars.samplesToMixFnp_, "--sampleToMixture", "Sample To Mixture, 2 columns 1)sample, 2)MixName", true);
	setUp.setOption(conBenchPars.mixSetUpFnp_, "--mixtureSetUp", "Mixture Set Up, 3 columns 1)MixName, 2)strain, 3)relative_abundance", true);

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
	if (!missingSamples.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error "
				<< "missing the following samples in " << processClustersDir << ": "
				<< njh::conToStrEndSpecial(missingSamples, ", ", " and ") << "\n";
		throw std::runtime_error { ss.str() };
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


	OutputStream haplotypesClassified(njh::files::make_path(setUp.pars_.directoryName_, "classifiedHaplotypes.tab.txt"));
	haplotypesClassified << "AnalysisName\tsample\tmix\tseqName\tfrac\tmatchExpcted\texpectedRef\texpectedFrac";
	OutputStream performanceOut(njh::files::make_path(setUp.pars_.directoryName_, "performancePerTarget.tab.txt"));
	performanceOut << "AnalysisName\tsample\tmix\trecoveredHaps\tfalseHaps\ttotalHaps\ttotalExpectedHaps\thapRecovery\tfalseHapsRate\tRMSE";
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
		//read in result sequences
		auto resultsSeqsFnp = analysisMaster.getSampleFinalHapsPath(sname);
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

		for(const auto & seq : resultSeqs){
			haplotypesClassified << name
					<< "\t" << sname
					<< "\t" << bencher.samplesToMix_[sname]
					<< "\t" << seq.name_
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

