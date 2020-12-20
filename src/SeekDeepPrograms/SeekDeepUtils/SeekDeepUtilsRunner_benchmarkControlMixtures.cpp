/*
 * SeekDeepUtilsRunner_benchmarkControlMixtures.cpp
 *
 *  Created on: May 6, 2019
 *      Author: nicholashathaway
 */


#include "SeekDeepUtilsRunner.hpp"
#include "SeekDeep/objects/ControlBenchmarking.h"
#include <njhseq/objects/dataContainers/tables/TableReader.hpp>
namespace njhseq {


int SeekDeepUtilsRunner::benchmarkControlMixtures(
		const njh::progutils::CmdArgs &inputCommands) {
	ControlBencher::ControlBencherPars conBenchPars;
	bfs::path processClustersDir = "";
	bfs::path expectedSeqsFnp = "";
	std::string name = "";
	bfs::path metaFnp = "";
	bool skipMissingSamples = false;
	bool fillInMissingSamples = false;
	bfs::path popSeqsFnp = "";
	comparison allowableError;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(processClustersDir, "--processClustersDir", "Process Clusters Dir", true);
	setUp.setOption(expectedSeqsFnp, "--expectedSeqsFnp", "Expected Seqs fasta file", true);
	setUp.setOption(name, "--name", "Name to give the current analysis", true);
	setUp.setOption(metaFnp, "--metaFnp", "meta data for the control samples");
	setUp.setOption(popSeqsFnp, "--popSeqsFnp", "Population Sequences");
	setUp.setOption(conBenchPars.samplesToMixFnp_, "--sampleToMixture", "Sample To Mixture, 2 columns 1)sample, 2)MixName", true);
	setUp.setOption(conBenchPars.mixSetUpFnp_, "--mixtureSetUp", "Mixture Set Up, 3 columns 1)MixName, 2)strain, 3)relative_abundance", true);
	setUp.setOption(skipMissingSamples, "--skipMissingSamples", "Skip Samples if they are missing");
	setUp.setOption(fillInMissingSamples, "--fillInMissingSamples", "Fill in missing Samples with placeholders");

	setUp.processComparison(allowableError);
	setUp.processAlignerDefualts();
	setUp.processDirectoryOutputName(name + "_benchmarkControlMixtures_TODAY", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	auto coreJsonFnp = njh::files::make_path(processClustersDir,"coreInfo.json");

	njh::files::checkExistenceThrow( { processClustersDir, coreJsonFnp, expectedSeqsFnp,
			conBenchPars.samplesToMixFnp_, conBenchPars.mixSetUpFnp_ },
			__PRETTY_FUNCTION__);

	Json::Value coreJson = njh::json::parseFile(coreJsonFnp.string());

	collapse::SampleCollapseCollection analysisMaster(coreJson);

	std::unordered_map<std::string, std::unordered_map<std::string, double>> readCountsPerHapPerSample;
	std::unordered_map<std::string, std::string> cNameToPopUID;
	std::unordered_map<std::string, uint32_t> cNamePopSampCnt;
	std::unordered_map<std::string, std::string> hPopUID_to_hConsensus;
	{
		auto popInfoFnp = analysisMaster.getPopInfoPath();
		TableReader popInfoReader(TableIOOpts::genTabFileIn(popInfoFnp, true));
		VecStr row;
		while(popInfoReader.getNextRow(row)){
			auto hapName = row[popInfoReader.header_.getColPos("h_popUID")];
			auto h_Consensus = row[popInfoReader.header_.getColPos("h_Consensus")];
			hPopUID_to_hConsensus[hapName] = h_Consensus;
		}
	}
	uint64_t maxLen = 0;
	auto sampInfoFnp = analysisMaster.getSampInfoPath();
	TableReader sampInfoReader(TableIOOpts::genTabFileIn(sampInfoFnp, true));
	VecStr row;
	std::unordered_set<std::string> allSamplesInOutput;

	std::unordered_map<std::string, std::vector<seqInfo>> allResultSeqs;
	while(sampInfoReader.getNextRow(row)){
		auto sample = row[sampInfoReader.header_.getColPos("s_Name")];
		auto hapName = row[sampInfoReader.header_.getColPos("c_name")];
		auto c_ReadCnt = njh::StrToNumConverter::stoToNum<double>(row[sampInfoReader.header_.getColPos("c_ReadCnt")]);
		auto c_AveragedFrac = njh::StrToNumConverter::stoToNum<double>(row[sampInfoReader.header_.getColPos("c_AveragedFrac")]);
		auto readCnt = njh::StrToNumConverter::stoToNum<double>(row[sampInfoReader.header_.getColPos("c_ReadCnt")]);
		readCountsPerHapPerSample[sample][hapName] = readCnt;
		auto h_popUID = row[sampInfoReader.header_.getColPos("h_popUID")];
		auto h_SampCnt = row[sampInfoReader.header_.getColPos("h_SampCnt")];
		seqInfo clus(hapName, njh::mapAt(hPopUID_to_hConsensus, h_popUID));
		clus.cnt_ = c_ReadCnt;
		clus.frac_ = c_AveragedFrac;
		bool add = true;
		for(auto & seq : allResultSeqs[sample]){
			if(seq.seq_ == clus.seq_){
				add = false;
				seq.cnt_ += c_ReadCnt;
				seq.frac_ += c_AveragedFrac;
				break;
			}
		}
		if(add){
			allResultSeqs[sample].emplace_back(clus);
		}

		readVec::getMaxLength(clus, maxLen);
		cNameToPopUID[hapName] = h_popUID;
		cNamePopSampCnt[hapName] = njh::StrToNumConverter::stoToNum<uint32_t>(h_SampCnt);
		allSamplesInOutput.emplace(sample);
	}

	ControlBencher bencher(conBenchPars);
	//check for samples in analysis
	auto controlSamples = bencher.getSamples();
	VecStr missingSamples;
	for (const auto & sname : controlSamples) {
		//if (!analysisMaster.hasSample(sname)) {
		if (!njh::in(sname, allSamplesInOutput)) {
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

	//population seqs;
	std::vector<seqInfo> popSeqs;
	if("" != popSeqsFnp){
		popSeqs = SeqInput::getSeqVec<seqInfo>(SeqIOOptions::genFastaIn(popSeqsFnp));
	}
	//check the needed expected names
	std::set<std::string> expNames;
	std::unordered_map<std::string, uint32_t> initialExpSeqsPositions;
	std::unordered_map<std::string, std::string> expSeqsKey;
	std::unordered_map<std::string, VecStr> matchingExpectedSeqs;

	std::vector<std::shared_ptr<seqInfo>> expSeqs;

	VecStr removeSeqs;
	for(const auto & expSeq : initialExpSeqs){
		if (std::all_of(expSeq->seq_.begin(), expSeq->seq_.end(),
				[](const char base) {
					return 'N' == base || 'n' == base;
				})) {
			removeSeqs.emplace_back(expSeq->name_);
			continue;
		}
		if(njh::in(expSeq->name_, expNames)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << " already have expected sequence named " << expSeq->name_<< "\n";
			throw std::runtime_error{ss.str()};
		}
		expNames.emplace(expSeq->name_);
		bool found = false;
		for(const auto otherSeqPos : iter::range(expSeqs.size())){
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
	//remove ref seqs that were zero;
	bencher.removeStrains(removeSeqs);


	std::unordered_map<std::string, uint32_t> finalExpSeqsPositions;

	//

	for(const auto expSeqPos : iter::range(expSeqs.size())){
		readVec::getMaxLength(expSeqs[expSeqPos], maxLen);
		expSeqsKey[expSeqs[expSeqPos]->name_] = expSeqs[expSeqPos]->name_;
		finalExpSeqsPositions[expSeqs[expSeqPos]->name_] = expSeqPos;
	}

	bencher.checkForStrainsThrow(expNames, __PRETTY_FUNCTION__);


	OutputStream falseHaplotypesToExpClassified(njh::files::make_path(setUp.pars_.directoryName_, "falseHaplotypesComparedToExpected.tab.txt"));
	falseHaplotypesToExpClassified << "AnalysisName\tsample\tmix\tc_name\treadCnt\tfrac\tRefName\tExpectedRefFreq\tExpectedMajor\tscore\tbestMatchScore\tmismatches\toneBaseIndels\ttwoBaseIndels\tlargeIndels\ttotalErrors";
	OutputStream falseHaplotypesToOtherResultsClassified(njh::files::make_path(setUp.pars_.directoryName_, "falseHaplotypesComparedToOthers.tab.txt"));
	falseHaplotypesToOtherResultsClassified << "AnalysisName\tsample\tmix\tc_name\treadCnt\tfrac\tOtherName\tOtherReadCnt\tOtherFrac\tratio\tOtherMajor\tOtherMatchesExpected\tOtherExpectedMatchName\tscore\tbestMatchScore\tmismatches\toneBaseIndels\ttwoBaseIndels\tlargeIndels\ttotalErrors";


	OutputStream haplotypesClassified(njh::files::make_path(setUp.pars_.directoryName_, "classifiedHaplotypes.tab.txt"));
	haplotypesClassified << "AnalysisName\tsample\tmix\tc_name\th_popUID\th_SampCnt\treadCnt\tfrac\tmatchExpcted\texpectedRef\texpectedFrac\tMajorOrMinor\tmatchingPopulation\tPopName";
	OutputStream performanceOut(njh::files::make_path(setUp.pars_.directoryName_, "performancePerTarget.tab.txt"));
	performanceOut << "AnalysisName\tsample\tmix\ttotalReads\trecoveredHaps\tfalseHaps\ttotalHaps\ttotalExpectedHaps\thapRecovery\tfalseHapsRate\tmissingExpecteds\tRMSE";
	VecStr metalevels;
	if(nullptr != analysisMaster.groupMetaData_){
		metalevels = getVectorOfMapKeys(analysisMaster.groupMetaData_->groupData_);
		njh::sort(metalevels);
		for(const auto & meta : metalevels){
			haplotypesClassified << "\t" << meta;
			performanceOut << "\t" << meta;
			falseHaplotypesToExpClassified << "\t" << meta;
			falseHaplotypesToOtherResultsClassified << "\t" << meta;
		}
	}

	haplotypesClassified << std::endl;
	performanceOut << std::endl;
	falseHaplotypesToExpClassified << std::endl;
	falseHaplotypesToOtherResultsClassified << std::endl;
	aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0));
	alignerObj.weighHomopolymers_ = false;
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, false);
	for(const auto & sname : controlSamples){
		//skip completely missing
		if(skipMissingSamples && njh::in(sname, missingSamples)){
			continue;
		}
		//read in result sequences
		//auto resultsSeqsFnp = analysisMaster.getSampleFinalHapsPath(sname);
		if(!njh::in(sname, allSamplesInOutput) && skipMissingSamples){
			OutOptions outOptsMissing(njh::files::make_path(setUp.pars_.directoryName_, "missingSamples.txt"));
			outOptsMissing.append_ = true;
			OutputStream outMissing(outOptsMissing);
			outMissing << sname << std::endl;
			continue;
		} else if(!njh::in(sname, allSamplesInOutput)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "missing results for the following sample: " << sname << "\n";
			throw std::runtime_error{ss.str()};
		}


		std::vector<seqInfo> resultSeqs = allResultSeqs[sname];

		std::unordered_map<std::string, uint32_t> resSeqToPos;
		double maxResFrac = 0;
		for(const auto pos : iter::range(resultSeqs.size())){
			resSeqToPos[resultSeqs[pos].name_] = pos;
			if(resultSeqs[pos].frac_ > maxResFrac){
				maxResFrac = resultSeqs[pos].frac_;
			}
		}
		//get current expected seqs
		std::unordered_map<std::string, double> currentExpectedSeqsFrac;
		double maxExpFrac = 0;
		std::unordered_map<std::string, VecStr> expectedSeqNameToCurrentSeqsKey;

		for(const auto & expSeqFrac : bencher.mixSetups_.at(bencher.samplesToMix_.at(sname)).relativeAbundances_){
			currentExpectedSeqsFrac[expSeqs[initialExpSeqsPositions[expSeqFrac.first]]->name_] += expSeqFrac.second;
			expectedSeqNameToCurrentSeqsKey[expSeqs[initialExpSeqsPositions[expSeqFrac.first]]->name_].emplace_back(expSeqFrac.first);
		}
		for(auto & key : expectedSeqNameToCurrentSeqsKey){
//			std::cout << "key: " << key.first << std::endl;
//			std::cout << "\t" << njh::conToStr(key.second, ",") << std::endl;
			njh::sort(key.second);
		}
		for(const auto & expFrac : currentExpectedSeqsFrac){
			if(expFrac.second > maxExpFrac){
				maxExpFrac = expFrac.second;
			}
		}

		std::unordered_map<std::string, std::string> expectedToMajorClass;
		for(const auto & expFrac : currentExpectedSeqsFrac){
			if(expFrac.second == maxExpFrac){
				expectedToMajorClass[expFrac.first] = "major";
			}else{
				expectedToMajorClass[expFrac.first] = "minor";
			}
		}
		std::unordered_map<std::string, std::string> resultsToMajorClass;
		for(const auto & res : resultSeqs){
			if(res.frac_ == maxResFrac){
				resultsToMajorClass[res.name_] = "major";
			}else{
				resultsToMajorClass[res.name_] = "minor";
			}
		}


		std::vector<std::shared_ptr<seqInfo>> currentExpectedSeqs;
		for(const auto & finalExp : currentExpectedSeqsFrac){
			currentExpectedSeqs.push_back(expSeqs[finalExpSeqsPositions[finalExp.first]]);
		}

		auto res = bencher.benchmark(resultSeqs, currentExpectedSeqs,
				currentExpectedSeqsFrac,
				expSeqsKey, alignerObj, allowableError);
		double total = 0;


		//haplotype classification
		for(const auto & seq : resultSeqs){
			total += readCountsPerHapPerSample[sname][seq.name_];
			haplotypesClassified << name
					<< "\t" << sname
					<< "\t" << bencher.samplesToMix_[sname]
					<< "\t" << seq.name_
					<< "\t" << cNameToPopUID[seq.name_]
					<< "\t" << cNamePopSampCnt[seq.name_]
					<< "\t" << readCountsPerHapPerSample[sname][seq.name_]
					<< "\t" << seq.frac_
					<< "\t" << ("" == res.resSeqToExpSeq_[seq.name_] ? "FALSE": "TRUE")
					<< "\t" << njh::conToStr(expectedSeqNameToCurrentSeqsKey[res.resSeqToExpSeq_[seq.name_]], ",")
					<< "\t" << currentExpectedSeqsFrac[res.resSeqToExpSeq_[seq.name_]]
					<< "\t" << expectedToMajorClass[res.resSeqToExpSeq_[seq.name_]];

			std::string matchingPop = "";
			for(const auto & popSeq : popSeqs){
				if(popSeq.seq_ == seq.seq_){
					matchingPop = popSeq.name_;
					break;
				}
			}
			haplotypesClassified
			    << "\t" << ("" == matchingPop ? "TRUE" : "FALSE")
					<< "\t" << matchingPop;
			if(nullptr != analysisMaster.groupMetaData_){
				for(const auto & meta : metalevels){
					haplotypesClassified << "\t" << analysisMaster.groupMetaData_->groupData_[meta]->getGroupForSample(sname);
				}
			}
			haplotypesClassified << std::endl;
		}
		for(const auto & missing : res.missingExpecteds_){
			haplotypesClassified << name
					<< "\t" << sname
					<< "\t" << bencher.samplesToMix_[sname]
					<< "\t" << "NA"
					<< "\t" << "NA"
					<< "\t" << "NA"
					<< "\t" << "NA"
					<< "\t" << "NA"
					<< "\t" << "NA"
					<< "\t" << njh::conToStr(expectedSeqNameToCurrentSeqsKey[missing], ",")
					<< "\t" << currentExpectedSeqsFrac[missing]
					<< "\t" << expectedToMajorClass[missing];

			haplotypesClassified
					<< "\t" << "NA"
					<< "\t" << "NA";
			if(nullptr != analysisMaster.groupMetaData_){
				for(const auto & meta : metalevels){
					haplotypesClassified << "\t" << analysisMaster.groupMetaData_->groupData_[meta]->getGroupForSample(sname);
				}
			}
			haplotypesClassified << std::endl;
		}
		//performance
		VecStr missingExpectedDecoded;
		for(const auto & missing: res.missingExpecteds_){
			missingExpectedDecoded.emplace_back(njh::conToStr(expectedSeqNameToCurrentSeqsKey[missing], ","));
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
				<< "\t" << njh::conToStr(missingExpectedDecoded, ";")
				<< "\t" << res.RMSE();
		if(nullptr != analysisMaster.groupMetaData_){
			for(const auto & meta : metalevels){
				performanceOut << "\t" << analysisMaster.groupMetaData_->groupData_[meta]->getGroupForSample(sname);
			}
		}
		performanceOut << std::endl;


		//further classification of false haplotypes to expected
		for(const auto & falseHap : res.falseHapsCompsToExpected){
			double bestScore = 0;
			for(const auto & exp : falseHap.second){
				if(exp.second.distances_.eventBasedIdentityHq_ > bestScore){
					bestScore = exp.second.distances_.eventBasedIdentityHq_;
				}
			}
			for(const auto & exp : falseHap.second){
				falseHaplotypesToExpClassified<< name
						<< "\t" << sname
						<< "\t" << bencher.samplesToMix_[sname]
						<< "\t" << falseHap.first
						<< "\t" << readCountsPerHapPerSample[sname][falseHap.first]
						<< "\t" << resultSeqs[resSeqToPos[falseHap.first]].frac_
						<< "\t" << exp.first
						<< "\t" << currentExpectedSeqsFrac[exp.first]
						<< "\t" << expectedToMajorClass[exp.first]
						<< "\t" << exp.second.distances_.eventBasedIdentityHq_
						<< "\t" << (exp.second.distances_.eventBasedIdentityHq_ == bestScore ? "TRUE" : "FALSE")
						<< "\t" << exp.second.hqMismatches_ + exp.second.lqMismatches_ + exp.second.lowKmerMismatches_
						<< "\t" << exp.second.oneBaseIndel_
						<< "\t" << exp.second.twoBaseIndel_
						<< "\t" << exp.second.largeBaseIndel_
						<< "\t" << exp.second.distances_.getNumOfEvents(true);
				if(nullptr != analysisMaster.groupMetaData_){
					for(const auto & meta : metalevels){
						falseHaplotypesToExpClassified << "\t" << analysisMaster.groupMetaData_->groupData_[meta]->getGroupForSample(sname);
					}
				}
				falseHaplotypesToExpClassified << std::endl;
			}
		}

		//further classification of false haplotypes to other sequences
		for(const auto & falseHap : res.falseHapsCompsToOthers){
			double bestScore = 0;
			for(const auto & other : falseHap.second){
				if(other.second.distances_.eventBasedIdentityHq_ > bestScore){
					bestScore = other.second.distances_.eventBasedIdentityHq_;
				}
			}
			for(const auto & other : falseHap.second){
				falseHaplotypesToOtherResultsClassified<< name
						<< "\t" << sname
						<< "\t" << bencher.samplesToMix_[sname]
						<< "\t" << falseHap.first
						<< "\t" << readCountsPerHapPerSample[sname][falseHap.first]
						<< "\t" << resultSeqs[resSeqToPos[falseHap.first]].frac_
						<< "\t" << other.first
						<< "\t" << readCountsPerHapPerSample[sname][other.first]
						<< "\t" << resultSeqs[resSeqToPos[other.first]].frac_
						<< "\t" << resultSeqs[resSeqToPos[other.first]].frac_/resultSeqs[resSeqToPos[falseHap.first]].frac_
						<< "\t" << resultsToMajorClass[other.first]
						<< "\t" << ("" == res.resSeqToExpSeq_[other.first] ? "FALSE" : "TRUE")
						<< "\t" << res.resSeqToExpSeq_[other.first]
						<< "\t" << other.second.distances_.eventBasedIdentityHq_
						<< "\t" << (other.second.distances_.eventBasedIdentityHq_ == bestScore ? "TRUE" : "FALSE")
						<< "\t" << other.second.hqMismatches_ + other.second.lqMismatches_ + other.second.lowKmerMismatches_
						<< "\t" << other.second.oneBaseIndel_
						<< "\t" << other.second.twoBaseIndel_
						<< "\t" << other.second.largeBaseIndel_
						<< "\t" << other.second.distances_.getNumOfEvents(true);
				if(nullptr != analysisMaster.groupMetaData_){
					for(const auto & meta : metalevels){
						falseHaplotypesToOtherResultsClassified << "\t" << analysisMaster.groupMetaData_->groupData_[meta]->getGroupForSample(sname);
					}
				}
				falseHaplotypesToOtherResultsClassified << std::endl;
			}
		}
	}

	if(fillInMissingSamples){
		for(const auto & sname : missingSamples){

			std::vector<seqInfo> resultSeqs = allResultSeqs[sname];

			std::unordered_map<std::string, uint32_t> resSeqToPos;
			double maxResFrac = 0;
			for(const auto pos : iter::range(resultSeqs.size())){
				resSeqToPos[resultSeqs[pos].name_] = pos;
				if(resultSeqs[pos].frac_ > maxResFrac){
					maxResFrac = resultSeqs[pos].frac_;
				}
			}
			//get current expected seqs
			std::unordered_map<std::string, double> currentExpectedSeqsFrac;
			double maxExpFrac = 0;
			std::unordered_map<std::string, VecStr> expectedSeqNameToCurrentSeqsKey;


			for(const auto & expSeqFrac : bencher.mixSetups_.at(bencher.samplesToMix_.at(sname)).relativeAbundances_){
				currentExpectedSeqsFrac[expSeqs[initialExpSeqsPositions[expSeqFrac.first]]->name_] += expSeqFrac.second;
				expectedSeqNameToCurrentSeqsKey[expSeqs[initialExpSeqsPositions[expSeqFrac.first]]->name_].emplace_back(expSeqFrac.first);
			}
			for(auto & key : expectedSeqNameToCurrentSeqsKey){
				njh::sort(key.second);
			}
			for(const auto & expFrac : currentExpectedSeqsFrac){
				if(expFrac.second > maxExpFrac){
					maxExpFrac = expFrac.second;
				}
			}

			std::unordered_map<std::string, std::string> expectedToMajorClass;
			for(const auto & expFrac : currentExpectedSeqsFrac){
				if(expFrac.second == maxExpFrac){
					expectedToMajorClass[expFrac.first] = "major";
				}else{
					expectedToMajorClass[expFrac.first] = "minor";
				}
			}
			std::unordered_map<std::string, std::string> resultsToMajorClass;
			for(const auto & res : resultSeqs){
				if(res.frac_ == maxResFrac){
					resultsToMajorClass[res.name_] = "major";
				}else{
					resultsToMajorClass[res.name_] = "minor";
				}
			}


			std::vector<std::shared_ptr<seqInfo>> currentExpectedSeqs;
			for(const auto & finalExp : currentExpectedSeqsFrac){
				currentExpectedSeqs.push_back(expSeqs[finalExpSeqsPositions[finalExp.first]]);
			}
			VecStr missingExpectedDecoded;
			for(const auto & missing: currentExpectedSeqs){
				missingExpectedDecoded.emplace_back(njh::conToStr(expectedSeqNameToCurrentSeqsKey[missing->name_], ","));
			}
			//performance
			performanceOut  << name
					<< "\t" << sname
					<< "\t" << bencher.samplesToMix_[sname]
					<< "\t" << 0
					<< "\t" << 0
					<< "\t" << 0
					<< "\t" << 0
					<< "\t" << currentExpectedSeqs.size()
					<< "\t" << 0
					<< "\t" << 0
					<< "\t" << njh::conToStr(missingExpectedDecoded, ";")
					<< "\t" << 0;
			if(nullptr != analysisMaster.groupMetaData_){
				for(const auto & meta : metalevels){
					performanceOut << "\t" << analysisMaster.groupMetaData_->groupData_[meta]->getGroupForSample(sname);
				}
			}
			performanceOut << std::endl;
			for(const auto & missing : readVec::getNames(currentExpectedSeqs)){
				haplotypesClassified << name
						<< "\t" << sname
						<< "\t" << bencher.samplesToMix_[sname]
						<< "\t" << "NA"
						<< "\t" << "NA"
						<< "\t" << "NA"
						<< "\t" << "NA"
						<< "\t" << "NA"
						<< "\t" << "NA"
						<< "\t" << njh::conToStr(expectedSeqNameToCurrentSeqsKey[missing], ",")
						<< "\t" << currentExpectedSeqsFrac[missing]
						<< "\t" << expectedToMajorClass[missing];

				haplotypesClassified
						<< "\t" << "NA"
						<< "\t" << "NA";
				if(nullptr != analysisMaster.groupMetaData_){
					for(const auto & meta : metalevels){
						haplotypesClassified << "\t" << analysisMaster.groupMetaData_->groupData_[meta]->getGroupForSample(sname);
					}
				}
				haplotypesClassified << std::endl;
			}
		}
	}

	alignerObj.processAlnInfoOutput(setUp.pars_.outAlnInfoDirName_, false);

	return 0;
}

}  // namespace njhseq

