/*
 * gatherInfoOnTargetedAmpliconSeqFile.cpp
 *
 *  Created on: Oct 16, 2019
 *      Author: nicholashathaway
 */


#include "SeekDeepUtilsRunner.hpp"

#include "SeekDeep/objects.h"
#include "SeekDeep/parameters.h"


namespace njhseq {
int SeekDeepUtilsRunner::gatherInfoOnTargetedAmpliconSeqFile(
		const njh::progutils::CmdArgs & inputCommands) {

	uint32_t unrecogBaseSampling = 20;
	uint32_t precdingBaseFreqCutOff = 5;
	bfs::path idFnp = "";
	bool dontCollapsePossibleMIDs = false;
	ExtractorPairedEndPars pars;
	uint32_t testNumber = std::numeric_limits<uint32_t>::max();
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(testNumber, "--testNumber", "Just use this number of reads of the top of the file");
	setUp.setOption(dontCollapsePossibleMIDs, "--dontCollapsePossibleMIDs",
			"Don't Collapse Possible MIDs", false);
	setUp.setOption(unrecogBaseSampling, "--unrecogBaseSampling",
			"Number of bases to sample from file for unrecognized sequences", false);

	setUp.setOption(precdingBaseFreqCutOff, "--precdingBaseFreqCutOff", "Preceding Base Freq Cut Off", false);

	pars.corePars_.pDetPars.primerWithin_ = 30;
	setUp.setOption(pars.corePars_.pDetPars.primerWithin_, "--primerWithin", "Primer Within bases search", false, "Primer");


	bool primerToUpperCase = false;
	setUp.setOption(primerToUpperCase, "--primerUpper",
			"Leave primers in upper case", false, "Primer");
	pars.corePars_.pDetPars.primerToLowerCase_ = !primerToUpperCase;
	setUp.setOption(pars.corePars_.pDetPars.allowable_.distances_.query_.coverage_, "--primerCoverage",
			"Amount of primers found", false, "Primer");
	setUp.setOption(pars.corePars_.pDetPars.allowable_.hqMismatches_, "--primerNumOfMismatches",
			"Number of Mismatches to allow in primers", false, "Primer");
	setUp.setOption(pars.corePars_.pDetPars.allowable_.oneBaseIndel_, "--primerOneBaseIndels",
			"Number Of One base indels to allow in primers", false, "Primer");
	setUp.setOption(pars.corePars_.pDetPars.allowable_.twoBaseIndel_, "--primerTwoBaseIndels",
			"Number Of Two base indels to allow in primers", false, "Primer");

	setUp.pars_.gapInfo_.gapOpen_ = 5;
	setUp.pars_.gapInfo_.gapExtend_ = 1;
	setUp.pars_.gap_ = "5,1";
	setUp.pars_.gapInfo_.gapRightQueryOpen_ = 0;
	setUp.pars_.gapInfo_.gapRightQueryExtend_ = 0;
	setUp.pars_.gapInfo_.gapRightRefOpen_ = 0;
	setUp.pars_.gapInfo_.gapRightRefExtend_ = 0;
	setUp.pars_.gapRight_ = "0,0";
	setUp.pars_.gapInfo_.gapLeftQueryOpen_ = 0;
	setUp.pars_.gapInfo_.gapLeftQueryExtend_ = 0;
	setUp.pars_.gapInfo_.gapLeftRefOpen_ = 0;
	setUp.pars_.gapInfo_.gapLeftRefExtend_ = 0;
	setUp.pars_.gapLeft_ = "0,0";
	setUp.processGap();

	setUp.setOption(idFnp, "--id", "SeekDeep primers file", true);
	setUp.processReadInNames(VecStr{"--fastq1", "--fastq", "--fasta", "--fastq1gz", "--fastqgz", "--fastagz"});
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	PrimersAndMids ids(idFnp);
	ids.initPrimerDeterminator();

	uint64_t maxReadSize = 0;
	uint32_t readCount = 0;
	if (setUp.pars_.ioOptions_.isPairedIn()) {
		PairedRead seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		while(reader.readNextRead(seq)){
			++readCount;
			readVec::getMaxLength(seq.seqBase_, maxReadSize);
			readVec::getMaxLength(seq.mateSeqBase_, maxReadSize);
			if(setUp.pars_.verbose_ && readCount % 10000 == 0){
				std::cout << "\r" << readCount;
				std::cout.flush();
			}
			if(readCount >=testNumber){
				break;
			}
		}
		if(setUp.pars_.verbose_){
			std::cout << std::endl;
		}
	} else {
		seqInfo seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		while(reader.readNextRead(seq)){
			++readCount;
			readVec::getMaxLength(seq, maxReadSize);
			if(setUp.pars_.verbose_ && readCount % 10000 == 0){
				std::cout << "\r" << readCount;
				std::cout.flush();
			}
			if(readCount >=testNumber){
				break;
			}
		}
		if(setUp.pars_.verbose_){
			std::cout << std::endl;
		}
	}

	// create aligner for primer identification
	auto scoreMatrix = substituteMatrix::createDegenScoreMatrixNoNInRef(
			setUp.pars_.generalMatch_, setUp.pars_.generalMismatch_);
	gapScoringParameters gapPars(setUp.pars_.gapInfo_);
	KmerMaps emptyMaps;
	bool countEndGaps = false;
	//to avoid allocating an extremely large aligner matrix;

	aligner alignObj(maxReadSize, gapPars, scoreMatrix, emptyMaps,
			setUp.pars_.qScorePars_, countEndGaps, false);


	if (setUp.pars_.ioOptions_.isPairedIn()) {
		PairedRead seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();

		//key = forward primer name, reverse primer name, count
		std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> primerPairCountsTot;
		std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> primerPairCountsFor;
		std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> primerPairCountsComp;

		//key = primer name, forward primer bases, reverse primer bases, counts
		std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>>> precedingBasesCounts;
		std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>>> precedingBasesCountsComp;

		//key
		std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> unrecognizedCounts;
		njh::ProgressBar pBar(readCount);
		uint32_t newReadCount = 0;
		while(reader.readNextRead(seq)){
			++newReadCount;
			if(setUp.pars_.verbose_){
				pBar.outputProgAdd(std::cout, 1, true);
			}
			std::string forwardPrimerName = "";
			std::string forwardPrimerPrecedingBases = "";

			std::string reversePrimerName = "";
			std::string reversePrimerPrcedingBases = "";

			bool complement = false;
			auto forPos = ids.pDeterminator_->determineBestForwardPrimerPosFront(seq.seqBase_, pars.corePars_.pDetPars, alignObj);
			if("" == forPos.primerName_){
				auto forMate = ids.pDeterminator_->determineBestForwardPrimerPosFront(seq.mateSeqBase_, pars.corePars_.pDetPars, alignObj);
				if("" != forMate.primerName_){
					forwardPrimerName = forMate.primerName_;
					if(0 != forMate.start_ ){
						forwardPrimerPrecedingBases = seq.mateSeqBase_.seq_.substr(0, forMate.start_);
					}
					complement = true;
				}else{
					forwardPrimerName = "unrecognized";
				}
			}else{
				forwardPrimerName = forPos.primerName_;
				if(0 != forPos.start_ ){
					forwardPrimerPrecedingBases = seq.seqBase_.seq_.substr(0, forPos.start_);
				}
			}

			if("unrecognized" == forwardPrimerName){
				auto revPos = ids.pDeterminator_->determineBestReversePrimerPosFront(seq.mateSeqBase_, pars.corePars_.pDetPars, alignObj);
				if("" == revPos.primerName_){
					auto revPosOther = ids.pDeterminator_->determineBestReversePrimerPosFront(seq.seqBase_, pars.corePars_.pDetPars, alignObj);
					if("" != revPosOther.primerName_){
						//matched
						complement = true;
						reversePrimerName = revPosOther.primerName_;
						if(0 != revPos.start_ ){
							reversePrimerPrcedingBases = seq.seqBase_.seq_.substr(0, revPos.start_);
						}
					}else{
						//unmatched
						reversePrimerName = "unrecognized";
					}
				}else{
					//matched
					reversePrimerName = revPos.primerName_;
					if(0 != revPos.start_ ){
						reversePrimerPrcedingBases = seq.mateSeqBase_.seq_.substr(0, revPos.start_);
					}
				}
			}else{
				if(complement){
					auto revPos = ids.pDeterminator_->determineBestReversePrimerPosFront(seq.seqBase_, pars.corePars_.pDetPars, alignObj);
					if("" != revPos.primerName_){
						//matched
						reversePrimerName = revPos.primerName_;
						if(0 != revPos.start_ ){
							reversePrimerPrcedingBases = seq.seqBase_.seq_.substr(0, revPos.start_);
						}
					}else{
						reversePrimerName = "unrecognized";
					}
				}else{
					auto revPos = ids.pDeterminator_->determineBestReversePrimerPosFront(seq.mateSeqBase_, pars.corePars_.pDetPars, alignObj);
					if("" != revPos.primerName_){
						//matched
						reversePrimerName = revPos.primerName_;
						if(0 != revPos.start_ ){
							reversePrimerPrcedingBases = seq.mateSeqBase_.seq_.substr(0, revPos.start_);
						}
					}else{
						reversePrimerName = "unrecognized";
					}
				}
			}
			if(complement){
				primerPairCountsComp[forwardPrimerName][reversePrimerName]+= seq.seqBase_.cnt_;
			} else {
				primerPairCountsFor[forwardPrimerName][reversePrimerName]+= seq.seqBase_.cnt_;
			}
			primerPairCountsTot[forwardPrimerName][reversePrimerName]+= seq.seqBase_.cnt_;
			if(forwardPrimerName == reversePrimerName && "unrecognized" != forwardPrimerName){

			}
			if(forwardPrimerName == reversePrimerName && "unrecognized" != forwardPrimerName){
				if(complement){
					precedingBasesCountsComp[forwardPrimerName][forwardPrimerPrecedingBases][reversePrimerPrcedingBases] += seq.seqBase_.cnt_;
				}else{
					precedingBasesCounts[forwardPrimerName][forwardPrimerPrecedingBases][reversePrimerPrcedingBases] += seq.seqBase_.cnt_;
				}
			}else if(forwardPrimerName == reversePrimerName && "unrecognized" == forwardPrimerName){
				unrecognizedCounts[seq.seqBase_.seq_.substr(0, unrecogBaseSampling)][seq.mateSeqBase_.seq_.substr(0, unrecogBaseSampling)] += seq.seqBase_.cnt_;
			}
			if(newReadCount >=testNumber){
				break;
			}
		}

		//primer pairing counts
		table primerCountsTab(VecStr{"ForwardPrimer", "ReversePrimer", "ForwardCount", "ReverseCount", "ForwardFraction", "Total", "Fraction"});
		for(const auto & forPrim : primerPairCountsTot){
			for(const auto & revPrim : forPrim.second){
				primerCountsTab.addRow(forPrim.first, revPrim.first,
						primerPairCountsFor[forPrim.first][revPrim.first],
						primerPairCountsComp[forPrim.first][revPrim.first],
						primerPairCountsFor[forPrim.first][revPrim.first]/static_cast<double>(revPrim.second),
						revPrim.second,
						revPrim.second/static_cast<double>(readCount));
			}
		}
		auto primerTabOpts = TableIOOpts::genTabFileOut(njh::files::make_path(setUp.pars_.directoryName_, "primerCounts.tab.txt"));
		primerCountsTab.sortTable("ForwardPrimer", "ReversePrimer", false);
		primerCountsTab.outPutContents(primerTabOpts);

		//preceding number of bases counts
		table precedingBasesCountsTab(VecStr{"PrimerPair", "BasesPrecedingPrimerCase", "NumOfBases", "Count"});
		std::unordered_map<std::string, uint32_t> mostCommonForwardPreBases;
		std::unordered_map<std::string, uint32_t> mostCommonReversePreBases;

		for(const auto & primerPair : ids.targets_){
			std::unordered_map<std::string,std::unordered_map<uint32_t, uint32_t>> counts;
			std::unordered_map<uint32_t, uint32_t> forwardTotalCounts;
			std::unordered_map<uint32_t, uint32_t> reverseTotalCounts;

			for(const auto & forward : precedingBasesCounts[primerPair.first]){
				for(const auto & reverse : forward.second){
					forwardTotalCounts[forward.first.size()] += reverse.second;
					reverseTotalCounts[reverse.first.size()] += reverse.second;
					counts["ForwardPrimerPrecedingForward"][forward.first.size()] += reverse.second;
					counts["ReversePrimerPrecedingForward"][reverse.first.size()] += reverse.second;
				}
			}
			for(const auto & forward : precedingBasesCountsComp[primerPair.first]){
				for(const auto & reverse : forward.second){
					forwardTotalCounts[forward.first.size()] += reverse.second;
					reverseTotalCounts[reverse.first.size()] += reverse.second;
					counts["ForwardPrimerPrecedingReverse"][forward.first.size()] += reverse.second;
					counts["ReversePrimerPrecedingReverse"][reverse.first.size()] += reverse.second;
				}
			}
			{
				uint32_t highestCount = 0;
				for(const auto & forw : forwardTotalCounts){
					if(forw.second > highestCount){
						highestCount = forw.second;
						mostCommonForwardPreBases[primerPair.first] = forw.first;
					}
				}
			}
			{
				uint32_t highestCount = 0;
				for(const auto & rev : reverseTotalCounts){
					if(rev.second > highestCount){
						highestCount = rev.second;
						mostCommonReversePreBases[primerPair.first] = rev.first;
					}
				}
			}
			for (const auto & count : counts) {
				for (const auto & numBases : count.second) {
					precedingBasesCountsTab.addRow(primerPair.first, count.first, numBases.first, numBases.second);
				}
			}
		}
		auto precedingBasesCountsTabOpts = TableIOOpts::genTabFileOut(njh::files::make_path(setUp.pars_.directoryName_, "precedingBasesCounts.tab.txt"));
		precedingBasesCountsTab.sortTable("PrimerPair", "BasesPrecedingPrimerCase", "NumOfBases", false);
		precedingBasesCountsTab.outPutContents(precedingBasesCountsTabOpts);

		//preceding bases counts, to determine possible mids
		table possibleMidCounts(VecStr{"PrimerPair", "ForwardMID", "ReverseMID", "InForDirCount", "InRevDirCount", "TotalCount"});
		for(const auto & primerPair : ids.targets_){
			std::string primerPairName = primerPair.first;
			std::unordered_map<std::string, std::unordered_map<std::string, double>> totalCounts;
			for(const auto & forwMid : precedingBasesCounts[primerPairName]){
				for(const auto & revMid : forwMid.second){
					totalCounts[forwMid.first][revMid.first] += revMid.second;
				}
			}
			for(const auto & forwMid : precedingBasesCountsComp[primerPairName]){
				for(const auto & revMid : forwMid.second){
					totalCounts[forwMid.first][revMid.first] += revMid.second;
				}
			}
			if(dontCollapsePossibleMIDs){
				for(const auto & forwMid : totalCounts){
					for(const auto & revMid : forwMid.second){
						if(revMid.second > precdingBaseFreqCutOff){
							possibleMidCounts.addRow(primerPairName,
									forwMid.first, revMid.first,
									precedingBasesCounts[primerPairName][forwMid.first][revMid.first],
									precedingBasesCountsComp[primerPairName][forwMid.first][revMid.first],
									revMid.second);
						}
					}
				}
			} else {
				std::vector<std::pair<std::string, std::string>> midsPairs;
				for (const auto & forwMid : totalCounts) {
					for (const auto & revMid : forwMid.second) {
						midsPairs.emplace_back(std::make_pair(forwMid.first, revMid.first));
					}
				}
				njh::sort(midsPairs,[&totalCounts](const std::pair<std::string, std::string> & p1,
						const std::pair<std::string, std::string> & p2){
					return totalCounts[p1.first][p1.second] > totalCounts[p2.first][p2.second];
				});
				std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> alreadyAdded;
				std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> precedingBasesCountsAgain;
				std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> precedingBasesCountsCompAgain;


				for(const auto & p : midsPairs){
					bool add = true;
					for( auto & alreadyP1 : alreadyAdded){
						for( auto & alreadyP2 : alreadyP1.second){
							if(p.first.size() == alreadyP1.first.size() &&
								 p.second.size() == alreadyP2.first.size() &&
									totalCounts[p.first][p.second]/static_cast<double>(alreadyP2.second) < 0.1){
								if(totalCounts[p.first][p.second]/static_cast<double>(alreadyP2.second) < 0.001){
									if(numberOfMismatches(p.first,  alreadyP1.first) <= 3 &&
										 numberOfMismatches(p.second, alreadyP2.first) <= 3){
										add = false;
										alreadyP2.second += totalCounts[p.first][p.second];
										precedingBasesCountsAgain[alreadyP1.first][alreadyP2.first] += precedingBasesCounts[primerPairName][p.first][p.second];
										precedingBasesCountsCompAgain[alreadyP1.first][alreadyP2.first] += precedingBasesCountsComp[primerPairName][p.first][p.second];
									}
								}else if(numberOfMismatches(p.first,  alreadyP1.first) <= 1 &&
									 numberOfMismatches(p.second, alreadyP2.first) <= 1){
									add = false;
									alreadyP2.second += totalCounts[p.first][p.second];
									precedingBasesCountsAgain[alreadyP1.first][alreadyP2.first] += precedingBasesCounts[primerPairName][p.first][p.second];
									precedingBasesCountsCompAgain[alreadyP1.first][alreadyP2.first] += precedingBasesCountsComp[primerPairName][p.first][p.second];
								}
							}
						}
					}
					if(add){
						alreadyAdded[p.first][p.second] = totalCounts[p.first][p.second];
						precedingBasesCountsAgain[p.first][p.second] = precedingBasesCounts[primerPairName][p.first][p.second];
						precedingBasesCountsCompAgain[p.first][p.second]  = precedingBasesCountsComp[primerPairName][p.first][p.second];
					}
				}

				for(const auto & forwMid : alreadyAdded){
					for(const auto & revMid : forwMid.second){
						if(revMid.second > precdingBaseFreqCutOff){
							possibleMidCounts.addRow(primerPairName,
									forwMid.first, revMid.first,
									precedingBasesCountsAgain[forwMid.first][revMid.first],
									precedingBasesCountsCompAgain[forwMid.first][revMid.first],
									revMid.second);
						}
					}
				}
			}
		}
		auto possibleMidCountsTabOpts = TableIOOpts::genTabFileOut(njh::files::make_path(setUp.pars_.directoryName_, "possibleMidCounts.tab.txt"));
		possibleMidCounts.sortTable("PrimerPair", "ForwardMID", "ReverseMID", false);
		possibleMidCounts.outPutContents(possibleMidCountsTabOpts);

		auto possibleMidCountsMostCommonTabOpts = TableIOOpts::genTabFileOut(njh::files::make_path(setUp.pars_.directoryName_, "possibleMidCountsMostCommonSize.tab.txt"));
		table possibleMidCountsMostCommonTab(VecStr{"PrimerPair", "ForwardMID", "ReverseMID", "InForDirCount", "InRevDirCount", "TotalCount"});
		for(const auto & row : possibleMidCounts){
			if(mostCommonForwardPreBases[row[possibleMidCounts.getColPos("PrimerPair")]] == row[possibleMidCounts.getColPos("ForwardMID")].size() &&
				 mostCommonReversePreBases[row[possibleMidCounts.getColPos("PrimerPair")]] == row[possibleMidCounts.getColPos("ReverseMID")].size()){
				possibleMidCountsMostCommonTab.addRow(row);
			}
		}
		possibleMidCountsMostCommonTab.sortTable("PrimerPair", "ForwardMID", "ReverseMID", false);
		possibleMidCountsMostCommonTab.outPutContents(possibleMidCountsMostCommonTabOpts);

		//unrecognized counts
		table unrecoginzedCountsTab(VecStr{"forward", "reverse", "count", "fraction"});
		for(const auto & forw : unrecognizedCounts){
			for(const auto & rev : forw.second){
				unrecoginzedCountsTab.addRow(forw.first, rev.first, rev.second, rev.second/static_cast<double>(primerPairCountsTot["unrecognized"]["unrecognized"]));
			}
		}
		auto unrecoginzedCountsTabOpts = TableIOOpts::genTabFileOut(njh::files::make_path(setUp.pars_.directoryName_, "unrecoginzedCounts.tab.txt"));
		unrecoginzedCountsTab.sortTable("count", false);
		unrecoginzedCountsTab.outPutContents(unrecoginzedCountsTabOpts);
	} else {
		seqInfo seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		while(reader.readNextRead(seq)){

		}
	}


	return 0;

}




}  // namespace njhseq

