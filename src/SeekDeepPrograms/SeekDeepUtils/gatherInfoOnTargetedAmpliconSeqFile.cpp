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




class TarAmpSeqInvestigator{
public:

	struct TarAmpSeqInvestigatorPars{
		uint32_t unrecogBaseSampling = 20;
		uint32_t precdingBaseFreqCutOff = 5;
		bfs::path idFnp = "";
		bool dontCollapsePossibleMIDs = false;
		ExtractorPairedEndPars pars;
		uint32_t testNumber = std::numeric_limits<uint32_t>::max();

		gapScoringParameters gapInfo_;

		TarAmpSeqInvestigatorPars(){
			gapInfo_.gapOpen_ = 5;
			gapInfo_.gapExtend_ = 1;
			gapInfo_.gapRightQueryOpen_ = 0;
			gapInfo_.gapRightQueryExtend_ = 0;
			gapInfo_.gapRightRefOpen_ = 0;
			gapInfo_.gapRightRefExtend_ = 0;
			gapInfo_.gapLeftQueryOpen_ = 0;
			gapInfo_.gapLeftQueryExtend_ = 0;
			gapInfo_.gapLeftRefOpen_ = 0;
			gapInfo_.gapLeftRefExtend_ = 0;
		}

	};

	TarAmpSeqInvestigator(const TarAmpSeqInvestigatorPars & pars):pars_(pars), ids_(pars.idFnp){
		ids_.initPrimerDeterminator();
	}

	TarAmpSeqInvestigatorPars pars_;
	PrimersAndMids ids_;

	//key = forward primer name, reverse primer name, count
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> primerPairCountsTot_;
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> primerPairCountsFor_;
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> primerPairCountsComp_;

	//key = primer name, forward primer bases, reverse primer bases, counts
	std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>>> precedingBasesCounts_;
	std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>>> precedingBasesCountsComp_;
	//key1 = forward seq, key2 = reverse seq, value = count
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> unrecognizedCounts_;

	uint32_t totalReadCount_{0};

	//primer pairing counts
	table primerCountsTab_{VecStr{"ForwardPrimer", "ReversePrimer", "ForwardCount", "ReverseCount", "ForwardFraction", "Total", "Fraction"}};
	//preceding number of bases counts
	table precedingBasesCountsTab_{VecStr{"PrimerPair", "BasesPrecedingPrimerCase", "NumOfBases", "Count"}};
	//preceding bases counts, to determine possible mids
	table possibleMidCounts_{VecStr{"PrimerPair", "ForwardMID", "ReverseMID", "InForDirCount", "InRevDirCount", "TotalCount"}};
	table possibleMidCountsMostCommonTab_{VecStr{"PrimerPair", "ForwardMID", "ReverseMID", "InForDirCount", "InRevDirCount", "TotalCount"}};
	//unrecognized counts
	table unrecoginzedCountsTab_{VecStr{"forward", "reverse", "count", "fraction"}};

	void investigateSeq(const seqInfo & forwardSeq, const seqInfo & revCompSeq, aligner & alignObj){
		++totalReadCount_;
		std::string forwardPrimerName = "";
		std::string forwardPrimerPrecedingBases = "";

		std::string reversePrimerName = "";
		std::string reversePrimerPrcedingBases = "";

		bool complement = false;
		auto forPos = ids_.pDeterminator_->determineBestForwardPrimerPosFront(forwardSeq, pars_.pars.corePars_.pDetPars, alignObj);
		if("" == forPos.primerName_){
			auto forMate = ids_.pDeterminator_->determineBestForwardPrimerPosFront(revCompSeq, pars_.pars.corePars_.pDetPars, alignObj);
			if("" != forMate.primerName_){
				forwardPrimerName = forMate.primerName_;
				if(0 != forMate.start_ ){
					forwardPrimerPrecedingBases = revCompSeq.seq_.substr(0, forMate.start_);
				}
				complement = true;
			}else{
				forwardPrimerName = "unrecognized";
			}
		}else{
			forwardPrimerName = forPos.primerName_;
			if(0 != forPos.start_ ){
				forwardPrimerPrecedingBases = forwardSeq.seq_.substr(0, forPos.start_);
			}
		}

		if("unrecognized" == forwardPrimerName){
			auto revPos = ids_.pDeterminator_->determineBestReversePrimerPosFront(revCompSeq, pars_.pars.corePars_.pDetPars, alignObj);
			if("" == revPos.primerName_){
				auto revPosOther = ids_.pDeterminator_->determineBestReversePrimerPosFront(forwardSeq, pars_.pars.corePars_.pDetPars, alignObj);
				if("" != revPosOther.primerName_){
					//matched
					complement = true;
					reversePrimerName = revPosOther.primerName_;
					if(0 != revPos.start_ ){
						reversePrimerPrcedingBases = forwardSeq.seq_.substr(0, revPos.start_);
					}
				}else{
					//unmatched
					reversePrimerName = "unrecognized";
				}
			}else{
				//matched
				reversePrimerName = revPos.primerName_;
				if(0 != revPos.start_ ){
					reversePrimerPrcedingBases = revCompSeq.seq_.substr(0, revPos.start_);
				}
			}
		}else{
			if(complement){
				auto revPos = ids_.pDeterminator_->determineBestReversePrimerPosFront(forwardSeq, pars_.pars.corePars_.pDetPars, alignObj);
				if("" != revPos.primerName_){
					//matched
					reversePrimerName = revPos.primerName_;
					if(0 != revPos.start_ ){
						reversePrimerPrcedingBases = forwardSeq.seq_.substr(0, revPos.start_);
					}
				}else{
					reversePrimerName = "unrecognized";
				}
			}else{
				auto revPos = ids_.pDeterminator_->determineBestReversePrimerPosFront(revCompSeq, pars_.pars.corePars_.pDetPars, alignObj);
				if("" != revPos.primerName_){
					//matched
					reversePrimerName = revPos.primerName_;
					if(0 != revPos.start_ ){
						reversePrimerPrcedingBases = revCompSeq.seq_.substr(0, revPos.start_);
					}
				}else{
					reversePrimerName = "unrecognized";
				}
			}
		}
		if(complement){
			primerPairCountsComp_[forwardPrimerName][reversePrimerName]+= forwardSeq.cnt_;
		} else {
			primerPairCountsFor_[forwardPrimerName][reversePrimerName]+= forwardSeq.cnt_;
		}
		primerPairCountsTot_[forwardPrimerName][reversePrimerName]+= forwardSeq.cnt_;
		if(forwardPrimerName == reversePrimerName && "unrecognized" != forwardPrimerName){

		}
		if(forwardPrimerName == reversePrimerName && "unrecognized" != forwardPrimerName){
			if(complement){
				precedingBasesCountsComp_[forwardPrimerName][forwardPrimerPrecedingBases][reversePrimerPrcedingBases] += forwardSeq.cnt_;
			}else{
				precedingBasesCounts_[forwardPrimerName][forwardPrimerPrecedingBases][reversePrimerPrcedingBases] += forwardSeq.cnt_;
			}
		}else if(forwardPrimerName == reversePrimerName && "unrecognized" == forwardPrimerName){
			unrecognizedCounts_[forwardSeq.seq_.substr(0, pars_.unrecogBaseSampling)][revCompSeq.seq_.substr(0, pars_.unrecogBaseSampling)] += forwardSeq.cnt_;
		}
	}

	void processCounts(){




		for(const auto & forPrim : primerPairCountsTot_){
			for(const auto & revPrim : forPrim.second){
				primerCountsTab_.addRow(forPrim.first, revPrim.first,
						primerPairCountsFor_[forPrim.first][revPrim.first],
						primerPairCountsComp_[forPrim.first][revPrim.first],
						primerPairCountsFor_[forPrim.first][revPrim.first]/static_cast<double>(revPrim.second),
						revPrim.second,
						revPrim.second/static_cast<double>(totalReadCount_));
			}
		}



		primerCountsTab_.sortTable("ForwardPrimer", "ReversePrimer", false);

		std::unordered_map<std::string, uint32_t> mostCommonForwardPreBases;
		std::unordered_map<std::string, uint32_t> mostCommonReversePreBases;

		for(const auto & primerPair : ids_.targets_){
			std::unordered_map<std::string,std::unordered_map<uint32_t, uint32_t>> counts;
			std::unordered_map<uint32_t, uint32_t> forwardTotalCounts;
			std::unordered_map<uint32_t, uint32_t> reverseTotalCounts;

			for(const auto & forward : precedingBasesCounts_[primerPair.first]){
				for(const auto & reverse : forward.second){
					forwardTotalCounts[forward.first.size()] += reverse.second;
					reverseTotalCounts[reverse.first.size()] += reverse.second;
					counts["ForwardPrimerPrecedingForward"][forward.first.size()] += reverse.second;
					counts["ReversePrimerPrecedingForward"][reverse.first.size()] += reverse.second;
				}
			}
			for(const auto & forward : precedingBasesCountsComp_[primerPair.first]){
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
					precedingBasesCountsTab_.addRow(primerPair.first, count.first, numBases.first, numBases.second);
				}
			}
		}
		precedingBasesCountsTab_.sortTable("PrimerPair", "BasesPrecedingPrimerCase", "NumOfBases", false);


		for(const auto & primerPair : ids_.targets_){
			std::string primerPairName = primerPair.first;
			std::unordered_map<std::string, std::unordered_map<std::string, double>> totalCounts;
			for(const auto & forwMid : precedingBasesCounts_[primerPairName]){
				for(const auto & revMid : forwMid.second){
					totalCounts[forwMid.first][revMid.first] += revMid.second;
				}
			}
			for(const auto & forwMid : precedingBasesCountsComp_[primerPairName]){
				for(const auto & revMid : forwMid.second){
					totalCounts[forwMid.first][revMid.first] += revMid.second;
				}
			}
			if(pars_.dontCollapsePossibleMIDs){
				for(const auto & forwMid : totalCounts){
					for(const auto & revMid : forwMid.second){
						if(revMid.second > pars_.precdingBaseFreqCutOff){
							possibleMidCounts_.addRow(primerPairName,
									forwMid.first, revMid.first,
									precedingBasesCounts_[primerPairName][forwMid.first][revMid.first],
									precedingBasesCountsComp_[primerPairName][forwMid.first][revMid.first],
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
										precedingBasesCountsAgain[alreadyP1.first][alreadyP2.first] += precedingBasesCounts_[primerPairName][p.first][p.second];
										precedingBasesCountsCompAgain[alreadyP1.first][alreadyP2.first] += precedingBasesCountsComp_[primerPairName][p.first][p.second];
									}
								}else if(numberOfMismatches(p.first,  alreadyP1.first) <= 1 &&
									 numberOfMismatches(p.second, alreadyP2.first) <= 1){
									add = false;
									alreadyP2.second += totalCounts[p.first][p.second];
									precedingBasesCountsAgain[alreadyP1.first][alreadyP2.first] += precedingBasesCounts_[primerPairName][p.first][p.second];
									precedingBasesCountsCompAgain[alreadyP1.first][alreadyP2.first] += precedingBasesCountsComp_[primerPairName][p.first][p.second];
								}
							}
						}
					}
					if(add){
						alreadyAdded[p.first][p.second] = totalCounts[p.first][p.second];
						precedingBasesCountsAgain[p.first][p.second] = precedingBasesCounts_[primerPairName][p.first][p.second];
						precedingBasesCountsCompAgain[p.first][p.second]  = precedingBasesCountsComp_[primerPairName][p.first][p.second];
					}
				}

				for(const auto & forwMid : alreadyAdded){
					for(const auto & revMid : forwMid.second){
						if(revMid.second > pars_.precdingBaseFreqCutOff){
							possibleMidCounts_.addRow(primerPairName,
									forwMid.first, revMid.first,
									precedingBasesCountsAgain[forwMid.first][revMid.first],
									precedingBasesCountsCompAgain[forwMid.first][revMid.first],
									revMid.second);
						}
					}
				}
			}
		}
		possibleMidCounts_.sortTable("PrimerPair", "ForwardMID", "ReverseMID", false);




		for(const auto & row : possibleMidCounts_){
			if(mostCommonForwardPreBases[row[possibleMidCounts_.getColPos("PrimerPair")]] == row[possibleMidCounts_.getColPos("ForwardMID")].size() &&
				 mostCommonReversePreBases[row[possibleMidCounts_.getColPos("PrimerPair")]] == row[possibleMidCounts_.getColPos("ReverseMID")].size()){
				possibleMidCountsMostCommonTab_.addRow(row);
			}
		}
		possibleMidCountsMostCommonTab_.sortTable("PrimerPair", "ForwardMID", "ReverseMID", false);


		for(const auto & forw : unrecognizedCounts_){
			for(const auto & rev : forw.second){
				unrecoginzedCountsTab_.addRow(forw.first, rev.first, rev.second, rev.second/static_cast<double>(primerPairCountsTot_["unrecognized"]["unrecognized"]));
			}
		}
		unrecoginzedCountsTab_.sortTable("count", false);


	}

	void writeOutTables(const bfs::path & directory, bool overWrite){
		auto primerTabOpts = TableIOOpts::genTabFileOut(njh::files::make_path(directory, "primerCounts.tab.txt"));
		primerTabOpts.out_.overWriteFile_ = overWrite;
		primerCountsTab_.outPutContents(primerTabOpts);

		auto precedingBasesCountsTabOpts = TableIOOpts::genTabFileOut(njh::files::make_path(directory, "precedingBasesCounts.tab.txt"));
		precedingBasesCountsTabOpts.out_.overWriteFile_ = overWrite;
		precedingBasesCountsTab_.outPutContents(precedingBasesCountsTabOpts);

		auto possibleMidCountsTabOpts = TableIOOpts::genTabFileOut(njh::files::make_path(directory, "possibleMidCounts.tab.txt"));
		possibleMidCountsTabOpts.out_.overWriteFile_ = overWrite;
		possibleMidCounts_.outPutContents(possibleMidCountsTabOpts);

		auto possibleMidCountsMostCommonTabOpts = TableIOOpts::genTabFileOut(njh::files::make_path(directory, "possibleMidCountsMostCommonSize.tab.txt"));
		possibleMidCountsMostCommonTabOpts.out_.overWriteFile_ = overWrite;
		possibleMidCountsMostCommonTab_.outPutContents(possibleMidCountsMostCommonTabOpts);

		auto unrecoginzedCountsTabOpts = TableIOOpts::genTabFileOut(njh::files::make_path(directory, "unrecoginzedCounts.tab.txt"));
		unrecoginzedCountsTabOpts.out_.overWriteFile_ = overWrite;
		unrecoginzedCountsTab_.outPutContents(unrecoginzedCountsTabOpts);
	}

};


int SeekDeepUtilsRunner::gatherInfoOnTargetedAmpliconSeqFile(
		const njh::progutils::CmdArgs & inputCommands) {

	TarAmpSeqInvestigator::TarAmpSeqInvestigatorPars investPars;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(investPars.testNumber, "--testNumber", "Just use this number of reads of the top of the file");
	setUp.setOption(investPars.dontCollapsePossibleMIDs, "--dontCollapsePossibleMIDs",
			"Don't Collapse Possible MIDs", false);
	setUp.setOption(investPars.unrecogBaseSampling, "--unrecogBaseSampling",
			"Number of bases to sample from file for unrecognized sequences", false);

	setUp.setOption(investPars.precdingBaseFreqCutOff, "--precdingBaseFreqCutOff", "Preceding Base Freq Cut Off", false);

	investPars.pars.corePars_.pDetPars.primerWithin_ = 30;
	setUp.setOption(investPars.pars.corePars_.pDetPars.primerWithin_, "--primerWithin", "Primer Within bases search", false, "Primer");


	bool primerToUpperCase = false;
	setUp.setOption(primerToUpperCase, "--primerUpper",
			"Leave primers in upper case", false, "Primer");
	investPars.pars.corePars_.pDetPars.primerToLowerCase_ = !primerToUpperCase;
	setUp.setOption(investPars.pars.corePars_.pDetPars.allowable_.distances_.query_.coverage_, "--primerCoverage",
			"Amount of primers found", false, "Primer");
	setUp.setOption(investPars.pars.corePars_.pDetPars.allowable_.hqMismatches_, "--primerNumOfMismatches",
			"Number of Mismatches to allow in primers", false, "Primer");
	setUp.setOption(investPars.pars.corePars_.pDetPars.allowable_.oneBaseIndel_, "--primerOneBaseIndels",
			"Number Of One base indels to allow in primers", false, "Primer");
	setUp.setOption(investPars.pars.corePars_.pDetPars.allowable_.twoBaseIndel_, "--primerTwoBaseIndels",
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
	investPars.gapInfo_ = setUp.pars_.gapInfo_;
	setUp.setOption(investPars.idFnp, "--id", "SeekDeep primers file", true);
	setUp.processReadInNames(VecStr{"--fastq1", "--fastq", "--fasta", "--fastq1gz", "--fastqgz", "--fastagz"});
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);



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
			if(readCount >=investPars.testNumber){
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
			if(readCount >=investPars.testNumber){
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
//	KmerMaps emptyMaps;
//	bool countEndGaps = false;
	//to avoid allocating an extremely large aligner matrix;
	aligner alignObj(maxReadSize, investPars.gapInfo_, scoreMatrix);


	TarAmpSeqInvestigator investigator(investPars);

	if (setUp.pars_.ioOptions_.isPairedIn()) {
		PairedRead seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		njh::ProgressBar pBar(readCount);
		uint32_t newReadCount = 0;
		while(reader.readNextRead(seq)){
			++newReadCount;
			if(setUp.pars_.verbose_){
				pBar.outputProgAdd(std::cout, 1, true);
			}
			investigator.investigateSeq(seq.seqBase_, seq.mateSeqBase_, alignObj);
			if(newReadCount >=investPars.testNumber){
				break;
			}
		}
	} else {
		seqInfo seq;
		SeqInput reader(setUp.pars_.ioOptions_);
		reader.openIn();
		njh::ProgressBar pBar(readCount);
		uint32_t newReadCount = 0;
		while(reader.readNextRead(seq)){
			auto revCompSeq = seq;
			revCompSeq.reverseComplementRead(false, true);
			++newReadCount;
			if(setUp.pars_.verbose_){
				pBar.outputProgAdd(std::cout, 1, true);
			}
			investigator.investigateSeq(seq, revCompSeq, alignObj);
			if(newReadCount >=investPars.testNumber){
				break;
			}
		}
	}

	investigator.processCounts();
	investigator.writeOutTables(setUp.pars_.directoryName_, true);

	return 0;

}




}  // namespace njhseq

