/*
 * TarAmpSeqInvestigator.cpp
 *
 *  Created on: Jan 31, 2020
 *      Author: nicholashathaway
 */


#include "TarAmpSeqInvestigator.hpp"

namespace njhseq {


TarAmpSeqInvestigator::TarAmpSeqInvestigatorPars::TarAmpSeqInvestigatorPars(){
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

	midPars.searchStop_ = 40;
	pars.corePars_.pDetPars.primerWithin_ = 40;
}






TarAmpSeqInvestigator::TarAmpSeqInvestigator(const TarAmpSeqInvestigatorPars & pars) :
		pars_(pars), ids_(pars.idFnp) {
	ids_.initPrimerDeterminator();
	if(ids_.containsMids()){
		ids_.initMidDeterminator(pars.midPars);
	}
}



void TarAmpSeqInvestigator::investigateSeq(const seqInfo & forwardSeq, const seqInfo & revCompSeq, aligner & alignObj){
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

void TarAmpSeqInvestigator::addOtherCounts(const TarAmpSeqInvestigator & other){

	totalReadCount_ += other.totalReadCount_;
	for(const auto & fwd : other.primerPairCountsTot_){
		for(const auto & rev : fwd.second){
			primerPairCountsTot_[fwd.first][rev.first] += rev.second;
		}
	}
	for(const auto & fwd : other.primerPairCountsFor_){
		for(const auto & rev : fwd.second){
			primerPairCountsFor_[fwd.first][rev.first] += rev.second;
		}
	}
	for(const auto & fwd : other.primerPairCountsComp_){
		for(const auto & rev : fwd.second){
			primerPairCountsComp_[fwd.first][rev.first] += rev.second;
		}
	}

	for(const auto & pair : other.precedingBasesCounts_){
		for(const auto & fwd : pair.second){
			for(const auto & rev : fwd.second){
				precedingBasesCounts_[pair.first][fwd.first][rev.first] += rev.second;
			}
		}
	}
	for(const auto & pair : other.precedingBasesCountsComp_){
		for(const auto & fwd : pair.second){
			for(const auto & rev : fwd.second){
				precedingBasesCountsComp_[pair.first][fwd.first][rev.first] += rev.second;
			}
		}
	}

	for(const auto & fwd : other.unrecognizedCounts_){
		for(const auto & rev : fwd.second){
			unrecognizedCounts_[fwd.first][rev.first] += rev.second;
		}
	}
}

void TarAmpSeqInvestigator::processCounts(){




	for(const auto & forPrim : primerPairCountsTot_){
		for(const auto & revPrim : forPrim.second){
			primerCountsTab_.addRow(forPrim.first, revPrim.first,
					primerPairCountsFor_[forPrim.first][revPrim.first],
					primerPairCountsComp_[forPrim.first][revPrim.first],
					primerPairCountsComp_[forPrim.first][revPrim.first]/static_cast<double>(revPrim.second),
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

void TarAmpSeqInvestigator::writeOutTables(const bfs::path & directory, bool overWrite){
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


	if(ids_.containsMids()){
		auto precedingBasesCountsAboveMidTabOpts = TableIOOpts::genTabFileOut(njh::files::make_path(directory, "precedingBasesCountsAboveMIDSize.tab.txt"));
		precedingBasesCountsAboveMidTabOpts.out_.overWriteFile_ = overWrite;

		table precedingBasesCountsAboveMIDTab{VecStr{"PrimerPair", "BasesPrecedingPrimerCase", "AboveMID", "Fraction", "Total"}};
		std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> countsAbove;
		std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> countsBelow;

		auto midSize = ids_.getMaxMIDSize();
		for(const auto & row : precedingBasesCountsTab_){
			auto target = row[precedingBasesCountsTab_.getColPos("PrimerPair")];
			auto orientation = row[precedingBasesCountsTab_.getColPos("BasesPrecedingPrimerCase")];
			auto count = njh::StrToNumConverter::stoToNum<uint32_t>(row[precedingBasesCountsTab_.getColPos("Count")]);
			auto NumOfBases = njh::StrToNumConverter::stoToNum<uint32_t>(row[precedingBasesCountsTab_.getColPos("NumOfBases")]);
			if(NumOfBases <= midSize){
				countsBelow[target][orientation] += count;
			}else{
				countsAbove[target][orientation] += count;
			}
		}
		for(const auto & tar : countsAbove){
			for(const auto & orientation : tar.second){
				double total = countsAbove[tar.first][orientation.first] + countsBelow[tar.first][orientation.first];
				auto above = countsAbove[tar.first][orientation.first] ;
				auto frac = above/total;
				precedingBasesCountsAboveMIDTab.addRow(tar.first, orientation.first, above, frac, total);
			}
		}
		precedingBasesCountsAboveMIDTab.sortTable("PrimerPair", "BasesPrecedingPrimerCase", false);
		precedingBasesCountsAboveMIDTab.outPutContents(precedingBasesCountsAboveMidTabOpts);
	}
}


TarAmpSeqInvestigator::prepareForInvestiagteFileRes TarAmpSeqInvestigator::prepareForInvestiagteFile(const SeqIOOptions & opts, bool verbose){
	prepareForInvestiagteFileRes ret;

	if (opts.isPairedIn()) {
		PairedRead seq;
		SeqInput reader(opts);
		reader.openIn();
		while(reader.readNextRead(seq)){
			++ret.readCount;
			readVec::getMaxLength(seq.seqBase_, ret.maxReadSize);
			readVec::getMaxLength(seq.mateSeqBase_, ret.maxReadSize);
			if(verbose && ret.readCount % 10000 == 0){
				std::cout << "\r" << ret.readCount;
				std::cout.flush();
			}
			if(ret.readCount >=pars_.testNumber){
				break;
			}
		}
		if(verbose){
			std::cout << std::endl;
		}
	} else {
		seqInfo seq;
		SeqInput reader(opts);
		reader.openIn();
		while(reader.readNextRead(seq)){
			++ret.readCount;
			readVec::getMaxLength(seq, ret.maxReadSize);
			if(verbose && ret.readCount % 10000 == 0){
				std::cout << "\r" << ret.readCount;
				std::cout.flush();
			}
			if(ret.readCount >=pars_.testNumber){
				break;
			}
		}
		if(verbose){
			std::cout << std::endl;
		}
	}
	return ret;
}


void TarAmpSeqInvestigator::investigateFile(const SeqIOOptions & opts, const TarAmpSeqInvestigator::prepareForInvestiagteFileRes & counts, bool verbose){
	// create aligner for primer identification
	auto scoreMatrix = substituteMatrix::createDegenScoreMatrixNoNInRef(2, -2);
	//to avoid allocating an extremely large aligner matrix;
	aligner alignObj(counts.maxReadSize, pars_.gapInfo_, scoreMatrix);

	if (opts.isPairedIn()) {
		PairedRead seq;
		SeqInput reader(opts);
		reader.openIn();
		njh::ProgressBar pBar(counts.readCount);
		uint32_t newReadCount = 0;
		while(reader.readNextRead(seq)){
			++newReadCount;
			if(verbose){
				pBar.outputProgAdd(std::cout, 1, true);
			}
			investigateSeq(seq.seqBase_, seq.mateSeqBase_, alignObj);
			if(newReadCount >=pars_.testNumber){
				break;
			}
		}
	} else {
		seqInfo seq;
		SeqInput reader(opts);
		reader.openIn();
		njh::ProgressBar pBar(counts.readCount);
		uint32_t newReadCount = 0;
		while(reader.readNextRead(seq)){
			auto revCompSeq = seq;
			revCompSeq.reverseComplementRead(false, true);
			++newReadCount;
			if(verbose){
				pBar.outputProgAdd(std::cout, 1, true);
			}
			investigateSeq(seq, revCompSeq, alignObj);
			if(newReadCount >=pars_.testNumber){
				break;
			}
		}
	}
}
void TarAmpSeqInvestigator::investigateFile(const SeqIOOptions & opts, bool verbose){
	auto counts = prepareForInvestiagteFile(opts, verbose);
	investigateFile(opts, counts, verbose);
}


bool TarAmpSeqInvestigator::reverseComplementLikely(uint32_t minReadAmount, double cutOff) const {
	// bool ret = false;
	double totalReads = 0;
	double totalReverseCount = 0;
	for(const auto & row : primerCountsTab_){
		auto forwardPrimer = row[primerCountsTab_.getColPos("ForwardPrimer")];
		auto reversePrimer = row[primerCountsTab_.getColPos("ReversePrimer")];
		auto ReverseCount = njh::StrToNumConverter::stoToNum<double>(row[primerCountsTab_.getColPos("ReverseCount")]);
		auto reverseFraction = njh::StrToNumConverter::stoToNum<double>(row[primerCountsTab_.getColPos("ReverseFraction")]);
		auto Total = njh::StrToNumConverter::stoToNum<double>(row[primerCountsTab_.getColPos("Total")]);
		if("unrecognized" != forwardPrimer && "unrecognized" != reversePrimer && forwardPrimer == reversePrimer){
			totalReads+= Total;
			totalReverseCount += ReverseCount;
			if(Total >=minReadAmount && reverseFraction >=cutOff){
				//forward fraction should be less than 80%
				// ret = true;
				break;
			}
		}
	}
//	std::cout << "totalReverseCount/totalReads: " << totalReverseCount/totalReads << std::endl;
//	std::cout << "totalReads: " << totalReads << std::endl;
//	std::cout << "totalReverseCount: " << totalReverseCount << std::endl;

	return totalReads >=minReadAmount && (totalReverseCount/totalReads) >=cutOff;
}

bool TarAmpSeqInvestigator::hasPossibleRandomPrecedingBases(uint32_t midSize, uint32_t minReadAmount, double cutOff) const{

	//key1 = target, key2 = primerOrientation, , value = count
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> countsAbove;
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> countsBelow;

	for(const auto & row : precedingBasesCountsTab_){
		auto target = row[precedingBasesCountsTab_.getColPos("PrimerPair")];
		auto orientation = row[precedingBasesCountsTab_.getColPos("BasesPrecedingPrimerCase")];
		auto count = njh::StrToNumConverter::stoToNum<uint32_t>(row[precedingBasesCountsTab_.getColPos("Count")]);
		auto NumOfBases = njh::StrToNumConverter::stoToNum<uint32_t>(row[precedingBasesCountsTab_.getColPos("NumOfBases")]);
		if(NumOfBases <= midSize){
			countsBelow[target][orientation] += count;
		}else{
			countsAbove[target][orientation] += count;
		}
	}
	for(const auto & tar : countsAbove){
		for(const auto & orientation : tar.second){
			if(orientation.second > minReadAmount){
				double total = orientation.second + countsBelow[tar.first][orientation.first];
				if(orientation.second/total > cutOff){
					return true;
				}
			}
		}
	}
	return false;
}


bool TarAmpSeqInvestigator::hasPossibleRandomPrecedingBasesForwardPrimer(uint32_t midSize, uint32_t minReadAmount, double cutOff) const{

	//key1 = target, key2 = primerOrientation, , value = count
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> countsAbove;
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> countsBelow;

	for(const auto & row : precedingBasesCountsTab_){
		auto target = row[precedingBasesCountsTab_.getColPos("PrimerPair")];
		auto orientation = row[precedingBasesCountsTab_.getColPos("BasesPrecedingPrimerCase")];
		if(njh::beginsWith(orientation, "ForwardPrimer")){
			auto count = njh::StrToNumConverter::stoToNum<uint32_t>(row[precedingBasesCountsTab_.getColPos("Count")]);
			auto NumOfBases = njh::StrToNumConverter::stoToNum<uint32_t>(row[precedingBasesCountsTab_.getColPos("NumOfBases")]);
			if(NumOfBases <= midSize){
				countsBelow[target][orientation] += count;
			}else{
				countsAbove[target][orientation] += count;
			}
		}
	}
	for(const auto & tar : countsAbove){
		for(const auto & orientation : tar.second){
			if(orientation.second > minReadAmount){
				double total = orientation.second + countsBelow[tar.first][orientation.first];
				if(orientation.second/total > cutOff){
					return true;
				}
			}
		}
	}
	return false;
}

bool TarAmpSeqInvestigator::hasPossibleRandomPrecedingBasesReversePrimer(uint32_t midSize, uint32_t minReadAmount, double cutOff) const{

	//key1 = target, key2 = primerOrientation, , value = count
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> countsAbove;
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> countsBelow;

	for(const auto & row : precedingBasesCountsTab_){
		auto target = row[precedingBasesCountsTab_.getColPos("PrimerPair")];
		auto orientation = row[precedingBasesCountsTab_.getColPos("BasesPrecedingPrimerCase")];
		if(njh::beginsWith(orientation, "ReversePrimer")){
			auto count = njh::StrToNumConverter::stoToNum<uint32_t>(row[precedingBasesCountsTab_.getColPos("Count")]);
			auto NumOfBases = njh::StrToNumConverter::stoToNum<uint32_t>(row[precedingBasesCountsTab_.getColPos("NumOfBases")]);
			if(NumOfBases <= midSize){
				countsBelow[target][orientation] += count;
			}else{
				countsAbove[target][orientation] += count;
			}
		}
	}

	for(const auto & tar : countsAbove){
		for(const auto & orientation : tar.second){
			if(orientation.second > minReadAmount){
				double total = orientation.second + countsBelow[tar.first][orientation.first];
				std::cout << "orientation.first: " << orientation.first << std::endl;
				std::cout << "orientation.second/total: " << orientation.second/total << std::endl;
				std::cout << "orientation.second: " << orientation.second << std::endl;
				std::cout << "total: " << total << std::endl;

				if(orientation.second/total > cutOff){
					return true;
				}
			}
		}
	}
	return false;
}

uint32_t TarAmpSeqInvestigator::maxPrecedingBases() const{
	uint32_t ret = 0;
	for(const auto & row : precedingBasesCountsTab_){
		auto NumOfBases = njh::StrToNumConverter::stoToNum<uint32_t>(row[precedingBasesCountsTab_.getColPos("NumOfBases")]);
		// std::cout << "NumOfBases: " << NumOfBases << std::endl;
		if(NumOfBases > ret){
			ret = NumOfBases;
		}
	}
	return ret;
}

uint32_t TarAmpSeqInvestigator::maxPrecedingReversePrimerBases() const{
	uint32_t ret = 0;
	for(const auto & row : precedingBasesCountsTab_){
		auto BasesPrecedingPrimerCase  = row[precedingBasesCountsTab_.getColPos("BasesPrecedingPrimerCase")];
		if(njh::beginsWith(BasesPrecedingPrimerCase, "ReversePrimer")){
			auto NumOfBases = njh::StrToNumConverter::stoToNum<uint32_t>(row[precedingBasesCountsTab_.getColPos("NumOfBases")]);
			if(NumOfBases > ret){
				ret = NumOfBases;
			}
		}
	}
	return ret;
}

uint32_t TarAmpSeqInvestigator::maxPrecedingForwardPrimerBases() const{
	uint32_t ret = 0;
	for(const auto & row : precedingBasesCountsTab_){
		auto BasesPrecedingPrimerCase  = row[precedingBasesCountsTab_.getColPos("BasesPrecedingPrimerCase")];
		if(njh::beginsWith(BasesPrecedingPrimerCase, "ForwardPrimer")){
			auto NumOfBases = njh::StrToNumConverter::stoToNum<uint32_t>(row[precedingBasesCountsTab_.getColPos("NumOfBases")]);
			if(NumOfBases > ret){
				ret = NumOfBases;
			}
		}
	}
	return ret;
}




VecStr TarAmpSeqInvestigator::recommendSeekDeepExtractorFlags() const {

	auto possibleRevComp = reverseComplementLikely();
	auto possiblePrecedingRandomeBases = hasPossibleRandomPrecedingBases(ids_.getMaxMIDSize());
	auto maxPre = maxPrecedingBases();
	VecStr recFlags{};
	if(ids_.containsMids()){
		// bool containsDualBarcode = false;
		bool containsSingleBarcode = false;
		for(const auto & id : ids_.mids_){
			// if(nullptr != id.second.forwardBar_ && nullptr != id.second.reverseBar_){
			// 	containsDualBarcode = true;
			// }
			if(nullptr != id.second.reverseBar_){
				containsSingleBarcode = true;
			}
			if(nullptr != id.second.forwardBar_){
				containsSingleBarcode = true;
			}
		}
		if(possibleRevComp){
			recFlags.emplace_back("--checkRevComplementForMids");
		}
		if(possiblePrecedingRandomeBases){
			auto pre = maxPre - ids_.getMaxMIDSize();
			recFlags.emplace_back(njh::pasteAsStr("--midWithinStart ", pre));
		}
		if(containsSingleBarcode){
			auto possiblePrecedingRandomeBasesRp = hasPossibleRandomPrecedingBasesReversePrimer(ids_.getMaxMIDSize());
			if(possiblePrecedingRandomeBasesRp){
				auto maxPreRp = maxPrecedingReversePrimerBases();
				recFlags.emplace_back(njh::pasteAsStr("--primerWithinStart ", maxPreRp));
			}
		}
	} else {
		if(possibleRevComp){
			recFlags.emplace_back("--checkRevComplementForPrimers");
		}
		if(possiblePrecedingRandomeBases){
			recFlags.emplace_back(njh::pasteAsStr("--primerWithinStart ", maxPre));
		}
	}
	return recFlags;
}



}  // namespace njhseq

