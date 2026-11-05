/*
 * PrimerDimerDetector.cpp
 *
 *  Created on: Nov 4, 2026
 *      Author: nicholashathaway
 */

// SeekDeep - A program for analyzing targeted amplicon sequencing
//
// This file is part of SeekDeep.
//
// SeekDeep is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// SeekDeep is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with SeekDeep.  If not, see <http://www.gnu.org/licenses/>.
//

#include "PrimerDimerDetector.hpp"

namespace njhseq {

PrimerDimerDetector::DetectorPars::DetectorPars(){
	allowableErrors_.hqMismatches_ = 2;
	allowableErrors_.lqMismatches_ = 3;
	allowableErrors_.oneBaseIndel_ = 2;
}

std::string PrimerDimerDetector::getReadCategoryStr(ReadCategory category){
	switch (category) {
		case ReadCategory::DIMER:
			return "dimer";
		case ReadCategory::PRIMERS_FAR_APART:
			return "primersFarApart";
		case ReadCategory::UNRESOLVED_DIMER:
			return "unresolvedDimer";
		case ReadCategory::ONE_PRIMER:
			return "onePrimer";
		case ReadCategory::NO_PRIMER:
			return "noPrimer";
	}
	return "notHandled";//shouldn't be getting here
}

uint32_t PrimerDimerDetector::calcMinPrimerMatchBases(uint64_t numPrimers){
	const double chanceAllowed = 0.001;
	const auto basesForChance = static_cast<uint32_t>(std::ceil(std::log(std::max<double>(numPrimers, 1) / chanceAllowed) / std::log(4.0)));
	return std::max<uint32_t>(8, basesForChance);
}

uint32_t PrimerDimerDetector::PrimerMatch::primerBasesMatched() const{
	return primerEnd_ - primerStart_;
}

double PrimerDimerDetector::PrimerMatch::primerFractionMatched() const{
	return 0 == primerLen_ ? 0.0 : static_cast<double>(primerBasesMatched()) / primerLen_;
}

uint32_t PrimerDimerDetector::PrimerMatch::getGap(const PrimerMatch & other) const{
	if(seqEnd_ < other.seqStart_){
		return other.seqStart_ - seqEnd_;
	}
	if(other.seqEnd_ < seqStart_){
		return seqStart_ - other.seqEnd_;
	}
	return 0;
}

uint32_t PrimerDimerDetector::PrimerMatch::getOverlap(const PrimerMatch & other) const{
	const auto overlapStart = std::max(seqStart_, other.seqStart_);
	const auto overlapEnd = std::min(seqEnd_, other.seqEnd_);
	return overlapEnd > overlapStart ? overlapEnd - overlapStart : 0;
}

PrimerDimerDetector::PrimerPair::PrimerPair(PrimerMatch front, PrimerMatch end) :
		front_(std::move(front)), end_(std::move(end)){
}

uint32_t PrimerDimerDetector::PrimerPair::getGap() const{
	return front_.getGap(end_);
}

uint32_t PrimerDimerDetector::PrimerPair::getOverlap() const{
	return front_.getOverlap(end_);
}

double PrimerDimerDetector::PrimerPair::getAlnScoreSum() const{
	return front_.comp_.alnScore_ + end_.comp_.alnScore_;
}

double PrimerDimerDetector::PrimerPair::getIdentitySum() const{
	return front_.comp_.distances_.eventBasedIdentity_ + end_.comp_.distances_.eventBasedIdentity_;
}

PrimerDimerDetector::PrimerDimerDetector(const PrimerDeterminator & primers, DetectorPars pars) :
		pars_(std::move(pars)){
	maxPrimerSize_ = primers.getMaxPrimerSize();
	const auto addPrimers = [this](const std::string & name, const std::vector<PrimerDeterminator::primerInfo::PrimerSeq> & seqs){
		for(const auto & primer : seqs){
			primers_.emplace_back(NamedPrimer{name, primer.info_, primer.infoRC_});
		}
	};
	for(const auto & primer : primers.primers_){
		addPrimers(primer.first + "_F", primer.second.fwds_);
		addPrimers(primer.first + "_R", primer.second.revs_);
	}
	if(0 == pars_.minPrimerMatchBases_){
		pars_.minPrimerMatchBases_ = calcMinPrimerMatchBases(primers_.size());
	}
}

aligner PrimerDimerDetector::genAligner(uint64_t maxSeqSize) const{
	return aligner(std::max(maxSeqSize, maxPrimerSize_), gapScoringParameters::genSemiGlobalQueryOnly(5, 1), substituteMatrix::createDegenScoreMatrix(2, -4));
}

uint64_t PrimerDimerDetector::getMaxPrimerSize() const{
	return maxPrimerSize_;
}

PrimerDimerDetector::PrimerMatch PrimerDimerDetector::genMatch(const std::string & name, const seqInfo & primer, bool revComp, aligner & alignerObj){
	PrimerMatch match;
	match.primerName_ = name;
	match.comp_ = alignerObj.comp_;
	match.alignedPrimer_ = alignerObj.alignObjectB_.seqBase_;
	match.seqStart_ = alignerObj.parts_.lHolder_.localAStart_;
	match.seqEnd_ = alignerObj.parts_.lHolder_.localAStart_ + alignerObj.parts_.lHolder_.localASize_;
	match.primerStart_ = alignerObj.parts_.lHolder_.localBStart_;
	match.primerEnd_ = alignerObj.parts_.lHolder_.localBStart_ + alignerObj.parts_.lHolder_.localBSize_;
	match.primerLen_ = len(primer);
	match.revComp_ = revComp;
	return match;
}

PrimerDimerDetector::DetectorResults PrimerDimerDetector::detect(const seqInfo & seq, aligner & alignerObj) const{
	DetectorResults ret;
	const auto alignPrimer = [&seq, &alignerObj](const seqInfo & primer){
		alignerObj.alignCacheLocal(seq, primer);
		alignerObj.rearrangeObjsLocal(seq, primer);
		alignerObj.profileAlignment(seq, primer, false, true, false);
	};
	std::vector<PrimerMatch> frontCandidates;
	std::vector<PrimerMatch> endCandidates;
	for(const auto & primer : primers_){
		alignPrimer(primer.info_);
		if(pars_.allowableErrors_.passErrorProfile(alignerObj.comp_) &&
				0 == alignerObj.parts_.lHolder_.localAStart_ &&
				0 == alignerObj.parts_.lHolder_.localBStart_){
			if(alignerObj.parts_.lHolder_.localBSize_ >= pars_.minPrimerMatchBases_){
				ret.frontMatches_.emplace_back(genMatch(primer.name_, primer.info_, false, alignerObj));
			}else if(alignerObj.parts_.lHolder_.localBSize_ >= pars_.minCandidateMatchBases_){
				frontCandidates.emplace_back(genMatch(primer.name_, primer.info_, false, alignerObj));
			}
		}
		alignPrimer(primer.infoRC_);
		if(pars_.allowableErrors_.passErrorProfile(alignerObj.comp_) &&
				len(seq) == alignerObj.parts_.lHolder_.localAStart_ + alignerObj.parts_.lHolder_.localASize_ &&
				len(primer.infoRC_) == alignerObj.parts_.lHolder_.localBStart_ + alignerObj.parts_.lHolder_.localBSize_){
			if(alignerObj.parts_.lHolder_.localBSize_ >= pars_.minPrimerMatchBases_){
				ret.endMatches_.emplace_back(genMatch(primer.name_, primer.infoRC_, true, alignerObj));
			}else if(alignerObj.parts_.lHolder_.localBSize_ >= pars_.minCandidateMatchBases_){
				endCandidates.emplace_back(genMatch(primer.name_, primer.infoRC_, true, alignerObj));
			}
		}
	}
	//a primer trusted at only one end, keep what could be its partner at the other end
	const auto addCandidates = [this, &ret](const std::vector<PrimerMatch> & trusted, const std::vector<PrimerMatch> & candidates){
		for(const auto & candidate : candidates){
			if(std::any_of(trusted.begin(), trusted.end(), [this, &candidate](const PrimerMatch & match){
				return candidate.getGap(match) <= pars_.maxGapBetweenPrimers_;
			})){
				ret.candidatePartners_.emplace_back(candidate);
			}
		}
	};
	if(ret.endMatches_.empty()){
		addCandidates(ret.frontMatches_, endCandidates);
	}else if(ret.frontMatches_.empty()){
		addCandidates(ret.endMatches_, frontCandidates);
	}
	double bestAlnScore = std::numeric_limits<double>::lowest();
	double bestIdentity = std::numeric_limits<double>::lowest();
	for(const auto & front : ret.frontMatches_){
		for(const auto & end : ret.endMatches_){
			if(front.getGap(end) > pars_.maxGapBetweenPrimers_){
				continue;
			}
			PrimerPair pair(front, end);
			const auto alnScore = pair.getAlnScoreSum();
			const auto identity = pair.getIdentitySum();
			if(alnScore > bestAlnScore || (alnScore == bestAlnScore && identity > bestIdentity)){
				bestAlnScore = alnScore;
				bestIdentity = identity;
				ret.dimers_.clear();
			}
			if(alnScore == bestAlnScore && identity == bestIdentity){
				ret.dimers_.emplace_back(std::move(pair));
			}
		}
	}
	if(!ret.dimers_.empty()){
		ret.category_ = ReadCategory::DIMER;
	}else if(!ret.frontMatches_.empty() && !ret.endMatches_.empty()){
		ret.category_ = ReadCategory::PRIMERS_FAR_APART;
	}else if(!ret.candidatePartners_.empty()){
		ret.category_ = ReadCategory::UNRESOLVED_DIMER;
	}else if(!ret.frontMatches_.empty() || !ret.endMatches_.empty()){
		ret.category_ = ReadCategory::ONE_PRIMER;
	}
	return ret;
}

}  // namespace njhseq
