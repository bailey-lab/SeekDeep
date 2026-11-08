/*
 * NanoporeAmpliconReadQC.cpp
 *
 *  Created on: Nov 6, 2026
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

#include "NanoporeAmpliconReadQC.hpp"
#include <njhseq/helpers/seqUtil.hpp>

namespace njhseq {

namespace {

constexpr uint32_t kSeedLen = 8;
constexpr uint32_t kRefKmerLen = 11;
constexpr uint32_t kBarcodeKmerLen = 15;
constexpr uint32_t kMinBarcodeKmers = 3;
constexpr uint32_t kNoKmer = std::numeric_limits<uint32_t>::max();

inline int32_t baseCode(char c){
	switch(c){
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
		default: return -1;
	}
}

/**@brief The packed kmer starting at each position of a sequence, kNoKmer where it has anything other than ACGT
 */
std::vector<uint32_t> packKmers(const std::string & seq, uint32_t kmerLen){
	std::vector<uint32_t> ret;
	if(seq.size() < kmerLen){
		return ret;
	}
	ret.assign(seq.size() - kmerLen + 1, kNoKmer);
	const uint32_t mask = (static_cast<uint64_t>(1) << (2 * kmerLen)) - 1;
	uint32_t current = 0;
	uint32_t valid = 0;
	for(size_t pos = 0; pos < seq.size(); ++pos){
		const auto code = baseCode(seq[pos]);
		if(code < 0){
			valid = 0;
			current = 0;
			continue;
		}
		current = ((current << 2) | static_cast<uint32_t>(code)) & mask;
		++valid;
		if(valid >= kmerLen){
			ret[pos + 1 - kmerLen] = current;
		}
	}
	return ret;
}

/**@brief The best place for a primer in a window with free ends in the window
 *
 * @return edits, start, end
 */
std::tuple<uint32_t, uint32_t, uint32_t> bestPlace(const std::string & primer, const std::string & seq, size_t winStart, size_t winEnd){
	const size_t n = winEnd - winStart;
	std::vector<uint32_t> prev(n + 1, 0);
	std::vector<uint32_t> prevStart(n + 1);
	std::iota(prevStart.begin(), prevStart.end(), 0);
	std::vector<uint32_t> cur(n + 1);
	std::vector<uint32_t> curStart(n + 1);
	for(size_t pos = 0; pos < primer.size(); ++pos){
		cur[0] = pos + 1;
		curStart[0] = 0;
		const char base = primer[pos];
		for(size_t col = 1; col <= n; ++col){
			const uint32_t diag = prev[col - 1] + (base != seq[winStart + col - 1] ? 1 : 0);
			const uint32_t up = prev[col] + 1;
			const uint32_t left = cur[col - 1] + 1;
			if(diag <= up && diag <= left){
				cur[col] = diag;
				curStart[col] = prevStart[col - 1];
			}else if(up <= left){
				cur[col] = up;
				curStart[col] = prevStart[col];
			}else{
				cur[col] = left;
				curStart[col] = curStart[col - 1];
			}
		}
		prev.swap(cur);
		prevStart.swap(curStart);
	}
	size_t best = 0;
	for(size_t col = 1; col <= n; ++col){
		if(prev[col] < prev[best]){
			best = col;
		}
	}
	return {prev[best], prevStart[best], best};
}

std::string gapBin(uint32_t gap){
	if(gap <= 70){
		return "gap<=70";
	}
	return gap <= 200 ? "gap71-200" : "gap>200";
}

/**@brief Shuffle the same way on every platform so the decoys, and what's reported about them, don't depend on the library
 */
std::string shuffled(const std::string & seq, std::mt19937 & gen){
	std::string ret = seq;
	for(size_t pos = ret.size(); pos > 1; --pos){
		std::swap(ret[pos - 1], ret[gen() % pos]);
	}
	return ret;
}

}  // namespace

const std::vector<NanoporeAmpliconReadQC::ReadCategory> NanoporeAmpliconReadQC::allCategories_{
	ReadCategory::AMPLICON,
	ReadCategory::AMPLICON_PLUS_OTHER,
	ReadCategory::AMPLICON_WITH_FOLD,
	ReadCategory::SAME_TARGET_INVERTED,
	ReadCategory::JOIN_SAME_TARGET,
	ReadCategory::JOIN_DIFFERENT_TARGETS,
	ReadCategory::FOLD_ON_TARGET,
	ReadCategory::FOLD_OFF_TARGET,
	ReadCategory::ON_TARGET_UNEXPECTED_LENGTH,
	ReadCategory::ON_TARGET_ONE_END,
	ReadCategory::TARGET_NO_PRIMERS,
	ReadCategory::PARTIAL_MULTI_TARGET,
	ReadCategory::PRIMER_DIMER,
	ReadCategory::TWO_PRIMER_OFF_TARGET,
	ReadCategory::ONE_PRIMER_OFF_TARGET,
	ReadCategory::NO_PRIMER
};

std::string NanoporeAmpliconReadQC::getCategoryStr(ReadCategory category){
	switch(category){
		case ReadCategory::AMPLICON: return "amplicon";
		case ReadCategory::AMPLICON_PLUS_OTHER: return "ampliconPlusOther";
		case ReadCategory::AMPLICON_WITH_FOLD: return "ampliconWithFold";
		case ReadCategory::SAME_TARGET_INVERTED: return "sameTargetInverted";
		case ReadCategory::JOIN_SAME_TARGET: return "joinSameTarget";
		case ReadCategory::JOIN_DIFFERENT_TARGETS: return "joinDifferentTargets";
		case ReadCategory::FOLD_ON_TARGET: return "foldOnTarget";
		case ReadCategory::FOLD_OFF_TARGET: return "foldOffTarget";
		case ReadCategory::ON_TARGET_UNEXPECTED_LENGTH: return "onTargetUnexpectedLength";
		case ReadCategory::ON_TARGET_ONE_END: return "onTargetOneEnd";
		case ReadCategory::TARGET_NO_PRIMERS: return "targetNoPrimers";
		case ReadCategory::PARTIAL_MULTI_TARGET: return "partialMultiTarget";
		case ReadCategory::PRIMER_DIMER: return "primerDimer";
		case ReadCategory::TWO_PRIMER_OFF_TARGET: return "twoPrimerOffTarget";
		case ReadCategory::ONE_PRIMER_OFF_TARGET: return "onePrimerOffTarget";
		case ReadCategory::NO_PRIMER: return "noPrimer";
	}
	std::stringstream ss;
	ss << __PRETTY_FUNCTION__ << ", error unknown category" << "\n";
	throw std::runtime_error{ss.str()};
}

NanoporeAmpliconReadQC::Panel NanoporeAmpliconReadQC::buildPanel(std::vector<Probe> probes){
	Panel ret;
	ret.probes_ = std::move(probes);
	ret.seedIndex_.resize(static_cast<size_t>(1) << (2 * kSeedLen));
	for(uint32_t probePos = 0; probePos < ret.probes_.size(); ++probePos){
		const auto seeds = packKmers(ret.probes_[probePos].seq_, kSeedLen);
		for(uint32_t pos = 0; pos < seeds.size(); ++pos){
			if(kNoKmer != seeds[pos]){
				ret.seedIndex_[seeds[pos]].emplace_back(probePos, pos);
			}
		}
	}
	return ret;
}

NanoporeAmpliconReadQC::NanoporeAmpliconReadQC(const PrimersAndMids & ids, Pars pars) : pars_(std::move(pars)), detector_(pars_.foldPars_){
	std::vector<Probe> probes;
	std::vector<Probe> decoys;
	//fixed so the decoys are the same every time
	std::mt19937 gen(7);
	const auto addPrimers = [&](const std::string & primers, uint32_t target, bool forward, uint32_t & firstLen){
		for(const auto & primer : tokenizeString(primers, ",")){
			const auto upperPrimer = njh::strToUpperRet(primer);
			if(0 == firstLen){
				firstLen = upperPrimer.size();
			}
			for(const auto & seq : createDegenStrs(upperPrimer)){
				probes.emplace_back(Probe{seq, target, forward, false});
				probes.emplace_back(Probe{seqUtil::reverseComplement(seq, "DNA"), target, forward, true});
				const auto decoy = shuffled(seq, gen);
				decoys.emplace_back(Probe{decoy, target, forward, false});
				decoys.emplace_back(Probe{seqUtil::reverseComplement(decoy, "DNA"), target, forward, true});
			}
		}
	};
	for(const auto & tar : ids.targets_){
		const uint32_t target = targetNames_.size();
		targetNames_.emplace_back(tar.first);
		uint32_t forwardLen = 0;
		uint32_t reverseLen = 0;
		addPrimers(tar.second.info_.forwardPrimerRaw_, target, true, forwardLen);
		addPrimers(tar.second.info_.reversePrimerRaw_, target, false, reverseLen);
		std::unordered_set<uint32_t> kmers;
		std::vector<uint32_t> refLens;
		for(const auto & ref : tar.second.refs_){
			auto seq = ref.seq_;
			njh::strToUpper(seq);
			refLens.emplace_back(seq.size());
			for(const auto & strand : {seq, seqUtil::reverseComplement(seq, "DNA")}){
				for(const auto kmer : packKmers(strand, kRefKmerLen)){
					if(kNoKmer != kmer){
						kmers.emplace(kmer);
					}
				}
			}
		}
		//the reference sequences don't have the primers on them
		expectedLens_.emplace_back(refLens.empty() ? 0 : static_cast<uint32_t>(vectorMedianRef(refLens)) + forwardLen + reverseLen);
		for(const auto kmer : kmers){
			kmerToTargets_[kmer].emplace_back(target);
		}
		targetKmers_.emplace_back(std::move(kmers));
	}
	panel_ = buildPanel(std::move(probes));
	decoyPanel_ = buildPanel(std::move(decoys));
}

void NanoporeAmpliconReadQC::setBarcodes(const bfs::path & barcodesFnp, uint32_t barcodeStart, uint32_t barcodeLen){
	barcodeKmers_.clear();
	SeqInput reader(SeqIOOptions::genFastaIn(barcodesFnp));
	reader.openIn();
	seqInfo seq;
	while(reader.readNextRead(seq)){
		std::string digits;
		for(const auto c : seq.name_){
			if(std::isdigit(static_cast<unsigned char>(c))){
				digits.push_back(c);
			}
		}
		if(digits.empty()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error no number in barcode name " << seq.name_ << " in " << barcodesFnp << "\n";
			throw std::runtime_error{ss.str()};
		}
		if(seq.seq_.size() < barcodeStart + barcodeLen){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << seq.name_ << " is shorter than a barcode start of " << barcodeStart << " and length of " << barcodeLen << "\n";
			throw std::runtime_error{ss.str()};
		}
		const auto barcode = njh::strToUpperRet(seq.seq_.substr(barcodeStart, barcodeLen));
		for(const auto & strand : {barcode, seqUtil::reverseComplement(barcode, "DNA")}){
			for(const auto kmer : packKmers(strand, kBarcodeKmerLen)){
				if(kNoKmer != kmer){
					auto & found = barcodeKmers_[kmer];
					const int32_t number = std::stoi(digits);
					if(!njh::in(number, found)){
						found.emplace_back(number);
					}
				}
			}
		}
	}
}

std::string NanoporeAmpliconReadQC::junctionBarcodeStatus(const ReadResult & result, uint32_t junction, const std::set<int32_t> & ownBarcodes){
	if(result.junctionBarcodes_.empty()){
		return "notChecked";
	}
	const auto barcode = result.junctionBarcodes_[junction];
	if(barcode < 0){
		return "noBarcode";
	}
	if(ownBarcodes.empty()){
		return "barcodeUnknownSample";
	}
	return njh::in(barcode, ownBarcodes) ? "own" : "foreign";
}

const std::vector<std::string> & NanoporeAmpliconReadQC::targetNames() const{
	return targetNames_;
}

uint32_t NanoporeAmpliconReadQC::expectedLen(uint32_t target) const{
	return expectedLens_[target];
}

const NanoporeAmpliconReadQC::Probe & NanoporeAmpliconReadQC::probe(uint32_t pos) const{
	return panel_.probes_[pos];
}

std::vector<NanoporeAmpliconReadQC::PrimerHit> NanoporeAmpliconReadQC::findHits(const Panel & panel, const std::string & seq) const{
	std::vector<PrimerHit> ret;
	const auto seeds = packKmers(seq, kSeedLen);
	//seeds that say about the same place for a probe are only looked at once, each has a window with room for indels
	std::vector<std::pair<uint32_t, int64_t>> candidates;
	const auto bucketDown = [](int64_t val){
		return val >= 0 ? val / 6 : -((-val + 5) / 6);
	};
	for(size_t pos = 0; pos < seeds.size(); ++pos){
		if(kNoKmer == seeds[pos]){
			continue;
		}
		for(const auto & found : panel.seedIndex_[seeds[pos]]){
			candidates.emplace_back(found.first, bucketDown(static_cast<int64_t>(pos) - static_cast<int64_t>(found.second) + 3));
		}
	}
	std::sort(candidates.begin(), candidates.end());
	candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
	std::vector<PrimerHit> forProbe;
	const auto addBest = [&ret, &forProbe](){
		std::stable_sort(forProbe.begin(), forProbe.end(), [](const PrimerHit & h1, const PrimerHit & h2){
			return std::tie(h1.edits_, h1.start_) < std::tie(h2.edits_, h2.start_);
		});
		const size_t alreadyFrom = ret.size();
		for(const auto & hit : forProbe){
			bool overlaps = false;
			for(size_t pos = alreadyFrom; pos < ret.size(); ++pos){
				if(hit.end_ > ret[pos].start_ && hit.start_ < ret[pos].end_){
					overlaps = true;
					break;
				}
			}
			if(!overlaps){
				ret.emplace_back(hit);
			}
		}
		forProbe.clear();
	};
	int64_t lastProbe = -1;
	for(const auto & candidate : candidates){
		if(static_cast<int64_t>(candidate.first) != lastProbe){
			addBest();
			lastProbe = candidate.first;
		}
		const auto & primer = panel.probes_[candidate.first].seq_;
		const size_t winStart = candidate.second * 6 - 6 > 0 ? candidate.second * 6 - 6 : 0;
		if(winStart >= seq.size()){
			continue;
		}
		const size_t winEnd = std::min(winStart + primer.size() + 12, seq.size());
		if(winEnd - winStart < 0.7 * primer.size()){
			continue;
		}
		const auto [edits, start, end] = bestPlace(primer, seq, winStart, winEnd);
		if(edits <= pars_.maxPrimerError_ * primer.size()){
			forProbe.emplace_back(PrimerHit{static_cast<uint32_t>(winStart + start), static_cast<uint32_t>(winStart + end), candidate.first, edits});
		}
	}
	addBest();
	return ret;
}

std::vector<NanoporeAmpliconReadQC::Unit> NanoporeAmpliconReadQC::findUnits(const Panel & panel, const std::vector<PrimerHit> & hits, const std::vector<uint32_t> & kmers, double primerError) const{
	struct ByTarget{
		std::vector<const PrimerHit *> forwardAsWritten;
		std::vector<const PrimerHit *> forwardRevComp;
		std::vector<const PrimerHit *> reverseAsWritten;
		std::vector<const PrimerHit *> reverseRevComp;
	};
	std::map<uint32_t, ByTarget> byTarget;
	for(const auto & hit : hits){
		const auto & hitProbe = panel.probes_[hit.probe_];
		if(hit.edits_ > primerError * hitProbe.seq_.size()){
			continue;
		}
		auto & tar = byTarget[hitProbe.target_];
		if(hitProbe.forwardPrimer_){
			(hitProbe.revComp_ ? tar.forwardRevComp : tar.forwardAsWritten).emplace_back(&hit);
		}else{
			(hitProbe.revComp_ ? tar.reverseRevComp : tar.reverseAsWritten).emplace_back(&hit);
		}
	}
	std::vector<Unit> candidates;
	const auto addCandidates = [&](uint32_t target, const std::vector<const PrimerHit *> & firsts, const std::vector<const PrimerHit *> & lasts, bool forward){
		const auto expected = expectedLens_[target];
		if(0 == expected){
			return;
		}
		const double tolerance = std::max<double>(pars_.minLenTolerance_, pars_.lenToleranceFrac_ * expected);
		for(const auto * first : firsts){
			for(const auto * last : lasts){
				if(last->start_ + 2 < first->end_ || last->end_ < first->start_){
					continue;
				}
				const double length = static_cast<double>(last->end_) - first->start_;
				if(std::abs(length - expected) > tolerance){
					continue;
				}
				uint32_t innerStart = first->end_;
				uint32_t innerEnd = last->start_ > first->end_ ? last->start_ : first->end_;
				if(innerEnd - innerStart < kRefKmerLen + 5){
					innerStart = first->start_;
					innerEnd = last->end_;
				}
				const uint32_t total = innerEnd - innerStart >= kRefKmerLen ? innerEnd - innerStart - kRefKmerLen + 1 : 1;
				uint32_t matched = 0;
				if(innerEnd - innerStart >= kRefKmerLen){
					const auto & refKmers = targetKmers_[target];
					for(uint32_t pos = innerStart; pos < innerStart + total; ++pos){
						if(kNoKmer != kmers[pos] && njh::in(kmers[pos], refKmers)){
							++matched;
						}
					}
				}
				const double frac = static_cast<double>(matched) / total;
				if(frac >= pars_.minInnerKmerFrac_){
					candidates.emplace_back(Unit{first->start_, last->end_, target, forward, first->edits_, last->edits_, frac});
				}
			}
		}
	};
	for(const auto & tar : byTarget){
		addCandidates(tar.first, tar.second.forwardAsWritten, tar.second.reverseRevComp, true);
		addCandidates(tar.first, tar.second.reverseAsWritten, tar.second.forwardRevComp, false);
	}
	std::stable_sort(candidates.begin(), candidates.end(), [](const Unit & u1, const Unit & u2){
		const auto e1 = u1.startEdits_ + u1.endEdits_;
		const auto e2 = u2.startEdits_ + u2.endEdits_;
		return std::make_tuple(u1.end_, e1, -u1.innerFrac_) < std::make_tuple(u2.end_, e2, -u2.innerFrac_);
	});
	std::vector<Unit> ret;
	int64_t lastEnd = -1;
	for(const auto & unit : candidates){
		//primers of neighbours can overlap by a few bases
		if(static_cast<int64_t>(unit.start_) >= lastEnd - 3){
			ret.emplace_back(unit);
			lastEnd = unit.end_;
		}
	}
	return ret;
}

std::tuple<int64_t, double, double> NanoporeAmpliconReadQC::bestTarget(const std::vector<uint32_t> & kmers, uint32_t start, uint32_t end) const{
	if(kmers.empty() || end <= start || start >= kmers.size()){
		return {-1, 0, 0};
	}
	const uint32_t stop = std::min<uint64_t>(end >= kRefKmerLen ? end - kRefKmerLen + 1 : 0, kmers.size());
	if(stop <= start){
		return {-1, 0, 0};
	}
	std::unordered_map<uint32_t, uint32_t> counts;
	for(uint32_t pos = start; pos < stop; ++pos){
		if(kNoKmer == kmers[pos]){
			continue;
		}
		const auto found = kmerToTargets_.find(kmers[pos]);
		if(found != kmerToTargets_.end()){
			for(const auto target : found->second){
				++counts[target];
			}
		}
	}
	int64_t best = -1;
	uint32_t bestCount = 0;
	uint32_t nextCount = 0;
	for(const auto & count : counts){
		//the lower target wins a tie so that this doesn't depend on the hash's order
		if(count.second > bestCount || (count.second == bestCount && static_cast<int64_t>(count.first) < best)){
			nextCount = bestCount;
			bestCount = count.second;
			best = count.first;
		}else if(count.second > nextCount){
			nextCount = count.second;
		}
	}
	const double total = stop - start;
	return {best, bestCount / total, nextCount / total};
}

NanoporeAmpliconReadQC::ReadResult NanoporeAmpliconReadQC::processRead(const std::string & seq, bool checkDecoys) const{
	ReadResult ret;
	ret.len_ = seq.size();
	auto upperSeq = seq;
	njh::strToUpper(upperSeq);
	const auto refKmers = packKmers(upperSeq, kRefKmerLen);
	const auto hits = findHits(panel_, upperSeq);
	ret.primerHits_ = hits.size();
	ret.units_ = findUnits(panel_, hits, refKmers, pars_.maxPrimerError_);
	for(const auto & primerError : stricterPrimerErrors_){
		ret.unitsAtStricterErrors_.emplace_back(findUnits(panel_, hits, refKmers, primerError).size());
	}
	if(checkDecoys){
		ret.decoyChecked_ = true;
		const auto decoyHits = findHits(decoyPanel_, upperSeq);
		ret.decoyHits_ = decoyHits.size();
		ret.decoyUnitsAtErrors_.emplace_back(findUnits(decoyPanel_, decoyHits, refKmers, pars_.maxPrimerError_).size());
		ret.unitsAtErrors_.emplace_back(ret.units_.size());
		for(const auto & primerError : stricterPrimerErrors_){
			ret.decoyUnitsAtErrors_.emplace_back(findUnits(decoyPanel_, decoyHits, refKmers, primerError).size());
		}
		ret.unitsAtErrors_.insert(ret.unitsAtErrors_.end(), ret.unitsAtStricterErrors_.begin(), ret.unitsAtStricterErrors_.end());
	}
	ret.fold_ = detector_.detect(upperSeq);
	if(!barcodeKmers_.empty() && ret.units_.size() >= 2){
		ret.junctionBarcodes_.assign(ret.units_.size() - 1, -1);
		for(size_t pos = 0; pos + 1 < ret.units_.size(); ++pos){
			//a little of each amplicon so a barcode right against the primer still has all its kmers
			const size_t start = ret.units_[pos].end_ > 3 ? ret.units_[pos].end_ - 3 : 0;
			const size_t end = std::min<size_t>(ret.units_[pos + 1].start_ + 3, upperSeq.size());
			if(end < start + kBarcodeKmerLen){
				continue;
			}
			std::map<int32_t, uint32_t> counts;
			for(const auto kmer : packKmers(upperSeq.substr(start, end - start), kBarcodeKmerLen)){
				if(kNoKmer == kmer){
					continue;
				}
				const auto found = barcodeKmers_.find(kmer);
				if(found != barcodeKmers_.end()){
					for(const auto barcode : found->second){
						++counts[barcode];
					}
				}
			}
			uint32_t bestCount = 0;
			for(const auto & count : counts){
				if(count.second > bestCount){
					bestCount = count.second;
					ret.junctionBarcodes_[pos] = count.first;
				}
			}
			if(bestCount < kMinBarcodeKmers){
				ret.junctionBarcodes_[pos] = -1;
			}
		}
	}

	for(const auto & hit : hits){
		const auto & hitProbe = panel_.probes_[hit.probe_];
		if(!hitProbe.revComp_ && hit.start_ <= pars_.primerAtStartWithin_
				&& (ret.startProbe_ < 0 || hit.start_ < ret.startProbePos_)){
			ret.startProbe_ = hit.probe_;
			ret.startProbePos_ = hit.start_;
		}
	}

	ret.keepOnce_.assign(ret.units_.size(), true);
	for(size_t pos = 0; pos + 1 < ret.units_.size(); ++pos){
		if(ret.units_[pos].target_ == ret.units_[pos + 1].target_
				&& ret.units_[pos].startsWithForwardPrimer_ != ret.units_[pos + 1].startsWithForwardPrimer_
				&& ret.keepOnce_[pos]){
			ret.keepOnce_[pos + 1] = false;
		}
	}

	if(!ret.units_.empty()){
		ret.target_ = ret.units_.front().target_;
		if(1 == ret.units_.size()){
			if(ret.fold_.called_){
				ret.category_ = ReadCategory::AMPLICON_WITH_FOLD;
			}else{
				const uint64_t outside = seq.size() - (ret.units_.front().end_ - ret.units_.front().start_);
				ret.category_ = outside > pars_.otherBasesCutOff_ ? ReadCategory::AMPLICON_PLUS_OTHER : ReadCategory::AMPLICON;
			}
			return ret;
		}
		bool sameTarget = true;
		bool allInverted = true;
		for(size_t pos = 0; pos + 1 < ret.units_.size(); ++pos){
			sameTarget = sameTarget && ret.units_[pos].target_ == ret.units_[pos + 1].target_;
			allInverted = allInverted && ret.units_[pos].startsWithForwardPrimer_ != ret.units_[pos + 1].startsWithForwardPrimer_;
		}
		if(!sameTarget){
			ret.category_ = ReadCategory::JOIN_DIFFERENT_TARGETS;
		}else{
			ret.category_ = allInverted ? ReadCategory::SAME_TARGET_INVERTED : ReadCategory::JOIN_SAME_TARGET;
		}
		return ret;
	}

	if(ret.fold_.called_){
		const auto [armTarget, armFrac, armOtherFrac] = bestTarget(refKmers, ret.fold_.leftStart_, ret.fold_.leftEnd_);
		if(armTarget >= 0 && armFrac >= pars_.minTargetKmerFrac_){
			ret.category_ = ReadCategory::FOLD_ON_TARGET;
			ret.target_ = armTarget;
		}else{
			ret.category_ = ReadCategory::FOLD_OFF_TARGET;
		}
		return ret;
	}

	const auto [readTarget, readFrac, readOtherFrac] = bestTarget(refKmers, 0, seq.size());
	if(readTarget >= 0 && readFrac >= pars_.minTargetKmerFrac_){
		ret.target_ = readTarget;
		if(readOtherFrac >= pars_.otherTargetKmerFrac_){
			ret.category_ = ReadCategory::PARTIAL_MULTI_TARGET;
			return ret;
		}
		bool asWritten = false;
		bool revComp = false;
		bool anyPrimer = false;
		for(const auto & hit : hits){
			const auto & hitProbe = panel_.probes_[hit.probe_];
			if(static_cast<int64_t>(hitProbe.target_) == readTarget){
				anyPrimer = true;
				asWritten = asWritten || !hitProbe.revComp_;
				revComp = revComp || hitProbe.revComp_;
			}
		}
		if(asWritten && revComp){
			ret.category_ = ReadCategory::ON_TARGET_UNEXPECTED_LENGTH;
		}else{
			ret.category_ = anyPrimer ? ReadCategory::ON_TARGET_ONE_END : ReadCategory::TARGET_NO_PRIMERS;
		}
		return ret;
	}

	//off target, what primers are on it
	if(hits.empty()){
		ret.category_ = ReadCategory::NO_PRIMER;
		return ret;
	}
	{
		std::set<uint32_t> primers;
		for(const auto & hit : hits){
			primers.emplace(primerId(panel_.probes_[hit.probe_]));
		}
		ret.offTargetPrimers_.assign(primers.begin(), primers.end());
	}
	bool twoPrimers = false;
	const PrimerHit * outerFirst = nullptr;
	const PrimerHit * outerLast = nullptr;
	for(const auto & first : hits){
		if(panel_.probes_[first.probe_].revComp_){
			continue;
		}
		for(const auto & last : hits){
			if(!panel_.probes_[last.probe_].revComp_ || last.start_ + 2 < first.end_ || last.end_ <= first.start_){
				continue;
			}
			twoPrimers = true;
			//the pair furthest apart is the one that made the product, the others are inside it
			if(nullptr == outerFirst || last.end_ - first.start_ > outerLast->end_ - outerFirst->start_){
				outerFirst = &first;
				outerLast = &last;
			}
			//a dimer is the whole read, two primers close together on a long read are something else with primers in it
			const int64_t gap = static_cast<int64_t>(last.start_) - first.end_;
			const uint64_t outside = seq.size() - (last.end_ - first.start_);
			if(gap <= pars_.maxDimerInsert_ && outside <= pars_.otherBasesCutOff_){
				ret.category_ = ReadCategory::PRIMER_DIMER;
				ret.dimerProbe1_ = first.probe_;
				ret.dimerProbe2_ = last.probe_;
				return ret;
			}
		}
	}
	ret.category_ = twoPrimers ? ReadCategory::TWO_PRIMER_OFF_TARGET : ReadCategory::ONE_PRIMER_OFF_TARGET;
	if(twoPrimers){
		//the read can be either strand, F...rc(R) on one is R...rc(F) on the other, so the pair is kept in one order
		const auto firstPrimer = primerId(panel_.probes_[outerFirst->probe_]);
		const auto lastPrimer = primerId(panel_.probes_[outerLast->probe_]);
		ret.offTargetPrimer1_ = std::min(firstPrimer, lastPrimer);
		ret.offTargetPrimer2_ = std::max(firstPrimer, lastPrimer);
		ret.offTargetInsert_ = outerLast->start_ > outerFirst->end_ ? outerLast->start_ - outerFirst->end_ : 0;
		ret.offTargetAtEnds_ = outerFirst->start_ <= pars_.primerAtStartWithin_ && seq.size() - outerLast->end_ <= pars_.primerAtStartWithin_;
	}
	return ret;
}

uint32_t NanoporeAmpliconReadQC::primerId(const Probe & probe){
	return probe.target_ * 2 + (probe.forwardPrimer_ ? 0 : 1);
}

std::string NanoporeAmpliconReadQC::primerName(uint32_t primer) const{
	return targetNames_[primer / 2] + (0 == primer % 2 ? "_F" : "_R");
}

std::string NanoporeAmpliconReadQC::primerPairType(uint32_t primer1, uint32_t primer2){
	if(primer1 == primer2){
		return "samePrimer";
	}
	return primer1 / 2 == primer2 / 2 ? "ownPair" : "diffTargets";
}

bool NanoporeAmpliconReadQC::foldTurnPos(const ReadResult & result, uint32_t & pos) const{
	if(!result.fold_.called_ || result.target_ < 0 || result.startProbe_ < 0){
		return false;
	}
	const auto & startProbe = panel_.probes_[result.startProbe_];
	if(static_cast<int64_t>(startProbe.target_) != result.target_ || result.fold_.rightStart_ < result.startProbePos_){
		return false;
	}
	const uint32_t fromStart = result.fold_.rightStart_ - result.startProbePos_;
	const auto expected = expectedLens_[result.target_];
	if(startProbe.forwardPrimer_){
		pos = fromStart;
	}else{
		pos = expected > fromStart ? expected - fromStart : 0;
	}
	return true;
}

void NanoporeAmpliconReadQC::addToTallies(const ReadResult & result, SampleTallies & tallies, const std::set<int32_t> & ownBarcodes) const{
	++tallies.reads_;
	tallies.bases_ += result.len_;
	++tallies.categoryReads_[result.category_];
	tallies.categoryBases_[result.category_] += result.len_;
	std::set<uint32_t> readTargets;
	for(const auto & unit : result.units_){
		readTargets.emplace(unit.target_);
	}
	if(result.target_ >= 0){
		readTargets.emplace(result.target_);
	}
	for(const auto target : readTargets){
		++tallies.targetCategoryReads_[{target, result.category_}];
		++tallies.targetReads_[target];
	}
	if(1 == result.units_.size()){
		++tallies.singleAmpliconReads_[result.units_.front().target_];
	}else if(result.units_.size() >= 2){
		++tallies.readsWithUnits_[std::min<uint32_t>(result.units_.size(), 5)];
		if(tallies.joinedReadsAtStricterErrors_.empty()){
			tallies.joinedReadsAtStricterErrors_.assign(stricterPrimerErrors_.size(), 0);
		}
		for(size_t pos = 0; pos < stricterPrimerErrors_.size(); ++pos){
			if(result.unitsAtStricterErrors_[pos] >= 2){
				++tallies.joinedReadsAtStricterErrors_[pos];
			}
		}
		bool foreign = false;
		for(size_t pos = 0; pos + 1 < result.units_.size(); ++pos){
			const auto & first = result.units_[pos];
			const auto & second = result.units_[pos + 1];
			const uint32_t gap = second.start_ > first.end_ ? second.start_ - first.end_ : 0;
			const auto status = junctionBarcodeStatus(result, pos, ownBarcodes);
			foreign = foreign || "foreign" == status;
			const auto key = std::make_tuple(first.target_ == second.target_, first.startsWithForwardPrimer_ != second.startsWithForwardPrimer_, gapBin(gap), status);
			++tallies.junctions_[key];
			tallies.junctionGaps_[key].emplace_back(gap);
		}
		if(foreign){
			++tallies.readsWithForeignBarcode_;
		}else{
			//amplicons from a read with another sample's barcode in it can't be said to belong to this sample
			for(size_t pos = 0; pos < result.units_.size(); ++pos){
				++tallies.unitsInJoins_[result.units_[pos].target_];
				if(result.keepOnce_[pos]){
					++tallies.unitsCountedOnce_[result.units_[pos].target_];
				}
			}
		}
	}
	if(result.fold_.called_){
		++tallies.foldsCalled_;
		if(result.target_ >= 0 && (ReadCategory::FOLD_ON_TARGET == result.category_ || ReadCategory::AMPLICON_WITH_FOLD == result.category_)){
			const auto target = static_cast<uint32_t>(result.target_);
			++tallies.targetFolds_[target];
			tallies.targetFoldLoops_[target].emplace_back(result.fold_.loopLen());
			tallies.targetFoldArms_[target].emplace_back(result.fold_.shorterArmLen());
			uint32_t pos = 0;
			if(foldTurnPos(result, pos)){
				const bool startedForward = panel_.probes_[result.startProbe_].forwardPrimer_;
				++tallies.targetFoldStrand_[{target, startedForward}];
				if(result.fold_.loopLen() >= pars_.minLoopForPileup_){
					++tallies.turnPileup_[std::make_tuple(target, startedForward, pos / pars_.turnBinSize_ * pars_.turnBinSize_)];
				}
			}
		}
	}
	if(ReadCategory::PRIMER_DIMER == result.category_){
		const auto primerName = [this](int64_t probePos){
			const auto & dimerProbe = panel_.probes_[probePos];
			return targetNames_[dimerProbe.target_] + (dimerProbe.forwardPrimer_ ? "_F" : "_R");
		};
		++tallies.dimers_[{primerName(result.dimerProbe1_), primerName(result.dimerProbe2_)}];
	}
	if(ReadCategory::TWO_PRIMER_OFF_TARGET == result.category_ || ReadCategory::ONE_PRIMER_OFF_TARGET == result.category_){
		const size_t col = ReadCategory::TWO_PRIMER_OFF_TARGET == result.category_ ? 0 : 1;
		for(const auto primer : result.offTargetPrimers_){
			++tallies.offTargetPrimerReads_[primer][col];
		}
	}
	if(ReadCategory::TWO_PRIMER_OFF_TARGET == result.category_){
		auto & pair = tallies.offTargetPairs_[{static_cast<uint32_t>(result.offTargetPrimer1_), static_cast<uint32_t>(result.offTargetPrimer2_)}];
		++pair.reads_;
		pair.bases_ += result.len_;
		pair.atEnds_ += result.offTargetAtEnds_ ? 1 : 0;
		pair.lens_.emplace_back(result.len_);
		pair.inserts_.emplace_back(result.offTargetInsert_);
	}
	if(result.decoyChecked_){
		++tallies.decoyReads_;
		tallies.decoyHits_ += result.decoyHits_;
		tallies.realHitsOnDecoyReads_ += result.primerHits_;
		std::vector<double> errors{pars_.maxPrimerError_};
		errors.insert(errors.end(), stricterPrimerErrors_.begin(), stricterPrimerErrors_.end());
		for(size_t pos = 0; pos < errors.size(); ++pos){
			auto & real = tallies.calibration_[{false, errors[pos]}];
			auto & decoy = tallies.calibration_[{true, errors[pos]}];
			real[0] += result.unitsAtErrors_[pos] >= 1 ? 1 : 0;
			real[1] += result.unitsAtErrors_[pos] >= 2 ? 1 : 0;
			decoy[0] += result.decoyUnitsAtErrors_[pos] >= 1 ? 1 : 0;
			decoy[1] += result.decoyUnitsAtErrors_[pos] >= 2 ? 1 : 0;
		}
	}
}

std::string NanoporeAmpliconReadQC::unitsStr(const ReadResult & result) const{
	std::stringstream out;
	for(size_t pos = 0; pos < result.units_.size(); ++pos){
		const auto & unit = result.units_[pos];
		out << (pos > 0 ? ";" : "") << targetNames_[unit.target_]
				<< "|" << (unit.startsWithForwardPrimer_ ? "F" : "R")
				<< "|" << unit.start_ << "|" << unit.end_
				<< "|" << unit.startEdits_ << "|" << unit.endEdits_
				<< "|" << unit.innerFrac_;
	}
	return out.str();
}

std::string NanoporeAmpliconReadQC::junctionsStr(const ReadResult & result, const std::set<int32_t> & ownBarcodes) const{
	std::stringstream out;
	for(size_t pos = 0; pos + 1 < result.units_.size(); ++pos){
		const auto & first = result.units_[pos];
		const auto & second = result.units_[pos + 1];
		out << (pos > 0 ? "," : "") << (first.target_ == second.target_ ? "sameTarget" : "diffTarget")
				<< "|" << (first.startsWithForwardPrimer_ != second.startsWithForwardPrimer_ ? "inverted" : "headToTail")
				<< "|" << (second.start_ > first.end_ ? second.start_ - first.end_ : 0)
				<< "|" << junctionBarcodeStatus(result, pos, ownBarcodes);
	}
	return out.str();
}

}  // namespace njhseq
