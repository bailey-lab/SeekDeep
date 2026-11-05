/*
 * AmpliconReadQC.cpp
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

#include "AmpliconReadQC.hpp"
#include <njhseq/IO/SeqIO/SeqInput.hpp>
#include <njhseq/IO/SeqIO/SeqOutput.hpp>
#include <njhseq/objects/BioDataObject/BioDataFileIO.hpp>
#include <njhseq/objects/BioDataObject/BLASTHitTabular.hpp>
#include <njhseq/objects/BioDataObject/reading.hpp>
#include <njhseq/objects/kmer/kmerInfo.hpp>
#include <njhseq/system/BioCmdsUtils.hpp>

namespace njhseq {

const std::vector<AmpliconReadQC::ReadPairCategory> AmpliconReadQC::allCategories_{
	ReadPairCategory::ON_TARGET,
	ReadPairCategory::UNEXPECTED_STITCH,
	ReadPairCategory::KMER_MISMATCH,
	ReadPairCategory::PRIMER_DIMER,
	ReadPairCategory::UNRESOLVED_DIMER,
	ReadPairCategory::SHORT_PRODUCT,
	ReadPairCategory::SHORT_ONE_PRIMER,
	ReadPairCategory::SHORT_NO_PRIMER,
	ReadPairCategory::MISMATCHED_PRIMERS,
	ReadPairCategory::FAILED_PRIMERS
};

std::string AmpliconReadQC::getCategoryStr(ReadPairCategory category){
	switch (category) {
		case ReadPairCategory::ON_TARGET:
			return "onTarget";
		case ReadPairCategory::UNEXPECTED_STITCH:
			return "unexpectedStitch";
		case ReadPairCategory::KMER_MISMATCH:
			return "kmerMismatch";
		case ReadPairCategory::PRIMER_DIMER:
			return "primerDimer";
		case ReadPairCategory::UNRESOLVED_DIMER:
			return "unresolvedDimer";
		case ReadPairCategory::SHORT_PRODUCT:
			return "shortProduct";
		case ReadPairCategory::SHORT_ONE_PRIMER:
			return "shortOnePrimer";
		case ReadPairCategory::SHORT_NO_PRIMER:
			return "shortNoPrimer";
		case ReadPairCategory::MISMATCHED_PRIMERS:
			return "mismatchedPrimers";
		case ReadPairCategory::FAILED_PRIMERS:
			return "failedPrimers";
	}
	return "notHandled";//shouldn't be getting here
}

AmpliconReadQC::QcPars::QcPars(){
	//the same as what extractorPairedEnd defaults to
	primerPars_.allowable_.hqMismatches_ = 3;
	primerPars_.allowable_.lqMismatches_ = 5;
	primerPars_.allowable_.distances_.query_.coverage_ = .90;
	primerPars_.allowable_.largeBaseIndel_ = 0.99;
	primerPars_.allowable_.oneBaseIndel_ = 0.5;
	primerPars_.allowable_.twoBaseIndel_ = 0.5;
	primerPars_.useMotif_ = true;
	stitchPars_.r1Trim_ = 1;
	stitchPars_.r2Trim_ = 1;
	stitchPars_.minOverlap_ = 8;
	stitchPars_.errorAllowed_ = 0.02;
	stitchPars_.stitchCacheMaxEntries_ = 500000;
}

AmpliconReadQC::AmpliconReadQC(const PrimersAndMids & ids, QcPars pars) :
		ids_(ids), pars_(std::move(pars)), detector_(*ids.pDeterminator_, pars_.dimerPars_){
}

bool AmpliconReadQC::expectedStatus(const std::string & target, PairedReadProcessor::ReadPairOverLapStatus status) const{
	using OverlapStatus = PairedReadProcessor::ReadPairOverLapStatus;
	const auto & expected = ids_.targets_.at(target).overlapStatuses_;
	if(njh::in(status, expected)){
		return true;
	}
	//the mates being the same size as the fragment can happen to a target expected to overlap in either way
	return OverlapStatus::PERFECTOVERLAP == status &&
			(njh::in(OverlapStatus::R1BEGINSINR2, expected) || njh::in(OverlapStatus::R1ENDSINR2, expected));
}

AmpliconReadQC::SampleResults AmpliconReadQC::processSample(const SeqIOOptions & readOpts, const std::string & sampleName) const{
	using OverlapStatus = PairedReadProcessor::ReadPairOverLapStatus;
	const std::string unrecognized = "unrecognized";
	SampleResults ret;
	ret.sampleName_ = sampleName;

	PairedReadProcessor stitcher(pars_.stitchPars_);
	auto retryPars = pars_.stitchPars_;
	retryPars.errorAllowed_ = pars_.dimerOverLapErrorAllowed_;
	PairedReadProcessor retryStitcher(retryPars);
	PairedReadProcessor::ProcessedResultsCounts stitchCounts;
	//the alignment of the mates doesn't depend on the errors allowed so both attempts can share the aligner and its cache
	const uint64_t startingMateSize = 350;
	aligner stitchingAligner = stitcher.genStitchingAligner(startingMateSize);
	stitchingAligner.qScorePars_.qualThresWindow_ = 0;

	aligner dimerAligner = detector_.genAligner(detector_.getMaxDimerSize());

	const auto & primerDeterminator = *ids_.pDeterminator_;
	const size_t primerWindowSize = pars_.primerPars_.primerWithin_ + primerDeterminator.getMaxPrimerSize() + 5;
	aligner primerAligner(primerWindowSize + 10, gapScoringParameters(5, 1, 0, 0, 0, 0), substituteMatrix::createDegenScoreMatrixLessN(2, -2));
	//the search is done on a copy of the read and only what it found is remembered, the primer found gets lower cased in the copy
	struct PrimerCall{
		std::string name_;
		uint32_t primerLen_{0};
	};
	std::array<std::unordered_map<std::string, PrimerCall>, 2> primerCallCaches;
	const auto callPrimer = [&](size_t cacheIdx, const seqInfo & read) -> PrimerCall {
		const auto determine = [&](){
			seqInfo readCopy = read;
			PrimerCall call;
			call.name_ = 0 == cacheIdx ?
					primerDeterminator.determineForwardPrimer(readCopy, pars_.primerPars_, primerAligner) :
					primerDeterminator.determineWithReversePrimer(readCopy, pars_.primerPars_, primerAligner);
			call.primerLen_ = std::count_if(readCopy.seq_.begin(), readCopy.seq_.end(), [](char base){
				return std::islower(static_cast<unsigned char>(base));
			});
			return call;
		};
		if(0 == pars_.primerCallCacheMaxEntries_ || read.seq_.size() < pars_.primerPars_.primerStart_){
			return determine();
		}
		const std::string window = read.seq_.substr(pars_.primerPars_.primerStart_, primerWindowSize);
		auto & cache = primerCallCaches[cacheIdx];
		const auto found = cache.find(window);
		if(cache.end() != found){
			return found->second;
		}
		if(cache.size() >= pars_.primerCallCacheMaxEntries_){
			cache.clear();
		}
		auto call = determine();
		cache.emplace(window, call);
		return call;
	};

	//the same comparison as filtering with --compareSeq in extractorPairedEnd, primers are left out since the references don't have them
	std::unordered_map<std::string, bool> refMatchCache;
	const auto compareToRefs = [&](const std::string & target, const std::string & cacheKey, const std::function<bool(const kmerInfo &)> & matchesRef){
		const auto found = refMatchCache.find(cacheKey);
		if(refMatchCache.end() != found){
			return found->second;
		}
		bool matches = false;
		for(const auto & refKInfo : ids_.targets_.at(target).refKInfos_){
			if(matchesRef(refKInfo)){
				matches = true;
				break;
			}
		}
		if(refMatchCache.size() >= pars_.stitchPars_.stitchCacheMaxEntries_){
			refMatchCache.clear();
		}
		refMatchCache.emplace(cacheKey, matches);
		return matches;
	};
	//returns true if it matches the references, or there isn't enough to compare, query is set to what was compared
	const auto stitchedMatchesRefs = [&](const std::string & target, const std::string & combined, uint32_t forwardLen, uint32_t reverseLen, std::string & query){
		if(ids_.targets_.at(target).refKInfos_.empty() || combined.size() <= forwardLen + reverseLen + pars_.compKmerLen_){
			return true;
		}
		query = njh::strToUpperRet(combined.substr(forwardLen, combined.size() - forwardLen - reverseLen));
		const kmerInfo seqKInfo(query, pars_.compKmerLen_, false);
		return compareToRefs(target, target + "\t" + query, [&](const kmerInfo & refKInfo){
			return refKInfo.compareKmers(seqKInfo).second >= pars_.compKmerSimCutOff_;
		});
	};
	const auto matesMatchRefs = [&](const std::string & target, const seqInfo & mate1, const seqInfo & mate2, uint32_t forwardLen, uint32_t reverseLen, std::string & query){
		const uint64_t mate1Trim = forwardLen + pars_.stitchPars_.r1Trim_;
		const uint64_t mate2Trim = reverseLen + pars_.stitchPars_.r2Trim_;
		if(ids_.targets_.at(target).refKInfos_.empty() || len(mate1) <= mate1Trim + pars_.compKmerLen_ || len(mate2) <= mate2Trim + pars_.compKmerLen_){
			return true;
		}
		query = njh::strToUpperRet(mate1.seq_.substr(forwardLen, len(mate1) - mate1Trim));
		const auto mate2Query = njh::strToUpperRet(mate2.seq_.substr(reverseLen, len(mate2) - mate2Trim));
		const kmerInfo seqKInfo(query, pars_.compKmerLen_, false);
		const kmerInfo mateKInfo(mate2Query, pars_.compKmerLen_, true);
		return compareToRefs(target, target + "\t" + query + "\t" + mate2Query, [&](const kmerInfo & refKInfo){
			return refKInfo.compareKmers(seqKInfo).second >= pars_.compKmerSimCutOff_ &&
					refKInfo.compareKmersRevComp(mateKInfo).second >= pars_.compKmerSimCutOff_;
		});
	};

	const auto countCategory = [&ret](ReadPairCategory category){
		++ret.categoryCounts_[category];
	};
	const auto addShortProduct = [&ret, &countCategory](ReadPairCategory category, const std::string & seq,
			const std::string & frontPrimer, const std::string & endPrimer, const std::string & gapSeq, uint32_t gap){
		countCategory(category);
		auto & product = ret.shortProducts_[seq];
		product.category_ = category;
		product.frontPrimer_ = frontPrimer;
		product.endPrimer_ = endPrimer;
		product.gapSeq_ = gapSeq;
		product.gap_ = gap;
		++product.count_;
	};
	//returns true if this was a dimer of either kind and has been counted
	const auto detectDimer = [&](seqInfo & combined, const std::string & forwardName, const std::string & reverseName,
			PrimerDimerDetector::DetectorResults & res){
		njh::strToUpper(combined.seq_);
		res = unrecognized == forwardName || unrecognized == reverseName ?
				detector_.detect(combined, dimerAligner) :
				detector_.detect(combined, forwardName + "_F", reverseName + "_R", dimerAligner);
		if(PrimerDimerDetector::ReadCategory::DIMER == res.category_ || PrimerDimerDetector::ReadCategory::UNRESOLVED_DIMER == res.category_){
			countCategory(PrimerDimerDetector::ReadCategory::DIMER == res.category_ ? ReadPairCategory::PRIMER_DIMER : ReadPairCategory::UNRESOLVED_DIMER);
			ret.dimers_.add("all", res, len(combined));
			return true;
		}
		return false;
	};
	const auto addLength = [](LengthTally & tally, const PairedReadProcessor::ProcessedPairRes & pairRes){
		++tally.count_;
		if(pairRes.combinedSeq_ != nullptr){
			tally.lengths_.emplace_back(len(*pairRes.combinedSeq_));
		}
	};

	PairedRead seq;
	SeqInput reader(readOpts);
	reader.openIn();
	while(reader.readNextRead(seq)){
		++ret.total_;
		stitchingAligner.parts_.setMaxSize(std::max(len(seq.seqBase_), len(seq.mateSeqBase_)));
		auto forwardCall = callPrimer(0, seq.seqBase_);
		if(pars_.primerPars_.checkComplement_ && unrecognized == forwardCall.name_){
			//the same as extractorPairedEnd, the whole pair could be in the other direction
			auto mateForwardCall = callPrimer(0, seq.mateSeqBase_);
			if(unrecognized != mateForwardCall.name_){
				std::swap(seq.seqBase_, seq.mateSeqBase_);
				forwardCall = std::move(mateForwardCall);
			}
		}
		const auto reverseCall = callPrimer(1, seq.mateSeqBase_);
		const auto & forwardName = forwardCall.name_;
		const auto & reverseName = reverseCall.name_;
		const bool bothFound = unrecognized != forwardName && unrecognized != reverseName;
		const bool primersMatch = bothFound && forwardName == reverseName;
		if(primersMatch){
			++ret.primersMatched_;
		}else if(bothFound){
			++ret.primersMismatched_;
		}else{
			++ret.primersFailed_;
		}

		const seqInfo firstMate = seq.seqBase_;
		const seqInfo secondMate = seq.mateSeqBase_;
		seq.mateSeqBase_.reverseComplementRead(false, true);
		seq.mateRComplemented_ = true;
		//stitching trims the mates so give it a copy in case there's a second attempt
		auto seqToStitch = seq;
		const auto pairRes = stitcher.processPairedEnd(seqToStitch, stitchCounts, stitchingAligner);
		const auto statusStr = PairedReadProcessor::getOverlapStatusStr(pairRes.status_);
		//the primers are from the same target, count it for the target
		const auto countForTarget = [&](){
			std::string query;
			const bool matchesRefs = pairRes.combinedSeq_ != nullptr ?
					stitchedMatchesRefs(forwardName, pairRes.combinedSeq_->seq_, forwardCall.primerLen_, reverseCall.primerLen_, query) :
					matesMatchRefs(forwardName, firstMate, secondMate, forwardCall.primerLen_, reverseCall.primerLen_, query);
			auto category = ReadPairCategory::KMER_MISMATCH;
			if(!matchesRefs){
				++ret.kmerMismatchSeqs_[std::make_pair(forwardName, query)];
			}else{
				category = expectedStatus(forwardName, pairRes.status_) ? ReadPairCategory::ON_TARGET : ReadPairCategory::UNEXPECTED_STITCH;
			}
			countCategory(category);
			addLength(ret.targetCounts_[std::make_tuple(forwardName, statusStr, category)], pairRes);
		};

		if(OverlapStatus::R1BEGINSINR2 == pairRes.status_){
			seqInfo combined = *pairRes.combinedSeq_;
			PrimerDimerDetector::DetectorResults res;
			const bool couldBeDimer = len(combined) <= detector_.getMaxDimerSize();
			if(couldBeDimer && detectDimer(combined, forwardName, reverseName, res)){
				continue;
			}
			if(primersMatch && expectedStatus(forwardName, pairRes.status_)){
				countForTarget();
				continue;
			}
			if(!couldBeDimer){
				//too long to be a dimer but not what's expected, the primers can be any of them and on either mate so search all of them
				dimerAligner.parts_.setMaxSize(len(combined));
				if(detectDimer(combined, unrecognized, unrecognized, res)){
					continue;
				}
			}
			if(PrimerDimerDetector::ReadCategory::PRIMERS_FAR_APART == res.category_){
				const auto & front = res.frontMatches_.front();
				const auto & end = res.endMatches_.front();
				addShortProduct(ReadPairCategory::SHORT_PRODUCT, combined.seq_, front.primerName_, end.primerName_,
						combined.seq_.substr(front.seqEnd_, end.seqStart_ - front.seqEnd_), front.getGap(end));
			}else if(PrimerDimerDetector::ReadCategory::ONE_PRIMER == res.category_){
				addShortProduct(ReadPairCategory::SHORT_ONE_PRIMER, combined.seq_,
						res.frontMatches_.empty() ? std::string("") : res.frontMatches_.front().primerName_,
						res.endMatches_.empty() ? std::string("") : res.endMatches_.front().primerName_, "", 0);
			}else{
				addShortProduct(ReadPairCategory::SHORT_NO_PRIMER, combined.seq_, "", "", "", 0);
			}
			continue;
		}

		if(OverlapStatus::NOOVERLAP == pairRes.status_){
			PairedReadProcessor::ProcessedResultsCounts retryCounts;
			const auto retryRes = retryStitcher.processPairedEnd(seq, retryCounts, stitchingAligner);
			if(OverlapStatus::R1BEGINSINR2 == retryRes.status_ && len(*retryRes.combinedSeq_) <= detector_.getMaxDimerSize()){
				seqInfo combined = *retryRes.combinedSeq_;
				PrimerDimerDetector::DetectorResults res;
				if(detectDimer(combined, forwardName, reverseName, res)){
					++ret.dimersFromRetryStitch_;
					continue;
				}
			}
		}

		if(primersMatch){
			countForTarget();
		}else{
			const auto category = bothFound ? ReadPairCategory::MISMATCHED_PRIMERS : ReadPairCategory::FAILED_PRIMERS;
			countCategory(category);
			if(bothFound){
				//primers of different targets at each end, where what's between them came from is what's of interest
				const std::string & seqWithPrimers = pairRes.combinedSeq_ != nullptr ? pairRes.combinedSeq_->seq_ : firstMate.seq_;
				const uint64_t endBases = pairRes.combinedSeq_ != nullptr ? reverseCall.primerLen_ : pars_.stitchPars_.r1Trim_;
				if(seqWithPrimers.size() > forwardCall.primerLen_ + endBases){
					++ret.mismatchedPrimerSeqs_[std::make_tuple(forwardName, reverseName,
							njh::strToUpperRet(seqWithPrimers.substr(forwardCall.primerLen_, seqWithPrimers.size() - forwardCall.primerLen_ - endBases)))];
				}
			}
			addLength(ret.otherPrimerCounts_[std::make_tuple(forwardName, reverseName, statusStr, category)], pairRes);
		}
	}
	return ret;
}

AmpliconProductPlacer::PlacerPars::PlacerPars(){
	//what's left of the TruSeq adapters when reads run off of the end of a fragment and the tails that primers commonly carry
	adapters_ = VecStr{
		"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
		"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
		"ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
	};
}

std::string AmpliconProductPlacer::Placement::getLabel(const VecStr & primerTargets) const{
	if("target" != placement_){
		return placement_;
	}
	for(const auto & target : overlappingTargets_){
		if(njh::in(target, primerTargets)){
			return "inPrimerTarget";
		}
	}
	return overlappingTargets_.empty() ? "targetGenomeOffTarget" : "inOtherPanelTarget";
}

AmpliconProductPlacer::AmpliconProductPlacer(PlacerPars pars) : pars_(std::move(pars)){
	const auto addKmers = [this](const std::string & adapter){
		for(size_t pos = 0; pos + pars_.adapterKmerLen_ <= adapter.size(); ++pos){
			adapterKmers_.emplace(adapter.substr(pos, pars_.adapterKmerLen_));
		}
	};
	for(const auto & adapter : pars_.adapters_){
		const auto upperAdapter = njh::strToUpperRet(adapter);
		addKmers(upperAdapter);
		addKmers(seqUtil::reverseComplement(upperAdapter, "DNA"));
	}
}

bool AmpliconProductPlacer::containsAdapter(const std::string & seq) const{
	for(size_t pos = 0; pos + pars_.adapterKmerLen_ <= seq.size(); ++pos){
		if(njh::in(seq.substr(pos, pars_.adapterKmerLen_), adapterKmers_)){
			return true;
		}
	}
	return false;
}

std::unordered_map<std::string, AmpliconProductPlacer::Placement> AmpliconProductPlacer::placeSeqs(
		const std::unordered_map<std::string, uint64_t> & seqCounts, const bfs::path & workingDir) const{
	std::unordered_map<std::string, Placement> ret;
	std::vector<std::pair<std::string, uint64_t>> toBlast;
	for(const auto & seqCount : seqCounts){
		auto & placement = ret[seqCount.first];
		if(containsAdapter(seqCount.first)){
			placement.placement_ = "adapter";
		}else if(seqCount.first.size() < pars_.minBlastLen_){
			placement.placement_ = "tooShort";
		}else{
			toBlast.emplace_back(seqCount);
		}
	}
	if(toBlast.empty() || (pars_.targetGenomes_.empty() && pars_.hostGenomeFnp_.empty())){
		return ret;
	}
	njh::sort(toBlast, [](const std::pair<std::string, uint64_t> & seq1, const std::pair<std::string, uint64_t> & seq2){
		return seq1.second == seq2.second ? seq1.first < seq2.first : seq1.second > seq2.second;
	});
	if(toBlast.size() > pars_.maxSeqsToBlast_){
		toBlast.resize(pars_.maxSeqsToBlast_);
	}

	const auto shortQueryFnp = njh::files::make_path(workingDir, "blastQueriesShort.fasta");
	const auto longQueryFnp = njh::files::make_path(workingDir, "blastQueriesLong.fasta");
	bool anyShort = false;
	bool anyLong = false;
	{
		SeqOutput shortWriter(SeqIOOptions::genFastaOut(shortQueryFnp));
		SeqOutput longWriter(SeqIOOptions::genFastaOut(longQueryFnp));
		shortWriter.openOut();
		longWriter.openOut();
		for(size_t pos = 0; pos < toBlast.size(); ++pos){
			const bool isLong = toBlast[pos].first.size() >= pars_.longSeqLen_;
			anyLong = anyLong || isLong;
			anyShort = anyShort || !isLong;
			(isLong ? longWriter : shortWriter).write(seqInfo(njh::pasteAsStr("seq", pos), toBlast[pos].first));
		}
	}

	struct BestHit{
		double bitScore_{-1};
		std::string genome_;
		BLASTHitTab hit_;
	};
	//blast needs these as plain files but they aren't needed that way once it's done with them
	const auto gzAndRemove = [](const bfs::path & fnp){
		gzZipFile(IoOptions(InOptions(fnp), OutOptions(bfs::path(fnp.string() + ".gz"))));
		bfs::remove(fnp);
	};
	BioCmdsUtils bioRunner;
	njh::sys::requireExternalProgramThrow("blastn");
	//the best hit of each query, index is the position in toBlast
	const auto blastGenome = [&](const std::string & genomeName, const bfs::path & genomeFnp, std::vector<BestHit> & bestHits){
		bioRunner.RunMakeblastdb(genomeFnp);
		const auto blastDb = bfs::path(genomeFnp).replace_extension("").string() + "_blastdb";
		//hits that wouldn't count are dropped as they're found as repeats can otherwise be millions of them. A sequence's hits are only compared to hits from the same search
		const auto runSearch = [&](const bfs::path & queryFnp, const std::string & searchName, const std::string & taskArgs){
			const auto hitsFnp = njh::files::make_path(workingDir, genomeName + "_" + searchName + "_blastHits.tsv");
			bioRunner.runCmdCheck(njh::pasteAsStr("blastn -query ", queryFnp.string(), " -db ", blastDb, " ", taskArgs,
					" -perc_identity ", pars_.minPercentId_, " -qcov_hsp_perc ", 100 * pars_.minCoverage_,
					" -evalue 1 -max_target_seqs 20 -outfmt 6 -num_threads ", pars_.numThreads_, " -out ", hitsFnp.string()),
					queryFnp, hitsFnp);
			{
				BioDataFileIO<BLASTHitTab> reader{IoOptions(InOptions(hitsFnp))};
				reader.openIn();
				BLASTHitTab hit;
				while(reader.readNextRecord(hit)){
					const auto queryPos = njh::StrToNumConverter::stoToNum<size_t>(hit.queryName_.substr(3));
					const double coverage = static_cast<double>(hit.alignLen_) / toBlast[queryPos].first.size();
					if(hit.perId_ >= pars_.minPercentId_ && coverage >= pars_.minCoverage_ && hit.bitScore_ > bestHits[queryPos].bitScore_){
						bestHits[queryPos].bitScore_ = hit.bitScore_;
						bestHits[queryPos].genome_ = genomeName;
						bestHits[queryPos].hit_ = hit;
					}
				}
			}
			gzAndRemove(hitsFnp);
		};
		if(anyShort){
			runSearch(shortQueryFnp, "short", "-task blastn-short -word_size 11 -dust no");
		}
		if(anyLong){
			runSearch(longQueryFnp, "long", "-task megablast -word_size 28 -dust '" + pars_.longSeqDust_ + "'");
		}
	};
	std::vector<BestHit> bestTargetHits(toBlast.size());
	for(const auto & genome : pars_.targetGenomes_){
		blastGenome(genome.name_, genome.genomeFnp_, bestTargetHits);
	}
	std::vector<BestHit> bestHostHits(toBlast.size());
	if(!pars_.hostGenomeFnp_.empty()){
		blastGenome("host", pars_.hostGenomeFnp_, bestHostHits);
	}
	gzAndRemove(shortQueryFnp);
	gzAndRemove(longQueryFnp);

	//the targets are looked up by genome and chromosome
	std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>>> targetLocs;
	for(const auto & genome : pars_.targetGenomes_){
		for(const auto & loc : getBeds(genome.targetBedFnp_)){
			//the mates of targets that don't overlap also get their own lines
			if(!njh::endsWith(loc->name_, "-R1") && !njh::endsWith(loc->name_, "-R2")){
				targetLocs[genome.name_][loc->chrom_].emplace_back(loc);
			}
		}
	}
	for(size_t pos = 0; pos < toBlast.size(); ++pos){
		auto & placement = ret[toBlast[pos].first];
		const auto & targetHit = bestTargetHits[pos];
		const auto & hostHit = bestHostHits[pos];
		if(targetHit.bitScore_ < 0 && hostHit.bitScore_ < 0){
			placement.placement_ = "unplaced";
			continue;
		}
		//ties go to the target genomes
		const bool inTarget = targetHit.bitScore_ >= hostHit.bitScore_;
		const auto & best = inTarget ? targetHit : hostHit;
		placement.placement_ = inTarget ? "target" : "host";
		placement.genome_ = best.genome_;
		placement.chrom_ = best.hit_.subjectName_;
		placement.start_ = std::min(best.hit_.sStart_, best.hit_.sEnd_) - 1;
		placement.end_ = std::max(best.hit_.sStart_, best.hit_.sEnd_);
		placement.percentId_ = best.hit_.perId_;
		placement.coverage_ = static_cast<double>(best.hit_.alignLen_) / toBlast[pos].first.size();
		if(inTarget){
			for(const auto & loc : targetLocs[best.genome_][placement.chrom_]){
				const uint64_t regionStart = loc->chromStart_ > pars_.targetRegionFlank_ ? loc->chromStart_ - pars_.targetRegionFlank_ : 0;
				if(placement.start_ < loc->chromEnd_ + pars_.targetRegionFlank_ && placement.end_ > regionStart && !njh::in(loc->name_, placement.overlappingTargets_)){
					placement.overlappingTargets_.emplace_back(loc->name_);
				}
			}
		}
	}
	return ret;
}

}  // namespace njhseq
