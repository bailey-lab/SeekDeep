/*
 * NanoporeReadPlacer.cpp
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

#include "NanoporeReadPlacer.hpp"
#include <cstdio>

namespace njhseq {

namespace {

std::string shellQuote(const std::string & str){
	std::string ret = "'";
	for(const auto c : str){
		if('\'' == c){
			ret += "'\\''";
		}else{
			ret.push_back(c);
		}
	}
	ret.push_back('\'');
	return ret;
}

int64_t overlap(int64_t start1, int64_t end1, int64_t start2, int64_t end2){
	return std::max<int64_t>(0, std::min(end1, end2) - std::max(start1, start2));
}

/**@brief A line of a PAF file, returns false for a secondary alignment
 */
bool parsePafLine(char * line, std::string & name, NanoporeReadPlacer::Alignment & aln){
	char * fields[12];
	char * pos = line;
	for(uint32_t field = 0; field < 12; ++field){
		fields[field] = pos;
		pos = std::strchr(pos, '\t');
		if(nullptr == pos){
			if(11 != field){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error too few columns in PAF line" << "\n";
				throw std::runtime_error{ss.str()};
			}
			break;
		}
		*pos = '\0';
		++pos;
	}
	const char * tags = pos;
	name.assign(fields[0], std::strcspn(fields[0], " \t["));
	aln.qStart_ = std::strtoul(fields[2], nullptr, 10);
	aln.qEnd_ = std::strtoul(fields[3], nullptr, 10);
	aln.strand_ = fields[4][0];
	aln.chrom_.assign(fields[5]);
	aln.rStart_ = std::strtoull(fields[7], nullptr, 10);
	aln.rEnd_ = std::strtoull(fields[8], nullptr, 10);
	aln.matches_ = std::strtoul(fields[9], nullptr, 10);
	aln.mapq_ = std::strtoul(fields[11], nullptr, 10);
	aln.cigar_.clear();
	if(nullptr != tags){
		if(nullptr != std::strstr(tags, "tp:A:S")){
			return false;
		}
		const char * cigar = std::strstr(tags, "cg:Z:");
		if(nullptr != cigar){
			cigar += 5;
			while(std::isdigit(static_cast<unsigned char>(*cigar))){
				char * after = nullptr;
				const auto opLen = std::strtoul(cigar, &after, 10);
				aln.cigar_.emplace_back(static_cast<uint32_t>(opLen), *after);
				cigar = after + 1;
			}
		}
	}
	return true;
}

/**@brief The reference position at a read position, walking the alignment along the read
 */
int64_t refAt(const NanoporeReadPlacer::Alignment & aln, int64_t queryPos){
	queryPos = std::min<int64_t>(std::max<int64_t>(queryPos, aln.qStart_), aln.qEnd_);
	const bool forward = '+' == aln.strand_;
	int64_t query = aln.qStart_;
	int64_t ref = forward ? aln.rStart_ : aln.rEnd_;
	const auto step = [&](const std::pair<uint32_t, char> & op, bool & done){
		if('M' == op.second || '=' == op.second || 'X' == op.second){
			if(query + op.first >= queryPos){
				const auto dist = queryPos - query;
				ref = forward ? ref + dist : ref - dist;
				done = true;
				return;
			}
			query += op.first;
			ref = forward ? ref + op.first : ref - op.first;
		}else if('I' == op.second){
			if(query + op.first >= queryPos){
				done = true;
				return;
			}
			query += op.first;
		}else if('D' == op.second || 'N' == op.second){
			ref = forward ? ref + op.first : ref - op.first;
		}
	};
	bool done = false;
	//for the reverse strand the cigar runs along the reference, so backwards along the read
	if(forward){
		for(auto op = aln.cigar_.begin(); op != aln.cigar_.end() && !done; ++op){
			step(*op, done);
		}
	}else{
		for(auto op = aln.cigar_.rbegin(); op != aln.cigar_.rend() && !done; ++op){
			step(*op, done);
		}
	}
	return ref;
}

/**@brief The alignment covering the most of a stretch of a read
 */
const NanoporeReadPlacer::Alignment * bestForArm(const std::vector<NanoporeReadPlacer::Alignment> & alns, int64_t armStart, int64_t armEnd, int64_t & covered){
	const NanoporeReadPlacer::Alignment * best = nullptr;
	covered = 0;
	for(const auto & aln : alns){
		const auto cov = overlap(aln.qStart_, aln.qEnd_, armStart, armEnd);
		if(cov > covered){
			best = &aln;
			covered = cov;
		}
	}
	return best;
}

}  // namespace

NanoporeReadPlacer::NanoporeReadPlacer(Pars pars) : pars_(std::move(pars)){
	njh::sys::requireExternalProgramThrow("minimap2");
	InputStream bed{InOptions{pars_.targetBedFnp_}};
	std::string line;
	while(njh::files::crossPlatGetline(bed, line)){
		if(line.empty() || '#' == line[0]){
			continue;
		}
		const auto toks = tokenizeString(line, "\t");
		if(toks.size() < 4){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << pars_.targetBedFnp_ << " needs at least 4 columns, line: " << line << "\n";
			throw std::runtime_error{ss.str()};
		}
		panel_[toks[0]].emplace_back(PanelTarget{njh::StrToNumConverter::stoToNum<int64_t>(toks[1]), njh::StrToNumConverter::stoToNum<int64_t>(toks[2]), toks[3]});
	}
	InputStream genome{InOptions{pars_.targetGenomeFnp_}};
	while(njh::files::crossPlatGetline(genome, line)){
		if(!line.empty() && '>' == line[0]){
			targetContigs_.emplace(line.substr(1, line.find_first_of(" \t") - 1));
		}
	}
	if(targetContigs_.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error found no contigs in " << pars_.targetGenomeFnp_ << "\n";
		throw std::runtime_error{ss.str()};
	}
}

const VecStr & NanoporeReadPlacer::allPlacements(){
	static const VecStr placements{"panelTarget", "otherTargetGenome", "host", "unplaced"};
	return placements;
}

void NanoporeReadPlacer::mapSeqs(const std::vector<bfs::path> & seqsFnps, const bfs::path & logFnp, bool forSegments,
		const std::function<void(const std::string &, const std::vector<Alignment> &)> & onSeq) const{
	std::stringstream cmd;
	//-c for the real number of matching bases and the cigar, a fold's two sides are short with most of their sequence being all that's there to map
	cmd << "minimap2 -x map-ont -c --secondary=no"
			<< (forSegments ? " -m 20 -s 40" : "")
			<< " -t " << std::max<uint32_t>(pars_.numThreads_, 1)
			<< " " << shellQuote(pars_.mappingIndex_.string());
	for(const auto & seqsFnp : seqsFnps){
		cmd << " " << shellQuote(seqsFnp.string());
	}
	cmd << " 2> " << shellQuote(logFnp.string());
	FILE * pipe = popen(cmd.str().c_str(), "r");
	if(nullptr == pipe){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error couldn't run " << cmd.str() << "\n";
		throw std::runtime_error{ss.str()};
	}
	char * buffer = nullptr;
	size_t capacity = 0;
	ssize_t lineLen = 0;
	std::string currentName;
	std::string name;
	std::vector<Alignment> alns;
	Alignment aln;
	while((lineLen = getline(&buffer, &capacity, pipe)) != -1){
		if(lineLen > 0 && '\n' == buffer[lineLen - 1]){
			buffer[lineLen - 1] = '\0';
		}
		const bool use = parsePafLine(buffer, name, aln);
		if(name != currentName){
			if(!currentName.empty()){
				onSeq(currentName, alns);
			}
			currentName = name;
			alns.clear();
		}
		if(use){
			alns.emplace_back(aln);
		}
	}
	if(!currentName.empty()){
		onSeq(currentName, alns);
	}
	free(buffer);
	const int status = pclose(pipe);
	if(0 != status){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << cmd.str() << " exited with " << status << ", check " << logFnp << "\n";
		throw std::runtime_error{ss.str()};
	}
}

NanoporeReadPlacer::ReadPlacement NanoporeReadPlacer::placeRead(const std::vector<Alignment> & alns) const{
	ReadPlacement ret;
	int64_t panelBases = 0;
	int64_t otherBases = 0;
	int64_t hostBases = 0;
	struct Locus{
		std::string chrom;
		int64_t start;
		int64_t end;
	};
	std::vector<Locus> loci;
	for(const auto & aln : alns){
		if(aln.matches_ < pars_.minAlnMatches_){
			continue;
		}
		const int64_t queryBases = aln.qEnd_ - aln.qStart_;
		const int64_t start = aln.rStart_;
		const int64_t end = aln.rEnd_;
		if(njh::in(aln.chrom_, targetContigs_)){
			bool inPanel = false;
			const auto targets = panel_.find(aln.chrom_);
			if(targets != panel_.end()){
				for(const auto & target : targets->second){
					const auto bases = overlap(start, end, target.start_, target.end_);
					if(bases >= pars_.minPanelOverlap_ * (end - start)){
						inPanel = true;
					}
					if(bases >= pars_.minWholeUnitFrac_ * (target.end_ - target.start_) && bases >= pars_.minWholeUnitFrac_ * (end - start)){
						++ret.wholeUnits_;
						break;
					}
				}
			}
			(inPanel ? panelBases : otherBases) += queryBases;
		}else{
			hostBases += queryBases;
		}
		if(queryBases >= pars_.minLocusAlnLen_){
			bool found = false;
			for(auto & locus : loci){
				if(locus.chrom == aln.chrom_ && start <= locus.end + pars_.maxLocusGap_ && end + pars_.maxLocusGap_ >= locus.start){
					locus.start = std::min(locus.start, start);
					locus.end = std::max(locus.end, end);
					found = true;
					break;
				}
			}
			if(!found){
				loci.emplace_back(Locus{aln.chrom_, start, end});
			}
		}
	}
	ret.twoLoci_ = loci.size() >= 2;
	if(panelBases + otherBases + hostBases > 0){
		if(panelBases >= otherBases && panelBases >= hostBases){
			ret.placement_ = "panelTarget";
		}else{
			ret.placement_ = otherBases >= hostBases ? "otherTargetGenome" : "host";
		}
	}
	return ret;
}

NanoporeReadPlacer::FoldPlacement NanoporeReadPlacer::placeFold(const FoldInput & input, const std::vector<Alignment> & leftAlns, const std::vector<Alignment> & rightAlns) const{
	FoldPlacement ret;
	const auto & fold = input.fold_;
	const int64_t leftEnd = fold.leftEnd_;
	const int64_t rightStart = fold.rightStart_;
	int64_t leftCovered = 0;
	int64_t rightCovered = 0;
	const auto * left = bestForArm(leftAlns, fold.leftStart_, leftEnd, leftCovered);
	//the right side starts at leftEnd
	const auto * right = bestForArm(rightAlns, rightStart - leftEnd, static_cast<int64_t>(fold.rightEnd_) - leftEnd, rightCovered);
	const bool leftOk = nullptr != left && leftCovered >= 0.5 * fold.leftArmLen();
	const bool rightOk = nullptr != right && rightCovered >= 0.5 * fold.rightArmLen();
	if(!leftOk && !rightOk){
		ret.status_ = "unplaced";
		return ret;
	}
	const auto * one = leftOk ? left : right;
	ret.organism_ = njh::in(one->chrom_, targetContigs_) ? "target" : "host";
	ret.chrom_ = one->chrom_;
	if(!leftOk || !rightOk){
		ret.status_ = "oneArmPlaced";
		return ret;
	}
	if(left->chrom_ != right->chrom_){
		ret.status_ = "armsOnDifferentContigs";
		return ret;
	}
	if(left->strand_ == right->strand_){
		ret.status_ = "armsSameStrand";
		return ret;
	}
	const int64_t leftRStart = left->rStart_;
	const int64_t leftREnd = left->rEnd_;
	const int64_t rightRStart = right->rStart_;
	const int64_t rightREnd = right->rEnd_;
	if(overlap(leftRStart, leftREnd, rightRStart, rightREnd) < 0.5 * std::min(leftREnd - leftRStart, rightREnd - rightRStart)){
		ret.status_ = "armsDontOverlap";
		return ret;
	}
	ret.status_ = "placed";
	ret.placed_ = true;
	const bool forward = '+' == left->strand_;
	ret.extendStrand_ = left->strand_;
	const int64_t edgeLeft = forward ? leftREnd : leftRStart;
	const int64_t edgeRight = forward ? rightREnd : rightRStart;
	ret.far_ = forward ? std::max(edgeLeft, edgeRight) : std::min(edgeLeft, edgeRight);
	ret.entryLeft_ = refAt(*left, leftEnd);
	ret.entryRight_ = refAt(*right, rightStart - leftEnd);
	ret.queryGap_ = static_cast<int64_t>(right->qStart_) + leftEnd - static_cast<int64_t>(left->qEnd_);
	ret.entry_ = (ret.entryLeft_ + ret.entryRight_) / 2;
	ret.reach_ = forward ? ret.far_ - ret.entry_ : ret.entry_ - ret.far_;
	if(edgeLeft == edgeRight){
		ret.loopSide_ = "neither";
	}else{
		ret.loopSide_ = std::abs(edgeLeft - ret.far_) <= std::abs(edgeRight - ret.far_) ? "left" : "right";
	}
	ret.spanLow_ = std::min(leftRStart, rightRStart);
	ret.spanHigh_ = std::max(leftREnd, rightREnd);
	ret.mapqLeft_ = left->mapq_;
	ret.mapqRight_ = right->mapq_;
	const PanelTarget * bestTarget = nullptr;
	int64_t bestOverlap = 0;
	const auto targets = panel_.find(left->chrom_);
	if(targets != panel_.end()){
		for(const auto & target : targets->second){
			const auto bases = overlap(ret.spanLow_, ret.spanHigh_, target.start_, target.end_);
			if(bases > bestOverlap){
				bestTarget = &target;
				bestOverlap = bases;
			}
		}
	}
	if(nullptr != bestTarget && bestOverlap >= 0.5 * (ret.spanHigh_ - ret.spanLow_)){
		const int64_t targetEdge = forward ? bestTarget->end_ : bestTarget->start_;
		ret.target_ = bestTarget->name_;
		ret.targetStart_ = bestTarget->start_;
		ret.targetEnd_ = bestTarget->end_;
		ret.farToEnd_ = forward ? targetEdge - ret.far_ : ret.far_ - targetEdge;
		ret.entryToEnd_ = forward ? targetEdge - ret.entry_ : ret.entry_ - targetEdge;
		if(ret.queryGap_ > 10){
			ret.turnClass_ = "loopNotPlaced";
		}else if(std::abs(ret.farToEnd_) <= 15){
			ret.turnClass_ = ret.reach_ <= 10 ? "terminalJoin" : "reachedEnd";
		}else{
			ret.turnClass_ = "internal";
		}
	}else{
		ret.turnClass_ = ret.queryGap_ > 10 ? "loopNotPlaced" : "offPanel";
	}
	return ret;
}

}  // namespace njhseq
