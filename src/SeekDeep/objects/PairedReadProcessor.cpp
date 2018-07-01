/*
 * PairedReadProcessor.cpp
 *
 *  Created on: Jan 14, 2018
 *      Author: nick
 */
//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
// Jeffrey Bailey <Jeffrey.Bailey@umassmed.edu>
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

#include "PairedReadProcessor.hpp"

namespace bibseq {


PairedReadProcessor::PairedReadProcessor(ProcessParams params):params_(params){
	setDefaultConsensusBuilderFunc();

}

void PairedReadProcessor::setDefaultConsensusBuilderFunc(){
	addToConsensus = [](
			uint32_t pos,
			const seqInfo & r1,
			const seqInfo & r2,
			std::string & cseq,
			std::vector<uint32_t> & cquals,
			aligner & alignerObj){
		if('-' == r1.seq_[pos]){
			//gap in r1
			if(r2.qual_[pos] >= alignerObj.qScorePars_.primaryQual_){

				cseq.push_back(r2.seq_[pos]);
				cquals.push_back(r2.qual_[pos]);
			}
		}else if('-' == r2.seq_[pos]){
			//gap in r2
			if(r1.qual_[pos] >= alignerObj.qScorePars_.primaryQual_){
				cseq.push_back(r1.seq_[pos]);
				cquals.push_back(r1.qual_[pos]);
			}
		}else if(alignerObj.parts_.scoring_.mat_[r1.seq_[pos]][r2.seq_[pos]] > 0){
			//match (use the higher quality score)
			cseq.push_back(r1.seq_[pos]);
			if(islower(r1.seq_[pos]) || islower(r2.seq_[pos])){
				cseq.back() = tolower(cseq.back());
			}
			cquals.push_back(r1.qual_[pos] >= r2.qual_[pos] ? r1.qual_[pos] : r2.qual_[pos]);
		}else{
			//mismatch (take the higher quality giving preference to r1)
			cseq.push_back(r1.qual_[pos] >= r2.qual_[pos] ? r1.seq_[pos] : r2.seq_[pos]);
			cquals.push_back(r1.qual_[pos] >= r2.qual_[pos] ? r1.qual_[pos] : r2.qual_[pos]);
			if(islower(r1.seq_[pos]) || islower(r2.seq_[pos])){
				cseq.back() = tolower(cseq.back());
			}
		}
	};
}

uint64_t PairedReadProcessor::guessMaxReadLenFromFile(const SeqIOOptions & inputOpts){
	PairedRead seq;
	SeqInput reader(inputOpts);
	reader.openIn();
	uint64_t maxsize = 0;
	uint32_t seqCount = 0;
	while(reader.readNextRead(seq)){
		readVec::getMaxLength(seq.seqBase_, maxsize);
		readVec::getMaxLength(seq.mateSeqBase_, maxsize);
		++seqCount;
		if(seqCount > params_.checkAmount_){
			break;
		}
	}
	return maxsize;
}


PairedReadProcessor::ProcessorOutWriters::ProcessorOutWriters(){
};
PairedReadProcessor::ProcessorOutWriters::ProcessorOutWriters(const OutOptions & outOpts){
	auto perfectOverlapCombinedOpts    = SeqIOOptions::genFastqOut (outOpts.outFilename_.string() + "_perfectOverlap.fastq");
	auto r1EndsInR2CombinedOpts    = SeqIOOptions::genFastqOut (outOpts.outFilename_.string() + "_r1EndsInR2.fastq");
	auto r1BeginsInR2CombinedOpts    = SeqIOOptions::genFastqOut (outOpts.outFilename_.string() + "_r1BegingsInR2.fastq");
	auto notCombinedOpts = SeqIOOptions::genPairedOut(outOpts.outFilename_.string() + "_notCombined");
	auto overhangsOpts =   SeqIOOptions::genPairedOut(outOpts.outFilename_.string() + "_overhangs");
	perfectOverlapCombinedOpts.out_.transferOverwriteOpts(outOpts);
	r1EndsInR2CombinedOpts.out_.transferOverwriteOpts(outOpts);
	r1BeginsInR2CombinedOpts.out_.transferOverwriteOpts(outOpts);
	notCombinedOpts.out_.transferOverwriteOpts(outOpts);
	overhangsOpts.out_.transferOverwriteOpts(outOpts);
	perfectOverlapCombinedWriter = std::make_unique<SeqOutput>(perfectOverlapCombinedOpts);
	r1EndsInR2CombinedWriter= std::make_unique<SeqOutput>(r1EndsInR2CombinedOpts);
	r1BeginsInR2CombinedWriter= std::make_unique<SeqOutput>(r1BeginsInR2CombinedOpts);
	notCombinedWriter= std::make_unique<SeqOutput>(notCombinedOpts);
	overhangsWriter= std::make_unique<SeqOutput>(overhangsOpts);
}



void PairedReadProcessor::ProcessorOutWriters::checkWritersSet(const std::string & funcName){
	bool failed = false;
	std::stringstream ss;
	ss << funcName << ", error the following writers aren't set " << "\n";
	auto checkWriter = [&failed,&ss](const std::unique_ptr<SeqOutput> & writer, const std::string & writerName){
		if(nullptr == writer){
			failed = true;
			ss << writerName << " is not set" << "\n";
		}
	};
	checkWriter(perfectOverlapCombinedWriter, "perfectOverlapCombinedWriter");
	checkWriter(r1EndsInR2CombinedWriter, "r1EndsInR2CombinedWriter");
	checkWriter(r1BeginsInR2CombinedWriter, "r1BegingsInR2CombinedWriter");
	checkWriter(notCombinedWriter, "notCombinedWriter");
	checkWriter(overhangsWriter, "overhangsWriter");

	if(failed){
		throw std::runtime_error{ss.str()};
	}
}

void PairedReadProcessor::ProcessorOutWriters::unsetWriters(){
	perfectOverlapCombinedWriter = nullptr;//(perfectOverlapCombinedOpts);
	r1EndsInR2CombinedWriter = nullptr;//(r1EndsInR2CombinedOpts);
	r1BeginsInR2CombinedWriter = nullptr;//(r1BeginsInR2CombinedOpts);
	notCombinedWriter = nullptr;//(notCombinedOpts);
	overhangsWriter = nullptr;//(overhangsOpts);
}


Json::Value PairedReadProcessor::ProcessedResults::toJson() const{
	Json::Value outVal;
	outVal["overlapFail"] = overlapFail;
	outVal["overlapFailPerc"] = (overlapFail/static_cast<double>(total)) * 100;
	outVal["overhangFail"] = overhangFail;
	outVal["overhangFailPerc"] = (overhangFail/static_cast<double>(total)) * 100;
	outVal["perfectOverlapCombined"] = perfectOverlapCombined;
	outVal["perfectOverlapCombinedPerc"] = (perfectOverlapCombined/static_cast<double>(total)) * 100;
	outVal["r1EndsInR2Combined"] = r1EndsInR2Combined;
	outVal["r1EndsInR2CombinedPerc"] = (r1EndsInR2Combined/static_cast<double>(total)) *100;
	outVal["r1BeginsInR2Combined"] = r1BeginsInR2Combined;
	outVal["r1BeginsInR2CombinedPerc"] = (r1BeginsInR2Combined/static_cast<double>(total)) *100;
	outVal["total"] = total;
	if(nullptr != perfectOverlapCombinedOpts){
		outVal["perfectOverlapCombinedOpts"] = bib::json::toJson(perfectOverlapCombinedOpts);
	}
	if(nullptr != r1EndsInR2CombinedOpts){
		outVal["r1EndsInR2CombinedOpts"] = bib::json::toJson(r1EndsInR2CombinedOpts);
	}
	if(nullptr != r1BeginsInR2CombinedOpts){
		outVal["r1BeginsInR2CombinedOpts"] = bib::json::toJson(r1BeginsInR2CombinedOpts);
	}
	if(nullptr != notCombinedOpts){
		outVal["notCombinedOpts"] = bib::json::toJson(notCombinedOpts);
	}
	if(nullptr != overhangsOpts){
		outVal["overhangsOpts"] = bib::json::toJson(overhangsOpts);
	}
	return outVal;
}

Json::Value PairedReadProcessor::ProcessedResults::toJsonCounts() const{
	Json::Value outVal;
	outVal["overlapFail"] = overlapFail;
	outVal["overlapFailPerc"] = (overlapFail/static_cast<double>(total)) * 100;
	outVal["overhangFail"] = overhangFail;
	outVal["overhangFailPerc"] = (overhangFail/static_cast<double>(total)) * 100;
	outVal["perfectOverlapCombined"] = perfectOverlapCombined;
	outVal["perfectOverlapCombinedPerc"] = (perfectOverlapCombined/static_cast<double>(total)) * 100;
	outVal["r1EndsInR2Combined"] = r1EndsInR2Combined;
	outVal["r1EndsInR2CombinedPerc"] = (r1EndsInR2Combined/static_cast<double>(total)) *100;
	outVal["r1BeginsInR2Combined"] = r1BeginsInR2Combined;
	outVal["r1BeginsInR2CombinedPerc"] = (r1BeginsInR2Combined/static_cast<double>(total)) *100;
	outVal["total"] = total;
	return outVal;
}

PairedReadProcessor::ProcessedResults PairedReadProcessor::processPairedEnd(
		SeqInput & reader,
		ProcessorOutWriters & writers,
		aligner & alignerObj){
	writers.checkWritersSet(__PRETTY_FUNCTION__);
	ProcessedResults res;
	if(!reader.inOpen()){
		reader.openIn();
	}

	PairedRead seq;

	while(processPairedEnd(reader, seq, writers, alignerObj, res)){
		if(res.total >= params_.testNumber_){
			break;
		}
	}

	if(writers.perfectOverlapCombinedWriter->outOpen()){
		res.perfectOverlapCombinedOpts = std::make_shared<SeqIOOptions>(SeqIOOptions::genFastqIn(writers.perfectOverlapCombinedWriter->getPrimaryOutFnp()));
	}
	if(writers.r1EndsInR2CombinedWriter->outOpen()){
		res.r1EndsInR2CombinedOpts = std::make_shared<SeqIOOptions>(SeqIOOptions::genFastqIn(writers.r1EndsInR2CombinedWriter->getPrimaryOutFnp()));
	}
	if(writers.r1BeginsInR2CombinedWriter->outOpen()){
		res.r1BeginsInR2CombinedOpts = std::make_shared<SeqIOOptions>(SeqIOOptions::genFastqIn(writers.r1BeginsInR2CombinedWriter->getPrimaryOutFnp()));
	}
	if(writers.notCombinedWriter->outOpen()){
		res.notCombinedOpts = std::make_shared<SeqIOOptions>(SeqIOOptions::genPairedIn(
				writers.notCombinedWriter->getPrimaryOutFnp(),
				writers.notCombinedWriter->getSecondaryOutFnp()));
	}
	if(writers.overhangsWriter->outOpen()){
		res.overhangsOpts = std::make_shared<SeqIOOptions>(SeqIOOptions::genPairedIn(
				writers.overhangsWriter->getPrimaryOutFnp(),
				writers.overhangsWriter->getSecondaryOutFnp()));
	}
	return res;
}

bool PairedReadProcessor::processPairedEnd(
		SeqInput & reader,
		PairedRead & seq,
		ProcessorOutWriters & writers,
		aligner & alignerObj,
		ProcessedResults & res){

//	auto seqOutTemp = SeqIOOptions::genFastqOut("test.fastq");
//	seqOutTemp.out_.append_ = true;
//	SeqOutput tempWriter(seqOutTemp);
//	tempWriter.openOut();

	if(reader.readNextRead(seq)){
		double percentId = 1 - params_.errorAllowed_ ;
		++res.total;
//		std::cout << '\r' << res.total;
//		std::cout.flush();
		if(res.total % 25000 == 0 && params_.verbose_){
			std::cout << res.total << std::endl;
		}
		//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ <<  std::endl;
		alignerObj.alignRegGlobalNoInternalGaps(seq.seqBase_, seq.mateSeqBase_);
		//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ <<  std::endl;
		alignerObj.profileAlignment(seq.seqBase_, seq.mateSeqBase_, false, true, true);
		//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ <<  std::endl;

		if( alignerObj.comp_.distances_.eventBasedIdentityHq_ >= percentId &&
				alignerObj.comp_.distances_.basesInAln_ >= params_.minOverlap_ &&
				alignerObj.comp_.hqMismatches_ + alignerObj.comp_.lqMismatches_ <= params_.hardMismatchCutOff_){
			AlignOverlapEnd frontCase = AlignOverlapEnd::UNHANDLEED;
			AlignOverlapEnd backCase  = AlignOverlapEnd::UNHANDLEED;
			if( '-' != alignerObj.alignObjectA_.seqBase_.seq_.front() &&
					'-' != alignerObj.alignObjectB_.seqBase_.seq_.front()){
				frontCase = AlignOverlapEnd::NOOVERHANG;
			}else if('-' != alignerObj.alignObjectA_.seqBase_.seq_.front() &&
							 '-' == alignerObj.alignObjectB_.seqBase_.seq_.front()){
				frontCase = AlignOverlapEnd::R1OVERHANG;
			}else if('-' == alignerObj.alignObjectA_.seqBase_.seq_.front() &&
							 '-' != alignerObj.alignObjectB_.seqBase_.seq_.front()){
				frontCase = AlignOverlapEnd::R2OVERHANG;
			}else{
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error not handled case for " << seq.seqBase_.name_ << "\n";
				ss << "R1.front(): " << alignerObj.alignObjectA_.seqBase_.seq_.front() << ", R2.front(): " << alignerObj.alignObjectB_.seqBase_.seq_.front();
				throw std::runtime_error{ss.str()};
			}
			if( '-' != alignerObj.alignObjectA_.seqBase_.seq_.back() &&
					'-' != alignerObj.alignObjectB_.seqBase_.seq_.back()){
				backCase = AlignOverlapEnd::NOOVERHANG;
			}else if('-' != alignerObj.alignObjectA_.seqBase_.seq_.back() &&
							 '-' == alignerObj.alignObjectB_.seqBase_.seq_.back()){
				backCase = AlignOverlapEnd::R1OVERHANG;
			}else if('-' == alignerObj.alignObjectA_.seqBase_.seq_.back() &&
							 '-' != alignerObj.alignObjectB_.seqBase_.seq_.back()){
				backCase = AlignOverlapEnd::R2OVERHANG;
			}else{
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error not handled case for " << seq.seqBase_.name_ << "\n";
				ss << "R1.front(): " << alignerObj.alignObjectA_.seqBase_.seq_.front() << ", R2.front(): " << alignerObj.alignObjectB_.seqBase_.seq_.front();
				throw std::runtime_error{ss.str()};
			}

			if(AlignOverlapEnd::NOOVERHANG == frontCase && AlignOverlapEnd::NOOVERHANG == backCase){
				//no over hangs, perfect overlap
				std::string cseq;
				cseq.reserve(alignerObj.comp_.distances_.basesInAln_);
				std::vector<uint32_t> quals;
				quals.reserve(alignerObj.comp_.distances_.basesInAln_);
				for(const auto pos : iter::range(len(alignerObj.alignObjectA_))){
					addToConsensus(pos,
							alignerObj.alignObjectA_.seqBase_,
							alignerObj.alignObjectB_.seqBase_,
							cseq,
							quals,
							alignerObj);
				}

				seqInfo combinedSeq(seq.seqBase_.name_, cseq, quals);
				writers.perfectOverlapCombinedWriter->openWrite(combinedSeq);
				++res.perfectOverlapCombined;
			}else if((AlignOverlapEnd::NOOVERHANG == frontCase || AlignOverlapEnd::R1OVERHANG == frontCase) &&
							 (AlignOverlapEnd::NOOVERHANG == backCase  || AlignOverlapEnd::R2OVERHANG == backCase)){
				//ideal situation, R1 end overlaps R2 beg
				std::string cseq;
				cseq.reserve(alignerObj.alignObjectA_.seqBase_.seq_.size());
				std::vector<uint32_t> quals;
				quals.reserve(alignerObj.alignObjectA_.seqBase_.seq_.size());
				uint32_t r1End = alignerObj.alignObjectA_.seqBase_.seq_.find_last_not_of('-') + 1;
				uint32_t r2Start = alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-');
				//add r1 beginning
				if(0 != r2Start ){
					cseq.append(alignerObj.alignObjectA_.seqBase_.seq_.substr(0, r2Start));
					addOtherVec(quals, getSubVector(alignerObj.alignObjectA_.seqBase_.qual_, 0,r2Start));
				}
				//get consensus of middle
				for(const auto pos : iter::range(r2Start, r1End)){
					addToConsensus(pos,
							alignerObj.alignObjectA_.seqBase_,
							alignerObj.alignObjectB_.seqBase_,
							cseq,
							quals,
							alignerObj);
				}
				//add r2 ending
				if(alignerObj.alignObjectA_.seqBase_.seq_.size() != r1End){
					cseq.append(alignerObj.alignObjectB_.seqBase_.seq_.substr(r1End));
					addOtherVec(quals, getSubVector(alignerObj.alignObjectB_.seqBase_.qual_, r1End));
				}
				seqInfo combinedSeq(seq.seqBase_.name_, cseq, quals);
				writers.r1EndsInR2CombinedWriter->openWrite(combinedSeq);
				++res.r1EndsInR2Combined;

			}else if((AlignOverlapEnd::NOOVERHANG == frontCase || AlignOverlapEnd::R2OVERHANG == frontCase) &&
							 (AlignOverlapEnd::NOOVERHANG == backCase  || AlignOverlapEnd::R1OVERHANG == backCase)){
				//read through situation, R2 end overlaps R1 beg, overhang is likely illumina adaptor/primer
				uint32_t r1Start =
						alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-');
				uint32_t r2End =
						alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-') + 1;
				//write out overhangs
				seqInfo back = alignerObj.alignObjectB_.seqBase_.getSubRead(0, r1Start);
				seqInfo front = alignerObj.alignObjectA_.seqBase_.getSubRead(r2End);
				back.reverseComplementRead(false, true);
				std::shared_ptr<PairedRead> overhang;
				if(std::string::npos != back.name_.find("_Comp")){
					overhang = std::make_shared<PairedRead>(back, front);
				}else{
					overhang = std::make_shared<PairedRead>(front, back);
				}
				if(params_.writeOverHangs_){
					writers.overhangsWriter->openWrite(overhang);

				}
//				tempWriter.write(alignerObj.alignObjectA_);
//				tempWriter.write(alignerObj.alignObjectB_);
				std::string cseq;
				cseq.reserve(alignerObj.comp_.distances_.basesInAln_);
				std::vector<uint32_t> quals;
				quals.reserve(alignerObj.comp_.distances_.basesInAln_);
				//get consensus of middle
				for(const auto pos : iter::range(r1Start, r2End)){
					addToConsensus(pos,
							alignerObj.alignObjectA_.seqBase_,
							alignerObj.alignObjectB_.seqBase_,
							cseq,
							quals,
							alignerObj);
				}
				seqInfo combinedSeq(seq.seqBase_.name_, cseq, quals);
				writers.r1BeginsInR2CombinedWriter->openWrite(combinedSeq);
				++res.r1BeginsInR2Combined;
			}else{
				//failure
				writers.notCombinedWriter->openWrite(seq);
				++res.overhangFail;
			}
		}else{
			++res.overlapFail;
			writers.notCombinedWriter->openWrite(seq);
		}
		return true;
	}else{
		return false;
	}
}


}  // namespace bibseq
