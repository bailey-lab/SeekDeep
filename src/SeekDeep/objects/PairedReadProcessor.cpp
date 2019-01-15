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

namespace njhseq {


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
	auto r1AllInR2CombinedOpts    = SeqIOOptions::genFastqOut (outOpts.outFilename_.string() + "_r1AllInR2.fastq");
	auto r2AllInR1ombinedOpts    = SeqIOOptions::genFastqOut (outOpts.outFilename_.string() + "_r2AllInR1.fastq");
	auto notCombinedOpts = SeqIOOptions::genPairedOut(outOpts.outFilename_.string() + "_notCombined");
	auto overhangsOpts =   SeqIOOptions::genPairedOut(outOpts.outFilename_.string() + "_overhangs");
	perfectOverlapCombinedOpts.out_.transferOverwriteOpts(outOpts);
	r1EndsInR2CombinedOpts.out_.transferOverwriteOpts(outOpts);
	r1BeginsInR2CombinedOpts.out_.transferOverwriteOpts(outOpts);
	r1AllInR2CombinedOpts.out_.transferOverwriteOpts(outOpts);
	r2AllInR1ombinedOpts.out_.transferOverwriteOpts(outOpts);
	notCombinedOpts.out_.transferOverwriteOpts(outOpts);
	overhangsOpts.out_.transferOverwriteOpts(outOpts);
	perfectOverlapCombinedWriter = std::make_unique<SeqOutput>(perfectOverlapCombinedOpts);
	r1EndsInR2CombinedWriter= std::make_unique<SeqOutput>(r1EndsInR2CombinedOpts);
	r1BeginsInR2CombinedWriter= std::make_unique<SeqOutput>(r1BeginsInR2CombinedOpts);
	r1AllInR2CombinedWriter= std::make_unique<SeqOutput>(r1AllInR2CombinedOpts);
	r2AllInR1CombinedWriter= std::make_unique<SeqOutput>(r2AllInR1ombinedOpts);
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
	checkWriter(r1AllInR2CombinedWriter, "r1AllInR2CombinedWriter");
	checkWriter(r2AllInR1CombinedWriter, "r2AllInR1CombinedWriter");
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
	r1AllInR2CombinedWriter = nullptr;//();
	r2AllInR1CombinedWriter = nullptr;//();
	notCombinedWriter = nullptr;//(notCombinedOpts);
	overhangsWriter = nullptr;//(overhangsOpts);
}


Json::Value PairedReadProcessor::ProcessedResultsCounts::toJson() const{
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

	outVal["r1AllInR2Combined"] = r1AllInR2Combined;
	outVal["r1AllInR2CombinedPerc"] = (r1AllInR2Combined/static_cast<double>(total)) *100;
	outVal["r2AllInR1Combined"] = r2AllInR1Combined;
	outVal["r2AllInR1CombinedPerc"] = (r2AllInR1Combined/static_cast<double>(total)) *100;


	outVal["total"] = total;
	if(nullptr != perfectOverlapCombinedOpts){
		outVal["perfectOverlapCombinedOpts"] = njh::json::toJson(perfectOverlapCombinedOpts);
	}
	if(nullptr != r1EndsInR2CombinedOpts){
		outVal["r1EndsInR2CombinedOpts"] = njh::json::toJson(r1EndsInR2CombinedOpts);
	}
	if(nullptr != r1BeginsInR2CombinedOpts){
		outVal["r1BeginsInR2CombinedOpts"] = njh::json::toJson(r1BeginsInR2CombinedOpts);
	}
	if(nullptr != r1AllInR2CombinedOpts){
		outVal["r1AllInR2CombinedOpts"] = njh::json::toJson(r1AllInR2CombinedOpts);
	}
	if(nullptr != r2AllInR1CombinedOpts){
		outVal["r2AllInR1CombinedOpts"] = njh::json::toJson(r2AllInR1CombinedOpts);
	}
	if(nullptr != notCombinedOpts){
		outVal["notCombinedOpts"] = njh::json::toJson(notCombinedOpts);
	}
	if(nullptr != overhangsOpts){
		outVal["overhangsOpts"] = njh::json::toJson(overhangsOpts);
	}
	return outVal;
}

Json::Value PairedReadProcessor::ProcessedResultsCounts::toJsonCounts() const{
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

	outVal["r1AllInR2Combined"] = r1AllInR2Combined;
	outVal["r1AllInR2CombinedPerc"] = (r1AllInR2Combined/static_cast<double>(total)) *100;
	outVal["r2AllInR1Combined"] = r2AllInR1Combined;
	outVal["r2AllInR1CombinedPerc"] = (r2AllInR1Combined/static_cast<double>(total)) *100;

	outVal["total"] = total;
	return outVal;
}

PairedReadProcessor::ProcessedResultsCounts PairedReadProcessor::processPairedEnd(
		SeqInput & reader,
		ProcessorOutWriters & writers,
		aligner & alignerObj){
	writers.checkWritersSet(__PRETTY_FUNCTION__);
	ProcessedResultsCounts res;
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

	if(writers.r1AllInR2CombinedWriter->outOpen()){
		res.r1AllInR2CombinedOpts = std::make_shared<SeqIOOptions>(SeqIOOptions::genFastqIn(writers.r1AllInR2CombinedWriter->getPrimaryOutFnp()));
	}
	if(writers.r2AllInR1CombinedWriter->outOpen()){
		res.r2AllInR1CombinedOpts = std::make_shared<SeqIOOptions>(SeqIOOptions::genFastqIn(writers.r2AllInR1CombinedWriter->getPrimaryOutFnp()));
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

PairedReadProcessor::ProcessedPairRes PairedReadProcessor::processPairedEnd(
		PairedRead & seq,
		ProcessedResultsCounts & counts,
		aligner & alignerObj) const {


	ProcessedPairRes ret;
	double percentId = 1 - params_.errorAllowed_ ;
	if(0 != params_.r1Trim_){
		readVecTrimmer::trimOffEndBases(seq.seqBase_, params_.r1Trim_);
	}
	if(0 != params_.r2Trim_){
		if(seq.mateRComplemented_){
			readVecTrimmer::trimOffForwardBases(seq.mateSeqBase_, params_.r2Trim_);
		}else{
			readVecTrimmer::trimOffEndBases(seq.mateSeqBase_, params_.r2Trim_);
		}
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ <<  std::endl;
	alignerObj.alignRegGlobalNoInternalGaps(seq.seqBase_, seq.mateSeqBase_);
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ <<  std::endl;
	alignerObj.profileAlignment(seq.seqBase_, seq.mateSeqBase_, false, true, true);
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ <<  std::endl;
//	std::cout << "alignerObj.comp_.distances_.eventBasedIdentityHq_: " << alignerObj.comp_.distances_.eventBasedIdentityHq_<< std::endl;
//	std::cout << "alignerObj.comp_.distances_.basesInAln_: " << alignerObj.comp_.distances_.basesInAln_ << std::endl;
//	std::cout << "alignerObj.comp_.hqMismatches_: " << alignerObj.comp_.hqMismatches_ << std::endl;
//
//	OutOptions tempOutR1BEGINSINR2Opts(bfs::path("temp_failedOverLap.fastq"));
//	tempOutR1BEGINSINR2Opts.append_ = true;
//	OutputStream tempOutR1BEGINSINR2(tempOutR1BEGINSINR2Opts);
//

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
			ss << "R1.front(): " << alignerObj.alignObjectA_.seqBase_.seq_.front() << ", R2.front(): " << alignerObj.alignObjectB_.seqBase_.seq_.front() << "\n";
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
			ss << "R1.front(): " << alignerObj.alignObjectA_.seqBase_.seq_.front() << ", R2.front(): " << alignerObj.alignObjectB_.seqBase_.seq_.front() << "\n";
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

			ret.combinedSeq_ = std::make_shared<seqInfo>(seq.seqBase_.name_, cseq, quals);
			ret.status_ = ReadPairOverLapStatus::PERFECTOVERLAP;
			//writers.perfectOverlapCombinedWriter->openWrite(combinedSeq);
			++counts.perfectOverlapCombined;
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
			ret.combinedSeq_ = std::make_shared<seqInfo>(seq.seqBase_.name_, cseq, quals);
			ret.status_ = ReadPairOverLapStatus::R1ENDSINR2;
			++counts.r1EndsInR2Combined;
		}else if((AlignOverlapEnd::NOOVERHANG == frontCase || AlignOverlapEnd::R2OVERHANG == frontCase) &&
						 (AlignOverlapEnd::NOOVERHANG == backCase  || AlignOverlapEnd::R1OVERHANG == backCase)){
			//read through situation, R2 end overlaps R1 beg, overhang is likely illumina adaptor/primer
			uint32_t r1Start = 0;
			if(AlignOverlapEnd::NOOVERHANG != frontCase){
				r1Start = alignerObj.alignObjectA_.seqBase_.seq_.find_first_not_of('-');
			}
			uint32_t r2End = alignerObj.alignObjectB_.seqBase_.seq_.size();
			if(AlignOverlapEnd::NOOVERHANG != backCase){
				r2End =
						alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-') + 1;
			}
			if(params_.writeOverHangs_){
				if(AlignOverlapEnd::NOOVERHANG != frontCase){
					ret.r2Overhang_ = std::make_shared<seqInfo>(alignerObj.alignObjectB_.seqBase_.getSubRead(0, r1Start));
					ret.r2Overhang_->reverseComplementRead(false, true);
				}
				if(AlignOverlapEnd::NOOVERHANG != backCase){
					ret.r1Overhang_ = std::make_shared<seqInfo>(alignerObj.alignObjectA_.seqBase_.getSubRead(r2End));
				}
				//writers.overhangsWriter->openWrite(overhang);
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
			ret.combinedSeq_ = std::make_shared<seqInfo>(seq.seqBase_.name_, cseq, quals);
			ret.status_ = ReadPairOverLapStatus::R1BEGINSINR2;
			++counts.r1BeginsInR2Combined;
		} else if(AlignOverlapEnd::R2OVERHANG == frontCase && AlignOverlapEnd::R2OVERHANG == backCase){
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
			ret.combinedSeq_ = std::make_shared<seqInfo>(seq.seqBase_.name_, cseq, quals);
			ret.status_ = ReadPairOverLapStatus::R1ALLINR2;
			//writers.perfectOverlapCombinedWriter->openWrite(combinedSeq);
			++counts.r1AllInR2Combined;
		} else if(AlignOverlapEnd::R1OVERHANG == frontCase && AlignOverlapEnd::R1OVERHANG == backCase){
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
			ret.combinedSeq_ = std::make_shared<seqInfo>(seq.seqBase_.name_, cseq, quals);
			ret.status_ = ReadPairOverLapStatus::R2ALLINR1;
			//writers.perfectOverlapCombinedWriter->openWrite(combinedSeq);
			++counts.r2AllInR1Combined;
		} else {
			//failure
//			writers.notCombinedWriter->openWrite(seq);
//			++res.overhangFail;
			OutOptions tempFileOpts(njh::files::findNonexitantFile(bfs::path("temp.fastq")));
			OutputStream tempFile(tempFileOpts);

			alignerObj.alignObjectA_.seqBase_.outPutFastq(tempFile);
			alignerObj.alignObjectB_.seqBase_.outPutFastq(tempFile);


			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error not handled case for " << seq.seqBase_.name_ << "\n";
			ss << "R1 Case: " << getAlignOverlapEndStr(frontCase) << ", R2 Case: " << getAlignOverlapEndStr(backCase) << "\n";
			throw std::runtime_error{ss.str()};
		}
	}else{
		++counts.overlapFail;
		//writers.notCombinedWriter->openWrite(seq);
		ret.status_ = ReadPairOverLapStatus::NOOVERLAP;
		ret.combinedSeq_ = nullptr; //not really needed as it should default construct to nullptr anyways
//		alignerObj.alignObjectA_.seqBase_.outPutFastq(tempOutR1BEGINSINR2);
//		alignerObj.alignObjectB_.seqBase_.outPutFastq(tempOutR1BEGINSINR2);
	}
	return ret;
}


bool PairedReadProcessor::processPairedEnd(
		SeqInput & reader,
		PairedRead & seq,
		ProcessorOutWriters & writers,
		aligner & alignerObj,
		ProcessedResultsCounts & res){
//  uncomment for debugging R1ENDSINR2
//	OutOptions tempOutR1BEGINSINR2Opts(bfs::path("temp_R1BEGINSINR2.fastq"));
//	tempOutR1BEGINSINR2Opts.append_ = true;
//	OutputStream tempOutR1BEGINSINR2(tempOutR1BEGINSINR2Opts);
//
//	OutOptions tempOutR1ENDSINR2Opts(bfs::path("temp_R1ENDSINR2.fastq"));
//	tempOutR1ENDSINR2Opts.append_ = true;
//	OutputStream tempOutR1ENDSINR2(tempOutR1ENDSINR2Opts);

	if(reader.readNextRead(seq)){
		++res.total;
		if(res.total % 25000 == 0 && params_.verbose_){
			std::cout << res.total << std::endl;
		}
		auto processedPairResults = processPairedEnd(seq, res, alignerObj);
//    uncomment for debugging
//		if(ReadPairOverLapStatus::R1BEGINSINR2 == processedPairResults.status_ ){
//			alignerObj.alignObjectA_.seqBase_.outPutFastq(tempOutR1BEGINSINR2);
//			alignerObj.alignObjectB_.seqBase_.outPutFastq(tempOutR1BEGINSINR2);
//		}
//		if(ReadPairOverLapStatus::R1ENDSINR2 == processedPairResults.status_ ){
//			alignerObj.alignObjectA_.seqBase_.outPutFastq(tempOutR1ENDSINR2);
//			alignerObj.alignObjectB_.seqBase_.outPutFastq(tempOutR1ENDSINR2);
//		}


		std::stringstream errorStream;
		std::stringstream errorStream2;
		switch (processedPairResults.status_) {
			case ReadPairOverLapStatus::NONE:
				errorStream << __PRETTY_FUNCTION__ << ", error not handled case for " << seq.seqBase_.name_ << ", case: " <<getOverlapStatusStr(processedPairResults.status_) << "\n";
				throw std::runtime_error{errorStream.str()};
				break;
			case ReadPairOverLapStatus::NOOVERLAP:
				writers.notCombinedWriter->openWrite(seq);
				break;
			case ReadPairOverLapStatus::PERFECTOVERLAP:
				writers.perfectOverlapCombinedWriter->openWrite(processedPairResults.combinedSeq_);
				break;
			case ReadPairOverLapStatus::R1BEGINSINR2:
				writers.r1BeginsInR2CombinedWriter->openWrite(processedPairResults.combinedSeq_);
				if(params_.writeOverHangs_){
					writers.overhangsWriter->openWrite(PairedRead(*processedPairResults.r1Overhang_, *processedPairResults.r2Overhang_, false));
				}
				break;
			case ReadPairOverLapStatus::R1ENDSINR2:
				writers.r1EndsInR2CombinedWriter->openWrite(processedPairResults.combinedSeq_);
				break;
			case ReadPairOverLapStatus::R1ALLINR2:
				writers.r1AllInR2CombinedWriter->openWrite(processedPairResults.combinedSeq_);
				break;
			case ReadPairOverLapStatus::R2ALLINR1:
				writers.r2AllInR1CombinedWriter->openWrite(processedPairResults.combinedSeq_);
				break;
			default:
				errorStream2 << __PRETTY_FUNCTION__ << ", error not handled case for " << seq.seqBase_.name_ << ", case: " <<getOverlapStatusStr(processedPairResults.status_) << "\n";
				throw std::runtime_error{errorStream2.str()};
				break;
		}
		return true;
	}else{
		return false;
	}
}


}  // namespace njhseq
