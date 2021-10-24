#pragma once
/*
 * PairedReadProcessor.hpp
 *
 *  Created on: Jan 14, 2018
 *      Author: nick
 */
//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012-2019 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include <njhseq/common.h>
#include <njhseq/IO/SeqIO.h>
#include <njhseq/alignment/aligner/aligner.hpp>

namespace njhseq {

class PairedReadProcessor{
public:

	enum class AlignOverlapEnd{
		NOOVERHANG,
		R1OVERHANG,
		R2OVERHANG,
		UNHANDLEED
	};

	enum class ReadPairOverLapStatus{
		NOOVERLAP,
		R1BEGINSINR2,
		R1ENDSINR2,
		R1ALLINR2,
		R2ALLINR1,
		PERFECTOVERLAP,
		AUTO,
		NONE
	};

	static std::string getOverlapStatusStr(const ReadPairOverLapStatus  status){
		switch (status) {
			case ReadPairOverLapStatus::NOOVERLAP:
				return "NOOVERLAP";
				break;
			case ReadPairOverLapStatus::PERFECTOVERLAP:
				return "PERFECTOVERLAP";
				break;
			case ReadPairOverLapStatus::R1BEGINSINR2:
				return "R1BEGINSINR2";
				break;
			case ReadPairOverLapStatus::R1ENDSINR2:
				return "R1ENDSINR2";
				break;
			case ReadPairOverLapStatus::R1ALLINR2:
				return "R1ALLINR2";
				break;
			case ReadPairOverLapStatus::R2ALLINR1:
				return "R2ALLINR1";
				break;
			case ReadPairOverLapStatus::AUTO:
				return "AUTO";
				break;
			default:
				return "NOTHANDLED";//shouldn't be getting here
				break;
		}
	}

	static std::string getAlignOverlapEndStr(const AlignOverlapEnd  status){
		switch (status) {
			case AlignOverlapEnd::NOOVERHANG:
				return "NOOVERHANG";
				break;
			case AlignOverlapEnd::R1OVERHANG:
				return "R1OVERHANG";
				break;
			case AlignOverlapEnd::R2OVERHANG:
				return "R2OVERHANG";
				break;
			case AlignOverlapEnd::UNHANDLEED:
				return "UNHANDLEED";
				break;
			default:
				return "NOTHANDLED";//shouldn't be getting here
				break;
		}
	}

	struct ProcessParams{
		uint32_t minOverlap_ = 10;
		double errorAllowed_ = 0.01;
		uint32_t hardMismatchCutOff_ = 20;
		uint32_t lqMismatchCutOff = 10;
		uint32_t hqMismatchCutOff = 10;

		uint32_t checkAmount_ = 100;
		uint32_t testNumber_ = std::numeric_limits<uint32_t>::max();
		bool verbose_ = false;
		bool debug_ = false;

		bool writeOverHangs_ = false;

		uint32_t primerDimmerSize_{100};

		uint32_t r1Trim_ = 0;
		uint32_t r2Trim_ = 0;

		struct QualWindowTrimPars{
			uint32_t windowSize_ {5};
			uint32_t windowStep_ {1};
			double avgQualCutOff_{25};
		};

		QualWindowTrimPars qualWindowPar_;
		bool trimLowQaulWindows_{true};
		double percentAfterTrimCutOff_ {.70};

	};

	PairedReadProcessor(ProcessParams params);

	void setDefaultConsensusBuilderFunc();
	ProcessParams params_;

	uint64_t guessMaxReadLenFromFile(const SeqIOOptions & inputOpts);

	struct ProcessorOutWriters{
		ProcessorOutWriters();
		ProcessorOutWriters(const OutOptions & outOpts);

		std::unique_ptr<SeqOutput> perfectOverlapCombinedWriter;//(perfectOverlapCombinedOpts);
		std::unique_ptr<SeqOutput> r1EndsInR2CombinedWriter;//(r1EndsInR2CombinedOpts);
		std::unique_ptr<SeqOutput> r1BeginsInR2CombinedWriter;//(r1BeginsInR2CombinedOpts);
		std::unique_ptr<SeqOutput> r1AllInR2CombinedWriter;//();
		std::unique_ptr<SeqOutput> r2AllInR1CombinedWriter;//();
		std::unique_ptr<SeqOutput> notCombinedWriter;//(notCombinedOpts);
		std::unique_ptr<SeqOutput> overhangsWriter;//(overhangsOpts);

		void checkWritersSet(const std::string & funcName);
		void unsetWriters();

		void closeAllOpenWriters();

	};

	struct ProcessedResultsCounts {
		uint32_t overlapFail = 0;
		uint32_t overhangFail = 0;
		uint32_t perfectOverlapCombined = 0;
		uint32_t r1EndsInR2Combined = 0;
		uint32_t r1BeginsInR2Combined = 0;
		uint32_t r1BeginsInR2CombinedAboveCutOff = 0;
		uint32_t r1AllInR2Combined = 0;
		uint32_t r2AllInR1Combined = 0;
		uint32_t total = 0;

		std::shared_ptr<SeqIOOptions> perfectOverlapCombinedOpts;
		std::shared_ptr<SeqIOOptions> r1EndsInR2CombinedOpts;
		std::shared_ptr<SeqIOOptions> r1BeginsInR2CombinedOpts;
		std::shared_ptr<SeqIOOptions> r1AllInR2CombinedOpts;
		std::shared_ptr<SeqIOOptions> r2AllInR1CombinedOpts;
		std::shared_ptr<SeqIOOptions> notCombinedOpts;
		std::shared_ptr<SeqIOOptions> overhangsOpts;

		void addOther(const ProcessedResultsCounts & otherCounts);

		Json::Value toJson() const;
		Json::Value toJsonCounts() const;
	};

	std::function<void(uint32_t, const seqInfo&,const seqInfo&,std::string&,std::vector<uint8_t>&,aligner&)> addToConsensus;


	ProcessedResultsCounts processPairedEnd(
			SeqInput & reader,
			ProcessorOutWriters & writers,
			aligner & alignerObj);

	struct ProcessedPairRes{
		std::shared_ptr<seqInfo> combinedSeq_;
		std::shared_ptr<seqInfo> r1Overhang_;
		std::shared_ptr<seqInfo> r2Overhang_;
		ReadPairOverLapStatus status_ {ReadPairOverLapStatus::NONE};
	};

	ProcessedPairRes processPairedEnd(
			PairedRead & seq,
			ProcessedResultsCounts & counts,
			aligner & alignerObj) const;

	bool processPairedEnd(
			SeqInput & reader,
			PairedRead & seq,
			ProcessorOutWriters & writers,
			aligner & alignerObj,
			ProcessedResultsCounts & res);

};


}  // namespace njhseq




