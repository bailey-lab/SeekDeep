#pragma once
/*
 * PairedReadProcessor.hpp
 *
 *  Created on: Jan 14, 2018
 *      Author: nick
 */

#include <bibseq.h>
namespace bibseq {

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
		R1BEGOVERR2END,
		R1ENDOVERR2BEG
	};

	struct ProcessParams{
		uint32_t minOverlap_ = 10;
		double errorAllowed_ = 0.01;
		uint32_t hardMismatchCutOff_ = 5;
		uint32_t checkAmount_ = 100;
		uint32_t testNumber_ = std::numeric_limits<uint32_t>::max();
		bool verbose_ = false;
		bool writeOverHangs_ = false;

		uint32_t r1Trim_ = 0;
		uint32_t r2Trim_ = 0;

	};

	PairedReadProcessor(ProcessParams params);

	void setDefaultConsensusBuilderFunc();
	ProcessParams params_;

	uint64_t guessMaxReadLenFromFile(const SeqIOOptions & inputOpts);

	struct ProcessorOutWriters{
		ProcessorOutWriters();
		ProcessorOutWriters(const OutOptions & outOpts);

		std::unique_ptr<SeqOutput> perfectOverlapCombinedWriter;//(perfectOverlapCombinedOpts);
		std::unique_ptr<SeqOutput> r1EndOverR2BegCombinedWriter;//(r1EndOverR2BegCombinedOpts);
		std::unique_ptr<SeqOutput> r1BegOverR2EndCombinedWriter;//(r1BegOverR2EndCombinedOpts);
		std::unique_ptr<SeqOutput> notCombinedWriter;//(notCombinedOpts);
		std::unique_ptr<SeqOutput> overhangsWriter;//(overhangsOpts);

		void checkWritersSet(const std::string & funcName);
		void unsetWriters();

	};

	struct ProcessedResults {
		uint32_t overlapFail = 0;
		uint32_t overhangFail = 0;
		uint32_t perfectOverlapCombined = 0;
		uint32_t r1EndOverR2BegCombined = 0;
		uint32_t r1BegOverR2EndCombined = 0;
		uint32_t total = 0;

		std::shared_ptr<SeqIOOptions> perfectOverlapCombinedOpts;
		std::shared_ptr<SeqIOOptions> r1EndOverR2BegCombinedOpts;
		std::shared_ptr<SeqIOOptions> r1BegOverR2EndCombinedOpts;
		std::shared_ptr<SeqIOOptions> notCombinedOpts;
		std::shared_ptr<SeqIOOptions> overhangsOpts;

		Json::Value toJson() const;
		Json::Value toJsonCounts() const;
	};

	std::function<void(uint32_t, const seqInfo&,const seqInfo&,std::string&,std::vector<uint32_t>&,aligner&)> addToConsensus;


	ProcessedResults processPairedEnd(
			SeqInput & reader,
			ProcessorOutWriters & writers,
			aligner & alignerObj);

	bool processPairedEnd(
			SeqInput & reader,
			PairedRead & seq,
			ProcessorOutWriters & writers,
			aligner & alignerObj,
			ProcessedResults & res);

};


}  // namespace bibseq




