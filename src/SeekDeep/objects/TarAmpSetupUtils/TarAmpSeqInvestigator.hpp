#pragma once

/*
 * TarAmpSeqInvestigator.hpp
 *
 *  Created on: Jan 31, 2020
 *      Author: nicholashathaway
 */





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
		MidDeterminator::MidDeterminePars midPars;
		uint32_t testNumber = std::numeric_limits<uint32_t>::max();

		gapScoringParameters gapInfo_;

		TarAmpSeqInvestigatorPars();

	};

	TarAmpSeqInvestigator(const TarAmpSeqInvestigatorPars & pars);

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

	void addOtherCounts(const TarAmpSeqInvestigator & other);

	//primer pairing counts
	table primerCountsTab_{VecStr{"ForwardPrimer", "ReversePrimer", "ForwardCount", "ReverseCount", "ReverseFraction", "Total", "Fraction"}};
	//preceding number of bases counts
	table precedingBasesCountsTab_{VecStr{"PrimerPair", "BasesPrecedingPrimerCase", "NumOfBases", "Count"}};
	//preceding bases counts, to determine possible mids
	table possibleMidCounts_{VecStr{"PrimerPair", "ForwardMID", "ReverseMID", "InForDirCount", "InRevDirCount", "TotalCount"}};
	table possibleMidCountsMostCommonTab_{VecStr{"PrimerPair", "ForwardMID", "ReverseMID", "InForDirCount", "InRevDirCount", "TotalCount"}};
	//unrecognized counts
	table unrecoginzedCountsTab_{VecStr{"forward", "reverse", "count", "fraction"}};
//	table onePrimerUnrecoginzedCountsTab_{VecStr{"recognizedPrimer", "", "count", "fraction"}};

	void investigateSeq(const seqInfo & forwardSeq, const seqInfo & revCompSeq, aligner & alignObj);
	void processCounts();
	void writeOutTables(const bfs::path & directory, bool overWrite);



	struct prepareForInvestiagteFileRes{
		uint64_t maxReadSize = 0;
		uint32_t readCount = 0;
	};

	prepareForInvestiagteFileRes prepareForInvestiagteFile(const SeqIOOptions & opts, bool verbose);
	void investigateFile(const SeqIOOptions & opts, const prepareForInvestiagteFileRes & counts, bool verbose);
	void investigateFile(const SeqIOOptions & opts, bool verbose);

	bool reverseComplementLikely(uint32_t minReadAmount = 250,
			double cutOff = 0.2) const;
	bool hasPossibleRandomPrecedingBases(uint32_t midSize,
			uint32_t minReadAmount = 250, double cutOff = 0.1) const;

	bool hasPossibleRandomPrecedingBasesForwardPrimer(uint32_t midSize,
			uint32_t minReadAmount = 250, double cutOff = 0.1) const;

	bool hasPossibleRandomPrecedingBasesReversePrimer(uint32_t midSize,
			uint32_t minReadAmount = 250, double cutOff = 0.1) const;

	uint32_t maxPrecedingBases() const;
	uint32_t maxPrecedingReversePrimerBases() const;
	uint32_t maxPrecedingForwardPrimerBases() const;
	VecStr recommendSeekDeepExtractorFlags() const;

};




}  // namespace njhseq



