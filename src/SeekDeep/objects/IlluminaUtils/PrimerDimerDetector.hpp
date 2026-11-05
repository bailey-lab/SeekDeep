#pragma once
/*
 * PrimerDimerDetector.hpp
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

#include <njhseq/alignment/aligner/aligner.hpp>
#include <njhseq/seqToolsUtils/determinators/PrimerDeterminator.hpp>

namespace njhseq {

/**@brief Find primer dimers in reads where the mates have already been stitched and read through each other
 *
 * A primer dimer reads as primer (5' end at the start of the read) then a reverse complemented primer (5' end at the
 * end of the read), so both matches are anchored at a primer's 5' end and either can stop short of its 3' end. Any
 * two primers can pair up, including two forward or two reverse primers.
 *
 * Only that one orientation is searched for, which assumes the current panel design where forward primers carry one
 * sequencing adapter and reverse primers the other, so that's the only form that gets both adapters and is sequenced.
 * A panel with both adapters on both primers could produce other forms and this would need to search for them.
 *
 * A primer's 3' end annealing near the 5' end of another primer leaves a read with only those few bases of the second
 * primer, so short matches have to count. How short can be trusted depends on the panel, the more primers there are the
 * more likely one of them matches a read's end by chance, see calcMinPrimerMatchBases(). A read with a trusted primer at
 * one end and only shorter matches at the other is an unresolved dimer, what's at the other end could be chance, and
 * the shorter matches are kept as candidates rather than choosing one.
 */
class PrimerDimerDetector{
public:

	struct DetectorPars{
		DetectorPars();
		comparison allowableErrors_; /**< errors allowed in each primer's match */
		uint32_t maxGapBetweenPrimers_{20}; /**< most bases allowed between the two primers, the primers can also overlap or abut */
		uint32_t minPrimerMatchBases_{0}; /**< fewest primer bases a match has to cover to be trusted, 0 to use calcMinPrimerMatchBases() */
		uint32_t minCandidateMatchBases_{5}; /**< fewest primer bases a shorter match has to cover to be kept as a candidate partner */
	};

	enum class ReadCategory{
		DIMER, /**< trusted primers at both ends within maxGapBetweenPrimers_ of each other */
		PRIMERS_FAR_APART, /**< trusted primers at both ends but with more than maxGapBetweenPrimers_ between them */
		UNRESOLVED_DIMER, /**< a trusted primer at one end and only candidate partners at the other */
		ONE_PRIMER, /**< a trusted primer at one end and nothing at the other */
		NO_PRIMER /**< no trusted primer at either end */
	};

	static std::string getReadCategoryStr(ReadCategory category);

	/**@brief The fewest primer bases to trust for a panel of this many primers
	 *
	 * The bases needed for the chance of any one of the primers matching a read's end to be about 0.1%, with no less than 8
	 *
	 * @param numPrimers the number of primer sequences, forward and reverse both count
	 * @return the number of bases
	 */
	static uint32_t calcMinPrimerMatchBases(uint64_t numPrimers);

	/**@brief Where a primer matched in a read
	 *
	 */
	struct PrimerMatch{
		std::string primerName_;
		comparison comp_;
		seqInfo alignedPrimer_; /**< the primer as aligned, for debugging */
		uint32_t seqStart_{0};
		uint32_t seqEnd_{0};
		uint32_t primerStart_{0};
		uint32_t primerEnd_{0};
		uint32_t primerLen_{0};
		bool revComp_{false}; /**< whether it was the reverse complement of the primer that matched */

		[[nodiscard]] uint32_t primerBasesMatched() const;
		[[nodiscard]] double primerFractionMatched() const;
		/**@brief bases between this match and the other in the read, 0 if they overlap or abut */
		[[nodiscard]] uint32_t getGap(const PrimerMatch & other) const;
		/**@brief how many bases of the read both this match and the other cover, 0 if they don't overlap */
		[[nodiscard]] uint32_t getOverlap(const PrimerMatch & other) const;
	};

	/**@brief A primer at the start of the read and one at the end close enough to be called a dimer
	 *
	 */
	struct PrimerPair{
		PrimerPair(PrimerMatch front, PrimerMatch end);
		PrimerMatch front_;
		PrimerMatch end_;

		[[nodiscard]] uint32_t getGap() const;
		[[nodiscard]] uint32_t getOverlap() const;
		[[nodiscard]] double getAlnScoreSum() const;
		[[nodiscard]] double getIdentitySum() const;
	};

	struct DetectorResults{
		std::vector<PrimerMatch> frontMatches_; /**< trusted primers matching the start of the read */
		std::vector<PrimerMatch> endMatches_; /**< trusted primers matching the end of the read */
		std::vector<PrimerPair> dimers_; /**< the best scoring pairs within maxGapBetweenPrimers_, more than one only if they tie */
		std::vector<PrimerMatch> candidatePartners_; /**< for an unresolved dimer, the shorter matches at the other end from the trusted primer within maxGapBetweenPrimers_ of it */
		ReadCategory category_{ReadCategory::NO_PRIMER};
	};

	PrimerDimerDetector(const PrimerDeterminator & primers, DetectorPars pars);

	DetectorPars pars_;

	/**@brief Create an aligner to give to detect()
	 *
	 * @param maxSeqSize the longest read to be tested, aligners don't grow on their own
	 * @return the aligner
	 */
	[[nodiscard]] aligner genAligner(uint64_t maxSeqSize) const;

	/**@brief Find the dimers in a read
	 *
	 * @param seq the read, should be the stitched mates so that it starts and ends with what was sequenced off of the fragment
	 * @param alignerObj an aligner from genAligner() big enough for seq
	 * @return the primers found and the best pairs
	 */
	[[nodiscard]] DetectorResults detect(const seqInfo & seq, aligner & alignerObj) const;

	[[nodiscard]] uint64_t getMaxPrimerSize() const;

private:
	struct NamedPrimer{
		std::string name_;
		seqInfo info_;
		seqInfo infoRC_;
	};
	std::vector<NamedPrimer> primers_;
	uint64_t maxPrimerSize_{0};

	static PrimerMatch genMatch(const std::string & name, const seqInfo & primer, bool revComp, aligner & alignerObj);
};

}  // namespace njhseq
