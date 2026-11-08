#pragma once
/*
 * NanoporeAmpliconReadQC.hpp
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

#include <njhseq/seqToolsUtils/FoldbackDetector.hpp>
#include "SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp"

namespace njhseq {

/**@brief Say what each nanopore read of a targeted amplicon sample is, directly from the read without a genome
 *
 * A complete amplicon on a read is one primer of a target as written, then the reverse complement of the target's other
 * primer at about the expected distance, with sequence between them that shares kmers with the target's reference
 * sequences. A read can hold several, which is what separates the categories first:
 *
 * Whole amplicons joined end to end are reads made of two molecules, with different targets or a target head to tail,
 * and each amplicon can be cut out and used. A target inverted against itself is usually one molecule read through a
 * turn, so it counts as one amplicon. These have to be told apart from foldbacks, FoldbackDetector calls some of each.
 *
 * Reads without a complete amplicon are then checked for a fold, for the target's kmers, and finally for any primers.
 * Whether an off target fold is host needs a reference, which this doesn't use.
 */
class NanoporeAmpliconReadQC{
public:

	enum class ReadCategory{
		AMPLICON, /**< one complete amplicon */
		AMPLICON_PLUS_OTHER, /**< one complete amplicon with a lot of other sequence on the read, no fold */
		AMPLICON_WITH_FOLD, /**< one complete amplicon and a fold */
		SAME_TARGET_INVERTED, /**< 2 or more complete amplicons, all one target and each inverted against the next, usually one molecule read through a turn */
		JOIN_SAME_TARGET, /**< 2 or more complete amplicons, all one target, at least one junction head to tail */
		JOIN_DIFFERENT_TARGETS, /**< 2 or more complete amplicons of different targets */
		FOLD_ON_TARGET, /**< no complete amplicon, a fold with arms made of a target's sequence */
		FOLD_OFF_TARGET, /**< no complete amplicon, a fold with arms that aren't */
		ON_TARGET_UNEXPECTED_LENGTH, /**< no fold, the target's sequence and both of its primers but not at the expected distance */
		ON_TARGET_ONE_END, /**< no fold, the target's sequence and only one of its primers */
		TARGET_NO_PRIMERS, /**< no fold, the target's sequence and none of its primers */
		PARTIAL_MULTI_TARGET, /**< no fold, no complete amplicon, sequence of more than one target */
		PRIMER_DIMER, /**< not a target's sequence, primers next to each other */
		TWO_PRIMER_OFF_TARGET, /**< not a target's sequence, primers at both ends of something */
		ONE_PRIMER_OFF_TARGET, /**< not a target's sequence, some primer on it */
		NO_PRIMER /**< nothing from the panel */
	};

	static const std::vector<ReadCategory> allCategories_;

	static std::string getCategoryStr(ReadCategory category);

	struct Pars{
		FoldbackDetector::Pars foldPars_;
		double maxPrimerError_{0.20}; /**< the fraction of a primer that can be edits, nanopore errors plus the demultiplexing's trimming */
		double minInnerKmerFrac_{0.30}; /**< the fraction of kmers between the primers that have to be in the target's reference sequences */
		uint32_t minLenTolerance_{40}; /**< how far off the expected length a complete amplicon can be, at least */
		double lenToleranceFrac_{0.25}; /**< and as a fraction of the expected length, whichever is larger */
		uint32_t maxDimerInsert_{30}; /**< the most bases between two primers on a read for it to be a dimer */
		uint32_t otherBasesCutOff_{100}; /**< the most bases that can be outside a read's only complete amplicon before it's noted */
		double minTargetKmerFrac_{0.30}; /**< the fraction of a read or fold arm's kmers that have to be a target's to say it's that target */
		double otherTargetKmerFrac_{0.10}; /**< the same for a second target to also count */
		uint32_t turnBinSize_{20}; /**< the bin size for where folds turn along an amplicon */
		uint32_t primerAtStartWithin_{60}; /**< a read's first primer has to start within this many bases to place its fold on the amplicon */
		uint32_t minLoopForPileup_{11}; /**< folds with a shorter loop are left out of where folds turn, the ends of the arms are not well defined for these */
		uint32_t decoyEvery_{20}; /**< check shuffled primers on 1 in this many reads, 0 to skip */
	};

	/**@brief A primer found on a read */
	struct PrimerHit{
		uint32_t start_{0}; /**< 0 based, inclusive */
		uint32_t end_{0}; /**< exclusive */
		uint32_t probe_{0};
		uint32_t edits_{0};
	};

	/**@brief What's being searched for, a primer in one direction */
	struct Probe{
		std::string seq_;
		uint32_t target_{0};
		bool forwardPrimer_{true};
		bool revComp_{false}; /**< false for the primer as written */
	};

	/**@brief A complete amplicon on a read */
	struct Unit{
		uint32_t start_{0};
		uint32_t end_{0};
		uint32_t target_{0};
		bool startsWithForwardPrimer_{true};
		uint32_t startEdits_{0};
		uint32_t endEdits_{0};
		double innerFrac_{0};
	};

	/**@brief Everything worked out from one read */
	struct ReadResult{
		ReadCategory category_{ReadCategory::NO_PRIMER};
		int64_t target_{-1}; /**< the read's target if it has one, -1 otherwise, the first unit's for joins */
		uint64_t len_{0};
		std::vector<Unit> units_; /**< at Pars::maxPrimerError_ */
		std::vector<uint32_t> unitsAtStricterErrors_; /**< the number of units at each of stricterPrimerErrors_ */
		std::vector<bool> keepOnce_; /**< per unit, false if it's the second of an inverted same target pair */
		FoldbackDetector::Fold fold_;
		uint32_t primerHits_{0};
		int64_t startProbe_{-1}; /**< the primer as written nearest the start of the read if it's within Pars::primerAtStartWithin_ */
		uint32_t startProbePos_{0};
		int64_t dimerProbe1_{-1}; /**< the first and second primer of a dimer */
		int64_t dimerProbe2_{-1};
		int64_t offTargetPrimer1_{-1}; /**< for twoPrimerOffTarget, the pair of primers furthest apart, as primer ids, the lower first since the read can be either strand */
		int64_t offTargetPrimer2_{-1};
		uint32_t offTargetInsert_{0}; /**< the bases between that pair */
		bool offTargetAtEnds_{false}; /**< both of that pair are within Pars::primerAtStartWithin_ of the read's ends */
		std::vector<uint32_t> offTargetPrimers_; /**< for off target reads with primers, every primer on it once, as primer ids */

		// decoys, only for the reads that are checked
		bool decoyChecked_{false};
		uint32_t decoyHits_{0};
		std::vector<uint32_t> decoyUnitsAtErrors_;
		std::vector<uint32_t> unitsAtErrors_;

		std::vector<int32_t> junctionBarcodes_; /**< per junction, the barcode found there or -1 for none, empty if barcodes weren't given */
	};

	/**@brief The counts of a sample, not thread safe, add to it serially */
	struct SampleTallies{
		std::string sampleName_;
		uint64_t reads_{0};
		uint64_t bases_{0};
		std::map<ReadCategory, uint64_t> categoryReads_;
		std::map<ReadCategory, uint64_t> categoryBases_;
		std::map<std::pair<uint32_t, ReadCategory>, uint64_t> targetCategoryReads_;
		std::map<uint32_t, uint64_t> singleAmpliconReads_;
		std::map<uint32_t, uint64_t> unitsInJoins_;
		std::map<uint32_t, uint64_t> unitsCountedOnce_;
		std::map<uint32_t, uint64_t> readsWithUnits_; /**< how many reads have 2+ by the count, 5 for 5 or more */
		std::vector<uint64_t> joinedReadsAtStricterErrors_;
		uint64_t readsWithForeignBarcode_{0}; /**< 2+ amplicons with another sample's barcode at a junction, these aren't in the unit counts */
		std::map<std::tuple<bool, bool, std::string, std::string>, uint64_t> junctions_; /**< key is same target, inverted, the gap's bin, the barcode status */
		std::map<std::tuple<bool, bool, std::string, std::string>, std::vector<uint32_t>> junctionGaps_;
		uint64_t foldsCalled_{0};
		std::map<uint32_t, uint64_t> targetReads_; /**< reads of a target by any category */
		std::map<uint32_t, uint64_t> targetFolds_; /**< of those, a fold that isn't only a same target inverted join */
		std::map<std::pair<uint32_t, bool>, uint64_t> targetFoldStrand_; /**< key is target, read started with the forward primer */
		std::map<uint32_t, std::vector<uint32_t>> targetFoldLoops_;
		std::map<uint32_t, std::vector<uint32_t>> targetFoldArms_;
		std::map<std::tuple<uint32_t, bool, uint32_t>, uint64_t> turnPileup_; /**< key is target, read started with the forward primer, bin along the amplicon */
		std::map<std::pair<std::string, std::string>, uint64_t> dimers_;
		struct OffTargetPair{
			uint64_t reads_{0};
			uint64_t bases_{0};
			uint64_t atEnds_{0};
			std::vector<uint32_t> lens_;
			std::vector<uint32_t> inserts_;
		};
		std::map<std::pair<uint32_t, uint32_t>, OffTargetPair> offTargetPairs_; /**< key is the primer ids of twoPrimerOffTarget reads */
		std::map<uint32_t, std::array<uint64_t, 2>> offTargetPrimerReads_; /**< key is primer id, value is twoPrimerOffTarget and onePrimerOffTarget reads with it */
		uint64_t decoyReads_{0};
		uint64_t decoyHits_{0};
		uint64_t realHitsOnDecoyReads_{0};
		std::map<std::pair<bool, double>, std::array<uint64_t, 2>> calibration_; /**< key is decoy, primer error, value is reads with 1+ units, 2+ units */
	};

	NanoporeAmpliconReadQC(const PrimersAndMids & ids, Pars pars);

	/**@brief Give the sample barcodes to look for at the junctions of joined amplicons, a read joined after barcoding has one there
	 *
	 * @param barcodesFnp a fasta of the barcodes, the number in each's name is the barcode's number
	 * @param barcodeStart where the barcode starts in each sequence, they can have flanking sequence
	 * @param barcodeLen the length of the barcode
	 */
	void setBarcodes(const bfs::path & barcodesFnp, uint32_t barcodeStart, uint32_t barcodeLen);

	/**@brief What a junction's barcode means for a read
	 *
	 * @param result what came out of processRead
	 * @param junction the junction
	 * @param ownBarcodes the barcodes of the read's sample
	 * @return notChecked, noBarcode, barcodeUnknownSample, own, or foreign
	 */
	[[nodiscard]] static std::string junctionBarcodeStatus(const ReadResult & result, uint32_t junction, const std::set<int32_t> & ownBarcodes);

	Pars pars_;
	std::vector<double> stricterPrimerErrors_{0.15, 0.10};

	/**@brief Work out what a read is, safe to call from multiple threads
	 *
	 * @param seq the read, it should have adapters and barcodes trimmed
	 * @param checkDecoys whether to also search for the shuffled primers
	 * @return what it is
	 */
	[[nodiscard]] ReadResult processRead(const std::string & seq, bool checkDecoys) const;

	/**@brief Add a read to the counts
	 *
	 * @param result what came out of processRead
	 * @param tallies the counts to add to
	 * @param ownBarcodes the barcodes of the read's sample
	 */
	void addToTallies(const ReadResult & result, SampleTallies & tallies, const std::set<int32_t> & ownBarcodes) const;

	[[nodiscard]] const std::vector<std::string> & targetNames() const;
	[[nodiscard]] uint32_t expectedLen(uint32_t target) const;
	[[nodiscard]] const Probe & probe(uint32_t pos) const;

	/**@brief A primer's id, target * 2 for the forward primer, plus 1 for the reverse */
	[[nodiscard]] static uint32_t primerId(const Probe & probe);
	/**@brief target_F or target_R */
	[[nodiscard]] std::string primerName(uint32_t primer) const;
	/**@brief samePrimer, ownPair (a target's forward and reverse), or diffTargets */
	[[nodiscard]] static std::string primerPairType(uint32_t primer1, uint32_t primer2);

	/**@brief Where a read's fold finished its first arm along the amplicon, from the start of the forward primer
	 *
	 * @param result what came out of processRead
	 * @param pos the position if it can be worked out
	 * @return true if the read's fold can be placed on its target's amplicon
	 */
	[[nodiscard]] bool foldTurnPos(const ReadResult & result, uint32_t & pos) const;

	/**@brief The units of a read in the format target|F or R|start|end|edits at start|edits at end|inner kmer fraction;...
	 */
	[[nodiscard]] std::string unitsStr(const ReadResult & result) const;
	/**@brief The junctions of a read in the format sameTarget or diffTarget|headToTail or inverted|gap|barcode status,...
	 */
	[[nodiscard]] std::string junctionsStr(const ReadResult & result, const std::set<int32_t> & ownBarcodes) const;

private:

	struct Panel{
		std::vector<Probe> probes_;
		std::vector<std::vector<std::pair<uint32_t, uint32_t>>> seedIndex_; /**< from a packed seed to probe and where the seed is on it */
	};

	[[nodiscard]] static Panel buildPanel(std::vector<Probe> probes);
	[[nodiscard]] std::vector<PrimerHit> findHits(const Panel & panel, const std::string & seq) const;
	[[nodiscard]] std::vector<Unit> findUnits(const Panel & panel, const std::vector<PrimerHit> & hits, const std::vector<uint32_t> & kmers, double primerError) const;

	/**@brief The target with the most of a stretch's kmers
	 *
	 * @return the best target (-1 if none), the fraction of kmers of it, and the fraction of kmers of the next best
	 */
	[[nodiscard]] std::tuple<int64_t, double, double> bestTarget(const std::vector<uint32_t> & kmers, uint32_t start, uint32_t end) const;

	std::vector<std::string> targetNames_;
	std::vector<uint32_t> expectedLens_;
	Panel panel_;
	Panel decoyPanel_;
	std::vector<std::unordered_set<uint32_t>> targetKmers_;
	std::unordered_map<uint32_t, std::vector<uint32_t>> kmerToTargets_;
	FoldbackDetector detector_;
	std::unordered_map<uint32_t, std::vector<int32_t>> barcodeKmers_;
};

}  // namespace njhseq
