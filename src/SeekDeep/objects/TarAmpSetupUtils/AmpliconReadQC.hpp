#pragma once
/*
 * AmpliconReadQC.hpp
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
#include <njhseq/IO/SeqIO/SeqIOOptions.hpp>
#include "SeekDeep/objects/IlluminaUtils/PairedReadProcessor.hpp"
#include "SeekDeep/objects/IlluminaUtils/PrimerDimerDetector.hpp"
#include "SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp"

namespace njhseq {

/**@brief Put every read pair of a targeted amplicon sample into one category without any of the extraction or clustering
 *
 * The mates are stitched with the settings given in QcPars::stitchPars_, the same as the extraction should use, so what's
 * reported about the stitching is what extraction would see. Pairs that don't stitch with those settings get another
 * attempt at QcPars::dimerOverLapErrorAllowed_ that's only used to find dimers, primer dimers are short enough that a
 * single sequencing error can stop them from stitching, it's never used for the other categories since allowing that
 * many errors can force a bad stitch for other reads.
 */
class AmpliconReadQC{
public:

	enum class ReadPairCategory{
		ON_TARGET, /**< the primers are from the same target and the mates stitched how that target is expected to */
		UNEXPECTED_STITCH, /**< the primers are from the same target but the mates didn't stitch how that target is expected to */
		KMER_MISMATCH, /**< the primers are from the same target but the sequence between them doesn't share enough kmers with the target's reference sequences, only if the target has any */
		PRIMER_DIMER, /**< read through each other and made up of primers from each end */
		UNRESOLVED_DIMER, /**< read through each other, a trusted primer at one end and only short matches at the other */
		SHORT_PRODUCT, /**< read through each other with primers at both ends and too much between them to be a dimer, and not expected for the target */
		SHORT_ONE_PRIMER, /**< read through each other with a primer at only one end */
		SHORT_NO_PRIMER, /**< read through each other with no primers */
		MISMATCHED_PRIMERS, /**< didn't read through each other, the primers are from different targets */
		FAILED_PRIMERS /**< didn't read through each other, a primer wasn't found on at least one of the mates */
	};

	static const std::vector<ReadPairCategory> allCategories_;

	static std::string getCategoryStr(ReadPairCategory category);

	struct QcPars{
		QcPars();
		PairedReadProcessor::ProcessParams stitchPars_;
		double dimerOverLapErrorAllowed_{0.05}; /**< the fraction of errors allowed when stitching only to find dimers */
		PrimerDimerDetector::DetectorPars dimerPars_;
		PrimerDeterminator::PrimerDeterminatorPars primerPars_;
		uint32_t primerCallCacheMaxEntries_{200000}; /**< the most distinct read starts to remember the primer found for, per sample, 0 to not remember them */
		uint32_t compKmerLen_{5}; /**< the kmer size for comparing to the target's reference sequences, has to be what the references were set up with */
		double compKmerSimCutOff_{0.50}; /**< the fraction of kmers that have to be shared with a reference sequence */
	};

	/**@brief A read that read through its mate but isn't a dimer, collapsed on its sequence
	 *
	 */
	struct ShortProduct{
		ReadPairCategory category_{ReadPairCategory::SHORT_NO_PRIMER};
		std::string frontPrimer_;
		std::string endPrimer_;
		std::string gapSeq_; /**< the sequence between the primers, empty if they aren't both known */
		uint32_t gap_{0};
		uint32_t count_{0};
	};

	struct LengthTally{
		uint32_t count_{0};
		std::vector<uint32_t> lengths_; /**< only for pairs that were stitched */
	};

	struct SampleResults{
		std::string sampleName_;
		uint64_t total_{0};
		std::map<ReadPairCategory, uint64_t> categoryCounts_;
		uint64_t primersMatched_{0}; /**< pairs where both mates' primers are from the same target */
		uint64_t primersMismatched_{0}; /**< pairs where both primers were found but are from different targets */
		uint64_t primersFailed_{0}; /**< pairs where a primer wasn't found on at least one of the mates */
		uint64_t dimersFromRetryStitch_{0}; /**< pairs only found to be dimers by the second attempt at stitching */
		PrimerDimerTallies dimers_;
		std::map<std::tuple<std::string, std::string, ReadPairCategory>, LengthTally> targetCounts_; /**< key is target, stitch status, category */
		std::map<std::tuple<std::string, std::string, std::string, ReadPairCategory>, LengthTally> otherPrimerCounts_; /**< key is forward primer, reverse primer, stitch status, category */
		std::map<std::string, ShortProduct> shortProducts_; /**< key is the read's sequence */
		std::map<std::pair<std::string, std::string>, uint32_t> kmerMismatchSeqs_; /**< key is target and the sequence between the primers (or just from the first mate if not stitched), value is the count */
		std::map<std::tuple<std::string, std::string, std::string>, uint32_t> mismatchedPrimerSeqs_; /**< key is forward primer, reverse primer, and the sequence between the primers (or just from the first mate if not stitched), value is the count */
	};

	/**@brief Set up for the targets
	 *
	 * @param ids the primers with the primer determinator initialized and the expected overlap statuses added, has to outlive this object
	 * @param pars the settings
	 */
	AmpliconReadQC(const PrimersAndMids & ids, QcPars pars);

	const PrimersAndMids & ids_;
	QcPars pars_;
	PrimerDimerDetector detector_;

	/**@brief Go through all the read pairs of a sample, safe to call from multiple threads for different samples
	 *
	 * @param readOpts the input for the sample, needs to be paired and shouldn't reverse complement the mate
	 * @param sampleName the name of the sample
	 * @return the counts for the sample
	 */
	[[nodiscard]] SampleResults processSample(const SeqIOOptions & readOpts, const std::string & sampleName) const;

	/**@brief Whether a target is expected to stitch this way
	 *
	 * @param target the target
	 * @param status how the mates stitched
	 * @return true if it's expected
	 */
	[[nodiscard]] bool expectedStatus(const std::string & target, PairedReadProcessor::ReadPairOverLapStatus status) const;
};

/**@brief Work out where the products that aren't what was expected came from
 *
 * Sequences are checked for adapter first as many database records have adapter left in them, a BLAST hit to adapter says
 * nothing about where a sequence came from, then for being too short to place, and the rest are BLASTed against the
 * target genomes and the host genome.
 */
class AmpliconProductPlacer{
public:
	struct TargetGenome{
		std::string name_;
		bfs::path genomeFnp_;
		bfs::path targetBedFnp_; /**< the locations of the targets in this genome, from genTargetInfoFromGenomes */
	};

	struct PlacerPars{
		PlacerPars();
		std::vector<TargetGenome> targetGenomes_;
		bfs::path hostGenomeFnp_; /**< can be left empty */
		VecStr adapters_; /**< the sequences to check for, both directions are checked */
		uint32_t adapterKmerLen_{8}; /**< how many bases in a row have to match an adapter */
		uint32_t minBlastLen_{25}; /**< sequences shorter than this can't be placed with any confidence */
		double minPercentId_{95}; /**< the lowest percent identity for a hit to count */
		double minCoverage_{0.85}; /**< the lowest fraction of the sequence a hit has to cover to count */
		uint32_t targetRegionFlank_{100}; /**< how far from a target a hit can be and still count as being in that target */
		uint32_t maxSeqsToBlast_{20000}; /**< only the most common sequences are BLASTed if there's more than this */
		uint32_t longSeqLen_{50}; /**< sequences at least this long are searched with larger seeds, much faster against a large genome and at this identity they'll still have one */
		std::string longSeqDust_{"40 64 1"}; /**< the low complexity masking for those searches, milder than the blast default of "20 64 1" which hides too much of an AT rich genome like plasmodium but off altogether is orders of magnitude slower */
		uint32_t numThreads_{1};
	};

	struct Placement{
		std::string placement_{"notBlasted"}; /**< adapter, tooShort, target, host, unplaced, or notBlasted */
		std::string genome_;
		std::string chrom_;
		uint64_t start_{0};
		uint64_t end_{0};
		double percentId_{0};
		double coverage_{0};
		VecStr overlappingTargets_; /**< panel targets a hit to the target genomes falls in */

		/**@brief Say what the placement means for a read with the primers of these targets
		 *
		 * @param primerTargets the targets of the read's primers
		 * @return inPrimerTarget, inOtherPanelTarget, or targetGenomeOffTarget for a hit to the target genomes, otherwise placement_
		 */
		[[nodiscard]] std::string getLabel(const VecStr & primerTargets) const;
	};

	explicit AmpliconProductPlacer(PlacerPars pars);

	PlacerPars pars_;

	[[nodiscard]] bool containsAdapter(const std::string & seq) const;

	/**@brief Place the sequences
	 *
	 * @param seqCounts the sequences and how many reads have them, sequences should be unique across all the samples so each is only BLASTed once
	 * @param workingDir where to write the BLAST input and output
	 * @return the placement of each sequence
	 */
	[[nodiscard]] std::unordered_map<std::string, Placement> placeSeqs(const std::unordered_map<std::string, uint64_t> & seqCounts, const bfs::path & workingDir) const;

private:
	std::unordered_set<std::string> adapterKmers_;
};

}  // namespace njhseq
