#pragma once
/*
 * NanoporeReadPlacer.hpp
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
#include "SeekDeep/utils.h"

namespace njhseq {

/**@brief Say where nanopore reads and their folds are, by mapping them with minimap2 to the target genome and the host
 *
 * The genomes should be in one index so that the host competes with the target genome for each read.
 *
 * A fold is mapped as two pieces rather than as a read. A self primed fold is S[0:far] + rc(S[0:entry]), paired arm, loop,
 * paired arm, and each side of it, the left arm and loop and then the loop and right arm, is a single strand of the
 * genome. Mapping the read as a whole leaves soft clipping to decide where the turn is.
 */
class NanoporeReadPlacer{
public:

	struct Pars{
		bfs::path mappingIndex_; /**< a minimap2 index of the target genome and the host */
		bfs::path targetGenomeFnp_; /**< says which contigs are the target genome's, the rest of the index is host */
		bfs::path targetBedFnp_; /**< the panel, the columns chrom, start, end, name, with the primers in the coordinates */
		uint32_t numThreads_{1};
		uint32_t minAlnMatches_{30}; /**< alignments with fewer matching bases aren't used for placing a read */
		double minPanelOverlap_{0.50}; /**< the fraction of an alignment in a target to say it's in that target */
		double minWholeUnitFrac_{0.70}; /**< the fraction of a target an alignment has to cover, and of the alignment that has to be in the target, to say it's the whole amplicon */
		uint32_t maxLocusGap_{1000}; /**< alignments of a read this close on a contig are one locus */
		uint32_t minLocusAlnLen_{50}; /**< shorter alignments don't count towards a read having more than one locus */
	};

	struct Alignment{
		uint32_t qStart_{0};
		uint32_t qEnd_{0};
		char strand_{'+'};
		std::string chrom_;
		uint64_t rStart_{0};
		uint64_t rEnd_{0};
		uint32_t matches_{0};
		uint32_t mapq_{0};
		std::vector<std::pair<uint32_t, char>> cigar_;
	};

	struct ReadPlacement{
		std::string placement_{"unplaced"}; /**< panelTarget, otherTargetGenome, host, or unplaced, where most of its aligned bases are */
		bool twoLoci_{false}; /**< aligned to more than one place */
		uint32_t wholeUnits_{0}; /**< alignments that are a whole panel amplicon */
	};

	/**@brief A fold's place in the read, from FoldbackDetector
	 */
	struct FoldInput{
		uint32_t len_{0};
		FoldbackDetector::Fold fold_;
	};

	struct FoldPlacement{
		std::string status_{"unplaced"}; /**< placed, unplaced, oneArmPlaced, armsOnDifferentContigs, armsSameStrand, armsDontOverlap */
		std::string organism_; /**< target or host */
		std::string chrom_;
		bool placed_{false}; /**< the values from here on are only for status placed */
		uint32_t mapqLeft_{0};
		uint32_t mapqRight_{0};
		char extendStrand_{'+'}; /**< the direction the 3' end was extending */
		int64_t far_{0}; /**< the furthest the 3' end got */
		int64_t entry_{0}; /**< where it primed on itself */
		int64_t reach_{0}; /**< the loop on the template */
		int64_t entryLeft_{0};
		int64_t entryRight_{0};
		int64_t queryGap_{0}; /**< read bases between the alignments of the two sides, large means the loop isn't placed */
		std::string loopSide_; /**< left, right, or neither */
		int64_t spanLow_{0};
		int64_t spanHigh_{0};
		std::string target_;
		int64_t targetStart_{0};
		int64_t targetEnd_{0};
		int64_t farToEnd_{0};
		int64_t entryToEnd_{0};
		std::string turnClass_; /**< loopNotPlaced, terminalJoin, reachedEnd, internal, or offPanel */
	};

	explicit NanoporeReadPlacer(Pars pars);

	Pars pars_;

	/**@brief Map sequences and go through the alignments as they come, a sequence with no alignments isn't reported
	 *
	 * @param seqsFnps the files of sequences to map, the index is only loaded once for all of them
	 * @param logFnp where to put minimap2's messages
	 * @param forSegments map with the settings for short pieces of a read
	 * @param onSeq called with a sequence's name, up to any information the demultiplexing added, and its alignments
	 */
	void mapSeqs(const std::vector<bfs::path> & seqsFnps, const bfs::path & logFnp, bool forSegments,
			const std::function<void(const std::string &, const std::vector<Alignment> &)> & onSeq) const;

	[[nodiscard]] ReadPlacement placeRead(const std::vector<Alignment> & alns) const;

	/**@brief Place a fold
	 *
	 * @param input the fold
	 * @param leftAlns alignments of read[:rightStart]
	 * @param rightAlns alignments of read[leftEnd:]
	 * @return where it is
	 */
	[[nodiscard]] FoldPlacement placeFold(const FoldInput & input, const std::vector<Alignment> & leftAlns, const std::vector<Alignment> & rightAlns) const;

	[[nodiscard]] static const VecStr & allPlacements();

private:
	struct PanelTarget{
		int64_t start_{0};
		int64_t end_{0};
		std::string name_;
	};
	std::unordered_map<std::string, std::vector<PanelTarget>> panel_;
	std::unordered_set<std::string> targetContigs_;
};

}  // namespace njhseq
