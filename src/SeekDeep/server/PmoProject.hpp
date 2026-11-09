#pragma once
//
// SeekDeep - A library for analyzing amplicon sequence data
// Copyright (C) 2012-2026 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
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
#include <njhseq/objects/dataContainers/tables/table.hpp>
#include <njhseq/IO/OutputStream.hpp>
#include <mutex>

namespace njhseq {

/**@brief The data of a PMO arranged per target for serving, the tables use the same column names as the processClusters tables
 *
 */
class PmoProject {
public:
	/**@brief The tables for a set of samples within a target, population columns are calculated from just those samples
	 *
	 */
	struct TargetTables {
		table sampInfo_;
		table popInfo_;
		table hapIdTab_;
	};

	/**@brief Read and validate a PMO file
	 *
	 * @param pmoFnp the PMO file, can be gzipped
	 */
	explicit PmoProject(const bfs::path & pmoFnp);

	[[nodiscard]] VecStr targetNames() const;
	[[nodiscard]] bool hasTarget(const std::string & target) const;

	/**@brief Library samples with haplotypes detected for this target, sorted
	 *
	 */
	[[nodiscard]] VecStr sampleNames(const std::string & target) const;
	[[nodiscard]] bool hasSample(const std::string & target, const std::string & sample) const;

	/**@brief Meta fields that split this target's samples into groups
	 *
	 */
	[[nodiscard]] VecStr groupNames(const std::string & target) const;
	[[nodiscard]] bool hasGroup(const std::string & target, const std::string & group) const;
	[[nodiscard]] VecStr subGroupNames(const std::string & target, const std::string & group) const;
	[[nodiscard]] bool hasSubGroup(const std::string & target, const std::string & group, const std::string & subGroup) const;
	[[nodiscard]] VecStr subGroupSampleNames(const std::string & target, const std::string & group, const std::string & subGroup) const;

	const TargetTables & targetTables(const std::string & target);
	const TargetTables & subGroupTables(const std::string & target, const std::string & group, const std::string & subGroup);
	const table & groupInfo(const std::string & target, const std::string & group);

	/**@brief One row per target with its primers and a summary of what was detected
	 *
	 */
	[[nodiscard]] table targetsInfo();

	/**@brief One row per library sample, including those with nothing detected, with its meta
	 *
	 */
	[[nodiscard]] table librarySamplesInfo() const;

	/**@brief The haplotypes of a target as fasta, names are the haplotype name with its average fraction
	 *
	 */
	[[nodiscard]] std::string popSeqsFasta(const std::string & target);

	/**@brief The haplotypes of a sample in a target as fasta, names are the sample and cluster position with its fraction
	 *
	 */
	[[nodiscard]] std::string sampleSeqsFasta(const std::string & target, const std::string & sample) const;

	[[nodiscard]] VecStr librarySampleNames() const;
	[[nodiscard]] bool hasLibrarySample(const std::string & sample) const;

	/**@brief One sample across all targets, one row per haplotype with how common it is in that target's population
	 *
	 */
	[[nodiscard]] table sampleAcrossTargets(const std::string & sample);

	/**@brief The targets a sample has no haplotypes for
	 *
	 */
	[[nodiscard]] VecStr targetsNotDetected(const std::string & sample) const;

	/**@brief The fields of the specimen for a library sample, including its name
	 *
	 */
	[[nodiscard]] std::map<std::string, std::string> sampleMeta(const std::string & sample) const;

	/**@brief The haplotypes of a sample across all targets as fasta, names are the target and cluster position with its fraction
	 *
	 */
	[[nodiscard]] std::string sampleAcrossTargetsFasta(const std::string & sample) const;

	/**@brief Meta fields that split the library samples into groups, the same fields as groupNames() but across all targets
	 *
	 */
	[[nodiscard]] VecStr metaGroupNames() const;
	[[nodiscard]] VecStr metaSubGroupNames(const std::string & group) const;
	[[nodiscard]] bool hasMetaSubGroup(const std::string & group, const std::string & subGroup) const;
	[[nodiscard]] VecStr metaSubGroupSampleNames(const std::string & group, const std::string & subGroup) const;

	/**@brief The data for a haplotype rainbow of a set of samples across all targets
	 *
	 * Haplotypes are ranked within each target by how many of the samples have them, rank 0 being the most common.
	 * calls is flattened into groups of four: sample position, target position, haplotype rank, reads.
	 * Samples with nothing detected in any target are listed in samplesWithoutData instead of samples.
	 *
	 * @param group the meta field to take samples from, empty for all library samples
	 * @param subGroup the value of that field
	 * @return the rainbow data, computed once per sub group
	 */
	const Json::Value & haplotypeRainbow(const std::string & group = "", const std::string & subGroup = "");

	/**@brief Whether any haplotype has amino acid calls (associated_protein_variants)
	 *
	 */
	[[nodiscard]] bool hasAminoAcidCalls() const;

	/**@brief Every amino acid seen at each typed protein position in a set of samples, calculated like AAChangesInfo from variantCallOnSeqAndProtein
	 *
	 * Like combining the per target vcfs, a position uses the target typed in the most of these samples, samples without it use the target with the most reads for them.
	 * AlleleCount and AlleleTotal count sample haplotypes, AlleleFrequency is the mean over the samples typed at the position of the amino acid's within sample read fraction,
	 * SampleCount is the samples with the amino acid and SamplePrevalence is that over SampleTotal, the samples typed at the position.
	 * Every position has a row for the reference amino acid even when no sample has it.
	 * There's a Gene column with the alternative gene names when any variant has one.
	 *
	 * @param group the meta field to take samples from, empty for all library samples
	 * @param subGroup the value of that field
	 * @return the table, computed once per sub group
	 */
	const table & aminoAcidChanges(const std::string & group = "", const std::string & subGroup = "");

	Json::Value header_; /**< the pmo_header of the file */

private:
	struct DetectedHap {
		uint32_t hapId_;
		uint32_t reads_;
	};
	struct AminoAcidCall {
		std::string transcript_;
		std::string gene_;
		std::string geneName_; /**< alternative_gene_name, empty without one*/
		uint32_t pos_; /**< 1-based*/
		std::string refAA_;
		std::string aa_;
	};
	struct Target {
		std::string name_;
		std::string forwardPrimer_;
		std::string reversePrimer_;
		bool fakePrimers_{false};
		VecStr hapNames_;
		VecStr hapSeqs_;
		std::vector<std::vector<AminoAcidCall>> hapAminoAcids_; /**< by haplotype id*/
		std::map<std::string, std::vector<DetectedHap>> detected_; /**< by sample, most reads first*/
	};
	bool hasAminoAcidCalls_{false};
	bool hasGeneNames_{false};

	std::map<std::string, Target> targets_;
	VecStr librarySamples_;
	std::map<std::string, std::string> specimenNames_;
	std::map<std::string, std::map<std::string, std::string>> meta_; /**< by sample then by field*/
	VecStr metaGroups_;

	std::mutex tablesMut_;
	std::map<std::string, TargetTables> tablesByTarget_;
	std::map<std::string, std::map<std::string, std::map<std::string, TargetTables>>> tablesBySubGroup_;
	std::map<std::string, std::map<std::string, table>> groupInfoByTarget_;
	std::map<std::pair<std::string, std::string>, Json::Value> rainbows_; /**< by group and sub group*/
	std::map<std::pair<std::string, std::string>, table> aminoAcidChanges_; /**< by group and sub group*/

	[[nodiscard]] Json::Value genHaplotypeRainbow(const VecStr & samples) const;
	[[nodiscard]] table genAminoAcidChanges(const VecStr & samples) const;

	[[nodiscard]] std::string metaValue(const std::string & sample, const std::string & group) const;
	[[nodiscard]] TargetTables genTables(const Target & target, const VecStr & samples) const;
	void computeGroup(const std::string & target, const std::string & group);
};

}  // namespace njhseq
