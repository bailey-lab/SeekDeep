#pragma once

/*
 * ExpectedControls.hpp
 *
 *  What each control sample is expected to contain for each target, built either from the mixture tables
 *  (--sampleToMixture, --mixtureSetUp, --expectedSeqsDirFnp and optionally --strainAllelesFnp) or from an expected PMO.
 */

#include <njhseq/common.h>
#include <njhseq/objects/BioDataObject/pmo.h>
#include <njhseq/objects/dataContainers/tables/table.hpp>

namespace njhseq {

class ExpectedControls {
public:
	/**@brief one expected allele in a sample, two alleles can have the same sequence (e.g. two strains that are identical at a target)
	 */
	struct ExpectedAllele {
		std::string name_; /**< the allele name, the strain name unless the strain has several alleles at this target*/
		std::string strain_; /**< the strain the allele comes from*/
		std::string seq_;
		double frac_{0}; /**< within sample frequency*/
	};

	/**@brief the expected alleles with the same sequence combined
	 */
	struct ExpectedHap {
		std::string name_; /**< the names of the combined alleles, sorted and comma joined, e.g. FCR3,W2*/
		std::string seq_;
		double frac_{0};
	};

	struct ExpectedSample {
		std::string mix_; /**< the mixture name*/
		std::vector<ExpectedAllele> alleles_; /**< can be empty, meaning nothing is expected*/

		/**@brief the alleles with the same sequence combined, sorted by name
		 */
		[[nodiscard]] std::vector<ExpectedHap> collapsed() const;
	};

	std::map<std::string, std::map<std::string, ExpectedSample>> targets_; /**< target -> sample -> what's expected*/

	/**@brief one allele of a strain that has more than one allele at a target (e.g. a copy number variant with different sequences)
	 */
	struct StrainAllele {
		std::string allele_; /**< the name of the sequence in the target's expected seqs fasta*/
		double relativeCopies_{1}; /**< the strain's abundance is split between its alleles by this*/
	};
	using StrainAlleles = std::map<std::string, std::map<std::string, std::vector<StrainAllele>>>;/**< target -> strain -> alleles*/

	/**@brief read in the alleles for strains with more than one allele at a target
	 *
	 * @param fnp tab delimited with columns target, strain, allele and optionally relative_copies (defaults to 1)
	 */
	[[nodiscard]] static StrainAlleles readStrainAlleles(const bfs::path & fnp);

	struct FromMixturesPars {
		bfs::path samplesToMixFnp_; /**< sample and MixName columns, optionally a target column to have different samples per target*/
		bfs::path mixSetUpFnp_; /**< MixName, strain and relative_abundance columns, optionally a target column*/
		bfs::path expectedSeqsDirFnp_; /**< a fasta per target, <target>.fasta or <target>.fasta.gz, named by strain or by allele from strainAllelesFnp_*/
		bfs::path strainAllelesFnp_; /**< optional, see readStrainAlleles*/
		std::string targetNameColName_{"p_name"}; /**< the target column in the two mixture tables*/
		uint32_t trimFront_{0}; /**< bases removed from the front of every expected seq*/
		uint32_t trimBack_{0}; /**< bases removed from the back of every expected seq*/
	};

	/**@brief the targets in the mixture tables if they have a target column, otherwise empty
	 */
	[[nodiscard]] static VecStr mixtureTableTargets(const FromMixturesPars & pars);

	/**@brief the targets with a fasta file in the expected seqs directory
	 */
	[[nodiscard]] static VecStr expectedSeqsTargets(const bfs::path & expectedSeqsDirFnp);

	/**@brief build from the mixture tables
	 *
	 * Expected seqs are trimmed by trimFront_ and trimBack_, then the ones that are all Ns are dropped (along with their strain if it has no other alleles).
	 *
	 * @param pars the inputs
	 * @param targets the targets to set up, each needs a fasta in the expected seqs directory
	 */
	[[nodiscard]] static ExpectedControls fromMixtures(const FromMixturesPars & pars, const VecStr & targets);

	/**@brief build from an expected PMO
	 *
	 * Library samples are the control samples and their specimen names are the mixture names, the detected microhaplotypes are the
	 * expected alleles with reads as pseudo counts that become within sample frequencies, the allele name is the
	 * microhaplotype_name and the strain is the "strain" extra field (the allele name if missing). A sample with no results for a target
	 * isn't benchmarked for that target, one with an empty list of microhaplotypes is expected to have nothing.
	 */
	[[nodiscard]] static ExpectedControls fromPmo(const pmo::PortableMicrohaplotypeObject & pmoObj);

	struct ToPmoPars {
		std::string programName_;
		std::string programVersion_;
		std::string panelName_{"expected_controls"};
		uint32_t pseudoCountTotal_{1000000}; /**< the reads for a sample's target add up to about this*/
		std::map<std::string, std::pair<std::string, std::string>> primers_; /**< target -> forward and reverse primer, targets missing get NA and the fake_primers attribute*/
	};

	/**@brief write out as an expected PMO, one representative microhaplotype per allele name and target, one specimen per mixture
	 */
	[[nodiscard]] pmo::PortableMicrohaplotypeObject toPmo(const ToPmoPars & pars) const;

	[[nodiscard]] VecStr targetNames() const;

	/**@brief the samples that have expectations for this target, sorted
	 */
	[[nodiscard]] VecStr samples(const std::string & target) const;

	/**@brief MixName, strain and relative_abundance of the mixtures in use, from summing the alleles of each strain, with a target column (named targetNameColName) if a mixture differs between targets
	 */
	[[nodiscard]] table genMixSetUpsTab(const std::string & targetNameColName) const;

	/**@brief sample and MixName
	 */
	[[nodiscard]] table genSamplesToMixTab() const;

	/**@brief the longest expected sequence
	 */
	[[nodiscard]] uint64_t maxSeqLen() const;
};

}  // namespace njhseq
