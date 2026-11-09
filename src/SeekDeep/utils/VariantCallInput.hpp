#pragma once
//
// The haplotypes variant calling is run on per target, built from final results the same way by variantCallOnSeqAndProtein and the population clustering programs
//

#include <njhseq/common.h>
#include <njhseq/objects/Meta/MultipleGroupMetaData.hpp>
#include <njhseq/objects/dataContainers/tables/table.hpp>
#include <njhseq/objects/BioDataObject/BedRecordCore.hpp>
#include <njhseq/PopulationGeneticsUtils/PopGenUtils.hpp>
#include <njhseq/programUtils/seqSetUp.hpp>

namespace njhseq {

/**@brief The columns of a final results table, e.g. selectedClustersInfo.tab.txt.gz
 *
 */
struct VariantCallInputColumns {
	std::string sample_ = "s_Sample";
	std::string readCnt_ = "c_ReadCnt";
	std::string popHapId_ = "h_popUID";
	std::string popHapSeq_ = "h_Consensus";
	std::string target_ = "p_name";
};

/**@brief The final haplotypes of each sample per target
 *
 */
class VariantCallInput {
public:
	//key1 = target, key2 = sample, a sample's haplotypes with the same sequence are summed
	std::unordered_map<std::string, std::unordered_map<std::string, std::vector<seqInfo>>> seqs_;
	//key1 = target, key2 = haplotype sequence, value = population haplotype id
	std::unordered_map<std::string, std::unordered_map<std::string, std::string>> popHapIdBySeq_;
	std::set<std::string> targets_;
	std::set<std::string> samples_;
	uint64_t maxLen_{0};

	/**@brief Add a sample's haplotype
	 *
	 * @param target the target
	 * @param sample the sample
	 * @param readCnt the read count as it is in the table, throws if it isn't a number
	 * @param popHapId the population haplotype id
	 * @param seq the haplotype sequence
	 * @param row the whole row for the error message
	 */
	void addHap(const std::string & target, const std::string & sample, const std::string & readCnt,
			const std::string & popHapId, const std::string & seq, const VecStr & row);

	/**@brief Add every row of a table for one target
	 *
	 * @param tab the final results table
	 * @param target the target to take rows for
	 * @param cols the column names
	 */
	void addTargetRows(const table & tab, const std::string & target, const VariantCallInputColumns & cols = VariantCallInputColumns{});

	/**@brief Write a target's haplotypes, named with the sample's meta, sample, read count rounded (at least 1) and originalIdentifier (the population haplotype id)
	 *
	 * @param target the target
	 * @param fnp the gzipped fasta to write
	 * @param meta the sample meta, can be null
	 */
	void writeInputSeqs(const std::string & target, const bfs::path & fnp, const MultipleGroupMetaData * meta) const;
};

/**@brief Read the genomic locations of targets from a bed file, the name is the target
 *
 * @param bedFnp the bed file
 * @param changePeriodToDash change periods in the names to dashes
 * @return the locations by name, throws on a name being there twice
 */
std::unordered_map<std::string, std::shared_ptr<Bed6RecordCore>> readVariantCallLocations(const bfs::path & bedFnp, bool changePeriodToDash);

/**@brief What the population clustering programs need to call variants on their final haplotypes like variantCallOnSeqAndProtein
 *
 */
struct FinalHapsVariantCallPars {
	bfs::path genomicLocationsFnp; /**< targets not in here aren't variant called*/
	bool genomicLocationsChangePeriodToDash{false};

	/**@brief Set the options that the programs don't already have
	 *
	 * @param setUp the program set up
	 * @param varCallPars the variant calling parameters to set
	 */
	void setOptions(seqSetUp & setUp, CollapseAndCallVariantsPars & varCallPars);
};

/**@brief Call variants on a target's final haplotypes the same way variantCallOnSeqAndProtein does
 *
 * The haplotypes are the target's rows of the final results table, written to inputSeqsFnp and read back before calling.
 * The output goes to varCallPars.outputDirectory, identifier is set to the target and refSeqRegion to the target's location when there is one.
 *
 * @param finalResults the final results table, e.g. what's written to selectedClustersInfo.tab.txt.gz
 * @param target the target
 * @param varCallPars the variant calling parameters
 * @param location the target's genomic location, can be null
 * @param metaFnp the sample meta, can be empty
 * @param inputSeqsFnp where to write the input haplotypes
 * @return the variant calls
 */
TranslatorByAlignment::TranslatorByAlignmentResult callVariantsOnFinalHaps(const table & finalResults,
		const std::string & target,
		CollapseAndCallVariantsPars varCallPars,
		const std::shared_ptr<Bed6RecordCore> & location,
		const bfs::path & metaFnp,
		const bfs::path & inputSeqsFnp);

}  // namespace njhseq
