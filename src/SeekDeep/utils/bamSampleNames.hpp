#pragma once
/*
 * bamSampleNames.hpp
 *
 *  Created on: Oct 28, 2026
 *      Author: nick
 */

#include <njhseq/common.h>

namespace njhseq {

/**@brief Get the sample names (SM) of the read groups in a bam's header, e.g. from a nanopore dorado demux unaligned bam
 *
 * @param bamFnp the bam file
 * @return the unique SM values, empty if there are none or the file couldn't be opened
 */
std::set<std::string> getSampleNamesFromBamHeader(const bfs::path & bamFnp);

/**@brief If the bam has exactly one SM in its header set sampleName to it, warn to std::cerr if there are several
 *
 * @param bamFnp the bam file
 * @param sampleName the name to set
 * @param sampleFlag the program flag used to set the sample name, used in the warning
 * @return whether sampleName was set
 */
bool setSampleNameFromBamHeader(const bfs::path & bamFnp, std::string & sampleName, const std::string & sampleFlag);

}  // namespace njhseq
