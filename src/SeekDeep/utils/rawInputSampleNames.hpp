#pragma once
/*
 * rawInputSampleNames.hpp
 *
 *  Created on: Oct 28, 2026
 *      Author: nick
 */

#include <njhseq/common.h>

namespace njhseq {

struct RawInputSampleNamesPars {
	bfs::path inputDir;
	bool recursive = false;
	std::string inputFilePat = ".*\\.(fastq(\\.gz)?|bam)";
	bool doNotUseBamSM = false; //name bams by file name rather than the SM of their read group
	bool prefixWithTopDir = false; //prefix names with the top directory under inputDir, e.g. run1/sub/barcode01.bam -> run1_barcode01
	std::string prefixSep = "_";
	VecStr ignoreSamples{"unclassified", "Undetermined", "undetermined"};
	bool failOnDuplicates = false;
};

/**@brief Determine the sample name each single end raw input file (fastq(.gz) or unaligned bam) would get and write a per file report
 *
 * Only bam headers are read so this is quick. Columns are file (relative to inputDir), sample, baseName, nameSource, topDir, status (ok, ignored or duplicate)
 *
 * @param pars the naming parameters
 * @param out where to write the report
 * @return the number of files with a status of duplicate, throws instead if pars.failOnDuplicates
 */
uint32_t writeRawInputSampleNamesPerFile(const RawInputSampleNamesPars & pars, std::ostream & out);

}  // namespace njhseq
