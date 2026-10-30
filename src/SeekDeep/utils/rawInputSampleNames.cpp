/*
 * rawInputSampleNames.cpp
 *
 *  Created on: Oct 28, 2026
 *      Author: nick
 */

#include "rawInputSampleNames.hpp"
#include "bamSampleNames.hpp"

namespace njhseq {

uint32_t writeRawInputSampleNamesPerFile(const RawInputSampleNamesPars & pars, std::ostream & out) {
	struct InputFileName {
		bfs::path relPath_;
		std::string topDir_;
		std::string baseName_;
		std::string nameSource_;
		std::string sample_;
		std::string status_{"ok"};
	};

	std::vector<InputFileName> inputs;
	std::regex inputFilePatReg{pars.inputFilePat};
	auto files = njh::files::listAllFilesThrowOnDupSymlink(pars.inputDir, pars.recursive, {inputFilePatReg});
	for (const auto & f : files) {
		if (f.second) {
			continue;
		}
		InputFileName input;
		input.relPath_ = f.first.lexically_relative(pars.inputDir);
		if (std::distance(input.relPath_.begin(), input.relPath_.end()) > 1) {
			input.topDir_ = input.relPath_.begin()->string();
		}
		input.baseName_ = njh::files::removeExtension(f.first.filename());
		input.nameSource_ = "filename";
		if (njh::endsWith(f.first.string(), ".bam") && !pars.doNotUseBamSM) {
			auto samples = getSampleNamesFromBamHeader(f.first);
			if (1 == samples.size()) {
				input.baseName_ = *samples.begin();
				input.nameSource_ = "SM";
			} else if (samples.size() > 1) {
				input.nameSource_ = "filename(multipleSM:" + njh::conToStr(samples, ",") + ")";
			} else {
				input.nameSource_ = "filename(noSM)";
			}
		}
		input.sample_ = pars.prefixWithTopDir && !input.topDir_.empty() ? input.topDir_ + pars.prefixSep + input.baseName_ : input.baseName_;
		if (njh::in(input.baseName_, pars.ignoreSamples)) {
			input.status_ = "ignored";
		}
		// demultiplexers write unassigned reads to e.g. bam_pass/unclassified/ with no SM and a run specific file name, so check the directories too
		for (const auto & dir : input.relPath_.parent_path()) {
			if (njh::in(dir.string(), pars.ignoreSamples)) {
				input.status_ = "ignored";
			}
		}
		inputs.emplace_back(input);
	}

	std::map<std::string, VecStr> filesPerSample;
	for (const auto & input : inputs) {
		if ("ignored" != input.status_) {
			filesPerSample[input.sample_].emplace_back(input.relPath_.string());
		}
	}
	if (filesPerSample.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "no input files matching " << pars.inputFilePat << " found in " << pars.inputDir << (pars.recursive ? " (searched recursively)" : " (try searching recursively)") << "\n";
		throw std::runtime_error{ss.str()};
	}
	std::stringstream duplicatesMessage;
	for (const auto & sampFiles : filesPerSample) {
		if (sampFiles.second.size() > 1) {
			duplicatesMessage << "\t" << sampFiles.first << ": " << njh::conToStr(sampFiles.second, ", ") << "\n";
		}
	}
	uint32_t duplicateCount = 0;
	for (auto & input : inputs) {
		if ("ignored" != input.status_ && filesPerSample[input.sample_].size() > 1) {
			input.status_ = "duplicate";
			++duplicateCount;
		}
	}

	out << "file\tsample\tbaseName\tnameSource\ttopDir\tstatus" << "\n";
	for (const auto & input : inputs) {
		out << input.relPath_.string()
				<< "\t" << input.sample_
				<< "\t" << input.baseName_
				<< "\t" << input.nameSource_
				<< "\t" << input.topDir_
				<< "\t" << input.status_ << "\n";
	}

	if (duplicateCount > 0) {
		std::stringstream ss;
		ss << "Multiple input files would be given the same sample name:\n" << duplicatesMessage.str();
		ss << "merge files split from the same sample first";
		if (!pars.prefixWithTopDir) {
			ss << ", or if they are from different runs in their own top directories prefix the names with the top directory";
		}
		if (!pars.doNotUseBamSM) {
			ss << ", or name bams by their file names rather than their SM";
		}
		ss << "\n";
		if (pars.failOnDuplicates) {
			out.flush();
			throw std::runtime_error{ss.str()};
		}
		std::cerr << "Warning, " << ss.str();
	}
	return duplicateCount;
}

}  // namespace njhseq
