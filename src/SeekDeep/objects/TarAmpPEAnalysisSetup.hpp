#pragma once

/*
 * TarAmpPEAnalysisSetup.hpp
 *
 *  Created on: Nov 25, 2016
 *      Author: nick
 */

#include <bibseq.h>

namespace bibseq {

class TarAmpPEAnalysisSetup {
public:
	TarAmpPEAnalysisSetup(const bfs::path & dir);

	struct Sample {
		explicit Sample(const std::string & name);
		std::string name_;
		std::vector<std::string> reps_;

		void addRep(const std::string & rep);

		void addReps(const VecStr & reps);

		VecStr getReps() const;

	};

	struct Samples {

		Samples(const std::string & target);

		std::string target_;
		std::unordered_map<std::string, Sample> samples_;

		bool hasSample(const std::string & sample);

		void addSample(const std::string & sample);

		void addSample(const std::string & sample, const VecStr & reps);

		VecStr getSamples() const;

		std::vector<std::string> getReps() const;
	};

	bfs::path dir_;
	bfs::path infoDir_;
	bfs::path logsDir_;
	bfs::path idsDir_;
	bfs::path serverConfigsDir_;
	bfs::path reportsDir_;
	std::unordered_map<std::string, Samples> samplesForTargets_;

	std::unique_ptr<MultipleGroupMetaData> groupMetaData_;

	std::set<std::string> getSamples() const;

	std::vector<std::string> getReps() const;

	void addGroupingMetaData(const bfs::path & groupingsFileFnp);

	void addGroupingsFile()const;

	void writeSampleNamesFile() const;

	VecStr getTargets() const;

	void addSamplesNames(const bfs::path & samplesNamesFnp);

};


}  // namespace bibseq



