#pragma once

/*
 * TarAmpPEAnalysisSetup.hpp
 *
 *  Created on: Nov 25, 2016
 *      Author: nick
 */

#include <bibseq.h>
#include "SeekDeep/objects/PrimersAndMids.hpp"

namespace bibseq {

class TarAmpAnalysisSetup {
public:
	struct TarAmpPars{
		bfs::path samplesNamesFnp = "";
		bfs::path outDir = "";
		bfs::path inputDir = "";
		bfs::path groupMeta = "";
		bfs::path idFile = "";
		bfs::path lenCutOffsFnp = "";
		bfs::path refSeqsDir = "";
		bfs::path targetsToIndexFnp = "";
		bool byIndex = false;

		uint32_t numThreads = 1;

		std::string technology = "illumina";

		std::string inputFilePat = ".*.fastq.gz";

		//Illumina specific
		//paired end specific
		uint32_t maxOverlap = 250;
		uint32_t r1Trim = std::numeric_limits<uint32_t>::max();
		uint32_t r2Trim = std::numeric_limits<uint32_t>::max();
		//
		bool noQualTrim = false;

		std::string extraExtractorCmds = "";
		std::string extraQlusterCmds = "";
		std::string extraProcessClusterCmds = "";

		std::string stitcherCmd= "flash";
#if defined( __APPLE__ ) || defined( __APPLE_CC__ ) || defined( macintosh ) || defined( __MACH__ )
		//apple zcat is stupid and requires files end with .Z because why not
		//using brew install gnutls instead
		std::string zcatCmd = "gzcat";
#else
		std::string zcatCmd = "zcat";
#endif


		//checks
		bool checkForStitcher(VecStr & warnings) const;
		bool checkForZcat(VecStr & warnings) const;

		bool checkForOutDir(VecStr & warnings) const;

		bool checkForRequiredFnpPars(VecStr & warnings) const;
		bool checkForOptionalFnpPars(VecStr & warnings) const;

		bool allChecks(VecStr & warnings) const;

		static bool checkIfFnpExists(const bfs::path & fnp, VecStr & warnings);

		bool techIs454()const;
		bool techIsIllumina()const;
		bool techIsIonTorrent()const;
/*
 * pars.technology != "454" && pars.technology != "iontorrent" && pars.technology != "illumina"
 */
	};

	TarAmpAnalysisSetup(const TarAmpPars & pars);

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

	struct TargetsInfoAgreement{

		VecStr missing_; /**< no info provided */
		VecStr notMatching_; /**< info found but not in analysis*/
	};

	TarAmpPars pars_;

	//directories
	bfs::path dir_;
	bfs::path infoDir_;
	bfs::path logsDir_;
	bfs::path idsDir_;
	bfs::path serverConfigsDir_;
	bfs::path reportsDir_;

	std::unordered_map<std::string, Samples> samples_;

	std::unique_ptr<MultipleGroupMetaData> groupMetaData_;
	std::unique_ptr<PrimersAndMids> idsMids_;

	std::unordered_map<std::string, VecStr> indexToTars_;

	TargetsInfoAgreement forRefSeqs_;
	TargetsInfoAgreement forLenCutOffs_;

	std::set<std::string> getSamples() const;

	VecStr getReps() const;

	void addGroupingMetaData(const bfs::path & groupingsFileFnp);

	void addGroupingsFile() const;

	void writeSampleNamesFile() const;

	VecStr getTargets() const;
	VecStr getIndexes() const;

	void addSamplesNames(const bfs::path & samplesNamesFnp);

	void addIndexToTargetsNames(const bfs::path & targetsToIndexFnp);

	void addRefSeqs(const bfs::path & refSeqsDir);

	void addLenCutOffs(const bfs::path & lenCutOffsFnp);



	void writeOutIdFiles() const;

	std::vector<VecStr> getTarCombos() const;

	VecStr getExpectantInputNames() const;

	void setUpPopClusteringDirs(bool verbose = false) const;

};


}  // namespace bibseq



