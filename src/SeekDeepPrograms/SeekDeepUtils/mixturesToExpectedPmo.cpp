/*
 * mixturesToExpectedPmo.cpp
 *
 *  Turn the control mixture tables into an expected PMO for benchmarkMultiTarAmpControlMixtures --expectedPmoFnp
 */

#include "SeekDeepUtilsRunner.hpp"
#include "SeekDeep/objects/ControlBenchmarking/ExpectedControls.hpp"
#include "SeekDeep/objects/TarAmpSetupUtils/PrimersAndMids.hpp"

namespace njhseq {

int SeekDeepUtilsRunner::mixturesToExpectedPmo(const njh::progutils::CmdArgs & inputCommands) {
	ExpectedControls::FromMixturesPars mixPars;
	ExpectedControls::ToPmoPars pmoPars;
	pmoPars.programName_ = "SeekDeep mixturesToExpectedPmo";
	pmoPars.programVersion_ = processClustersPars{}.pmoProgramVersion;
	bfs::path primersFnp;
	std::set<std::string> selectTargets;
	std::set<std::string> excludeTargets;
	OutOptions outOpts(bfs::path("expectedControls"), ".json.gz");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.description_ = "Make an expected PMO of what each control sample should have for benchmarkMultiTarAmpControlMixtures --expectedPmoFnp. "
											 "Library samples are the control samples and their specimen is their mixture. Each strain's allele is a representative microhaplotype named by the allele with the strain as an extra field. "
											 "A sample's reads are pseudo counts that add up to about --pseudoCountTotal, split by the mixture's relative abundances and then by a strain's alleles from --strainAllelesFnp. "
											 "Identical alleles are kept separate here, the benchmark combines them within each sample.";
	setUp.setOption(mixPars.samplesToMixFnp_, "--sampleToMixture",
									"Sample To Mixture, 2 columns 1)sample, 2)MixName, can have a target column (--targetNameColName) along with --mixtureSetUp to have different samples and mixtures per target", true);
	setUp.setOption(mixPars.mixSetUpFnp_, "--mixtureSetUp",
									"Mixture Set Up, 3 columns 1)MixName, 2)strain, 3)relative_abundance", true);
	setUp.setOption(mixPars.expectedSeqsDirFnp_, "--expectedSeqsDirFnp",
									"Expected Seqs directory of fasta files, one per target named <target>.fasta or <target>.fasta.gz, sequences named by strain or by allele from --strainAllelesFnp", true);
	setUp.setOption(mixPars.strainAllelesFnp_, "--strainAllelesFnp",
									"Optional, for strains with more than one allele at a target (e.g. a copy number variant with different sequences), tab delimited with columns target, strain, allele and optionally relative_copies (default 1)");
	setUp.setOption(mixPars.targetNameColName_, "--targetNameColName", "the target column in the mixture tables if they're split by target");
	setUp.setOption(mixPars.trimFront_, "--trimFront", "trim this many bases from the front of every expected seq, e.g. to match a pipeline that trims the edges of the insert");
	setUp.setOption(mixPars.trimBack_, "--trimBack", "trim this many bases from the back of every expected seq, e.g. to match a pipeline that trims the edges of the insert");
	setUp.setOption(selectTargets, "--selectTargets", "only these targets");
	setUp.setOption(excludeTargets, "--excludeTargets", "exclude these targets");
	setUp.setOption(primersFnp, "--primers",
									"Optional, SeekDeep id file with the targets' primers, targets without primers get NA primers and the fake_primers attribute");
	setUp.setOption(pmoPars.pseudoCountTotal_, "--pseudoCountTotal", "the pseudo counts for a sample's target add up to about this, the benchmark turns them back into frequencies");
	setUp.setOption(pmoPars.panelName_, "--panelName", "panel name");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	njh::files::checkExistenceThrow({mixPars.samplesToMixFnp_, mixPars.mixSetUpFnp_, mixPars.expectedSeqsDirFnp_}, __PRETTY_FUNCTION__);
	if (!mixPars.strainAllelesFnp_.empty()) {
		njh::files::checkExistenceThrow(std::vector<bfs::path>{mixPars.strainAllelesFnp_}, __PRETTY_FUNCTION__);
	}

	auto passTargetFilter = [&selectTargets, &excludeTargets](const std::string & target) {
		return (selectTargets.empty() || njh::in(target, selectTargets)) && njh::notIn(target, excludeTargets);
	};
	auto candidateTargets = ExpectedControls::mixtureTableTargets(mixPars);
	if (candidateTargets.empty()) {
		candidateTargets = ExpectedControls::expectedSeqsTargets(mixPars.expectedSeqsDirFnp_);
	}
	VecStr targets;
	for (const auto & target : candidateTargets) {
		if (passTargetFilter(target)) {
			targets.emplace_back(target);
		}
	}
	if (targets.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "no targets found in " << mixPars.expectedSeqsDirFnp_
			 << " (or the mixture tables), or none left after --selectTargets/--excludeTargets" << "\n";
		throw std::runtime_error{ss.str()};
	}

	if (!primersFnp.empty()) {
		PrimersAndMids ids(primersFnp);
		VecStr noPrimers;
		for (const auto & target : targets) {
			if (!ids.hasTarget(target)) {
				noPrimers.emplace_back(target);
				continue;
			}
			const auto & info = ids.targets_.at(target).info_;
			pmoPars.primers_[target] = {info.forwardPrimerRaw_, info.reversePrimerRaw_};
		}
		if (!noPrimers.empty()) {
			std::cerr << __PRETTY_FUNCTION__ << ", warning, " << "no primers in " << primersFnp
								<< " for the following targets, they'll get NA primers and the fake_primers attribute: "
								<< njh::conToStr(noPrimers, ",") << std::endl;
		}
	}

	auto expected = ExpectedControls::fromMixtures(mixPars, targets);
	pmo::writePmo(expected.toPmo(pmoPars), outOpts);
	if (setUp.pars_.verbose_) {
		std::cout << "wrote " << outOpts.outName() << " with " << expected.targetNames().size() << " targets and "
							<< expected.genSamplesToMixTab().nRow() << " samples" << std::endl;
	}
	return 0;
}

}  // namespace njhseq
