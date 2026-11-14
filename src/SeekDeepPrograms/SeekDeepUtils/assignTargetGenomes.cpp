//
// Write which genome each target is variant called against
//

#include "SeekDeepUtilsRunner.hpp"
#include <njhseq/GenomeUtils/TargetGenomeAssignment.hpp>

namespace njhseq {

int SeekDeepUtilsRunner::assignTargetGenomes(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path locationsDir;
	TargetGenomeChoicePars pars;
	pars.onlyPreferred_ = true;
	njh::files::MkdirPar outDirPars("assignedTargetGenomes");
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.description_ = "Assign each target to the first of the given genomes SeekDeep genTargetInfoFromGenomes found it in, writes targetGenomes.tsv "
			"and a <genome>_variantCall.bed per genome with the inner regions of its targets (and their -R1/-R2 mates) for --genomicLocations, "
			"targets in none of the genomes aren't assigned";
	setUp.examples_.emplace_back("MASTERPROGRAM SUBPROGRAM --locations genome_extraction/locationsByGenome --genomes Pf3D7,PvP01,Poc221 --dout assignedTargetGenomes");
	setUp.setOption(locationsDir, "--locations", "The locationsByGenome directory of SeekDeep genTargetInfoFromGenomes", true);
	setUp.setOption(pars.preferredGenomes_, "--genomes", "The genomes to assign targets to, in order, a target goes to the first it's in", true);
	setUp.setOption(outDirPars.dirName_, "--dout", "Output directory");
	setUp.setOption(outDirPars.overWriteDir_, "--overWriteDir", "Overwrite the output directory");
	setUp.finishSetUp(std::cout);

	njh::files::makeDir(outDirPars);
	const auto locations = TargetGenomeLocations::read(locationsDir);
	const auto chosen = chooseTargetGenomes(locations, pars);

	OutputStream assignedOut(njh::files::make_path(outDirPars.dirName_, "targetGenomes.tsv"));
	assignedOut << "target\tgenome" << "\n";
	for (const auto & target : chosen) {
		assignedOut << target.first << "\t" << target.second.genome_ << "\n";
	}
	for (const auto & genome : pars.preferredGenomes_) {
		OutputStream bedOut(njh::files::make_path(outDirPars.dirName_, genome + "_variantCall.bed"));
		uint32_t regions = 0;
		auto regionsSearch = locations.innerRegionsByGenome_.find(genome);
		if (locations.innerRegionsByGenome_.end() != regionsSearch) {
			//the inner bed's order, the regions at each target's chosen hit
			for (const auto & region : regionsSearch->second) {
				auto chosenSearch = chosen.find(locations.baseTarget(region->name_));
				if (chosen.end() != chosenSearch && genome == chosenSearch->second.genome_
						&& region == locations.innerRegion(region->name_, chosenSearch->second)) {
					bedOut << region->toDelimStrWithExtra() << "\n";
					++regions;
				}
			}
		}
		if (setUp.pars_.verbose_) {
			std::cout << genome << "\t" << regions << " regions" << std::endl;
		}
	}
	return 0;
}

}  // namespace njhseq
