/*
 * ExpectedControls.cpp
 *
 */

#include "ExpectedControls.hpp"
#include "SeekDeep/objects/ControlBenchmarking/ControlBencher.hpp"

#include <njhseq/objects/dataContainers/tables.h>
#include <njhseq/IO/SeqIO/SeqInput.hpp>

#include <regex>

namespace njhseq {

namespace {
const std::string MIX_BY_TARGET_FIELD = "mix_by_target";
const std::string MIX_NAME_FIELD = "mix_name";
const std::string STRAIN_FIELD = "strain";

//the characters the PMO schema allows in each kind of name
const std::regex SPECIMEN_NAME_CHAR{R"([A-z\-._0-9(),/ ])"};
const std::regex LIBRARY_SAMPLE_NAME_PAT{R"(^[A-z\-._0-9 ]+$)"};
const std::regex TARGET_NAME_PAT{R"(^[A-z\-._0-9]+$)"};

//the name with every character the PMO schema doesn't allow in a specimen name replaced by _
std::string toSpecimenName(const std::string & name) {
	std::string ret = name;
	for (auto & letter : ret) {
		if (!std::regex_match(std::string(1, letter), SPECIMEN_NAME_CHAR)) {
			letter = '_';
		}
	}
	return ret;
}

//throw listing the names that don't match the PMO schema's pattern
void checkNamesThrow(const std::set<std::string> & names, const std::regex & pattern, const std::string & what,
										 const std::string & allowed, const std::string & funcName) {
	VecStr bad;
	for (const auto & name : names) {
		if (!std::regex_match(name, pattern)) {
			bad.emplace_back(name);
		}
	}
	if (!bad.empty()) {
		std::stringstream ss;
		ss << funcName << ", error " << "the following " << what << " can't be used in a PMO, they can only have "
			 << allowed << ": " << njh::conToStr(bad, ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
}

bfs::path targetFastaFnp(const bfs::path & dir, const std::string & target) {
	auto fnp = njh::files::make_path(dir, target + ".fasta");
	if (!exists(fnp)) {
		fnp = njh::files::make_path(dir, target + ".fasta.gz");
	}
	if (!exists(fnp)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "could not find expected seqs file for target: " << target
			 << ", should be " << njh::files::make_path(dir, target + ".fasta") << " or "
			 << njh::files::make_path(dir, target + ".fasta.gz") << "\n";
		throw std::runtime_error{ss.str()};
	}
	return fnp;
}

bool allNs(const std::string & seq) {
	return std::all_of(seq.begin(), seq.end(), [](const char base) {
		return 'N' == base || 'n' == base;
	});
}
}  // namespace

std::vector<ExpectedControls::ExpectedHap> ExpectedControls::ExpectedSample::collapsed() const {
	std::vector<ExpectedHap> ret;
	std::vector<VecStr> names;
	for (const auto & allele : alleles_) {
		bool found = false;
		for (const auto pos : iter::range(ret.size())) {
			if (ret[pos].seq_ == allele.seq_) {
				ret[pos].frac_ += allele.frac_;
				names[pos].emplace_back(allele.name_);
				found = true;
				break;
			}
		}
		if (!found) {
			ret.emplace_back(ExpectedHap{"", allele.seq_, allele.frac_});
			names.emplace_back(VecStr{allele.name_});
		}
	}
	for (const auto pos : iter::range(ret.size())) {
		njh::sort(names[pos]);
		ret[pos].name_ = njh::conToStr(names[pos], ",");
	}
	std::sort(ret.begin(), ret.end(), [](const ExpectedHap & first, const ExpectedHap & second) {
		return first.name_ < second.name_;
	});
	return ret;
}

ExpectedControls::StrainAlleles ExpectedControls::readStrainAlleles(const bfs::path & fnp) {
	table tab(fnp, "\t", true);
	tab.checkForColumnsThrow(VecStr{"target", "strain", "allele"}, __PRETTY_FUNCTION__);
	bool hasCopies = njh::in(std::string("copies"), tab.columnNames_);
	if (!hasCopies && njh::in(std::string("relative_copies"), tab.columnNames_)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << fnp << " has a relative_copies column, the column should be copies, "
			 << "the number of copies of each allele per genome" << "\n";
		throw std::runtime_error{ss.str()};
	}
	StrainAlleles ret;
	for (const auto & row : tab) {
		const auto & target = row[tab.getColPos("target")];
		const auto & strain = row[tab.getColPos("strain")];
		StrainAllele allele;
		allele.allele_ = row[tab.getColPos("allele")];
		if (hasCopies) {
			allele.copies_ = njh::StrToNumConverter::stoToNum<double>(row[tab.getColPos("copies")]);
		}
		if (allele.copies_ <= 0) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "copies needs to be greater than 0, " << target << " " << strain
				 << " " << allele.allele_ << ": " << allele.copies_ << "\n";
			throw std::runtime_error{ss.str()};
		}
		auto & alleles = ret[target][strain];
		for (const auto & other : alleles) {
			if (other.allele_ == allele.allele_) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "already have allele " << allele.allele_ << " for strain " << strain
					 << " for target " << target << " in " << fnp << "\n";
				throw std::runtime_error{ss.str()};
			}
		}
		alleles.emplace_back(allele);
	}
	return ret;
}

VecStr ExpectedControls::mixtureTableTargets(const FromMixturesPars & pars) {
	table sampleToMixTab(pars.samplesToMixFnp_, "\t", true);
	table mixSetUpTab(pars.mixSetUpFnp_, "\t", true);
	if (!njh::in(pars.targetNameColName_, sampleToMixTab.columnNames_) ||
			!njh::in(pars.targetNameColName_, mixSetUpTab.columnNames_)) {
		return {};
	}
	auto sampleToMixTargets = njh::getVecOfMapKeys(sampleToMixTab.splitTableOnColumn(pars.targetNameColName_));
	auto mixSetUpTargets = njh::getVecOfMapKeys(mixSetUpTab.splitTableOnColumn(pars.targetNameColName_));
	VecStr inSampToMix;
	VecStr inMixSetup;
	VecStr shared;
	njh::decompose_sets(sampleToMixTargets.begin(), sampleToMixTargets.end(),
											mixSetUpTargets.begin(), mixSetUpTargets.end(),
											std::back_inserter(inSampToMix),
											std::back_inserter(inMixSetup),
											std::back_inserter(shared));
	if (!inSampToMix.empty() || !inMixSetup.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << " found the following errors" << "\n";
		if (!inSampToMix.empty()) {
			ss << "found the following targets in " << pars.samplesToMixFnp_ << " but not in " << pars.mixSetUpFnp_ << "\n"
				 << njh::conToStr(inSampToMix, ",") << "\n";
		}
		if (!inMixSetup.empty()) {
			ss << "found the following targets in " << pars.mixSetUpFnp_ << " but not in " << pars.samplesToMixFnp_ << "\n"
				 << njh::conToStr(inMixSetup, ",") << "\n";
		}
		throw std::runtime_error{ss.str()};
	}
	njh::naturalSortNameSet(shared);
	return shared;
}

VecStr ExpectedControls::expectedSeqsTargets(const bfs::path & expectedSeqsDirFnp) {
	std::set<std::string> targets;
	for (const auto & entry : bfs::directory_iterator(expectedSeqsDirFnp)) {
		auto name = entry.path().filename().string();
		for (const auto & ext : VecStr{".fasta.gz", ".fasta"}) {
			if (njh::endsWith(name, ext)) {
				targets.emplace(name.substr(0, name.size() - ext.size()));
				break;
			}
		}
	}
	VecStr ret(targets.begin(), targets.end());
	njh::naturalSortNameSet(ret);
	return ret;
}

ExpectedControls ExpectedControls::fromMixtures(const FromMixturesPars & pars, const VecStr & targets) {
	njh::files::checkExistenceThrow({pars.samplesToMixFnp_, pars.mixSetUpFnp_, pars.expectedSeqsDirFnp_}, __PRETTY_FUNCTION__);
	StrainAlleles strainAlleles;
	if (!pars.strainAllelesFnp_.empty()) {
		strainAlleles = readStrainAlleles(pars.strainAllelesFnp_);
	}
	table sampleToMixTab(pars.samplesToMixFnp_, "\t", true);
	table mixSetUpTab(pars.mixSetUpFnp_, "\t", true);
	const bool perTarget = !mixtureTableTargets(pars).empty();
	std::map<std::string, table> sampleToMixSplitTab;
	std::map<std::string, table> mixSetUpSplitTab;
	std::unique_ptr<ControlBencher> allTargetsBencher;
	if (perTarget) {
		sampleToMixSplitTab = sampleToMixTab.splitTableOnColumn(pars.targetNameColName_);
		mixSetUpSplitTab = mixSetUpTab.splitTableOnColumn(pars.targetNameColName_);
	} else {
		allTargetsBencher = std::make_unique<ControlBencher>(mixSetUpTab, sampleToMixTab);
	}

	ExpectedControls ret;
	for (const auto & target : targets) {
		std::unique_ptr<ControlBencher> bencher;
		if (perTarget) {
			if (!njh::in(target, mixSetUpSplitTab)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "target " << target << " isn't in " << pars.mixSetUpFnp_ << " or "
					 << pars.samplesToMixFnp_ << "\n";
				throw std::runtime_error{ss.str()};
			}
			bencher = std::make_unique<ControlBencher>(mixSetUpSplitTab.at(target), sampleToMixSplitTab.at(target));
		} else {
			bencher = std::make_unique<ControlBencher>(*allTargetsBencher);
		}

		//read in expected seqs
		auto expectedSeqsFnp = targetFastaFnp(pars.expectedSeqsDirFnp_, target);
		SeqInput reader{SeqIOOptions(expectedSeqsFnp, SeqIOOptions::getInFormatFromFnp(expectedSeqsFnp))};
		reader.openIn();
		seqInfo seq;
		std::unordered_map<std::string, std::string> seqByName;
		std::set<std::string> allNNames;
		while (reader.readNextRead(seq)) {
			if (njh::in(seq.name_, seqByName) || njh::in(seq.name_, allNNames)) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << " already have expected sequence named " << seq.name_ << " in "
					 << expectedSeqsFnp << "\n";
				throw std::runtime_error{ss.str()};
			}
			if (allNs(seq.seq_)) {
				allNNames.emplace(seq.name_);
				continue;
			}
			if (static_cast<uint64_t>(pars.trimFront_) + pars.trimBack_ >= seq.seq_.size()) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << "can't trim " << pars.trimFront_ << " from the front and "
					 << pars.trimBack_ << " from the back of " << seq.name_ << " in " << expectedSeqsFnp
					 << ", it's only " << seq.seq_.size() << " long" << "\n";
				throw std::runtime_error{ss.str()};
			}
			auto trimmed = seq.seq_.substr(pars.trimFront_, seq.seq_.size() - pars.trimFront_ - pars.trimBack_);
			if (allNs(trimmed)) {
				allNNames.emplace(seq.name_);
			} else {
				seqByName[seq.name_] = trimmed;
			}
		}

		//work out the alleles of each strain, dropping all N seqs and strains left without any alleles
		const auto targetStrainAlleles = njh::in(target, strainAlleles) ? strainAlleles.at(target)
																																		: std::map<std::string, std::vector<StrainAllele>>{};
		std::map<std::string, std::vector<StrainAllele>> allelesByStrain;
		VecStr removeStrains;
		VecStr missing;
		for (const auto & strain : bencher->getAllStrains()) {
			auto strainSearch = targetStrainAlleles.find(strain);
			std::vector<StrainAllele> alleles = targetStrainAlleles.end() == strainSearch
																							? std::vector<StrainAllele>{StrainAllele{strain, 1}}
																							: strainSearch->second;
			std::vector<StrainAllele> kept;
			for (const auto & allele : alleles) {
				if (njh::in(allele.allele_, seqByName)) {
					kept.emplace_back(allele);
				} else if (!njh::in(allele.allele_, allNNames)) {
					missing.emplace_back(allele.allele_);
				}
			}
			if (kept.empty()) {
				removeStrains.emplace_back(strain);
			} else {
				allelesByStrain[strain] = kept;
			}
		}
		if (!missing.empty()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error for: " << target << "\n"
				 << "the following strains are missing from expected seqs: "
				 << njh::conToStrEndSpecial(missing, ", ", " and ") << "\n";
			throw std::runtime_error{ss.str()};
		}
		//warn about mixtures left with nothing expected once strains whose alleles are all Ns are removed
		for (const auto & strain : removeStrains) {
			for (const auto & mix : bencher->mixSetups_) {
				if (njh::in(strain, mix.second.relativeAbundances_) && 1 == mix.second.relativeAbundances_.size()) {
					std::cerr << __PRETTY_FUNCTION__ << ", warning, " << "for " << target << " mixture " << mix.first
										<< " only had " << strain << " which is all Ns, so nothing is expected for it" << std::endl;
				}
			}
		}
		bencher->removeStrains(removeStrains);

		for (const auto & sample : bencher->getSamples()) {
			ExpectedSample expected;
			expected.mix_ = bencher->samplesToMix_.at(sample);
			//each allele's template is its strain's relative abundance times its copies, then these are made to sum to 1
			double totalTemplate = 0;
			for (const auto & strainFrac : bencher->mixSetups_.at(expected.mix_).relativeAbundances_) {
				for (const auto & allele : allelesByStrain.at(strainFrac.first)) {
					expected.alleles_.emplace_back(ExpectedAllele{allele.allele_, strainFrac.first, seqByName.at(allele.allele_),
																												strainFrac.second * allele.copies_});
					totalTemplate += strainFrac.second * allele.copies_;
				}
			}
			for (auto & allele : expected.alleles_) {
				allele.frac_ /= totalTemplate;
			}
			std::sort(expected.alleles_.begin(), expected.alleles_.end(), [](const ExpectedAllele & first, const ExpectedAllele & second) {
				return first.name_ < second.name_;
			});
			ret.targets_[target][sample] = expected;
		}
	}
	return ret;
}

ExpectedControls ExpectedControls::fromPmo(const pmo::PortableMicrohaplotypeObject & pmoObj) {
	ExpectedControls ret;
	//a target can be in representative_microhaplotypes more than once, mhaps_target_id is the position in that list and mhap_id the position within that entry
	std::vector<std::string> targetNameByMhapsId;
	for (const auto & representatives : pmoObj.representative_microhaplotypes_.targets_) {
		targetNameByMhapsId.emplace_back(pmoObj.target_info_.at(representatives.target_id_).target_name_);
	}
	for (const auto & run : pmoObj.detected_microhaplotypes_) {
		for (const auto & sampleResults : run.library_samples_) {
			const auto & library = pmoObj.library_sample_info_.at(sampleResults.library_sample_id_);
			const auto & sample = library.library_sample_name_;
			//the mixture is the specimen's mix_name field if it has one, otherwise its name
			const auto & specimen = pmoObj.specimen_info_.at(library.specimen_id_);
			std::string specimenName = specimen.specimen_name_;
			auto mixName = specimen.extras_.find(MIX_NAME_FIELD);
			if (specimen.extras_.end() != mixName && mixName->second.is_string()) {
				specimenName = mixName->second.get<std::string>();
			}
			auto mixByTarget = library.extras_.find(MIX_BY_TARGET_FIELD);
			for (const auto & targetResults : sampleResults.target_results_) {
				const auto & target = targetNameByMhapsId.at(targetResults.mhaps_target_id_);
				if (njh::in(sample, ret.targets_[target])) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << sample << " has more than one set of expected results for " << target << "\n";
					throw std::runtime_error{ss.str()};
				}
				ExpectedSample expected;
				expected.mix_ = specimenName;
				if (library.extras_.end() != mixByTarget && mixByTarget->second.contains(target)) {
					expected.mix_ = mixByTarget->second.at(target).get<std::string>();
				}
				const auto & representatives = pmoObj.representative_microhaplotypes_.targets_.at(targetResults.mhaps_target_id_);
				double totalReads = 0;
				for (const auto & mhap : targetResults.mhaps_) {
					totalReads += mhap.reads_;
				}
				if (!targetResults.mhaps_.empty() && 0 == totalReads) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << sample << " has 0 reads for all its expected alleles for " << target << "\n";
					throw std::runtime_error{ss.str()};
				}
				std::map<std::string, ExpectedAllele> allelesByName;
				for (const auto & mhap : targetResults.mhaps_) {
					const auto & rep = representatives.microhaplotypes_.at(mhap.mhap_id_);
					ExpectedAllele allele;
					allele.name_ = rep.microhaplotype_name_.has_value()
												 ? rep.microhaplotype_name_.value()
												 : njh::pasteAsStr(target, ".", targetResults.mhaps_target_id_, ".", mhap.mhap_id_);
					allele.strain_ = allele.name_;
					auto strain = rep.extras_.find(STRAIN_FIELD);
					if (rep.extras_.end() != strain && strain->second.is_string()) {
						allele.strain_ = strain->second.get<std::string>();
					}
					allele.seq_ = rep.seq_;
					allele.frac_ = mhap.reads_ / totalReads;
					auto search = allelesByName.find(allele.name_);
					if (allelesByName.end() == search) {
						allelesByName.emplace(allele.name_, allele);
					} else if (search->second.seq_ != allele.seq_) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << sample << " has two different sequences named " << allele.name_
							 << " for " << target << "\n";
						throw std::runtime_error{ss.str()};
					} else {
						search->second.frac_ += allele.frac_;
					}
				}
				for (const auto & allele : allelesByName) {
					expected.alleles_.emplace_back(allele.second);
				}
				ret.targets_[target][sample] = expected;
			}
		}
	}
	return ret;
}

pmo::PortableMicrohaplotypeObject ExpectedControls::toPmo(const ToPmoPars & pars) const {
	pmo::PortableMicrohaplotypeObject ret;
	auto today = njh::getCurrentDate();
	ret.pmo_header_.pmo_version_ = pmo::PMO_SCHEMA_VERSION;
	ret.pmo_header_.creation_date_ = today.substr(0, today.find('_'));
	pmo::PmoGenerationMethod generationMethod;
	generationMethod.program_name_ = pars.programName_;
	generationMethod.program_version_ = pars.programVersion_;
	ret.pmo_header_.generation_method_ = generationMethod;

	auto targets = targetNames();
	std::set<std::string> allSamples;
	for (const auto & target : targets_) {
		for (const auto & sample : target.second) {
			allSamples.emplace(sample.first);
		}
	}
	checkNamesThrow(std::set<std::string>(targets.begin(), targets.end()), TARGET_NAME_PAT, "target names",
									"letters, digits, - . and _", __PRETTY_FUNCTION__);
	checkNamesThrow(std::set<std::string>{pars.panelName_}, TARGET_NAME_PAT, "panel names",
									"letters, digits, - . and _", __PRETTY_FUNCTION__);
	checkNamesThrow(allSamples, LIBRARY_SAMPLE_NAME_PAT, "control sample names",
									"letters, digits, - . _ and spaces", __PRETTY_FUNCTION__);

	pmo::ReactionInfo reaction;
	reaction.reaction_name_ = "full";
	std::map<std::string, uint32_t> targetIds;
	for (const auto & target : targets) {
		pmo::TargetInfo info;
		info.target_name_ = target;
		auto primers = pars.primers_.find(target);
		if (pars.primers_.end() == primers) {
			info.forward_primer_.seq_ = "NA";
			info.reverse_primer_.seq_ = "NA";
			info.target_attributes_ = VecStr{pmo::PMO_FAKE_PRIMERS_ATTRIBUTE};
		} else {
			info.forward_primer_.seq_ = primers->second.first;
			info.reverse_primer_.seq_ = primers->second.second;
		}
		targetIds[target] = ret.target_info_.size();
		reaction.panel_targets_.emplace_back(ret.target_info_.size());
		ret.target_info_.emplace_back(info);
	}
	pmo::PanelInfo panel;
	panel.panel_name_ = pars.panelName_;
	panel.reactions_.emplace_back(reaction);
	ret.panel_info_.emplace_back(panel);

	//one representative per allele name
	std::map<std::string, std::map<std::string, uint32_t>> mhapIdByName;
	for (const auto & target : targets) {
		pmo::RepresentativeMicrohaplotypesForTarget representatives;
		representatives.target_id_ = targetIds.at(target);
		auto & ids = mhapIdByName[target];
		for (const auto & sample : targets_.at(target)) {
			for (const auto & allele : sample.second.alleles_) {
				auto search = ids.find(allele.name_);
				if (ids.end() != search) {
					const auto & rep = representatives.microhaplotypes_.at(search->second);
					if (rep.seq_ != allele.seq_ || rep.extras_.at(STRAIN_FIELD).get<std::string>() != allele.strain_) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "allele " << allele.name_ << " for " << target
							 << " has different sequences or strains in different samples" << "\n";
						throw std::runtime_error{ss.str()};
					}
					continue;
				}
				ids[allele.name_] = representatives.microhaplotypes_.size();
				pmo::RepresentativeMicrohaplotype rep;
				rep.microhaplotype_name_ = allele.name_;
				rep.seq_ = allele.seq_;
				rep.extras_[STRAIN_FIELD] = allele.strain_;
				representatives.microhaplotypes_.emplace_back(rep);
			}
		}
		ret.representative_microhaplotypes_.targets_.emplace_back(representatives);
	}

	//library samples with their mixture as their specimen, unless the sample is a different mixture for different targets
	std::map<std::string, std::map<std::string, std::string>> mixByTargetBySample;
	for (const auto & target : targets_) {
		for (const auto & sample : target.second) {
			mixByTargetBySample[sample.first][target.first] = sample.second.mix_;
		}
	}
	//specimen names have the characters the PMO doesn't allow replaced by _ (with a numbered suffix if that makes two the same),
	//a specimen whose name had to change keeps the mixture's name in its mix_name field
	std::map<std::string, uint32_t> specimenIds;
	std::set<std::string> specimenNames;
	auto specimenId = [&ret, &specimenIds, &specimenNames](const std::string & name) {
		auto search = specimenIds.find(name);
		if (specimenIds.end() == search) {
			search = specimenIds.emplace(name, ret.specimen_info_.size()).first;
			pmo::SpecimenInfo specimen;
			specimen.specimen_name_ = toSpecimenName(name);
			for (uint32_t suffix = 2; njh::in(specimen.specimen_name_, specimenNames); ++suffix) {
				specimen.specimen_name_ = njh::pasteAsStr(toSpecimenName(name), "_", suffix);
			}
			specimenNames.emplace(specimen.specimen_name_);
			if (specimen.specimen_name_ != name) {
				specimen.extras_[MIX_NAME_FIELD] = name;
			}
			ret.specimen_info_.emplace_back(specimen);
		}
		return search->second;
	};
	auto samples = njh::getVecOfMapKeys(mixByTargetBySample);
	njh::naturalSortNameSet(samples);
	pmo::DetectedMicrohaplotypes detected;
	for (const auto & sample : samples) {
		const auto & mixByTarget = mixByTargetBySample.at(sample);
		std::set<std::string> mixes;
		for (const auto & mix : mixByTarget) {
			mixes.emplace(mix.second);
		}
		pmo::LibrarySampleInfo library;
		library.library_sample_name_ = sample;
		library.panel_id_ = 0;
		if (1 == mixes.size()) {
			library.specimen_id_ = specimenId(*mixes.begin());
		} else {
			library.specimen_id_ = specimenId(sample);
			nlohmann::json mixJson = nlohmann::json::object();
			for (const auto & mix : mixByTarget) {
				mixJson[mix.first] = mix.second;
			}
			library.extras_[MIX_BY_TARGET_FIELD] = mixJson;
		}
		const uint32_t sampleId = ret.library_sample_info_.size();
		ret.library_sample_info_.emplace_back(library);

		pmo::DetectedMicrohaplotypesForSample sampleResult;
		sampleResult.library_sample_id_ = sampleId;
		for (const auto & target : targets) {
			auto expected = targets_.at(target).find(sample);
			if (targets_.at(target).end() == expected) {
				continue;
			}
			std::map<uint32_t, uint32_t> readsByMhapId;
			for (const auto & allele : expected->second.alleles_) {
				auto reads = static_cast<uint32_t>(std::round(allele.frac_ * pars.pseudoCountTotal_));
				readsByMhapId[mhapIdByName.at(target).at(allele.name_)] += std::max<uint32_t>(1, reads);
			}
			pmo::DetectedMicrohaplotypesForTarget targetResult;
			targetResult.mhaps_target_id_ = targetIds.at(target);
			for (const auto & mhapReads : readsByMhapId) {
				pmo::MicrohaplotypeForTarget mhap;
				mhap.mhap_id_ = mhapReads.first;
				mhap.reads_ = mhapReads.second;
				targetResult.mhaps_.emplace_back(mhap);
			}
			sampleResult.target_results_.emplace_back(targetResult);
		}
		detected.library_samples_.emplace_back(sampleResult);
	}
	ret.detected_microhaplotypes_.emplace_back(detected);
	ret.validate();
	return ret;
}

VecStr ExpectedControls::targetNames() const {
	auto ret = njh::getVecOfMapKeys(targets_);
	njh::naturalSortNameSet(ret);
	return ret;
}

VecStr ExpectedControls::samples(const std::string & target) const {
	auto search = targets_.find(target);
	if (targets_.end() == search) {
		return {};
	}
	auto ret = njh::getVecOfMapKeys(search->second);
	njh::sort(ret);
	return ret;
}

table ExpectedControls::genMixSetUpsTab(const std::string & targetNameColName) const {
	//target -> mix -> strain -> relative abundance
	std::map<std::string, std::map<std::string, std::map<std::string, double>>> strainAbundances;
	for (const auto & target : targets_) {
		for (const auto & sample : target.second) {
			auto & mix = strainAbundances[target.first][sample.second.mix_];
			if (!mix.empty()) {
				continue;
			}
			for (const auto & allele : sample.second.alleles_) {
				mix[allele.strain_] += allele.frac_;
			}
		}
	}
	//check if the mixtures are the same for every target they're in
	std::map<std::string, std::map<std::string, double>> mixes;
	bool sameForAllTargets = true;
	for (const auto & target : strainAbundances) {
		for (const auto & mix : target.second) {
			auto search = mixes.find(mix.first);
			if (mixes.end() == search) {
				mixes.emplace(mix.first, mix.second);
				continue;
			}
			if (njh::getVecOfMapKeys(search->second) != njh::getVecOfMapKeys(mix.second)) {
				sameForAllTargets = false;
				continue;
			}
			for (const auto & strain : mix.second) {
				if (std::abs(search->second.at(strain.first) - strain.second) > 1e-9) {
					sameForAllTargets = false;
				}
			}
		}
	}
	if (sameForAllTargets) {
		table ret(VecStr{"MixName", "strain", "relative_abundance"});
		for (const auto & mix : mixes) {
			for (const auto & strain : mix.second) {
				ret.addRow(mix.first, strain.first, strain.second);
			}
		}
		return ret;
	}
	table ret(VecStr{"MixName", "strain", "relative_abundance", targetNameColName});
	auto targets = targetNames();
	for (const auto & target : targets) {
		for (const auto & mix : strainAbundances.at(target)) {
			for (const auto & strain : mix.second) {
				ret.addRow(mix.first, strain.first, strain.second, target);
			}
		}
	}
	return ret;
}

table ExpectedControls::genSamplesToMixTab() const {
	std::set<std::pair<std::string, std::string>> rows;
	for (const auto & target : targets_) {
		for (const auto & sample : target.second) {
			rows.emplace(sample.first, sample.second.mix_);
		}
	}
	table ret(VecStr{"sample", "MixName"});
	for (const auto & row : rows) {
		ret.addRow(row.first, row.second);
	}
	return ret;
}

uint64_t ExpectedControls::maxSeqLen() const {
	uint64_t ret = 0;
	for (const auto & target : targets_) {
		for (const auto & sample : target.second) {
			for (const auto & allele : sample.second.alleles_) {
				ret = std::max<uint64_t>(ret, allele.seq_.size());
			}
		}
	}
	return ret;
}

}  // namespace njhseq
