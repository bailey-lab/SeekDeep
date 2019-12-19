/*
 * extractTargetsFromGenomes.cpp
 *
 *  Created on: Aug 29, 2018
 *      Author: nick
 */


#include "extractTargetsFromGenomes.hpp"

namespace njhseq {

void extractBetweenSeqsPars::setUpCoreOptions(seqSetUp & setUp, bool needReadLength){

	setUp.setOption(writeOutAllSeqsFile, "--writeOutAllSeqsFile", "Write Out All Seqs File without collpasing to unique sequences");

	setUp.setOption(shortNames, "--shortNames", "Create short names for reference genomes extractions");
	setUp.setOption(lenCutOffSizeExpand, "--lenCutOffSizeExpand", "When creating length cut off file how much to expand the length of the found targets");
	setUp.setOption(pairedEndLength, "--pairedEndLength", "Paired End Read Length", needReadLength);
	setUp.setOption(barcodeSize, "--barcodeSize", "Barcode Size, if on both primer, the sum of the two barcodes");
	setUp.setOption(errors, "--errors", "Number of errors to allow in primers");
	setUp.setOption(sizeLimit, "--sizeLimit", "Output target extractions for only targets below this size");
	setUp.setOption(pars.numThreads_, "--numThreads", "Number of CPUs to utilize");
	setUp.setOption(removeRefAlignments, "--removeRefAlignments", "Remove Ref Alignments");
  setUp.setOption(pars.genomeDir_, "--genomeDir", "Name of the genome file fnp", true);
  setUp.setOption(pars.primaryGenome_, "--primaryGenome", "The primary reference genome");
  setUp.setOption(pars.gffDir_, "--gffDir", "A directory with a gff for the genomes in --genomeDir, should be named GENOME.gff (for GENOME.fasta)");
  setUp.setOption(gffExtraAttributesStr, "--gffExtraAttributes", "Extra attributes to add to genome that has an accompany gff");
  setUp.setOption(selectedGenomesStr, "--selectedGenomes", "Name of the other genomes in --genomeDir to be read in, leave blank to just do all fastas");
	setUp.setOption(outputDirPars.dirName_, "--dout", "Output directory", true);
	setUp.setOption(outputDirPars.overWriteDir_, "--overWriteDir", "Overwrite Output directory");
}

//
//inline std::vector<bfs::path> filesInFolderDev(bfs::path d) {
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	std::vector<bfs::path> ret;
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	if (bfs::is_directory(d)) {
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		njh::files::dir dir_iter(d);
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		for (const auto& e : dir_iter) {
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			ret.emplace_back(e);
//		}
//	}
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	return ret;
//}
//
//
//inline bool rmDirForceDev(const bfs::path & dirName) {
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	if (bfs::is_directory(dirName)) {
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		auto files = filesInFolderDev(dirName);
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		for (const auto & f : files) {
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			rmDirForceDev(f.string());
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		}
//	}
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	return bfs::remove(dirName);
//}
//
//inline int32_t makeDirDev(const njh::files::MkdirPar & pars) {
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	if (bfs::exists(pars.dirName_)) {
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		if (pars.overWriteDir_) {
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			bool removalStatus = rmDirForceDev(pars.dirName_);
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			if (!removalStatus) {
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
//				std::stringstream ss;
//				ss << "Error in: " << __PRETTY_FUNCTION__
//						<< ", when removing directory " << pars.dirName_ << std::endl;
//				throw std::runtime_error { ss.str() };
//			}
//		} else {
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			std::stringstream ss;
//			ss << "Error in: " << __PRETTY_FUNCTION__ << ", directory " << pars.dirName_
//					<< " already exists, use overWrite = true to overwrite contents "
//					<< std::endl;
//			throw std::runtime_error { ss.str() };
//		}
//	}
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	int32_t directoryStatus = mkdir(pars.dirName_.c_str(), pars.perms_);
//	if (directoryStatus != 0) {
//		std::stringstream ss;
//		ss << "Error in: " << __PRETTY_FUNCTION__ << ", in making directory "
//				<< pars.dirName_ << ", mkdir return stats: " << directoryStatus << std::endl;
//		throw std::runtime_error { ss.str() };
//	}
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	return directoryStatus;
//}


void extractBetweenSeqs(const PrimersAndMids & ids,
		const extractBetweenSeqsPars & extractPars){

	njh::concurrent::LockableQueue<std::string> targetsQueue(getVectorOfMapKeys(ids.targets_));

	std::unique_ptr<MultiGenomeMapper> gMapper;
	//check for existence of genome directory
	njh::files::checkExistenceThrow(extractPars.pars.genomeDir_, __PRETTY_FUNCTION__);

	//make output directory
	//std::cout << "extractPars.outputDirPars: " << extractPars.outputDirPars.dirName_ << std::endl;
	//int32_t directoryStatus = mkdir(extractPars.outputDirPars.dirName_.c_str(), extractPars.outputDirPars.perms_);
	njh::files::makeDir(extractPars.outputDirPars);
	//std::cout << "extractPars.outputDirPars: " << extractPars.outputDirPars.dirName_ << std::endl;

	bfs::path outputDir = extractPars.outputDirPars.dirName_;

	//set up genome mapper;
	gMapper = std::make_unique<MultiGenomeMapper>(extractPars.pars);

	//set up selected genomes
	gMapper->setSelectedGenomes(extractPars.selectedGenomesStr);
	//init is threaded
	gMapper->pars_.numThreads_ = extractPars.pars.numThreads_;
	gMapper->init();
	gMapper->pars_.numThreads_ = 1;

	//set threads;
	uint32_t numThreads = extractPars.pars.numThreads_;
	if(extractPars.pars.numThreads_ >= 4){
		gMapper->pars_.numThreads_ = 2;
		numThreads = extractPars.pars.numThreads_/2;
	}
	if(1 == ids.targets_.size()){
		gMapper->pars_.numThreads_ = extractPars.pars.numThreads_;
		numThreads = 1;
	}


	auto alignToGenome = [&extractPars](njh::concurrent::LockableQueue<bfs::path> & genomesQuque,
			const bfs::path & forwardFnp,
			const bfs::path & reverseFnp,
			const bfs::path & alignDir){
		bfs::path genomeFnp = "";
		while(genomesQuque.getVal(genomeFnp)){
			{
				BioCmdsUtils bioRunner;
				bioRunner.verbose_ = extractPars.verbose_;
				auto seqOpts = SeqIOOptions::genFastaIn(forwardFnp, false);
				seqOpts.out_.outFilename_ = njh::files::make_path(alignDir, bfs::basename(genomeFnp) + "_" + bfs::basename(forwardFnp) + ".sorted.bam");
				bioRunner.bowtie2Align(seqOpts, genomeFnp	, "-D 20 -R 3 -N 1 -L 15 -i S,1,0.5 -a --end-to-end");
			}
			{
				BioCmdsUtils bioRunner;
				bioRunner.verbose_ = extractPars.verbose_;
				auto seqOpts = SeqIOOptions::genFastaIn(reverseFnp, false);
				seqOpts.out_.outFilename_ = njh::files::make_path(alignDir, bfs::basename(genomeFnp) + "_" + bfs::basename(reverseFnp) + ".sorted.bam");
				bioRunner.bowtie2Align(seqOpts, genomeFnp, "-D 20 -R 3 -N 1 -L 15 -i S,1,0.5 -a --end-to-end"	);
			}
		}
	};

	comparison allowableErrors;
	allowableErrors.hqMismatches_ = extractPars.errors;


	auto extractPathway =
			[&targetsQueue,&extractPars,&outputDir,&ids,&alignToGenome,&allowableErrors](
					const std::unique_ptr<MultiGenomeMapper> & gMapper) {
				std::string target;
				while(targetsQueue.getVal(target)) {
					const auto & primerInfo = ids.pDeterminator_->primers_.at(target);
					auto primerDirectory = njh::files::makeDir(outputDir, njh::files::MkdirPar(primerInfo.primerPairName_));
					auto bedDirectory = njh::files::makeDir(primerDirectory, njh::files::MkdirPar("genomeLocations"));
					auto forwardOpts = SeqIOOptions::genFastaOut(njh::files::make_path(primerDirectory, "forwardPrimer"));
					auto reverseOpts = SeqIOOptions::genFastaOut(njh::files::make_path(primerDirectory, "reversePrimer"));
					VecStr forDegens;
					for(const auto & fwd : primerInfo.fwds_){
						addOtherVec(forDegens, createDegenStrs(fwd.info_.seq_));
					}
					VecStr revDegens;
					for(const auto & rev : primerInfo.revs_){
						addOtherVec(revDegens, createDegenStrs(rev.info_.seq_));
					}
					//
					std::vector<seqInfo> forSeqs;
					if(forDegens.size() == 1){
						forSeqs.emplace_back(primerInfo.fwds_.front().info_);
					}else{
						forSeqs = vecStrToReadObjs<seqInfo>(forDegens, primerInfo.primerPairName_);
					}
					std::vector<seqInfo> revSeqs;
					if(revDegens.size() == 1){
						revSeqs.emplace_back(primerInfo.revs_.front().info_);
					}else{
						revSeqs = vecStrToReadObjs<seqInfo>(revDegens, primerInfo.primerPairName_);
					}
					SeqOutput::write(forSeqs, forwardOpts);
					SeqOutput::write(revSeqs, reverseOpts);
					auto refAlignmentDir = njh::files::makeDir(primerDirectory,njh::files::MkdirPar{ "refAlignments"});
					std::vector<std::thread> threads;
					njh::concurrent::LockableQueue<bfs::path> genomesQueue(gMapper->getGenomeFnps());
					auto forOutFnp = forwardOpts.out_.outName();
					auto revOutFnp = reverseOpts.out_.outName() ;
					for(uint32_t t = 0; t < gMapper->pars_.numThreads_; ++t){
						threads.emplace_back(std::thread(alignToGenome,
								std::ref(genomesQueue),
								std::cref(forOutFnp),
								std::cref(revOutFnp),
								std::cref(refAlignmentDir)));
					}
					for(auto & t : threads){
						t.join();
					}
					struct GenExtracRes{
						uint32_t forwardHits_{0};
						uint32_t reverseHits_{0};
						uint32_t extractCounts_{0};
					};
					std::unordered_map<std::string, GenExtracRes> genomeExtractionsResults;
					for(const auto & genome : gMapper->genomes_){
						genomeExtractionsResults[genome.first] = GenExtracRes{};
					}
					std::unordered_map<std::string, std::vector<GenomeExtractResult>> genomeExtracts;
					for(const auto & genome : gMapper->genomes_) {
						auto forBamFnp = njh::files::make_path(refAlignmentDir,
								bfs::basename(genome.second->fnp_) + "_" + bfs::basename(forwardOpts.out_.outName()) + ".sorted.bam");
						auto revBamFnp = njh::files::make_path(refAlignmentDir,
														bfs::basename(genome.second->fnp_) + "_" + bfs::basename(reverseOpts.out_.outName()) + ".sorted.bam");
						auto forResults = gatherMapResults(forBamFnp, genome.second->fnpTwoBit_, allowableErrors);
						auto revResults = gatherMapResults(revBamFnp, genome.second->fnpTwoBit_, allowableErrors);
						genomeExtractionsResults[genome.first].forwardHits_ = forResults.size();
						genomeExtractionsResults[genome.first].reverseHits_ = revResults.size();
						if(!forResults.empty() && !revResults.empty()){
							auto uniForRes = getUniqueLocationResults(forResults);
							auto uniRevRes = getUniqueLocationResults(revResults);
							genomeExtracts[genome.first] = getPossibleGenomeExtracts(uniForRes, uniRevRes, extractPars.sizeLimit);
						}
					}
					std::vector<seqInfo> refSeqs;
					std::vector<seqInfo> refTrimmedSeqs;

					std::vector<seqInfo> allSeqs;
					std::vector<seqInfo> allSeqsTrimmedSeqs;

					for(auto & genome : genomeExtracts){
						if(genome.second.empty()){
							continue;
						}
						genomeExtractionsResults[genome.first].extractCounts_ = genome.second.size();
						OutputStream bedOut{OutOptions(njh::files::make_path(bedDirectory, genome.first + ".bed"))};
						OutputStream bedInnerOut{OutOptions(njh::files::make_path(bedDirectory, genome.first + "_inner.bed"))};
						uint32_t extractionCount = 0;
						for(auto & extract : genome.second){
							extract.setRegion();

							auto name = genome.first;
							if(0 != extractionCount){
								name.append("." + estd::to_string(extractionCount));
							}
							MetaDataInName meta;
							meta.addMeta("genome", genome.first);
							meta.addMeta("target", target);
							meta.addMeta("extractionCount", extractionCount);
							name += " " + meta.createMetaName();

							{
								auto bedRegion = extract.gRegion_->genBedRecordCore();
								bedRegion.name_ = name;
								bedOut << bedRegion.toDelimStr() << std::endl;
							}
							{
								auto bedRegion = extract.gRegionInner_->genBedRecordCore();
								bedRegion.name_ = name;
								bedInnerOut << bedRegion.toDelimStr() << std::endl;
							}

							TwoBit::TwoBitFile tReader(gMapper->genomes_.at(genome.first)->fnpTwoBit_);
							auto eSeq = extract.gRegion_->extractSeq(tReader);
							if(extractPars.shortNames){
								trimAtFirstWhitespace(name);
							}
							eSeq.name_ = name;
							bool refFound = false;
							if(extractPars.writeOutAllSeqsFile){
								allSeqs.emplace_back(eSeq);
							}
							for(auto & rSeq : refSeqs){
								if(rSeq.seq_ == eSeq.seq_){
									refFound = true;
									rSeq.name_.append("-" + eSeq.name_);
								}
							}
							if(!refFound){
								refSeqs.emplace_back(eSeq);
							}

							bool trimmed_refFound = false;
							auto innerSeq = extract.gRegionInner_->extractSeq(tReader);
							innerSeq.name_ = name;
							if(extractPars.writeOutAllSeqsFile){
								allSeqsTrimmedSeqs.emplace_back(innerSeq);
							}
							for(auto & rSeq : refTrimmedSeqs){
								if(rSeq.seq_ == innerSeq.seq_){
									trimmed_refFound = true;
									rSeq.name_.append("-" + innerSeq.name_);
								}
							}
							if(!trimmed_refFound){
								refTrimmedSeqs.emplace_back(innerSeq);
							}
							++extractionCount;
						}
					}

					if(!allSeqs.empty()){
						auto fullSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(primerDirectory, "all_" + primerInfo.primerPairName_ +".fasta"));
						SeqOutput::write(allSeqs, fullSeqOpts);
					}
					if(!allSeqsTrimmedSeqs.empty()){
						auto innerSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(primerDirectory, "all_" + primerInfo.primerPairName_ +"_primersRemoved.fasta"));
						SeqOutput::write(allSeqsTrimmedSeqs, innerSeqOpts);
					}

					if(!refSeqs.empty()){
						auto fullSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(primerDirectory, primerInfo.primerPairName_ +".fasta"));
						SeqOutput::write(refSeqs, fullSeqOpts);
					}
					if(!refTrimmedSeqs.empty()){
						auto innerSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(primerDirectory, primerInfo.primerPairName_ +"_primersRemoved.fasta"));
						SeqOutput::write(refTrimmedSeqs, innerSeqOpts);
					}

					table performanceTab(VecStr{"genome", "forwardPrimerHits", "reversePrimerHits", "extractionCounts", "target"});
					auto genomeKeys = getVectorOfMapKeys(genomeExtractionsResults);
					njh::sort(genomeKeys);
					for(const auto & genomeKey : genomeKeys){
						performanceTab.addRow(genomeKey,
								genomeExtractionsResults[genomeKey].forwardHits_,
								genomeExtractionsResults[genomeKey].reverseHits_,
								genomeExtractionsResults[genomeKey].extractCounts_,
								primerInfo.primerPairName_);
					}
					auto perTabOpts = TableIOOpts::genTabFileOut(njh::files::make_path(primerDirectory, "extractionCounts"),true);
					performanceTab.outPutContents(perTabOpts);
					if(extractPars.removeRefAlignments){
						njh::files::rmDirForce(refAlignmentDir);
					}
				}
			};
	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < numThreads; ++t){
		threads.emplace_back(std::thread(extractPathway,
				std::cref(gMapper)));
	}
	for(auto & t : threads){
		t.join();
	}

	auto locationsCombined = njh::files::make_path(outputDir, "locationsByGenome");
	njh::files::makeDir(njh::files::MkdirPar{locationsCombined});

	//outer with primer
	for(const auto & genome : gMapper->genomes_){
		auto genomeBedOpts = njh::files::make_path(locationsCombined, genome.first + ".bed");
		std::vector<std::shared_ptr<Bed6RecordCore>> allRegions;
		for(const auto & tar : ids.targets_){
			auto bedForTarFnp = njh::files::make_path(outputDir, tar.first, "genomeLocations", genome.first + ".bed");
			if(bfs::exists(bedForTarFnp)){
				auto locs = getBeds(bedForTarFnp);
				addOtherVec(allRegions, locs);
			}
		}
		if(!allRegions.empty()){
			if("" != genome.second->gffFnp_){
				intersectBedLocsWtihGffRecordsPars pars(genome.second->gffFnp_);
				pars.selectFeatures_ = VecStr{"gene"};
				pars.extraAttributes_ = tokenizeString(extractPars.gffExtraAttributesStr, ",");
				auto jsonValues = intersectBedLocsWtihGffRecords(allRegions, pars);
				auto geneIds = jsonValues.getMemberNames();
				std::set<std::string> geneIdsSet(geneIds.begin(), geneIds.end());
				auto genes = GeneFromGffs::getGenesFromGffForIds(genome.second->gffFnp_,geneIdsSet);
				for(const auto & reg : allRegions){
					if(reg->extraFields_.empty()){
						continue;
					}
					auto geneToks = njh::tokenizeString(reg->extraFields_.back(), "],[");
					if(geneToks.size() > 1){
						geneToks.front().append("]");
						for(const auto pos : iter::range<uint32_t>(1, geneToks.size() - 1)){
							geneToks[pos] = "[" + geneToks[pos] + "]";
						}
						geneToks.back() = "[" +geneToks.back();
					}
					std::string replacement = "";
					for(const auto & geneTok : geneToks){
						MetaDataInName geneMeta(geneTok);
						TwoBit::TwoBitFile tReader(genome.second->fnpTwoBit_);
						auto infos = njh::mapAt(genes, geneMeta.getMeta("ID"))->generateGeneSeqInfo(tReader, false);
						for(const auto & info : infos){
							auto posInfos = info.second->getInfosByGDNAPos();

							auto minPos = vectorMinimum(getVectorOfMapKeys(posInfos));
							auto maxPos = vectorMaximum(getVectorOfMapKeys(posInfos));
							auto startPos = std::max(minPos, reg->chromStart_);
							auto stopPos =  std::min(maxPos, reg->chromEnd_);
							auto aaStartPos = std::min(posInfos[startPos].aaPos_, posInfos[stopPos].aaPos_) + 1;
							auto aaStopPos =  std::max(posInfos[startPos].aaPos_, posInfos[stopPos].aaPos_) + 1;
							geneMeta.addMeta(info.first + "-AAStart", aaStartPos );
							geneMeta.addMeta(info.first + "-AAStop", aaStopPos );
						}
						if("" != replacement){
							replacement += ",";
						}
						replacement += geneMeta.createMetaName();
					}
					reg->extraFields_.back() = replacement;
				}
			}
			OutputStream genomeBedOut(genomeBedOpts);
			njh::sort(allRegions, [](
					const std::shared_ptr<Bed3RecordCore> & bed1,
					const std::shared_ptr<Bed3RecordCore> & bed2
					){
				if(bed1->chrom_ == bed2->chrom_){
					return bed1->chromStart_ < bed2->chromStart_;
				}else{
					return bed1->chrom_ < bed2->chrom_;
				}
			});
			for(const auto & loc : allRegions){
				genomeBedOut << loc->toDelimStrWithExtra() << std::endl;
			}
		}
	}

	//inner with primer
	for(const auto & genome : gMapper->genomes_){
		auto genomeBedOpts = njh::files::make_path(locationsCombined, genome.first + "_inner.bed");
		std::vector<std::shared_ptr<Bed6RecordCore>> allRegions;
		for(const auto & tar : ids.targets_){
			auto bedForTarFnp = njh::files::make_path(outputDir, tar.first, "genomeLocations", genome.first + "_inner.bed");
			if(bfs::exists(bedForTarFnp)){
				auto locs = getBeds(bedForTarFnp);
				addOtherVec(allRegions, locs);
			}
		}
		if(!allRegions.empty()){
			if("" != genome.second->gffFnp_){
				intersectBedLocsWtihGffRecordsPars pars(genome.second->gffFnp_);
				pars.selectFeatures_ = VecStr{"gene"};
				pars.extraAttributes_ = tokenizeString(extractPars.gffExtraAttributesStr, ",");
				auto jsonValues = intersectBedLocsWtihGffRecords(allRegions, pars);
				auto geneIds = jsonValues.getMemberNames();
				std::set<std::string> geneIdsSet(geneIds.begin(), geneIds.end());
				auto genes = GeneFromGffs::getGenesFromGffForIds(genome.second->gffFnp_,geneIdsSet);
				for(const auto & reg : allRegions){
					if(reg->extraFields_.empty()){
						continue;
					}
					auto geneToks = njh::tokenizeString(reg->extraFields_.back(), "],[");
					if(geneToks.size() > 1){
						geneToks.front().append("]");
						for(const auto pos : iter::range<uint32_t>(1, geneToks.size() - 1)){
							geneToks[pos] = "[" + geneToks[pos] + "]";
						}
						geneToks.back() = "[" +geneToks.back();
					}
					std::string replacement = "";
					for(const auto & geneTok : geneToks){
						MetaDataInName geneMeta(geneTok);
						TwoBit::TwoBitFile tReader(genome.second->fnpTwoBit_);
						auto infos = njh::mapAt(genes, geneMeta.getMeta("ID"))->generateGeneSeqInfo(tReader, false);
						for(const auto & info : infos){
							auto posInfos = info.second->getInfosByGDNAPos();

							auto minPos = vectorMinimum(getVectorOfMapKeys(posInfos));
							auto maxPos = vectorMaximum(getVectorOfMapKeys(posInfos));
							auto startPos = std::max(minPos, reg->chromStart_);
							auto stopPos =  std::min(maxPos, reg->chromEnd_);
							auto aaStartPos = std::min(posInfos[startPos].aaPos_, posInfos[stopPos].aaPos_) + 1;
							auto aaStopPos =  std::max(posInfos[startPos].aaPos_, posInfos[stopPos].aaPos_) + 1;
							geneMeta.addMeta(info.first + "-AAStart", aaStartPos );
							geneMeta.addMeta(info.first + "-AAStop", aaStopPos );
						}
						if("" != replacement){
							replacement += ",";
						}
						replacement += geneMeta.createMetaName();
					}
					reg->extraFields_.back() = replacement;
				}
			}
			OutputStream genomeBedOut(genomeBedOpts);
			njh::sort(allRegions, [](
					const std::shared_ptr<Bed3RecordCore> & bed1,
					const std::shared_ptr<Bed3RecordCore> & bed2
					){
				if(bed1->chrom_ == bed2->chrom_){
					return bed1->chromStart_ < bed2->chromStart_;
				}else{
					return bed1->chrom_ < bed2->chrom_;
				}
			});
			for(const auto & loc : allRegions){
				genomeBedOut << loc->toDelimStrWithExtra() << std::endl;
			}
		}
	}
}


}  // namespace njhseq
