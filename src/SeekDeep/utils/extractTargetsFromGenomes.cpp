/*
 * extractTargetsFromGenomes.cpp
 *
 *  Created on: Aug 29, 2018
 *      Author: nick
 */


#include "extractTargetsFromGenomes.hpp"

#include <njhseq/GenomeUtils/GenomeExtraction/ParsingAlignmentInfo/AlignmentResults.hpp>
#include <njhseq/GenomeUtils/GenomeExtraction/ParsingAlignmentInfo/GenomeExtractResult.hpp>
#include <njhseq/objects/Gene/GeneFromGffs.hpp>

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

	struct GenExtracRes{
		uint32_t forwardHits_{0};
		uint32_t reverseHits_{0};
		uint32_t extractCounts_{0};
	};

	std::unordered_map<std::string, GenExtracRes> genomeExtractionsResults;
	for(const auto & genome : gMapper->genomes_){
		genomeExtractionsResults[genome.first] = GenExtracRes{};
	}
	if(extractPars.errors > 0){
		auto maxPrimerSize = ids.pDeterminator_->getMaxPrimerSize();
		//setup directories
		for(const auto & target : ids.pDeterminator_->primers_){
			const auto & primerInfo = target.second;
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
		}
		struct PrimerPairSearchResults{
			std::vector<GenomicRegion> fPrimerPositions_;
			std::vector<GenomicRegion> rPrimerPositions_;


			class GenomeExtractResultByPrimerPair {
			public:
				GenomeExtractResultByPrimerPair(const GenomicRegion & fPrimerReg,
						const GenomicRegion & rPrimerReg): fPrimerReg_(fPrimerReg),rPrimerReg_(rPrimerReg){
					setRegion();
				}
				GenomicRegion fPrimerReg_;
				GenomicRegion rPrimerReg_;

				std::shared_ptr<GenomicRegion> gRegion_ ;//= nullptr; it's silly but this is due to eclipse for now, should be able to revert once cdt 9.3 is released june 28 2017
				std::shared_ptr<GenomicRegion> gRegionInner_ ;

				void setRegion(){
					if (fPrimerReg_.chrom_ != rPrimerReg_.chrom_) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error extention chrom, "
								<< fPrimerReg_.chrom_ << "doesn't equal ligation chrom "
								<< rPrimerReg_.chrom_ << "\n";
						throw std::runtime_error { ss.str() };
					}
					if (fPrimerReg_.reverseSrand_ == rPrimerReg_.reverseSrand_) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__
								<< ", error extention and ligation are on the same strand, should be mapping to opposite strands"
								<< "\n";
						throw std::runtime_error { ss.str() };
					}
					if (fPrimerReg_.reverseSrand_) {
						if (fPrimerReg_.start_ < rPrimerReg_.start_) {
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__
									<< ", error if extention is mapping to the reverse strand, it's start, "
									<< fPrimerReg_.start_
									<< ", should be greater than ligation start, "
									<< rPrimerReg_.start_ << "\n";
							throw std::runtime_error { ss.str() };
						}
					}else{
						if (fPrimerReg_.start_ > rPrimerReg_.start_) {
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__
									<< ", error if extention is mapping to the plus strand, it's start, "
									<< fPrimerReg_.start_
									<< ", should be less than than ligation start, "
									<< rPrimerReg_.start_ << "\n";
							throw std::runtime_error { ss.str() };
						}
					}

					size_t start = fPrimerReg_.start_;
					size_t end = rPrimerReg_.end_;
					size_t innerStart = fPrimerReg_.end_;
					size_t innereEnd = rPrimerReg_.start_;
					if(fPrimerReg_.reverseSrand_){
						start = rPrimerReg_.start_;
						end = fPrimerReg_.end_;
						innerStart = rPrimerReg_.end_;
						innereEnd = fPrimerReg_.start_;
					}

					gRegion_ = std::make_shared<GenomicRegion>(fPrimerReg_.uid_ + "-" + rPrimerReg_.uid_, fPrimerReg_.chrom_, start, end, fPrimerReg_.reverseSrand_);
					gRegionInner_ = std::make_shared<GenomicRegion>(fPrimerReg_.uid_ + "-" + rPrimerReg_.uid_, fPrimerReg_.chrom_, innerStart, innereEnd, fPrimerReg_.reverseSrand_);
				}

			};


			static std::vector<GenomeExtractResultByPrimerPair> getPossibleGenomeExtracts(const std::vector<GenomicRegion> &fPrimerPositions,
					const std::vector<GenomicRegion> & rPrimerPositions,
					const size_t insertSizeCutOff = std::numeric_limits<size_t>::max()){
				std::vector<GenomeExtractResultByPrimerPair> ret;
				//same chrom, opposite strands, less than the insert size
				for (const auto & fwd : fPrimerPositions) {
					for (const auto & rev : rPrimerPositions) {
						//need to be on the same chromosome
						//need to be on opposite strands (should both should be in 5'->3' direction
						//and they shouldn't overlap
						if (fwd.chrom_ == rev.chrom_
								&& fwd.reverseSrand_ != rev.reverseSrand_
								&& !fwd.overlaps(rev)
								&& fwd.start_ != rev.end_
								&& fwd.end_ != rev.start_ ) {

							if(fwd.reverseSrand_){
								if(fwd.start_ > rev.start_){
									GenomeExtractResultByPrimerPair extraction(fwd, rev);
									if (extraction.gRegion_->getLen() <= insertSizeCutOff) {
										ret.emplace_back(extraction);
									}
								}
							}else{
								if(fwd.start_ < rev.start_){
									GenomeExtractResultByPrimerPair extraction(fwd, rev);
									if (extraction.gRegion_->getLen() <= insertSizeCutOff) {
										ret.emplace_back(extraction);
									}
								}
							}
						}
					}
				}
				return ret;
			}

			std::vector<GenomeExtractResultByPrimerPair> regions_;


		};
		// target, genome, results for target for genome
		std::map<std::string, std::map<std::string, PrimerPairSearchResults>> resultsByTargetByGenome;
		std::mutex resultsByTargetByGenomeMut;
		for(const auto & genome : gMapper->genomes_){
//			SeqInput genomeReader(SeqIOOptions(genome.second->fnp_, SeqIOOptions::getInFormatFromFnp(genome.second->fnp_), false));
//			genomeReader.openIn();

			TwoBit::TwoBitFile globalTreader(genome.second->fnpTwoBit_);
			auto names = globalTreader.sequenceNames();
			njh::concurrent::LockableQueue<std::string> nameQueue(names);

			std::function<void()> extractGenomicPositions = [&maxPrimerSize,&ids,
																											 &resultsByTargetByGenomeMut,&resultsByTargetByGenome,
																											 &extractPars,&genome,&nameQueue](){
				//seqInfo seq;
				std::map<std::string, std::map<std::string, PrimerPairSearchResults>> resultsByTargetByGenomeCurrent;
				std::string seqName = "";
				TwoBit::TwoBitFile treader(genome.second->fnpTwoBit_);
				auto seqLens = treader.getSeqLens();
				seqInfo seq;
				while(nameQueue.getVal(seqName)){
				//while(genomeReader.readNextReadLock(seq)){
					treader[seqName]->getSequence(seq.seq_);
					seq.name_ = seqName;
//					auto chromName = seq.name_;
//					trimAtFirstWhitespace(chromName);
					auto chromName = seqName;

					kmerInfo info(seq.seq_, maxPrimerSize, false);
					for(const auto & primer : ids.pDeterminator_->primers_){
						std::vector<GenomicRegion> fPrimerPositions;
						std::vector<GenomicRegion> rPrimerPositions;

						//searching with `5-`3 for both primers so in order to match later the strands need to be opposite
						for(const auto & k : info.kmers_){
							for(const auto & fwd : primer.second.fwds_){
								if(fwd.mot_.frontPassNoCheck(k.first,extractPars.errors )){
									auto score = fwd.mot_.scoreMotif(k.first.begin(), k.first.begin() + fwd.mot_.size());
									auto errors = fwd.mot_.size() - score;
									for(const auto & pos : k.second.positions_){
										auto start = pos;
										auto end = pos + fwd.mot_.size();
										fPrimerPositions.emplace_back(GenomicRegion(njh::pasteAsStr(primer.second.primerPairName_, "-", "forwardPrimer[errors=", errors, ";]"),
												chromName,
												start, end, false));
									}
								}
							} //fwd primers
							for(const auto & rev : primer.second.revs_){
								if(rev.mot_.frontPassNoCheck(k.first,extractPars.errors )){
									auto score = rev.mot_.scoreMotif(k.first.begin(), k.first.begin() + rev.mot_.size());
									auto errors = rev.mot_.size() - score;
									for(const auto & pos : k.second.positions_){
										auto start = pos;
										auto end = pos + rev.mot_.size();
										rPrimerPositions.emplace_back(GenomicRegion(njh::pasteAsStr(primer.second.primerPairName_, "-", "reversePrimer[errors=", errors, ";]"),
												chromName,
												start, end, false));
									}
								}
							} //reverse primers
						} //forward strand search
						for(const auto & k : info.kmers_){
							for(const auto & fwd : primer.second.fwds_){
								if(fwd.motRC_.frontPassNoCheck(k.first,extractPars.errors )){
									auto score = fwd.motRC_.scoreMotif(k.first.begin(), k.first.begin() + fwd.motRC_.size());
									auto errors = fwd.motRC_.size() - score;
									for(const auto & pos : k.second.positions_){
										auto start = pos;
										auto end = pos + fwd.motRC_.size();
										fPrimerPositions.emplace_back(GenomicRegion(njh::pasteAsStr(primer.second.primerPairName_, "-", "forwardPrimer[errors=", errors, ";]"),
												chromName,
												start, end, true));
									}
								}
							} //fwd primers
							for(const auto & rev : primer.second.revs_){
								if(rev.motRC_.frontPassNoCheck(k.first,extractPars.errors )){
									auto score = rev.motRC_.scoreMotif(k.first.begin(), k.first.begin() + rev.motRC_.size());
									auto errors = rev.motRC_.size() - score;
									for(const auto & pos : k.second.positions_){
										auto start = pos;
										auto end = start + rev.motRC_.size();
										rPrimerPositions.emplace_back(GenomicRegion(njh::pasteAsStr(primer.second.primerPairName_, "-", "reversePrimer[errors=", errors, ";]"),
												chromName,
												start, end, true));
									}
								}
							} //reverse primers
						} //forward strand search
						{
							//add results
							std::vector<GenomicRegion> fPrimerPositionsUnique;
							njh::sort(fPrimerPositions, [](const GenomicRegion & reg1, const GenomicRegion & reg2){
								auto reg1Meta =MetaDataInName(reg1.uid_.substr(reg1.uid_.rfind("[")));
								auto reg2Meta =MetaDataInName(reg2.uid_.substr(reg2.uid_.rfind("[")));
								return reg1Meta.getMeta<uint32_t>("errors") < reg2Meta.getMeta<uint32_t>("errors");
							});
							for(const auto & fPrim : fPrimerPositions){
								bool add = true;
								for(const auto & alreadyAdded : fPrimerPositionsUnique){
									if(alreadyAdded.sameRegion(fPrim)){
										add = false;
										break;
									}
								}
								if(add){
									fPrimerPositionsUnique.emplace_back(fPrim);
								}
							}
							std::vector<GenomicRegion> rPrimerPositionsUnique;
							njh::sort(rPrimerPositions, [](const GenomicRegion & reg1, const GenomicRegion & reg2){
								auto reg1Meta =MetaDataInName(reg1.uid_.substr(reg1.uid_.rfind("[")));
								auto reg2Meta =MetaDataInName(reg2.uid_.substr(reg2.uid_.rfind("[")));
								return reg1Meta.getMeta<uint32_t>("errors") < reg2Meta.getMeta<uint32_t>("errors");
							});
							for(const auto & rPrim : rPrimerPositions){
								bool add = true;
								for(const auto & alreadyAdded : rPrimerPositionsUnique){
									if(alreadyAdded.sameRegion(rPrim)){
										add = false;
										break;
									}
								}
								if(add){
									rPrimerPositionsUnique.emplace_back(rPrim);
								}
							}

							addOtherVec(resultsByTargetByGenome[primer.second.primerPairName_][genome.first].fPrimerPositions_, fPrimerPositionsUnique);
							addOtherVec(resultsByTargetByGenome[primer.second.primerPairName_][genome.first].rPrimerPositions_, rPrimerPositionsUnique);
							addOtherVec(resultsByTargetByGenome[primer.second.primerPairName_][genome.first].regions_, PrimerPairSearchResults::getPossibleGenomeExtracts(fPrimerPositionsUnique, rPrimerPositionsUnique, extractPars.sizeLimit));
						}
					} //end of this primer's search
				} //end of chromosome search
				{
					//add results
					std::lock_guard<std::mutex> lock(resultsByTargetByGenomeMut);
					for(const auto & tar : resultsByTargetByGenomeCurrent){
						for(const auto & gen : tar.second){
							addOtherVec(resultsByTargetByGenome[tar.first][genome.first].fPrimerPositions_, gen.second.fPrimerPositions_);
							addOtherVec(resultsByTargetByGenome[tar.first][genome.first].rPrimerPositions_, gen.second.rPrimerPositions_);
							addOtherVec(resultsByTargetByGenome[tar.first][genome.first].regions_, gen.second.regions_);
						}
					}
				}
			};

			njh::concurrent::runVoidFunctionThreaded(extractGenomicPositions, extractPars.pars.numThreads_);
		} //end of genome search
		auto genRegSort = [](const GenomicRegion & reg1,const GenomicRegion & reg2){
			if(reg1.chrom_ == reg2.chrom_){
				if(reg1.start_ == reg2.start_){
					if(reg1.end_ == reg2.end_){
						return reg1.uid_ < reg2.uid_;
					}else{
						return reg1.end_ < reg2.end_;
					}
				}else{
					return reg1.start_ < reg2.start_;
				}
			}else{
				return reg1.chrom_ < reg2.chrom_;
			}
		};
		for(auto & target : resultsByTargetByGenome){
			for(auto & genome : target.second){

				njh::sort(genome.second.fPrimerPositions_,genRegSort);
				njh::sort(genome.second.rPrimerPositions_,genRegSort);
				njh::sort(genome.second.regions_,[](const PrimerPairSearchResults::GenomeExtractResultByPrimerPair & reg1,const PrimerPairSearchResults::GenomeExtractResultByPrimerPair & reg2){
					if(reg1.gRegion_->chrom_ == reg2.gRegion_->chrom_){
						if(reg1.gRegion_->start_ == reg2.gRegion_->start_){
							if(reg1.gRegion_->end_ == reg2.gRegion_->end_){
								return reg1.gRegion_->uid_ < reg2.gRegion_->uid_;
							}else{
								return reg1.gRegion_->end_ < reg2.gRegion_->end_;
							}
						}else{
							return reg1.gRegion_->start_ < reg2.gRegion_->start_;
						}
					}else{
						return reg1.gRegion_->chrom_ < reg2.gRegion_->chrom_;
					}
				});

			}
		}
		for(const auto & target : resultsByTargetByGenome){
			auto primerDirectory = njh::files::make_path(outputDir, target.first);
			auto bedDirectory = njh::files::make_path(primerDirectory, "genomeLocations");

			std::unordered_map<std::string, GenExtracRes> genomeExtractionsResults;
			for(const auto & genome : gMapper->genomes_){
				genomeExtractionsResults[genome.first] = GenExtracRes{};
			}
			for(const auto & genome : gMapper->genomes_) {

				auto primerLocationsFnp = njh::files::make_path(bedDirectory,bfs::basename(genome.second->fnp_) + "_primerLocs.bed");
				OutputStream primerLocationsOut(primerLocationsFnp);
				std::vector<GenomicRegion> allPrimerRegions;
				if(njh::in(genome.first, target.second)){
					addOtherVec(allPrimerRegions, target.second.at(genome.first).fPrimerPositions_);
					genomeExtractionsResults[genome.first].forwardHits_ = target.second.at(genome.first).fPrimerPositions_.size();
				}

				if(njh::in(genome.first, target.second)){
					addOtherVec(allPrimerRegions, target.second.at(genome.first).rPrimerPositions_);
					genomeExtractionsResults[genome.first].reverseHits_ = target.second.at(genome.first).rPrimerPositions_.size();
				}

				njh::sort(allPrimerRegions, genRegSort);
				for(const auto & primerRes : allPrimerRegions){
					auto bedOut = primerRes.genBedRecordCore();
					primerLocationsOut << bedOut.toDelimStrWithExtra() << std::endl;
				}
			}
			std::vector<seqInfo> refSeqs;
			std::vector<seqInfo> refTrimmedSeqs;

			std::vector<seqInfo> allSeqs;
			std::vector<seqInfo> allSeqsTrimmedSeqs;

			for(auto & genome : target.second){
				if(genome.second.regions_.empty()){
					continue;
				}
				genomeExtractionsResults[genome.first].extractCounts_ = genome.second.regions_.size();
				OutputStream bedOut{OutOptions(njh::files::make_path(bedDirectory, genome.first + ".bed"))};
				OutputStream bedInnerOut{OutOptions(njh::files::make_path(bedDirectory, genome.first + "_inner.bed"))};
				uint32_t extractionCount = 0;
				OutputStream regionInfoOut{OutOptions(njh::files::make_path(bedDirectory, genome.first + "_regionInfo.tab.txt"))};
				regionInfoOut << "#chrom\tfullStart\tfullStop\tname\tlength\tstrand\tgenome\ttarget\tfPrimerStart\tfPrimerStop\terrorsInFPrimer\tinsertStart\tinsertStop\trPrimerStart\trPrimerStop\terrorsInRPrimer" << std::endl;
				for(auto & extract : genome.second.regions_){
					auto name = genome.first;
					if(0 != extractionCount){
						name.append("." + estd::to_string(extractionCount));
					}
					MetaDataInName meta;
					meta.addMeta("genome", genome.first);
					meta.addMeta("target", target.first);
					meta.addMeta("extractionCount", extractionCount);

					name += " " + meta.createMetaName();
					auto fPrimerMeta = MetaDataInName(extract.fPrimerReg_.uid_.substr(extract.fPrimerReg_.uid_.rfind("[")));
					auto rPrimerMeta = MetaDataInName(extract.rPrimerReg_.uid_.substr(extract.rPrimerReg_.uid_.rfind("[")));

					regionInfoOut << extract.gRegion_->chrom_
							<< "\t" << extract.gRegion_->start_
							<< "\t" << extract.gRegion_->end_
							<< "\t" << name
							<< "\t" << extract.gRegion_->getLen()
							<< "\t" << (extract.gRegion_->reverseSrand_ ? '-' : '+')
							<< "\t" << genome.first
							<< "\t" << target.first
							<< "\t" << extract.fPrimerReg_.start_
							<< "\t" << extract.fPrimerReg_.end_
							<< "\t" << fPrimerMeta.getMeta("errors")
							<< "\t" << extract.gRegionInner_->start_
							<< "\t" << extract.gRegionInner_->end_
							<< "\t" << extract.rPrimerReg_.start_
							<< "\t" << extract.rPrimerReg_.end_
							<< "\t" << rPrimerMeta.getMeta("errors")
							<< std::endl;
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
				auto fullSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(primerDirectory, "all_" + target.first +".fasta"));
				SeqOutput::write(allSeqs, fullSeqOpts);
			}
			if(!allSeqsTrimmedSeqs.empty()){
				auto innerSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(primerDirectory, "all_" + target.first +"_primersRemoved.fasta"));
				SeqOutput::write(allSeqsTrimmedSeqs, innerSeqOpts);
			}

			if(!refSeqs.empty()){
				auto fullSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(primerDirectory, target.first +".fasta"));
				SeqOutput::write(refSeqs, fullSeqOpts);
			}
			if(!refTrimmedSeqs.empty()){
				auto innerSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(primerDirectory, target.first +"_primersRemoved.fasta"));
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
						target.first);
			}
			auto perTabOpts = TableIOOpts::genTabFileOut(njh::files::make_path(primerDirectory, "extractionCounts"),true);
			performanceTab.outPutContents(perTabOpts);
		}
	}else{
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
					bioRunner.bowtie2Align(seqOpts, genomeFnp	, "-D 20 -R 3 -N 1 -L 10 -i S,1,0.5 -a --end-to-end");
				}
				{
					BioCmdsUtils bioRunner;
					bioRunner.verbose_ = extractPars.verbose_;
					auto seqOpts = SeqIOOptions::genFastaIn(reverseFnp, false);
					seqOpts.out_.outFilename_ = njh::files::make_path(alignDir, bfs::basename(genomeFnp) + "_" + bfs::basename(reverseFnp) + ".sorted.bam");
					bioRunner.bowtie2Align(seqOpts, genomeFnp, "-D 20 -R 3 -N 1 -L 10 -i S,1,0.5 -a --end-to-end"	);
				}
			}
		};

		comparison allowableErrors;
		allowableErrors.hqMismatches_ = extractPars.errors;
		std::function<void()> extractPathway =
				[&targetsQueue,&extractPars,&outputDir,&ids,&alignToGenome,&allowableErrors,&gMapper]() {
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
						if(gMapper->pars_.numThreads_ <=1){
							alignToGenome(genomesQueue, forOutFnp, revOutFnp, refAlignmentDir);
						}else{
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
						}

						std::unordered_map<std::string, GenExtracRes> genomeExtractionsResults;
						for(const auto & genome : gMapper->genomes_){
							genomeExtractionsResults[genome.first] = GenExtracRes{};
						}
						std::unordered_map<std::string, std::vector<GenomeExtractResult>> genomeExtracts;
						for(const auto & genome : gMapper->genomes_) {
							auto forBamFnp = njh::files::make_path(refAlignmentDir,bfs::basename(genome.second->fnp_) + "_" + bfs::basename(forwardOpts.out_.outName()) + ".sorted.bam");
							auto revBamFnp = njh::files::make_path(refAlignmentDir,bfs::basename(genome.second->fnp_) + "_" + bfs::basename(reverseOpts.out_.outName()) + ".sorted.bam");
							auto forResults = gatherMapResults(forBamFnp, genome.second->fnpTwoBit_, allowableErrors);
							auto revResults = gatherMapResults(revBamFnp, genome.second->fnpTwoBit_, allowableErrors);
							auto primerLocationsFnp = njh::files::make_path(bedDirectory,bfs::basename(genome.second->fnp_) + "_primerLocs.bed");
							OutputStream primerLocationsOut(primerLocationsFnp);
							for(const auto & forRes : forResults){
								auto bedOut = forRes->gRegion_.genBedRecordCore();
								bedOut.name_ = njh::pasteAsStr(target, "-forwardPrimer[errors=", forRes->comp_.distances_.getNumOfEvents(true), ";]");
								primerLocationsOut << bedOut.toDelimStrWithExtra() << std::endl;
							}
							for(const auto & revRes : revResults){
								auto bedOut = revRes->gRegion_.genBedRecordCore();
								bedOut.name_ = njh::pasteAsStr(target, "-reversePrimer[errors=", revRes->comp_.distances_.getNumOfEvents(true), ";]");
								primerLocationsOut << bedOut.toDelimStrWithExtra() << std::endl;
							}
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
							OutputStream regionInfoOut{OutOptions(njh::files::make_path(bedDirectory, genome.first + "_regionInfo.tab.txt"))};
							regionInfoOut << "#chrom\tfullStart\tfullStop\tname\tlength\tstrand\tgenome\ttarget\tfPrimerStart\tfPrimerStop\terrorsInFPrimer\tinsertStart\tinsertStop\trPrimerStart\trPrimerStop\terrorsInRPrimer" << std::endl;
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
								regionInfoOut << extract.gRegion_->chrom_
										<< "\t" << extract.gRegion_->start_
										<< "\t" << extract.gRegion_->end_
										<< "\t" << name
										<< "\t" << extract.gRegion_->getLen()
										<< "\t" << (extract.gRegion_->reverseSrand_ ? '-' : '+')
										<< "\t" << genome.first
										<< "\t" << target
										<< "\t" << extract.ext_->gRegion_.start_
										<< "\t" << extract.ext_->gRegion_.end_
										<< "\t" << extract.ext_->comp_.distances_.getNumOfEvents(true)
										<< "\t" << extract.gRegionInner_->start_
										<< "\t" << extract.gRegionInner_->end_
										<< "\t" << extract.lig_->gRegion_.start_
										<< "\t" << extract.lig_->gRegion_.end_
										<< "\t" << extract.lig_->comp_.distances_.getNumOfEvents(true)
										<< std::endl;
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

		njh::concurrent::runVoidFunctionThreaded(extractPathway, numThreads);

	}


	auto locationsCombined = njh::files::make_path(outputDir, "locationsByGenome");
	njh::files::makeDir(njh::files::MkdirPar{locationsCombined});



	//primer location files
	for(const auto & genome : gMapper->genomes_){
		auto genomeBedOpts = njh::files::make_path(locationsCombined, genome.first + "_primersLocs.bed");
		std::vector<std::shared_ptr<Bed6RecordCore>> allPrimersRegions;
		for(const auto & tar : ids.targets_){
			auto bedForTarPrimersFnp = njh::files::make_path(outputDir, tar.first, "genomeLocations", genome.first + "_primerLocs.bed");
			if(bfs::exists(bedForTarPrimersFnp)){
				auto locs = getBeds(bedForTarPrimersFnp);
				addOtherVec(allPrimersRegions, locs);
			}
		}
		OutputStream genomeBedOut(genomeBedOpts);
		for(const auto & loc : allPrimersRegions){
			genomeBedOut << loc->toDelimStrWithExtra() << std::endl;
		}
	}

	njh::concurrent::LockableQueue<std::string> genomeQueue(getVectorOfMapKeys(gMapper->genomes_));
	std::function<void()> getOuterRegionInfos = [&genomeQueue,&gMapper,&extractPars,&locationsCombined,&outputDir,&ids](){
		std::string genome;
		while(genomeQueue.getVal(genome)){

			auto genomeBedOpts = njh::files::make_path(locationsCombined, genome + ".bed");
			std::vector<std::shared_ptr<Bed6RecordCore>> allRegions;
			for(const auto & tar : ids.targets_){
				auto bedForTarFnp = njh::files::make_path(outputDir, tar.first, "genomeLocations", genome + ".bed");
				if(bfs::exists(bedForTarFnp)){
					auto locs = getBeds(bedForTarFnp);
					addOtherVec(allRegions, locs);
				}
			}
			if(!allRegions.empty()){
				if("" != gMapper->genomes_.at(genome)->gffFnp_){
					intersectBedLocsWtihGffRecordsPars pars(gMapper->genomes_.at(genome)->gffFnp_);
					pars.selectFeatures_ = VecStr{"gene"};
					pars.extraAttributes_ = tokenizeString(extractPars.gffExtraAttributesStr, ",");
					auto jsonValues = intersectBedLocsWtihGffRecords(allRegions, pars);
					auto geneIds = jsonValues.getMemberNames();
					std::set<std::string> geneIdsSet(geneIds.begin(), geneIds.end());
					auto genes = GeneFromGffs::getGenesFromGffForIds(gMapper->genomes_.at(genome)->gffFnp_,geneIdsSet);
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
							TwoBit::TwoBitFile tReader(gMapper->genomes_.at(genome)->fnpTwoBit_);
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
	};
	njh::concurrent::runVoidFunctionThreaded(getOuterRegionInfos, extractPars.pars.numThreads_);

	struct InsertProteinInfo {
		InsertProteinInfo(const std::string & id, uint32_t aaStart, uint32_t aaStop,
				const std::string & desc) :
				id_(id), aaStart_(aaStart), aaStop_(aaStop), description_(desc) {
		}
		std::string id_;
		uint32_t aaStart_; //1-based
		uint32_t aaStop_; //1-based
		std::string description_;
	};
	std::unordered_map<std::string, std::vector<InsertProteinInfo>> proteinInsertInfoByName;
	std::mutex proteinInsertInfoByNameMut;
	njh::concurrent::LockableQueue<std::string> genomeQueueForInsert(getVectorOfMapKeys(gMapper->genomes_));
	std::function<void()> getInsertInfo = [&proteinInsertInfoByName,&genomeQueueForInsert,
																				 &gMapper,&extractPars,&locationsCombined,&outputDir,
																				 &proteinInsertInfoByNameMut,&ids](){
		std::string genome;
		std::unordered_map<std::string, std::vector<InsertProteinInfo>> proteinInsertInfoByNameCurrent;

		while(genomeQueueForInsert.getVal(genome)){
			auto genomeBedOpts = njh::files::make_path(locationsCombined, genome + "_inner.bed");
			std::vector<std::shared_ptr<Bed6RecordCore>> allRegions;
			for(const auto & tar : ids.targets_){
				auto bedForTarFnp = njh::files::make_path(outputDir, tar.first, "genomeLocations", genome + "_inner.bed");
				if(bfs::exists(bedForTarFnp)){
					auto locs = getBeds(bedForTarFnp);
					addOtherVec(allRegions, locs);
				}
			}
			if(!allRegions.empty()){
				if("" != gMapper->genomes_.at(genome)->gffFnp_){
					intersectBedLocsWtihGffRecordsPars pars(gMapper->genomes_.at(genome)->gffFnp_);
					pars.selectFeatures_ = VecStr{"gene"};
					pars.extraAttributes_ = tokenizeString(extractPars.gffExtraAttributesStr, ",");
					auto jsonValues = intersectBedLocsWtihGffRecords(allRegions, pars);
					auto geneIds = jsonValues.getMemberNames();
					std::set<std::string> geneIdsSet(geneIds.begin(), geneIds.end());
					auto genes = GeneFromGffs::getGenesFromGffForIds(gMapper->genomes_.at(genome)->gffFnp_,geneIdsSet);
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
							TwoBit::TwoBitFile tReader(gMapper->genomes_.at(genome)->fnpTwoBit_);
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
								proteinInsertInfoByNameCurrent[reg->name_].emplace_back(geneMeta.getMeta("ID"), aaStartPos, aaStopPos, geneMeta.getMeta("description"));
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
		{
			std::lock_guard<std::mutex> lock(proteinInsertInfoByNameMut);
			for(const auto & proteinInfo : proteinInsertInfoByNameCurrent){
				addOtherVec(proteinInsertInfoByName[proteinInfo.first],proteinInfo.second);
			}
		}
	};


	//inner with primer
	njh::concurrent::runVoidFunctionThreaded(getInsertInfo, extractPars.pars.numThreads_);

	//info files
	for(const auto & genome : gMapper->genomes_){
		auto targetGenomeInfoFnp = njh::files::make_path(locationsCombined, genome.first + "_infos.tab.txt");
		table outAllInfoTab;
		auto tarKeys = getVectorOfMapKeys(ids.targets_);
		njh::sort(tarKeys);
		for(const auto &  tar: tarKeys){
			auto infoFnp = njh::files::make_path(outputDir, tar, "genomeLocations", genome.first + "_regionInfo.tab.txt");
			if(bfs::exists(infoFnp)){
				table infoTab(infoFnp, "\t", true);
				if("" != genome.second->gffFnp_){
					table updatedInfoTab(toVecStr(infoTab.columnNames_, "insertGeneID", "insertGeneAAStart", "insertGeneAAStop", "insertGeneDescription"));
					for(auto & row : infoTab){
						if(proteinInsertInfoByName[row[infoTab.getColPos("name")]].empty()){
							updatedInfoTab.addRow(toVecStr(row, "", "", "", ""));
						}else{
							for(const auto & info : proteinInsertInfoByName[row[infoTab.getColPos("name")]]){
								updatedInfoTab.addRow(toVecStr(row, info.id_, info.aaStart_, info.aaStop_, info.description_));
							}
						}
					}
					infoTab = updatedInfoTab;
				}
				infoTab.addColumn({ids.targets_.at(tar).info_.forwardPrimerRaw_}, "Fwd_primer");
				infoTab.addColumn({ids.targets_.at(tar).info_.reversePrimerRaw_}, "Rev_primer");
				if(0 == outAllInfoTab.nRow()){
					outAllInfoTab = infoTab;
				}else{
					outAllInfoTab.rbind(infoTab, true);
				}
			}
		}
		OutputStream allInfoOut(targetGenomeInfoFnp);
		outAllInfoTab.outPutContents(allInfoOut, "\t");
	}

}


}  // namespace njhseq
