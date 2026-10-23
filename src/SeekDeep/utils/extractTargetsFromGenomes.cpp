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
#include <njhseq/objects/dataContainers/tables/TableReader.hpp>
#include <njhseq/objects/BioDataObject/BLASTHitTabular.hpp>
#include <njhseq/BamToolsUtils/ReAlignedSeq.hpp>
#include <njhseq/objects/BioDataObject/BioRecordsUtils/BedUtility.hpp>


namespace njhseq {


void extractBetweenSeqsPars::setUpCoreOptions(seqSetUp & setUp, bool needReadLength){

	setUp.setOption(longRangeAmplicon, "--longRangeAmplicon", "long Range Amplicon extraction");
	if(longRangeAmplicon){
		sizeLimit  = 8000;
		needReadLength = false;
	}

	setUp.setOption(doNotRenameBeds, "--doNotRenameBeds", "Do Not Rename Beds");

	setUp.setOption(writeOutAllSeqsFile, "--writeOutAllSeqsFile", "Write Out All Seqs File without collapsing to unique sequences");
	bool longNames = false;
	setUp.setOption(longNames, "--longNames", "Create long names for reference genomes extractions");
	shortNames = !longNames;
	setUp.setOption(shortNames, "--shortNames", "Create short names for reference genomes extractions");

	if(setUp.setOption(minLenCutOffSizeExpand, "--lenCutOffSizeExpand", "When creating length cut off file how much to expand the length of the found targets")){
		maxLenCutOffSizeExpand = minLenCutOffSizeExpand;
	}
	setUp.setOption(minLenCutOffSizeExpand, "--minLenCutOffSizeExpand", "When creating length cut off file how much to decrease the min length requirement of the found targets");
	setUp.setOption(maxLenCutOffSizeExpand, "--maxLenCutOffSizeExpand", "When creating length cut off file how much to increase the max length requirement the length of the found targets");

	setUp.setOption(pairedEndLength, "--pairedEndLength", "Paired End Read Length", needReadLength);
	setUp.setOption(barcodeSize, "--barcodeSize", "Barcode Size, if on both primer, the sum of the two barcodes");
	setUp.setOption(errors, "--errors", "Number of errors to allow in primers");
	setUp.setOption(sizeLimit, "--insertSizeLimit,--sizeLimit", "Output target extractions for only targets below this size");
	setUp.setOption(pars.numThreads_, "--numThreads", "Number of CPUs to utilize");
	setUp.setOption(removeRefAlignments, "--removeRefAlignments", "Remove Ref Alignments");
	setUp.setOption(pars.genomeDir_, "--genomeDir", "Name of the genome file fnp", true);
	setUp.setOption(pars.primaryGenome_, "--primaryGenome", "The primary reference genome");
	setUp.setOption(pars.gffDir_, "--gffDir", "A directory with a gff for the genomes in --genomeDir, should be named GENOME.gff (for GENOME.fasta)");
	setUp.setOption(gffExtraAttributesStr, "--gffExtraAttributes", "Extra attributes to add to genome that has an accompany gff");
	setUp.setOption(selectedGenomesStr, "--selectedGenomes", "Name of the other genomes in --genomeDir to be read in, leave blank to just do all fastas");
	setUp.setOption(outputDirPars.dirName_, "--dout", "Output directory", true);
	setUp.setOption(outputDirPars.overWriteDir_, "--overWriteDir", "Overwrite Output directory");

	setUp.setOption(useBlast, "--useBlast", "use blastn to search for primer locations");
	setUp.setOption(blastExpandSize, "--blastExpandSize", "when using blast, this is the expand size around the blast hit for the re-alignment");
}


namespace {

////////////////////////////////////////////////////////////////////////////////////////////////
// shared types
////////////////////////////////////////////////////////////////////////////////////////////////

struct GenExtracRes{
	uint32_t forwardHits_{0};
	uint32_t reverseHits_{0};
	uint32_t extractCounts_{0};
	uint32_t bestExtractionCounts_{0};
};

using GenExtracResByGenome = std::unordered_map<std::string, GenExtracRes>;

struct ExtractContext{
	const PrimersAndMids & ids_;
	const extractBetweenSeqsPars & extractPars_;
	MultiGenomeMapper & gMapper_;
	bfs::path outputDir_;
	uint32_t numThreads_; //threads used across targets, gMapper_.pars_.numThreads_ is used within a target
};

/**@brief The extraction of one primer pair from a genome, independent of how the primers were found
 */
struct ExtractionRecord{
	std::shared_ptr<GenomicRegion> gRegion_;
	std::shared_ptr<GenomicRegion> gRegionInner_;
	uint64_t fPrimerStart_{0};
	uint64_t fPrimerStop_{0};
	uint32_t fPrimerErrors_{0};
	uint64_t rPrimerStart_{0};
	uint64_t rPrimerStop_{0};
	uint32_t rPrimerErrors_{0};
};

//genome name, extractions found for that genome, in the order they should be written out
using GenomeRecords = std::vector<std::pair<std::string, std::vector<ExtractionRecord>>>;

ExtractionRecord toExtractionRecord(const GenomeExtractResult & extract){
	ExtractionRecord ret;
	ret.gRegion_ = extract.gRegion_;
	ret.gRegionInner_ = extract.gRegionInner_;
	ret.fPrimerStart_ = extract.extRegion_.start_;
	ret.fPrimerStop_ = extract.extRegion_.end_;
	ret.fPrimerErrors_ = extract.extComp_.distances_.getNumOfEvents(true);
	ret.rPrimerStart_ = extract.ligRegion_.start_;
	ret.rPrimerStop_ = extract.ligRegion_.end_;
	ret.rPrimerErrors_ = extract.ligComp_.distances_.getNumOfEvents(true);
	return ret;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// primer files
////////////////////////////////////////////////////////////////////////////////////////////////

struct PrimerFiles{
	bfs::path primerDirectory_;
	bfs::path bedDirectory_;
	bfs::path forwardFnp_;
	bfs::path reverseFnp_;
	std::vector<seqInfo> forSeqs_;
	std::vector<seqInfo> revSeqs_;
};

/**@brief Create the directories for a primer pair and write out all the (degenerate expanded) primer sequences
 */
PrimerFiles writePrimerFiles(const bfs::path & outputDir, const PrimerDeterminator::primerInfo & primerInfo){
	PrimerFiles ret;
	ret.primerDirectory_ = njh::files::makeDir(outputDir, njh::files::MkdirPar(primerInfo.primerPairName_));
	ret.bedDirectory_ = njh::files::makeDir(ret.primerDirectory_, njh::files::MkdirPar("genomeLocations"));
	auto forwardOpts = SeqIOOptions::genFastaOut(njh::files::make_path(ret.primerDirectory_, "forwardPrimer"));
	auto reverseOpts = SeqIOOptions::genFastaOut(njh::files::make_path(ret.primerDirectory_, "reversePrimer"));
	VecStr forDegens;
	for(const auto & fwd : primerInfo.fwds_){
		addOtherVec(forDegens, createDegenStrs(fwd.info_.seq_));
	}
	VecStr revDegens;
	for(const auto & rev : primerInfo.revs_){
		addOtherVec(revDegens, createDegenStrs(rev.info_.seq_));
	}
	if(forDegens.size() == 1){
		ret.forSeqs_.emplace_back(primerInfo.fwds_.front().info_);
	}else{
		ret.forSeqs_ = vecStrToReadObjs<seqInfo>(forDegens, primerInfo.primerPairName_);
	}
	if(revDegens.size() == 1){
		ret.revSeqs_.emplace_back(primerInfo.revs_.front().info_);
	}else{
		ret.revSeqs_ = vecStrToReadObjs<seqInfo>(revDegens, primerInfo.primerPairName_);
	}
	SeqOutput::write(ret.forSeqs_, forwardOpts);
	SeqOutput::write(ret.revSeqs_, reverseOpts);
	ret.forwardFnp_ = forwardOpts.out_.outName();
	ret.reverseFnp_ = reverseOpts.out_.outName();
	return ret;
}

/**@brief Run alignFunc(genomesQueue, forwardFnp, reverseFnp, alignDir) over all the genomes, threaded if gMapper has more than one thread
 */
template<typename ALIGNFUNC>
void alignPrimersToGenomes(const ExtractContext & ctx,
													 const PrimerFiles & primerFiles,
													 const bfs::path & alignDir,
													 const ALIGNFUNC & alignFunc){
	njh::concurrent::LockableQueue<bfs::path> genomesQueue(ctx.gMapper_.getGenomeFnps());
	if(ctx.gMapper_.pars_.numThreads_ <= 1){
		alignFunc(genomesQueue, primerFiles.forwardFnp_, primerFiles.reverseFnp_, alignDir);
	}else{
		std::vector<std::thread> threads;
		for(uint32_t t = 0; t < ctx.gMapper_.pars_.numThreads_; ++t){
			threads.emplace_back(alignFunc,
													 std::ref(genomesQueue),
													 std::cref(primerFiles.forwardFnp_),
													 std::cref(primerFiles.reverseFnp_),
													 std::cref(alignDir));
		}
		for(auto & t : threads){
			t.join();
		}
	}
}

/**@brief Write out the bed of primer locations for a genome
 */
template<typename RESULT>
void writePrimerLocations(const bfs::path & primerLocationsFnp,
													const std::string & target,
													const std::vector<std::shared_ptr<RESULT>> & forResults,
													const std::vector<std::shared_ptr<RESULT>> & revResults){
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
}

////////////////////////////////////////////////////////////////////////////////////////////////
// writing out the extractions for a target
////////////////////////////////////////////////////////////////////////////////////////////////

void writeCollapsedCounts(const bfs::path & outFnp,
													const std::vector<seqInfo> & allSeqs,
													const std::vector<seqInfo> & refSeqs){
	std::unordered_map<std::string, std::string> seqToNameKeys;
	std::unordered_map<std::string, uint32_t> seqCounts;
	std::unordered_map<std::string, std::set<std::string>> seqNames;
	for(const auto & seq : allSeqs){
		++seqCounts[seq.seq_];
		seqNames[seq.seq_].emplace(seq.name_);
	}
	for(const auto & refSeq : refSeqs){
		seqToNameKeys[refSeq.seq_] = refSeq.name_;
	}
	OutputStream collapsedSeqCountsOut(outFnp);
	VecStr refSeqNames = njh::getVecOfMapKeys(seqToNameKeys);
	//sort by counts
	njh::sort(refSeqNames,[&seqCounts](const std::string & n1, const std::string & n2){
		if(seqCounts[n1] == seqCounts[n2]){
			return n1 < n2;
		}else{
			return seqCounts[n1] > seqCounts[n2];
		}
	});
	collapsedSeqCountsOut << "collapsedName\tcount\textractedName" << std::endl;
	for(const auto & name : refSeqNames){
		for(const auto & eName : seqNames[name]){
			collapsedSeqCountsOut
					<< seqToNameKeys[name]
					<< "\t" << seqCounts[name]
					<< "\t" << eName << std::endl;
		}
	}
}

/**@brief The extracted sequences, both with and without the primers, for either all the extractions or just the best ones
 */
struct ExtractedSeqs{
	//sequences collapsed down to unique sequences, the names of the identical sequences get concatenated
	std::vector<seqInfo> refSeqs_;
	std::vector<seqInfo> refTrimmedSeqs_;
	//every sequence, only filled if --writeOutAllSeqsFile
	std::vector<seqInfo> allSeqs_;
	std::vector<seqInfo> allTrimmedSeqs_;

	static void collapseInto(std::vector<seqInfo> & refs, const seqInfo & seq){
		bool found = false;
		for(auto & ref : refs){
			if(ref.seq_ == seq.seq_){
				found = true;
				ref.name_.append("-" + seq.name_);
			}
		}
		if(!found){
			refs.emplace_back(seq);
		}
	}

	void add(const seqInfo & fullSeq, const seqInfo & innerSeq, bool keepAllSeqs){
		if(keepAllSeqs){
			allSeqs_.emplace_back(fullSeq);
			allTrimmedSeqs_.emplace_back(innerSeq);
		}
		collapseInto(refSeqs_, fullSeq);
		collapseInto(refTrimmedSeqs_, innerSeq);
	}

	/**@brief Write out the sequences
	 * @param suffix "_all" for the every extraction, "" for the best extractions
	 */
	void write(const bfs::path & dir, const std::string & target, const std::string & suffix) const {
		if(!allSeqs_.empty()){
			auto fullSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(dir, "separated_" + target + suffix + ".fasta"));
			SeqOutput::write(allSeqs_, fullSeqOpts);
			writeCollapsedCounts(njh::files::make_path(dir, target + "_collapsed_counts" + suffix + ".tab.txt"), allSeqs_, refSeqs_);
		}
		if(!allTrimmedSeqs_.empty()){
			auto innerSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(dir, "separated_" + target + "_primersRemoved" + suffix + ".fasta"));
			SeqOutput::write(allTrimmedSeqs_, innerSeqOpts);
			writeCollapsedCounts(njh::files::make_path(dir, target + "_primersRemoved_collapsed_counts" + suffix + ".tab.txt"), allTrimmedSeqs_, refTrimmedSeqs_);
		}
		if(!refSeqs_.empty()){
			auto fullSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(dir, target + suffix + ".fasta"));
			SeqOutput::write(refSeqs_, fullSeqOpts);
		}
		if(!refTrimmedSeqs_.empty()){
			auto innerSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(dir, target + "_primersRemoved" + suffix + ".fasta"));
			SeqOutput::write(refTrimmedSeqs_, innerSeqOpts);
		}
	}
};

const std::string regionInfoHeader = "#chrom\tfullStart\tfullStop\tname\tlength\tstrand\tgenome\ttarget\tfPrimerStart\tfPrimerStop\terrorsInFPrimer\tinsertStart\tinsertStop\trPrimerStart\trPrimerStop\terrorsInRPrimer";

void writeRegionInfoRow(OutputStream & out,
												const ExtractionRecord & extract,
												const std::string & name,
												const std::string & genomeName,
												const std::string & targetName){
	out << extract.gRegion_->chrom_
			<< "\t" << extract.gRegion_->start_
			<< "\t" << extract.gRegion_->end_
			<< "\t" << name
			<< "\t" << extract.gRegion_->getLen()
			<< "\t" << (extract.gRegion_->reverseSrand_ ? '-' : '+')
			<< "\t" << genomeName
			<< "\t" << targetName
			<< "\t" << extract.fPrimerStart_
			<< "\t" << extract.fPrimerStop_
			<< "\t" << extract.fPrimerErrors_
			<< "\t" << extract.gRegionInner_->start_
			<< "\t" << extract.gRegionInner_->end_
			<< "\t" << extract.rPrimerStart_
			<< "\t" << extract.rPrimerStop_
			<< "\t" << extract.rPrimerErrors_
			<< std::endl;
}

/**@brief Write the bed record of region to the all out and, if it is one of the best, to the best out too
 */
void writeBedRecord(const GenomicRegion & region,
										const std::string & allName, OutputStream & allOut,
										const std::string & bestName, OutputStream & bestOut,
										bool inBest){
	auto bedRegion = region.genBedRecordCore();
	bedRegion.name_ = allName;
	allOut << bedRegion.toDelimStr() << std::endl;
	if(inBest){
		bedRegion.name_ = bestName;
		bestOut << bedRegion.toDelimStr() << std::endl;
	}
}

/**@brief Find the extractions with the fewest primer errors, ties broken by the fewest errors in the worst primer
 * @return the coordinate uids of the best extractions
 */
std::unordered_set<std::string> findBestExtractions(const std::vector<ExtractionRecord> & records){
	std::unordered_set<std::string> bestHits;
	uint32_t bestHitPrimerErrorCount = std::numeric_limits<uint32_t>::max();
	uint32_t minPerPrimerHit = std::numeric_limits<uint32_t>::max();
	for(const auto & extract : records){
		uint32_t current_minPerPrimerHit = std::min(extract.fPrimerErrors_, extract.rPrimerErrors_);
		uint32_t current_primerErrorCount = extract.fPrimerErrors_ + extract.rPrimerErrors_;
		if (current_primerErrorCount < bestHitPrimerErrorCount ||
				(current_primerErrorCount == bestHitPrimerErrorCount && current_minPerPrimerHit < minPerPrimerHit)) {
			bestHits.clear();
			bestHits.emplace(extract.gRegion_->createUidFromCoordsStrand());
			bestHitPrimerErrorCount = current_primerErrorCount;
			minPerPrimerHit = current_minPerPrimerHit;
		} else if (current_primerErrorCount == bestHitPrimerErrorCount && current_minPerPrimerHit == minPerPrimerHit) {
			bestHits.emplace(extract.gRegion_->createUidFromCoordsStrand());
		}
	}
	return bestHits;
}

std::string makeExtractionName(const std::string & genomeName,
															 const std::string & metaTargetName,
															 uint32_t extractionCount){
	auto name = genomeName;
	if(0 != extractionCount){
		name.append("." + estd::to_string(extractionCount));
	}
	MetaDataInName meta;
	meta.addMeta("genome", genomeName);
	meta.addMeta("target", metaTargetName);
	meta.addMeta("extractionCount", extractionCount);
	name += " " + meta.createMetaName();
	return name;
}

/**@brief Write out everything for the extractions of one target: beds, region infos, fastas and the extraction counts table
 *
 * @param targetName the name used for output file names and the target column
 * @param metaTargetName the name used for the target in the extraction names
 * @param genomeRecords the extractions by genome, genomes without extractions are skipped
 * @param genomeExtractionsResults the per genome hit counts, extraction counts are filled in here
 */
void writeTargetExtractions(const ExtractContext & ctx,
														const std::string & targetName,
														const std::string & metaTargetName,
														const bfs::path & primerDirectory,
														const bfs::path & bedDirectory,
														const GenomeRecords & genomeRecords,
														GenExtracResByGenome & genomeExtractionsResults){
	const auto & extractPars = ctx.extractPars_;
	ExtractedSeqs bestSeqs;
	ExtractedSeqs allSeqs;

	for(const auto & genome : genomeRecords){
		const auto & genomeName = genome.first;
		const auto & records = genome.second;
		if(records.empty()){
			continue;
		}
		genomeExtractionsResults[genomeName].extractCounts_ = records.size();
		uint32_t all_extractionCount = 0;
		uint32_t best_extractionCount = 0;

		OutputStream all_bedOut{OutOptions(njh::files::make_path(bedDirectory, genomeName + "_all.bed"))};
		OutputStream all_bedInnerOut{OutOptions(njh::files::make_path(bedDirectory, genomeName + "_inner_all.bed"))};
		OutputStream all_regionInfoOut{OutOptions(njh::files::make_path(bedDirectory, genomeName + "_regionInfo_all.tab.txt"))};
		OutputStream best_bedOut{OutOptions(njh::files::make_path(bedDirectory, genomeName + ".bed"))};
		OutputStream best_bedInnerOut{OutOptions(njh::files::make_path(bedDirectory, genomeName + "_inner.bed"))};
		OutputStream best_regionInfoOut{OutOptions(njh::files::make_path(bedDirectory, genomeName + "_regionInfo.tab.txt"))};
		all_regionInfoOut << regionInfoHeader << std::endl;
		best_regionInfoOut << regionInfoHeader << std::endl;

		const auto bestHits = findBestExtractions(records);
		TwoBit::TwoBitFile tReader(ctx.gMapper_.genomes_.at(genomeName)->fnpTwoBit_);

		for(const auto & extract : records){
			bool inBest = njh::in(extract.gRegion_->createUidFromCoordsStrand(), bestHits);
			auto best_name = makeExtractionName(genomeName, metaTargetName, best_extractionCount);
			auto all_name = makeExtractionName(genomeName, metaTargetName, all_extractionCount);

			writeRegionInfoRow(all_regionInfoOut, extract, all_name, genomeName, targetName);
			if (inBest) {
				writeRegionInfoRow(best_regionInfoOut, extract, best_name, genomeName, targetName);
			}
			writeBedRecord(*extract.gRegion_, all_name, all_bedOut, best_name, best_bedOut, inBest);
			writeBedRecord(*extract.gRegionInner_, all_name, all_bedInnerOut, best_name, best_bedInnerOut, inBest);

			if(extractPars.shortNames){
				trimAtFirstWhitespace(all_name);
				trimAtFirstWhitespace(best_name);
			}
			auto fullSeq = extract.gRegion_->extractSeq(tReader);
			auto innerSeq = extract.gRegionInner_->extractSeq(tReader);
			fullSeq.name_ = all_name;
			innerSeq.name_ = all_name;
			allSeqs.add(fullSeq, innerSeq, extractPars.writeOutAllSeqsFile);
			if (inBest){
				fullSeq.name_ = best_name;
				innerSeq.name_ = best_name;
				bestSeqs.add(fullSeq, innerSeq, extractPars.writeOutAllSeqsFile);
			}

			++all_extractionCount;
			if (inBest) {
				++best_extractionCount;
			}
		}
		genomeExtractionsResults[genomeName].bestExtractionCounts_ = best_extractionCount;
	}
	allSeqs.write(primerDirectory, targetName, "_all");
	bestSeqs.write(primerDirectory, targetName, "");

	table performanceTab(VecStr{"genome", "forwardPrimerHits", "reversePrimerHits", "extractionCounts", "allExtractionCounts", "target"});
	auto genomeKeys = getVectorOfMapKeys(genomeExtractionsResults);
	njh::sort(genomeKeys);
	for(const auto & genomeKey : genomeKeys){
		performanceTab.addRow(genomeKey,
													genomeExtractionsResults[genomeKey].forwardHits_,
													genomeExtractionsResults[genomeKey].reverseHits_,
													genomeExtractionsResults[genomeKey].bestExtractionCounts_,
													genomeExtractionsResults[genomeKey].extractCounts_,
													targetName);
	}
	auto perTabOpts = TableIOOpts::genTabFileOut(njh::files::make_path(primerDirectory, "extractionCounts"),true);
	performanceTab.outPutContents(perTabOpts);
}

GenExtracResByGenome makeEmptyGenomeResults(const ExtractContext & ctx){
	GenExtracResByGenome ret;
	for(const auto & genome : ctx.gMapper_.genomes_){
		ret[genome.first] = GenExtracRes{};
	}
	return ret;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// finding the primers with blast
////////////////////////////////////////////////////////////////////////////////////////////////

/**@brief One primer sequence in the combined query file that all of the targets' primers are blasted with
 */
struct BlastQuery{
	std::string target_;
	bool forward_;
	std::string name_; //the name of the primer in the target's primer file
};

/**@brief Blast the primers of all the targets against each genome with a single blast run per genome
 *
 * The genome's blast database is only made and read once per genome rather than once per primer file. The hits are then
 * split back out into the per target hit files, with the original primer names, that the re-alignment reads in.
 *
 * @param primerFilesByTarget the primer files made for each target
 * @param refAlignmentDirByTarget where the hits for each target should be written
 */
void blastAllPrimersToGenomes(const ExtractContext & ctx,
															const std::map<std::string, PrimerFiles> & primerFilesByTarget,
															const std::map<std::string, bfs::path> & refAlignmentDirByTarget){
	//write out all the primers in one file, the names are replaced with a simple unique id
	std::vector<BlastQuery> queries;
	std::vector<seqInfo> querySeqs;
	for(const auto & primerFiles : primerFilesByTarget){
		for(const bool forward : {true, false}){
			for(const auto & seq : (forward ? primerFiles.second.forSeqs_ : primerFiles.second.revSeqs_)){
				queries.emplace_back(BlastQuery{primerFiles.first, forward, seq.name_});
				querySeqs.emplace_back(seq);
				querySeqs.back().name_ = "q" + estd::to_string(queries.size() - 1);
			}
		}
	}
	auto blastDir = njh::files::makeDir(ctx.outputDir_, njh::files::MkdirPar{"allPrimersBlast"});
	auto queryOpts = SeqIOOptions::genFastaOut(njh::files::make_path(blastDir, "allPrimers"));
	SeqOutput::write(querySeqs, queryOpts);
	const bfs::path queryFnp = queryOpts.out_.outName();

	//one blast per genome
	auto genomeFnps = ctx.gMapper_.getGenomeFnps();
	const uint32_t totalThreads = std::max<uint32_t>(1, ctx.extractPars_.pars.numThreads_);
	const uint32_t genomeThreads = std::max<uint32_t>(1, std::min<uint32_t>(genomeFnps.size(), totalThreads));
	const uint32_t blastThreads = std::max<uint32_t>(1, totalThreads/genomeThreads);
	njh::concurrent::LockableQueue<bfs::path> genomesQueue(genomeFnps);
	std::function<void()> blastGenomes = [&](){
		bfs::path genomeFnp = "";
		while(genomesQueue.getVal(genomeFnp)){
			BioCmdsUtils bioRunner;
			bioRunner.verbose_ = ctx.extractPars_.verbose_;
			auto blastDatabaseFnp = bfs::path(genomeFnp).replace_extension("").string() + "_blastdb";
			bioRunner.RunMakeblastdb(genomeFnp);
			bfs::path outputFnp = njh::files::make_path(blastDir, bfs::basename(genomeFnp) + "_allPrimers_blastHits.tsv");
			std::string cmd = "blastn -query " + queryFnp.string()  + " -perc_identity 80 -qcov_hsp_perc 80 -task blastn-short -word_size 6 -db " + blastDatabaseFnp + " -outfmt 6 -max_target_seqs 100000000 -num_threads " + estd::to_string(blastThreads) + " > " + outputFnp.string();
			bioRunner.runCmdCheck(cmd, queryFnp, outputFnp);
		}
	};
	njh::concurrent::runVoidFunctionThreaded(blastGenomes, genomeThreads);

	//split the hits back out by target and primer direction
	for(const auto & genomeFnp : genomeFnps){
		//target, forward -> the lines of the hits, a line for every query so that targets without hits still get a (empty) file
		std::map<std::pair<std::string, bool>, std::string> hitLinesByTargetByDirection;
		for(const auto & query : queries){
			hitLinesByTargetByDirection[{query.target_, query.forward_}];
		}
		std::ifstream hitsFile(njh::files::make_path(blastDir, bfs::basename(genomeFnp) + "_allPrimers_blastHits.tsv").string());
		std::string line;
		while(std::getline(hitsFile, line)){
			auto tabPos = line.find('\t');
			if(std::string::npos == tabPos){
				continue;
			}
			const auto & query = queries.at(std::stoul(line.substr(1, tabPos - 1)));
			auto & lines = hitLinesByTargetByDirection[{query.target_, query.forward_}];
			lines.append(query.name_);
			lines.append(line, tabPos, std::string::npos);
			lines.push_back('\n');
		}
		for(const auto & hitLines : hitLinesByTargetByDirection){
			const auto & primerFiles = primerFilesByTarget.at(hitLines.first.first);
			const auto & primerFnp = hitLines.first.second ? primerFiles.forwardFnp_ : primerFiles.reverseFnp_;
			OutputStream hitsOut(njh::files::make_path(refAlignmentDirByTarget.at(hitLines.first.first),
																								 bfs::basename(genomeFnp) + "_" + bfs::basename(primerFnp) + "_blastHits.tsv"));
			hitsOut << hitLines.second;
		}
	}
}

/**@brief Re-align the primers to the regions of the blast hits, keeping the full length ones with allowable errors
 */
std::vector<std::shared_ptr<ReAlignedSeq>> gatherBlastRealignments(const bfs::path & blastHitsFnp,
																																	 const std::unordered_map<std::string, std::string> & primersByName,
																																	 aligner & alignerObj,
																																	 TwoBit::TwoBitFile & tReader,
																																	 const ReAlignedSeq::genRealignmentPars & reAlnPars,
																																	 const comparison & allowableErrors){
	std::vector<std::shared_ptr<ReAlignedSeq>> ret;
	BioDataFileIO<BLASTHitTab> reader{IoOptions(InOptions(blastHitsFnp))};
	reader.openIn();
	BLASTHitTab hit;
	while(reader.readNextRecord(hit)){
		auto realignedSeq = ReAlignedSeq::genRealignment(hit, njh::mapAt(primersByName, hit.queryName_), alignerObj, tReader.getSeqLens(), tReader, reAlnPars);
		if(1 == realignedSeq.comp_.distances_.query_.coverage_ && allowableErrors.passErrorProfile(realignedSeq.comp_)){
			ret.emplace_back(std::make_shared<ReAlignedSeq>(realignedSeq));
		}
	}
	return ret;
}

std::unordered_map<std::string, std::string> primerSeqsByName(const std::vector<seqInfo> & seqs){
	std::unordered_map<std::string, std::string> ret;
	for(const auto & seq : seqs){
		ret[seq.name_] = seq.seq_;
	}
	return ret;
}

void extractTargetsWithBlast(const ExtractContext & ctx,
														 njh::concurrent::LockableQueue<std::string> & targetsQueue){
	const auto & extractPars = ctx.extractPars_;
	comparison allowableErrors;
	allowableErrors.hqMismatches_ = extractPars.errors;

	//set up the directories and primer files for all the targets and then blast all of the primers at once
	std::map<std::string, PrimerFiles> primerFilesByTarget;
	std::map<std::string, bfs::path> refAlignmentDirByTarget;
	for(const auto & tar : ctx.ids_.targets_){
		const auto & primerInfo = ctx.ids_.pDeterminator_->primers_.at(tar.first);
		auto primerFiles = writePrimerFiles(ctx.outputDir_, primerInfo);
		refAlignmentDirByTarget[tar.first] = njh::files::makeDir(primerFiles.primerDirectory_, njh::files::MkdirPar{"refAlignments"});
		primerFilesByTarget.emplace(tar.first, std::move(primerFiles));
	}
	blastAllPrimersToGenomes(ctx, primerFilesByTarget, refAlignmentDirByTarget);

	std::function<void()> extractPathway = [&]() {
		std::string target;
		while(targetsQueue.getVal(target)) {
			const auto & primerInfo = ctx.ids_.pDeterminator_->primers_.at(target);
			const auto & primerFiles = primerFilesByTarget.at(target);
			const auto & refAlignmentDir = refAlignmentDirByTarget.at(target);

			uint64_t maxLen = 0;
			readVec::getMaxLength(primerFiles.forSeqs_, maxLen);
			readVec::getMaxLength(primerFiles.revSeqs_, maxLen);
			aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,5,1,0,0,5,1), substituteMatrix(2, -2));
			auto forwardPrimersByName = primerSeqsByName(primerFiles.forSeqs_);
			auto reversePrimersByName = primerSeqsByName(primerFiles.revSeqs_);

			auto genomeExtractionsResults = makeEmptyGenomeResults(ctx);
			ReAlignedSeq::genRealignmentPars reAlnPars;
			reAlnPars.extendAmount = extractPars.blastExpandSize;
			reAlnPars.adjustForSoftClipping = true;
			std::unordered_map<std::string, std::vector<GenomeExtractResult>> genomeExtracts;
			for(const auto & genome : ctx.gMapper_.genomes_) {
				TwoBit::TwoBitFile tReader(genome.second->fnpTwoBit_);
				auto forBlastHitsFnp = njh::files::make_path(refAlignmentDir, bfs::basename(genome.second->fnp_) + "_" + bfs::basename(primerFiles.forwardFnp_) + "_blastHits.tsv");
				auto revBlastHitsFnp = njh::files::make_path(refAlignmentDir, bfs::basename(genome.second->fnp_) + "_" + bfs::basename(primerFiles.reverseFnp_) + "_blastHits.tsv");
				auto forRealignments = gatherBlastRealignments(forBlastHitsFnp, forwardPrimersByName, alignerObj, tReader, reAlnPars, allowableErrors);
				auto forResults = ReAlignedSeq::getUniqueLocationResults(forRealignments);
				auto revRealignments = gatherBlastRealignments(revBlastHitsFnp, reversePrimersByName, alignerObj, tReader, reAlnPars, allowableErrors);
				auto revResults = ReAlignedSeq::getUniqueLocationResults(revRealignments);

				writePrimerLocations(njh::files::make_path(primerFiles.bedDirectory_, bfs::basename(genome.second->fnp_) + "_primerLocs.bed"),
														 target, forResults, revResults);

				genomeExtractionsResults[genome.first].forwardHits_ = forResults.size();
				genomeExtractionsResults[genome.first].reverseHits_ = revResults.size();
				if(!forResults.empty() && !revResults.empty()){
					auto uniForRes = ReAlignedSeq::getUniqueLocationResults(forResults);
					auto uniRevRes = ReAlignedSeq::getUniqueLocationResults(revResults);
					genomeExtracts[genome.first] = getPossibleGenomeExtracts(uniForRes, uniRevRes, extractPars.sizeLimit);
				}
			}
			GenomeRecords genomeRecords;
			for(const auto & genome : genomeExtracts){
				std::vector<ExtractionRecord> records;
				for(const auto & extract : genome.second){
					records.emplace_back(toExtractionRecord(extract));
				}
				genomeRecords.emplace_back(genome.first, std::move(records));
			}
			writeTargetExtractions(ctx, target, primerInfo.primerPairName_,
														 primerFiles.primerDirectory_, primerFiles.bedDirectory_,
														 genomeRecords, genomeExtractionsResults);
			if(extractPars.removeRefAlignments){
				njh::files::rmDirForce(refAlignmentDir);
			}
		}
	};
	njh::concurrent::runVoidFunctionThreaded(extractPathway, ctx.numThreads_);
	if(extractPars.removeRefAlignments){
		njh::files::rmDirForce(ctx.outputDir_ / "allPrimersBlast");
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////
// finding the primers with bowtie2, exact matches only
////////////////////////////////////////////////////////////////////////////////////////////////

void bowtie2AlignPrimersToGenomes(const extractBetweenSeqsPars & extractPars,
																	njh::concurrent::LockableQueue<bfs::path> & genomesQueue,
																	const bfs::path & forwardFnp,
																	const bfs::path & reverseFnp,
																	const bfs::path & alignDir){
	bfs::path genomeFnp = "";
	while(genomesQueue.getVal(genomeFnp)){
		for(const auto & primerFnp : {forwardFnp, reverseFnp}){
			BioCmdsUtils bioRunner;
			bioRunner.verbose_ = extractPars.verbose_;
			auto seqOpts = SeqIOOptions::genFastaIn(primerFnp, false);
			seqOpts.out_.outFilename_ = njh::files::make_path(alignDir, bfs::basename(genomeFnp) + "_" + bfs::basename(primerFnp) + ".sorted.bam");
			bioRunner.bowtie2Align(seqOpts, genomeFnp, "-D 20 -R 3 -N 1 -L 10 -i S,1,0.5 -a --end-to-end");
		}
	}
}

void extractTargetsWithBowtie2(const ExtractContext & ctx,
															 njh::concurrent::LockableQueue<std::string> & targetsQueue){
	const auto & extractPars = ctx.extractPars_;
	comparison allowableErrors;
	allowableErrors.hqMismatches_ = extractPars.errors;
	auto bowtie2Align = [&extractPars](njh::concurrent::LockableQueue<bfs::path> & genomesQueue,
																		 const bfs::path & forwardFnp,
																		 const bfs::path & reverseFnp,
																		 const bfs::path & alignDir){
		bowtie2AlignPrimersToGenomes(extractPars, genomesQueue, forwardFnp, reverseFnp, alignDir);
	};

	std::function<void()> extractPathway = [&]() {
		std::string target;
		while(targetsQueue.getVal(target)) {
			const auto & primerInfo = ctx.ids_.pDeterminator_->primers_.at(target);
			auto primerFiles = writePrimerFiles(ctx.outputDir_, primerInfo);

			auto refAlignmentDir = njh::files::makeDir(primerFiles.primerDirectory_, njh::files::MkdirPar{"refAlignments"});
			alignPrimersToGenomes(ctx, primerFiles, refAlignmentDir, bowtie2Align);

			auto genomeExtractionsResults = makeEmptyGenomeResults(ctx);
			std::unordered_map<std::string, std::vector<GenomeExtractResult>> genomeExtracts;
			for(const auto & genome : ctx.gMapper_.genomes_) {
				auto forBamFnp = njh::files::make_path(refAlignmentDir, bfs::basename(genome.second->fnp_) + "_" + bfs::basename(primerFiles.forwardFnp_) + ".sorted.bam");
				auto revBamFnp = njh::files::make_path(refAlignmentDir, bfs::basename(genome.second->fnp_) + "_" + bfs::basename(primerFiles.reverseFnp_) + ".sorted.bam");
				auto forResults = gatherMapResults(forBamFnp, genome.second->fnpTwoBit_, allowableErrors);
				auto revResults = gatherMapResults(revBamFnp, genome.second->fnpTwoBit_, allowableErrors);
				writePrimerLocations(njh::files::make_path(primerFiles.bedDirectory_, bfs::basename(genome.second->fnp_) + "_primerLocs.bed"),
														 target, forResults, revResults);
				genomeExtractionsResults[genome.first].forwardHits_ = forResults.size();
				genomeExtractionsResults[genome.first].reverseHits_ = revResults.size();
				if(!forResults.empty() && !revResults.empty()){
					auto uniForRes = getUniqueLocationResults(forResults);
					auto uniRevRes = getUniqueLocationResults(revResults);
					genomeExtracts[genome.first] = getPossibleGenomeExtracts(uniForRes, uniRevRes, extractPars.sizeLimit);
				}
			}
			GenomeRecords genomeRecords;
			for(const auto & genome : genomeExtracts){
				std::vector<ExtractionRecord> records;
				for(const auto & extract : genome.second){
					records.emplace_back(toExtractionRecord(extract));
				}
				genomeRecords.emplace_back(genome.first, std::move(records));
			}
			writeTargetExtractions(ctx, target, primerInfo.primerPairName_,
														 primerFiles.primerDirectory_, primerFiles.bedDirectory_,
														 genomeRecords, genomeExtractionsResults);
			if(extractPars.removeRefAlignments){
				njh::files::rmDirForce(refAlignmentDir);
			}
		}
	};
	njh::concurrent::runVoidFunctionThreaded(extractPathway, ctx.numThreads_);
}

////////////////////////////////////////////////////////////////////////////////////////////////
// finding the primers by searching the genome's kmers, allows for errors in the primers
////////////////////////////////////////////////////////////////////////////////////////////////

struct PrimerPairSearchResults{
	std::vector<GenomicRegion> fPrimerPositions_;
	std::vector<GenomicRegion> rPrimerPositions_;

	class GenomeExtractResultByPrimerPair {
	public:
		GenomeExtractResultByPrimerPair(GenomicRegion fPrimerReg,
																		GenomicRegion rPrimerReg): fPrimerReg_(std::move(fPrimerReg)),rPrimerReg_(std::move(rPrimerReg)){
			setRegion();
		}
		GenomicRegion fPrimerReg_;
		GenomicRegion rPrimerReg_;

		std::shared_ptr<GenomicRegion> gRegion_ = nullptr;// it's silly but this is due to eclipse for now, should be able to revert once cdt 9.3 is released june 28 2017
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

	static std::vector<GenomeExtractResultByPrimerPair> getPossibleGenomeExtracts(const std::vector<GenomicRegion> & fPrimerPositions,
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
					//the forward primer has to be upstream of the reverse primer, in the direction of the strand
					if(fwd.reverseSrand_ ? fwd.start_ > rev.start_ : fwd.start_ < rev.start_){
						GenomeExtractResultByPrimerPair extraction(fwd, rev);
						if (extraction.gRegion_->getLen() <= insertSizeCutOff) {
							ret.emplace_back(extraction);
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
using SearchResultsByTargetByGenome = std::map<std::string, std::map<std::string, PrimerPairSearchResults>>;

uint32_t getErrorsFromPrimerUid(const GenomicRegion & region){
	return MetaDataInName(region.uid_.substr(region.uid_.rfind('['))).getMeta<uint32_t>("errors");
}

ExtractionRecord toExtractionRecord(const PrimerPairSearchResults::GenomeExtractResultByPrimerPair & extract){
	ExtractionRecord ret;
	ret.gRegion_ = extract.gRegion_;
	ret.gRegionInner_ = extract.gRegionInner_;
	ret.fPrimerStart_ = extract.fPrimerReg_.start_;
	ret.fPrimerStop_ = extract.fPrimerReg_.end_;
	ret.fPrimerErrors_ = getErrorsFromPrimerUid(extract.fPrimerReg_);
	ret.rPrimerStart_ = extract.rPrimerReg_.start_;
	ret.rPrimerStop_ = extract.rPrimerReg_.end_;
	ret.rPrimerErrors_ = getErrorsFromPrimerUid(extract.rPrimerReg_);
	return ret;
}

/**@brief Add the positions of all the kmers that match any of the primers
 *
 * @param reverseStrand search with the reverse complement of the primers
 * @param label forwardPrimer or reversePrimer, goes in the name of the regions
 */
template<typename PRIMERS>
void addKmerPrimerHits(std::vector<GenomicRegion> & positions,
											 const PRIMERS & primers,
											 bool reverseStrand,
											 const kmerInfo & info,
											 uint32_t allowedErrors,
											 const std::string & primerPairName,
											 const std::string & label,
											 const std::string & chromName){
	for(const auto & k : info.kmers_){
		for(const auto & primer : primers){
			const auto & mot = reverseStrand ? primer.motRC_ : primer.mot_;
			if(mot.frontPassNoCheck(k.first, allowedErrors)){
				auto score = mot.scoreMotif(k.first.begin(), k.first.begin() + mot.size());
				auto errors = mot.size() - score;
				for(const auto & pos : k.second.positions_){
					auto start = pos;
					auto end = pos + mot.size();
					positions.emplace_back(GenomicRegion(njh::pasteAsStr(primerPairName, "-", label, "[errors=", errors, ";]"),
																							 chromName,
																							 start, end, reverseStrand));
				}
			}
		}
	}
}

/**@brief Remove positions that are the same region as a position with fewer errors
 */
std::vector<GenomicRegion> uniquePositionsByFewestErrors(std::vector<GenomicRegion> positions){
	njh::sort(positions, [](const GenomicRegion & reg1, const GenomicRegion & reg2){
		return getErrorsFromPrimerUid(reg1) < getErrorsFromPrimerUid(reg2);
	});
	std::vector<GenomicRegion> ret;
	for(const auto & position : positions){
		bool add = true;
		for(const auto & alreadyAdded : ret){
			if(alreadyAdded.sameRegion(position)){
				add = false;
				break;
			}
		}
		if(add){
			ret.emplace_back(position);
		}
	}
	return ret;
}

bool genomicRegionCoordSort(const GenomicRegion & reg1, const GenomicRegion & reg2){
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
}

void addSearchResults(PrimerPairSearchResults & to, const PrimerPairSearchResults & from){
	addOtherVec(to.fPrimerPositions_, from.fPrimerPositions_);
	addOtherVec(to.rPrimerPositions_, from.rPrimerPositions_);
	addOtherVec(to.regions_, from.regions_);
}

/**@brief Search all of a genome's chromosomes for the primers of all the targets
 */
void searchGenomeForPrimers(const ExtractContext & ctx,
														const std::string & genomeName,
														const bfs::path & genomeTwoBitFnp,
														uint32_t maxPrimerSize,
														uint32_t minPrimerSize,
														SearchResultsByTargetByGenome & resultsByTargetByGenome,
														std::mutex & resultsMut){
	const auto & extractPars = ctx.extractPars_;
	TwoBit::TwoBitFile globalTreader(genomeTwoBitFnp);
	auto names = globalTreader.sequenceNames();
	njh::concurrent::LockableQueue<std::string> nameQueue(names);

	std::function<void()> extractGenomicPositions = [&](){
		SearchResultsByTargetByGenome resultsCurrent;
		std::string seqName = "";
		TwoBit::TwoBitFile treader(genomeTwoBitFnp);
		seqInfo seq;
		while(nameQueue.getVal(seqName)){
			treader[seqName]->getSequence(seq.seq_);
			seq.name_ = seqName;
			const auto & chromName = seqName;

			kmerInfo info(seq.seq_, maxPrimerSize, false);
			if(minPrimerSize < maxPrimerSize && seq.seq_.size() > maxPrimerSize){
				//add the kmers at the end of the chromosome that are shorter than the max primer size
				for(const auto pos : iter::range(minPrimerSize, maxPrimerSize)){
					auto currentK = seq.seq_.substr(len(seq) - pos);
					auto k = info.kmers_.find(currentK);
					if (k != info.kmers_.end()) {
						k->second.addPosition(len(seq) - pos);
					} else {
						info.kmers_.emplace(currentK, kmer(currentK, len(seq) - pos));
					}
				}
			}

			for(const auto & primer : ctx.ids_.pDeterminator_->primers_){
				const auto & pairName = primer.second.primerPairName_;
				std::vector<GenomicRegion> fPrimerPositions;
				std::vector<GenomicRegion> rPrimerPositions;
				//searching with `5-`3 for both primers so in order to match later the strands need to be opposite
				for(const bool reverseStrand : {false, true}){
					addKmerPrimerHits(fPrimerPositions, primer.second.fwds_, reverseStrand, info, extractPars.errors, pairName, "forwardPrimer", chromName);
					addKmerPrimerHits(rPrimerPositions, primer.second.revs_, reverseStrand, info, extractPars.errors, pairName, "reversePrimer", chromName);
				}
				auto fPrimerPositionsUnique = uniquePositionsByFewestErrors(fPrimerPositions);
				auto rPrimerPositionsUnique = uniquePositionsByFewestErrors(rPrimerPositions);
				auto & current = resultsCurrent[pairName][genomeName];
				addOtherVec(current.fPrimerPositions_, fPrimerPositionsUnique);
				addOtherVec(current.rPrimerPositions_, rPrimerPositionsUnique);
				addOtherVec(current.regions_, PrimerPairSearchResults::getPossibleGenomeExtracts(fPrimerPositionsUnique, rPrimerPositionsUnique, extractPars.sizeLimit));
			}
		}
		{
			std::lock_guard<std::mutex> lock(resultsMut);
			for(const auto & tar : resultsCurrent){
				for(const auto & gen : tar.second){
					addSearchResults(resultsByTargetByGenome[tar.first][genomeName], gen.second);
				}
			}
		}
	};
	njh::concurrent::runVoidFunctionThreaded(extractGenomicPositions, extractPars.pars.numThreads_);
}

void extractTargetsWithKmerSearch(const ExtractContext & ctx){
	const auto & primers = ctx.ids_.pDeterminator_->primers_;
	const auto maxPrimerSize = ctx.ids_.pDeterminator_->getMaxPrimerSize();
	const auto minPrimerSize = ctx.ids_.pDeterminator_->getMinPrimerSize();

	//setup directories and write out the primers
	for(const auto & target : primers){
		writePrimerFiles(ctx.outputDir_, target.second);
	}

	SearchResultsByTargetByGenome resultsByTargetByGenome;
	std::mutex resultsMut;
	for(const auto & genome : ctx.gMapper_.genomes_){
		searchGenomeForPrimers(ctx, genome.first, genome.second->fnpTwoBit_, maxPrimerSize, minPrimerSize, resultsByTargetByGenome, resultsMut);
	}

	for(auto & target : resultsByTargetByGenome){
		for(auto & genome : target.second){
			njh::sort(genome.second.fPrimerPositions_, genomicRegionCoordSort);
			njh::sort(genome.second.rPrimerPositions_, genomicRegionCoordSort);
			njh::sort(genome.second.regions_, [](const PrimerPairSearchResults::GenomeExtractResultByPrimerPair & reg1,
																					 const PrimerPairSearchResults::GenomeExtractResultByPrimerPair & reg2){
				return genomicRegionCoordSort(*reg1.gRegion_, *reg2.gRegion_);
			});
		}
	}

	for(const auto & target : resultsByTargetByGenome){
		auto primerDirectory = njh::files::make_path(ctx.outputDir_, target.first);
		auto bedDirectory = njh::files::make_path(primerDirectory, "genomeLocations");

		auto genomeExtractionsResults = makeEmptyGenomeResults(ctx);
		for(const auto & genome : ctx.gMapper_.genomes_) {
			auto primerLocationsFnp = njh::files::make_path(bedDirectory, bfs::basename(genome.second->fnp_) + "_primerLocs.bed");
			OutputStream primerLocationsOut(primerLocationsFnp);
			std::vector<GenomicRegion> allPrimerRegions;
			if(njh::in(genome.first, target.second)){
				const auto & results = target.second.at(genome.first);
				addOtherVec(allPrimerRegions, results.fPrimerPositions_);
				genomeExtractionsResults[genome.first].forwardHits_ = results.fPrimerPositions_.size();
				addOtherVec(allPrimerRegions, results.rPrimerPositions_);
				genomeExtractionsResults[genome.first].reverseHits_ = results.rPrimerPositions_.size();
			}
			njh::sort(allPrimerRegions, genomicRegionCoordSort);
			for(const auto & primerRes : allPrimerRegions){
				primerLocationsOut << primerRes.genBedRecordCore().toDelimStrWithExtra() << std::endl;
			}
		}
		GenomeRecords genomeRecords;
		for(const auto & genome : target.second){
			std::vector<ExtractionRecord> records;
			for(const auto & extract : genome.second.regions_){
				records.emplace_back(toExtractionRecord(extract));
			}
			genomeRecords.emplace_back(genome.first, std::move(records));
		}
		writeTargetExtractions(ctx, target.first, target.first, primerDirectory, bedDirectory, genomeRecords, genomeExtractionsResults);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////
// combining the per target outputs into the per genome outputs
////////////////////////////////////////////////////////////////////////////////////////////////

using BedRecords = std::vector<std::shared_ptr<Bed6RecordCore>>;

void combineExtractionCounts(const ExtractContext & ctx){
	auto extractionCountsCombinedOpts = njh::files::make_path(ctx.outputDir_, "allExtractionCounts.tab.txt");
	OutputStream allExtCountOut(extractionCountsCombinedOpts);
	uint32_t count = 0;
	for(const auto & tar : ctx.ids_.targets_){
		auto extractionCountsFnp = njh::files::make_path(ctx.outputDir_, tar.first, "extractionCounts.tab.txt");
		if(bfs::exists(extractionCountsFnp)){
			TableReader extReader(TableIOOpts::genTabFileIn(extractionCountsFnp, true));
			if(0 == count){
				extReader.header_.outPutContents(allExtCountOut, "\t");
			}
			VecStr row;
			while(extReader.getNextRow(row)){
				allExtCountOut << njh::conToStr(row, "\t") << "\n";
			}
			++count;
		}
	}
}

void combinePrimerLocations(const ExtractContext & ctx, const bfs::path & locationsCombined){
	for(const auto & genome : ctx.gMapper_.genomes_){
		auto genomeBedOpts = njh::files::make_path(locationsCombined, genome.first + "_primersLocs.bed");
		BedRecords allPrimersRegions;
		for(const auto & tar : ctx.ids_.targets_){
			auto bedForTarPrimersFnp = njh::files::make_path(ctx.outputDir_, tar.first, "genomeLocations", genome.first + "_primerLocs.bed");
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
}

/**@brief Read in the per target beds of a genome, the beds are named genomeLocations/GENOME + suffix
 *
 * Unless --doNotRenameBeds, the regions are renamed to the target name (target.1, target.2 etc. for the repeats)
 */
BedRecords readTargetBeds(const ExtractContext & ctx, const std::string & genome, const std::string & suffix){
	BedRecords ret;
	for(const auto & tar : ctx.ids_.targets_){
		auto bedForTarFnp = njh::files::make_path(ctx.outputDir_, tar.first, "genomeLocations", genome + suffix);
		if(bfs::exists(bedForTarFnp)){
			auto locs = getBeds(bedForTarFnp);
			if (!ctx.extractPars_.doNotRenameBeds) {
				uint32_t count = 0;
				for (auto & loc : locs) {
					if (0 == count) {
						loc->name_ = tar.first;
					} else {
						loc->name_ = njh::pasteAsStr(tar.first, ".", count);
					}
					++count;
				}
			}
			addOtherVec(ret, locs);
		}
	}
	return ret;
}

void sortBedsByChromStart(BedRecords & beds){
	njh::sort(beds, [](const std::shared_ptr<Bed3RecordCore> & bed1,
										 const std::shared_ptr<Bed3RecordCore> & bed2) {
		if (bed1->chrom_ == bed2->chrom_) {
			return bed1->chromStart_ < bed2->chromStart_;
		} else {
			return bed1->chrom_ < bed2->chrom_;
		}
	});
}

void writeBeds(const bfs::path & outFnp, const BedRecords & beds){
	OutputStream bedOut(outFnp);
	for(const auto & loc : beds){
		bedOut << loc->toDelimStrWithExtra() << std::endl;
	}
}

/**@brief Get the genes in the genome's gff that overlap the regions, ignoring the non-protein coding genes
 */
std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> getOverlappingProteinCodingGenes(const bfs::path & gffFnp,
																																															 BedRecords & regions,
																																															 const std::string & gffExtraAttributesStr){
	intersectBedLocsWtihGffRecordsPars pars(gffFnp);
	pars.selectFeatures_ = VecStr{"gene", "protein_coding_gene"};
	pars.extraAttributes_ = tokenizeString(gffExtraAttributesStr, ",");
	auto jsonValues = intersectBedLocsWtihGffRecords(regions, pars);
	auto rawGeneIds = jsonValues.getMemberNames();
	std::set<std::string> rawGeneIdsSet(rawGeneIds.begin(), rawGeneIds.end());

	auto rawGenes = GeneFromGffs::getGenesFromGffForIds(gffFnp, rawGeneIdsSet);
	std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> genes;
	for(const auto & gene : rawGenes){
		bool failFilter = false;
		for(const auto & transcript : gene.second->mRNAs_){
			if(njh::in(njh::strToLowerRet(transcript->type_), VecStr{"rrna", "trna", "snorna","snrna","ncrna"} ) ){
				failFilter = true;
				break;
			}
		}
		if(!failFilter){
			genes[gene.first] = gene.second;
		}
	}
	return genes;
}

/**@brief Split the last extra field of a region, which is a list of "[meta],[meta]", into the individual "[meta]"s
 */
VecStr splitGeneMetaToks(const std::string & extraField){
	auto geneToks = njh::tokenizeString(extraField, "],[");
	if(geneToks.size() > 1){
		geneToks.front().append("]");
		for(const auto pos : iter::range<uint32_t>(1, geneToks.size() - 1)){
			geneToks[pos] = "[" + geneToks[pos] + "]";
		}
		geneToks.back() = "[" +geneToks.back();
	}
	return geneToks;
}

/**@brief Get the 1 based amino acid start and stop of the part of the transcript that the region covers
 *
 * @return start and stop, both are max() if the region is completely within an intron
 */
template<typename POSINFOS>
std::pair<uint32_t, uint32_t> getAAStartStop(POSINFOS & posInfos, const Bed6RecordCore & reg){
	auto minPos = vectorMinimum(getVectorOfMapKeys(posInfos));
	auto maxPos = vectorMaximum(getVectorOfMapKeys(posInfos));
	auto startPos = std::max(minPos, reg.chromStart_);
	auto stopPos =  std::min(maxPos, reg.chromEnd_);
	if(stopPos == reg.chromEnd_){
		stopPos = reg.chromEnd_ - 1;
	}
	//stop position is inclusive
	uint32_t aaStartPos = std::numeric_limits<uint32_t>::max();
	uint32_t aaStopPos = 0;

	for(const auto posInGene : iter::range(startPos, stopPos + 1)){
		if(std::numeric_limits<uint32_t>::max() != posInfos[posInGene].aaPos_){
			if(posInfos[posInGene].aaPos_ < aaStartPos){
				aaStartPos = posInfos[posInGene].aaPos_;
			}
			if(posInfos[posInGene].aaPos_ > aaStopPos) {
				aaStopPos = posInfos[posInGene].aaPos_;
			}
		}
	}
	//if the region is completely within the intron the final stop and stops will be max() and 0;
	if(std::numeric_limits<uint32_t>::max() == aaStartPos){
		aaStopPos = std::numeric_limits<uint32_t>::max();
	}else{
		//make 1 based
		++aaStopPos;
		++aaStartPos;
	}
	return {aaStartPos, aaStopPos};
}

/**@brief Replace the gene info in the last extra field of the regions with the gene info plus the AA positions covered by the region
 */
void addAAPositionsToRegions(const ExtractContext & ctx, const std::string & genome, BedRecords & allRegions){
	const auto & genomeInfo = ctx.gMapper_.genomes_.at(genome);
	auto genes = getOverlappingProteinCodingGenes(genomeInfo->gffFnp_, allRegions, ctx.extractPars_.gffExtraAttributesStr);
	for(const auto & reg : allRegions){
		if(reg->extraFields_.empty()){
			continue;
		}
		std::string replacement;
		for(const auto & geneTok : splitGeneMetaToks(reg->extraFields_.back())){
			MetaDataInName geneMeta(geneTok);
			if(njh::in(geneMeta.getMeta("ID"), genes)){
				TwoBit::TwoBitFile tReader(genomeInfo->fnpTwoBit_);
				auto infos = njh::mapAt(genes, geneMeta.getMeta("ID"))->generateGeneSeqInfo(tReader, false);
				for(const auto & info : infos){
					auto posInfos = info.second->getInfosByGDNAPos();
					auto aaStartStop = getAAStartStop(posInfos, *reg);
					geneMeta.addMeta(info.first + "-AAStart", aaStartStop.first);
					geneMeta.addMeta(info.first + "-AAStop", aaStartStop.second);
				}
				if(!replacement.empty()){
					replacement += ",";
				}
				replacement += geneMeta.createMetaName();
			}
		}
		reg->extraFields_.back() = replacement;
	}
}

/**@brief Same as addAAPositionsToRegions but also add the gene descriptions and record the protein info by the region's coordinates
 */
void addAAPositionsAndDescriptionsToRegions(const ExtractContext & ctx,
																						const std::string & genome,
																						BedRecords & allRegions,
																						std::unordered_map<std::string, std::vector<MultiGenomeMapper::IntersectedProteinInfo>> & proteinInsertInfoByName){
	const auto & genomeInfo = ctx.gMapper_.genomes_.at(genome);
	auto genes = getOverlappingProteinCodingGenes(genomeInfo->gffFnp_, allRegions, ctx.extractPars_.gffExtraAttributesStr);
	for(const auto & reg : allRegions){
		if(reg->extraFields_.empty()){
			continue;
		}
		std::string replacement;
		for(const auto & geneTok : splitGeneMetaToks(reg->extraFields_.back())){
			MetaDataInName geneMeta(geneTok);
			if(njh::in(geneMeta.getMeta("ID"), genes)){
				TwoBit::TwoBitFile tReader(genomeInfo->fnpTwoBit_);
				const auto & gene = njh::mapAt(genes, geneMeta.getMeta("ID"));
				auto infos = gene->generateGeneSeqInfo(tReader, false);
				auto detailedName = gene->getGeneDetailedName();
				geneMeta.addMeta("detailedDescription", detailedName.begin()->second);
				for(const auto & info : infos){
					auto posInfos = info.second->getInfosByGDNAPos();
					auto aaStartStop = getAAStartStop(posInfos, *reg);
					geneMeta.addMeta(info.first + "-AAStart", aaStartStop.first);
					geneMeta.addMeta(info.first + "-AAStop", aaStartStop.second);
					std::string description = detailedName[info.first];
					std::string geneName = geneMeta.getMeta("ID");
					if (gene->gene_->hasAttr("Name") && "NA" != gene->gene_->getAttr("Name") && !gene->gene_->getAttr("Name").empty()) {
						geneName = gene->gene_->getAttr("Name");
					}
					proteinInsertInfoByName[reg->genUIDFromCoords()].emplace_back(
							geneMeta.getMeta("ID"),
							info.first,
							geneName, aaStartStop.first, aaStartStop.second, description);
				}
			}
			//unlike addAAPositionsToRegions, genes not in the passing gene set are kept as they were
			if(!replacement.empty()){
				replacement += ",";
			}
			replacement += geneMeta.createMetaName();
		}
		reg->extraFields_.back() = replacement;
	}
}

/**@brief Combine the per target beds of the full extractions (primers included) into per genome beds
 */
void combineOuterRegions(const ExtractContext & ctx, const bfs::path & locationsCombined, const bfs::path & locationsAllHitsCombined){
	njh::concurrent::LockableQueue<std::string> genomeQueue(getVectorOfMapKeys(ctx.gMapper_.genomes_));
	std::function<void()> combineOuter = [&](){
		std::string genome;
		while(genomeQueue.getVal(genome)){
			auto allRegions = readTargetBeds(ctx, genome, ".bed");
			{
				auto allHitsRegions = readTargetBeds(ctx, genome, "_all.bed");
				sortBedsByChromStart(allHitsRegions);
				writeBeds(njh::files::make_path(locationsAllHitsCombined, genome + ".bed"), allHitsRegions);
			}
			if(!allRegions.empty()){
				if(!ctx.gMapper_.genomes_.at(genome)->gffFnp_.empty()){
					addAAPositionsToRegions(ctx, genome, allRegions);
				}
				sortBedsByChromStart(allRegions);
				writeBeds(njh::files::make_path(locationsCombined, genome + ".bed"), allRegions);
			}
		}
	};
	njh::concurrent::runVoidFunctionThreaded(combineOuter, ctx.extractPars_.pars.numThreads_);
}

/**@brief Combine the per target beds of the inserts (primers removed) into per genome beds
 *
 * @return the protein info for the inserts by the insert's coordinates
 */
std::unordered_map<std::string, std::vector<MultiGenomeMapper::IntersectedProteinInfo>> combineInsertRegions(
		const ExtractContext & ctx, const bfs::path & locationsCombined, const bfs::path & locationsAllHitsCombined){
	std::unordered_map<std::string, std::vector<MultiGenomeMapper::IntersectedProteinInfo>> proteinInsertInfoByName;
	std::mutex proteinInsertInfoByNameMut;
	njh::concurrent::LockableQueue<std::string> genomeQueue(getVectorOfMapKeys(ctx.gMapper_.genomes_));
	std::function<void()> combineInserts = [&](){
		std::unordered_map<std::string, std::vector<MultiGenomeMapper::IntersectedProteinInfo>> proteinInsertInfoByNameCurrent;
		std::string genome;
		while(genomeQueue.getVal(genome)){
			auto allRegions = readTargetBeds(ctx, genome, "_inner.bed");
			{
				auto allHitsRegions = readTargetBeds(ctx, genome, "_inner_all.bed");
				sortBedsByChromStart(allHitsRegions);
				writeBeds(njh::files::make_path(locationsAllHitsCombined, genome + "_inner.bed"), allHitsRegions);
			}
			if(!allRegions.empty()){
				if(!ctx.gMapper_.genomes_.at(genome)->gffFnp_.empty()){
					addAAPositionsAndDescriptionsToRegions(ctx, genome, allRegions, proteinInsertInfoByNameCurrent);
				}
				BedUtility::coordSort(allRegions);
				writeBeds(njh::files::make_path(locationsCombined, genome + "_inner.bed"), allRegions);
			}
		}
		{
			std::lock_guard<std::mutex> lock(proteinInsertInfoByNameMut);
			for(const auto & proteinInfo : proteinInsertInfoByNameCurrent){
				addOtherVec(proteinInsertInfoByName[proteinInfo.first], proteinInfo.second);
			}
		}
	};
	njh::concurrent::runVoidFunctionThreaded(combineInserts, ctx.extractPars_.pars.numThreads_);
	return proteinInsertInfoByName;
}

/**@brief Add the info on the genes the insert is in to the region info table
 */
table addInsertGeneInfo(const table & infoTab,
												std::unordered_map<std::string, std::vector<MultiGenomeMapper::IntersectedProteinInfo>> & proteinInsertInfoByName){
	table updatedInfoTab(toVecStr(infoTab.columnNames_, "insertGeneID", "insertGeneName", "insertTranscriptID", "insertGeneAAStart", "insertGeneAAStop", "insertGeneDescription"));
	for(const auto & row : infoTab){
		//coord name
		std::string coordName = njh::pasteAsStr(row[infoTab.getColPos("#chrom")], "-", row[infoTab.getColPos("insertStart")], "-", row[infoTab.getColPos("insertStop")]);
		if(proteinInsertInfoByName[coordName].empty()){
			updatedInfoTab.addRow(toVecStr(row, "", "", "", "", "", ""));
		}else{
			std::set<std::string> insertGeneID;
			std::set<std::string> insertGeneName;
			std::set<std::string> insertGeneDescription;
			VecStr insertTranscriptID;
			VecStr insertGeneAAStart;
			VecStr insertGeneAAStop;
			for(const auto & info : proteinInsertInfoByName[coordName]) {
				insertGeneID.emplace(info.gene_id_);
				insertGeneName.emplace(info.Name_);
				insertGeneDescription.emplace(info.description_);
				insertTranscriptID.emplace_back(info.id_);
				insertGeneAAStart.emplace_back(std::numeric_limits<uint32_t>::max() == info.aaStart_ ? "NA" : njh::pasteAsStr(info.aaStart_));
				insertGeneAAStop.emplace_back(std::numeric_limits<uint32_t>::max() == info.aaStop_ ? "NA" : njh::pasteAsStr(info.aaStop_));
			}
			updatedInfoTab.addRow(toVecStr(row,
																		 njh::conToStr(insertGeneID, "::"),
																		 njh::conToStr(insertGeneName, "::"),
																		 njh::conToStr(insertTranscriptID, "::"),
																		 njh::conToStr(insertGeneAAStart, "::"),
																		 njh::conToStr(insertGeneAAStop, "::"),
																		 njh::conToStr(insertGeneDescription, "::")));
		}
	}
	return updatedInfoTab;
}

/**@brief Combine the per target region info tables into per genome tables
 */
void combineRegionInfos(const ExtractContext & ctx,
												const bfs::path & locationsCombined,
												const bfs::path & locationsAllHitsCombined,
												std::unordered_map<std::string, std::vector<MultiGenomeMapper::IntersectedProteinInfo>> & proteinInsertInfoByName){
	auto tarKeys = getVectorOfMapKeys(ctx.ids_.targets_);
	njh::sort(tarKeys);
	for(const auto & genome : ctx.gMapper_.genomes_){
		{
			//all hits
			auto targetGenomeInfoFnp = njh::files::make_path(locationsAllHitsCombined, genome.first + "_infos.tab.txt");
			table outAllInfoTab;
			for(const auto & tar : tarKeys) {
				auto infoFnp = njh::files::make_path(ctx.outputDir_, tar, "genomeLocations", genome.first + "_regionInfo_all.tab.txt");
				if(bfs::exists(infoFnp)) {
					table infoTab(infoFnp, "\t", true);
					if (0 == outAllInfoTab.nRow()) {
						outAllInfoTab = infoTab;
					} else {
						outAllInfoTab.rbind(infoTab, true);
					}
				}
			}
			outAllInfoTab.outPutContents(TableIOOpts::genTabFileOut(targetGenomeInfoFnp));
		}
		//best hits
		auto targetGenomeInfoFnp = njh::files::make_path(locationsCombined, genome.first + "_infos.tab.txt");
		table outAllInfoTab;
		for(const auto & tar : tarKeys){
			auto infoFnp = njh::files::make_path(ctx.outputDir_, tar, "genomeLocations", genome.first + "_regionInfo.tab.txt");
			if(bfs::exists(infoFnp)){
				table infoTab(infoFnp, "\t", true);
				if(!genome.second->gffFnp_.empty()){
					infoTab = addInsertGeneInfo(infoTab, proteinInsertInfoByName);
				}
				infoTab.addColumn({ctx.ids_.targets_.at(tar).info_.forwardPrimerRaw_}, "Fwd_primer");
				infoTab.addColumn({ctx.ids_.targets_.at(tar).info_.reversePrimerRaw_}, "Rev_primer");
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

}  // namespace


void extractBetweenSeqs(const PrimersAndMids & ids,
												const extractBetweenSeqsPars & extractPars){

	njh::concurrent::LockableQueue<std::string> targetsQueue(getVectorOfMapKeys(ids.targets_));

	//check for existence of genome directory
	njh::files::checkExistenceThrow(extractPars.pars.genomeDir_, __PRETTY_FUNCTION__);

	//make output directory
	njh::files::makeDir(extractPars.outputDirPars);
	bfs::path outputDir = extractPars.outputDirPars.dirName_;

	//set up genome mapper;
	MultiGenomeMapper gMapper(extractPars.pars);
	//set up selected genomes
	gMapper.setSelectedGenomes(extractPars.selectedGenomesStr);
	//init is threaded
	gMapper.pars_.numThreads_ = extractPars.pars.numThreads_;
	gMapper.init();
	gMapper.pars_.numThreads_ = 1;

	//set threads;
	uint32_t numThreads = extractPars.pars.numThreads_;
	if(extractPars.pars.numThreads_ >= 4){
		gMapper.pars_.numThreads_ = 2;
		numThreads = extractPars.pars.numThreads_/2;
	}
	if(1 == ids.targets_.size()){
		gMapper.pars_.numThreads_ = extractPars.pars.numThreads_;
		numThreads = 1;
	}

	ExtractContext ctx{ids, extractPars, gMapper, outputDir, numThreads};

	//find the primers in the genomes and write out the extractions for each target
	if(extractPars.useBlast){
		extractTargetsWithBlast(ctx, targetsQueue);
	}else if(extractPars.errors > 0){
		extractTargetsWithKmerSearch(ctx);
	}else{
		extractTargetsWithBowtie2(ctx, targetsQueue);
	}

	//combine the per target outputs into per genome outputs
	combineExtractionCounts(ctx);

	auto locationsCombined = njh::files::make_path(outputDir, "locationsByGenome");
	njh::files::makeDir(njh::files::MkdirPar{locationsCombined});
	auto locationsAllHitsCombined = njh::files::make_path(locationsCombined, "allHits");
	njh::files::makeDir(njh::files::MkdirPar{locationsAllHitsCombined});

	combinePrimerLocations(ctx, locationsCombined);
	combineOuterRegions(ctx, locationsCombined, locationsAllHitsCombined);
	auto proteinInsertInfoByName = combineInsertRegions(ctx, locationsCombined, locationsAllHitsCombined);
	combineRegionInfos(ctx, locationsCombined, locationsAllHitsCombined, proteinInsertInfoByName);
}


}  // namespace njhseq
