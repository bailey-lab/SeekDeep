//
// Created by Nicholas Hathaway on 3/6/23.
//
#include "SeekDeepPrograms/SeekDeepProgram/SeekDeepRunner.hpp"

#include <njhseq/objects/kmer/SimpleKmerHash.hpp>
#include <njhseq/objects/kmer/KmerSetInvertedIndex.hpp>
#include <njhseq/objects/kmer/KmerCollectionReadMatcher.hpp>

namespace njhseq {

namespace {

/**
 * @brief Where each target sits along a read that hit several targets, from the positions of k-mers unique to one target
 *
 * Reads that join several amplicons end to end (e.g. amplicons ligated to each other before barcoding) hold each amplicon's
 * unique k-mers in one stretch. Hits of the same target and orientation are chained while they are within maxGap_ of each
 * other and the chain stays no longer than the target can be, so two copies of one target in a row stay separate; chains with
 * fewer than minKmers_ hits are dropped as noise.
 */
struct MultihitSplitter {
	struct Pars {
		bool split_ = true;
		uint32_t minKmers_ = 20;      ///< unique k-mer hits a stretch needs to count as an amplicon
		uint32_t maxGap_ = 100;       ///< largest gap between unique k-mer hits within one stretch
		uint32_t primerSlop_ = 20;    ///< sequence kept past the longest primer on each side of a stretch
		std::string oneEnd_ = "all"; ///< reads demultiplexed on one end: split none of them, keep the amplicon nearest that barcode, or keep all
	};

	struct Segment {
		uint32_t setIdx_;
		bool revComp_;
		uint32_t start_;  ///< first base of the first k-mer
		uint32_t end_;    ///< one past the last base of the last k-mer
		uint32_t kmers_;
		uint32_t lastPos_;
	};

	/**
	 * @param allowed by index order, the targets a stretch may come from
	 * @param maxSpan by index order, the longest a stretch of that target can be
	 */
	static std::vector<Segment> findSegments(const std::string & seq, const KmerSetInvertedIndex & index, const SimpleKmerHash & hasher,
			uint32_t klen, const std::vector<char> & allowed, const std::vector<uint32_t> & maxSpan, const Pars & pars) {
		std::vector<Segment> done;
		if (seq.size() <= klen) {
			return done;
		}
		//the stretch currently open for each target and orientation, index*2 + revComp
		std::vector<int64_t> open(index.setNames_.size() * 2, -1);
		std::vector<Segment> all;
		for (uint32_t pos = 0; pos < seq.size() - klen + 1; ++pos) {
			for (const bool revComp : {false, true}) {
				const uint64_t hash = revComp ? hasher.revCompHash(seq, pos, klen) : hasher.hash(seq, pos, klen);
				const auto search = index.kmerToSet_.find(hash);
				if (search == index.kmerToSet_.end() || (search->second & KmerSetInvertedIndex::SHARED_FLAG)) {
					continue;
				}
				const uint32_t setIdx = search->second;
				if (!allowed[setIdx]) {
					continue;
				}
				auto & openIdx = open[setIdx * 2 + revComp];
				if (openIdx >= 0) {
					auto & seg = all[openIdx];
					if (pos <= seg.lastPos_ + pars.maxGap_ && pos + klen - seg.start_ <= maxSpan[setIdx]) {
						seg.lastPos_ = pos;
						seg.end_ = pos + klen;
						++seg.kmers_;
						continue;
					}
				}
				openIdx = static_cast<int64_t>(all.size());
				all.emplace_back(Segment{setIdx, revComp, pos, pos + klen, 1, pos});
			}
		}
		for (const auto & seg : all) {
			if (seg.kmers_ >= pars.minKmers_) {
				done.emplace_back(seg);
			}
		}
		std::sort(done.begin(), done.end(), [](const Segment & a, const Segment & b) { return a.start_ < b.start_; });
		return done;
	}
};

/**
 * @brief A region to trim the reads of an amplified target to, marked off by its own pair of primers
 *
 * The region is found in a read by its k-mers, on their own so regions can overlap, the trimming primers are then looked for
 * on either side of that stretch the same way the amplified primers are.
 */
struct TrimRegion {
	std::string name_;
	uint32_t maxSpan_{std::numeric_limits<uint32_t>::max()}; ///< the longest the region can be
	uint32_t forwardPrimerLen_{0};
	uint32_t reversePrimerLen_{0};
	bool revComp_{false}; ///< whether the region is in the opposite direction to the amplified target, e.g. a gene on the other strand
	std::shared_ptr<const KmerSetInvertedIndex> index_; ///< this region's k-mers only

	/**
	 * @brief The part of a read, in the region's forward direction, that holds this region
	 * @return the first base and one past the last base, 0,0 if the region isn't in the read
	 */
	std::pair<uint32_t, uint32_t> find(const std::string & seq, const SimpleKmerHash & hasher, uint32_t klen, const MultihitSplitter::Pars & pars) const {
		const auto segments = MultihitSplitter::findSegments(seq, *index_, hasher, klen, std::vector<char>{1}, std::vector<uint32_t>{maxSpan_}, pars);
		//the stretch with the most k-mers, with the stretches beside it that still fit within the longest the region can be (an error dense
		//patch breaks a region up), a second copy of the region further along doesn't fit
		std::vector<MultihitSplitter::Segment> forward;
		for(const auto & seg : segments){
			if(!seg.revComp_){
				forward.emplace_back(seg);
			}
		}
		if(forward.empty()){
			return {0, 0};
		}
		size_t best = 0;
		for(size_t idx = 1; idx < forward.size(); ++idx){
			if(forward[idx].kmers_ > forward[best].kmers_){
				best = idx;
			}
		}
		uint32_t start = forward[best].start_;
		uint32_t end = forward[best].end_;
		for(size_t idx = best; idx-- > 0;){
			if(std::max(end, forward[idx].end_) - forward[idx].start_ > maxSpan_){
				break;
			}
			start = std::min(start, forward[idx].start_);
			end = std::max(end, forward[idx].end_);
		}
		for(size_t idx = best + 1; idx < forward.size(); ++idx){
			if(std::max(end, forward[idx].end_) - start > maxSpan_){
				break;
			}
			end = std::max(end, forward[idx].end_);
		}
		return {start, end};
	}
};

/**
 * @brief Cuts sequencing adapters out of reads, found by exact k-mers of the adapters in either direction
 *
 * Each adapter k-mer in a read places the whole adapter on the read (from where in the adapter the k-mer sits), placements that
 * overlap are one adapter, and one with at least minKmerHits_ k-mers counts. The read is cut at every adapter found and the longest
 * piece left is kept, e.g. a nanopore rapid kit read whose transposase cut inside an amplicon starts with the rapid adapter followed
 * by part of the amplicon, which without the cut keeps the adapter in front of the amplicon's trimming region.
 */
struct AdapterTrimmer {
	struct Pars {
		bool trim_ = true;
		uint32_t kmerLen_ = 16;     ///< bases in a row that have to match an adapter, at most 32
		uint32_t minKmerHits_ = 2;  ///< adapter k-mers needed to call an adapter
		VecStr adapters_;

		Pars() {
			adapters_ = VecStr{
				//nanopore rapid adapter (rapid and rapid barcoding kits, e.g. SQK-RAD004, SQK-RBK004, SQK-RBK114), the part next to the insert
				"GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
				//nanopore ligation adapter, kit 14 (e.g. SQK-LSK114, SQK-NBD114)
				"CCTGTACTTCGTTCAGTTACGTATTGCT",
				//nanopore ligation adapter, kits 9 and 10 (e.g. SQK-LSK109, SQK-LSK110)
				"AATGTACTTCGTTCAGTTACGTATTGCT",
				//nanopore PCR barcoding tails (e.g. SQK-PBK004, EXP-PBC096), in front of the target primers
				"TTAACCTTTCTGTTGGTGCTGATATTGC",
				"TTAACCTACTTGCCTGTCGCTCTATCTTC",
				//Illumina TruSeq adapters, for reads that run off of the end of a fragment
				"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
				"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
				"ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
				"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
				//Illumina Nextera mosaic end
				"AGATGTGTATAAGAGACAG"
			};
		}
	};

	/**
	 * @brief The part of a read to keep
	 */
	struct Kept {
		bool adapterFound_{false};
		uint32_t start_{0};
		uint32_t end_{0};  ///< one past the last base
	};

	explicit AdapterTrimmer(Pars pars) : pars_(std::move(pars)) {
		if (0 == pars_.kmerLen_ || pars_.kmerLen_ > 32) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "adapter k-mer length should be 1 to 32, not " << pars_.kmerLen_ << "\n";
			throw std::runtime_error { ss.str() };
		}
		for (const auto & adapter : pars_.adapters_) {
			const auto upperAdapter = njh::strToUpperRet(adapter);
			for (const auto & seq : {upperAdapter, seqUtil::reverseComplement(upperAdapter, "DNA")}) {
				forEachKmer(seq, [this, &seq](uint64_t kmer, uint32_t pos) {
					kmers_[kmer].emplace_back(Placement{pos, static_cast<uint32_t>(seq.size())});
				});
			}
		}
	}

	[[nodiscard]] Kept find(const std::string & seq) const {
		Kept ret;
		ret.end_ = seq.size();
		if (kmers_.empty()) {
			return ret;
		}
		//where each adapter k-mer puts its adapter on the read
		std::vector<std::pair<int64_t, int64_t>> spans;
		forEachKmer(seq, [this, &spans](uint64_t kmer, uint32_t pos) {
			const auto search = kmers_.find(kmer);
			if (kmers_.end() != search) {
				for (const auto & placement : search->second) {
					const int64_t start = static_cast<int64_t>(pos) - static_cast<int64_t>(placement.offset_);
					spans.emplace_back(start, start + placement.adapterLen_);
				}
			}
		});
		if (spans.empty()) {
			return ret;
		}
		std::sort(spans.begin(), spans.end());
		//overlapping placements are one adapter
		std::vector<std::pair<uint32_t, uint32_t>> adapters;
		int64_t curStart = spans.front().first;
		int64_t curEnd = spans.front().second;
		uint32_t curHits = 0;
		const auto addCurrent = [&]() {
			if (curHits >= pars_.minKmerHits_) {
				adapters.emplace_back(static_cast<uint32_t>(std::max<int64_t>(0, curStart)),
						static_cast<uint32_t>(std::min<int64_t>(static_cast<int64_t>(seq.size()), curEnd)));
			}
		};
		for (const auto & span : spans) {
			if (span.first >= curEnd) {
				addCurrent();
				curStart = span.first;
				curEnd = span.second;
				curHits = 0;
			}
			curEnd = std::max(curEnd, span.second);
			++curHits;
		}
		addCurrent();
		if (adapters.empty()) {
			return ret;
		}
		//the longest piece between the adapters
		ret.adapterFound_ = true;
		ret.end_ = 0;
		uint32_t pieceStart = 0;
		for (size_t idx = 0; idx <= adapters.size(); ++idx) {
			const uint32_t pieceEnd = idx < adapters.size() ? adapters[idx].first : static_cast<uint32_t>(seq.size());
			if (pieceEnd > pieceStart && pieceEnd - pieceStart > ret.end_ - ret.start_) {
				ret.start_ = pieceStart;
				ret.end_ = pieceEnd;
			}
			if (idx < adapters.size()) {
				pieceStart = std::max(pieceStart, adapters[idx].second);
			}
		}
		return ret;
	}

	Pars pars_;

private:
	struct Placement {
		uint32_t offset_;     ///< where the k-mer starts in the adapter
		uint32_t adapterLen_;
	};
	std::unordered_map<uint64_t, std::vector<Placement>> kmers_;

	/**
	 * @brief Calls func(kmer, position) for each k-mer of only A, C, G and T, the k-mer packed 2 bits a base
	 */
	template<typename FUNC>
	void forEachKmer(const std::string & seq, FUNC func) const {
		const uint64_t mask = 32 == pars_.kmerLen_ ? std::numeric_limits<uint64_t>::max() : (uint64_t{1} << (2 * pars_.kmerLen_)) - 1;
		uint64_t kmer = 0;
		uint32_t basesInKmer = 0;
		for (uint32_t pos = 0; pos < seq.size(); ++pos) {
			uint64_t code = 0;
			switch (seq[pos]) {
				case 'A': case 'a': code = 0; break;
				case 'C': case 'c': code = 1; break;
				case 'G': case 'g': code = 2; break;
				case 'T': case 't': code = 3; break;
				default:
					basesInKmer = 0;
					kmer = 0;
					continue;
			}
			kmer = ((kmer << 2) | code) & mask;
			if (++basesInKmer >= pars_.kmerLen_) {
				func(kmer, pos + 1 - pars_.kmerLen_);
			}
		}
	}
};

}  // namespace

int SeekDeepRunner::extractorByKmerMatching(const njh::progutils::CmdArgs &inputCommands) {
  bfs::path idFnp;
  bfs::path lenCutOffsFnp;
  bfs::path uniqueKmersPerTargetFnp;
  bfs::path refSeqDir;
  CoreExtractorPars corePars;
  corePars.smallFragmentCutoff = 100;
  corePars.largeFragmentCutOff = 10000;
  uint32_t numThreads = 1;

  std::string sampleName = "sample";
  bool rename = false;
  bool ignoreMissingKmerSets = false;
  SeekDeepSetUp setUp(inputCommands);
  setUp.description_ = "Extract each target's reads from one sample's long reads (e.g. nanopore or pacbio) by the target k-mers they share, from genTargetInfoFromGenomes, "
      "then check primers, length and quality, reads of several targets are split into one read per target. "
      "Writes passing and failing reads per target plus extractionProfile.tab.txt and extractionStats.tab.txt to the output directory";
  setUp.setOption(ignoreMissingKmerSets, "--ignoreMissingKmerSets", "ignore Missing Kmer Sets");

  //id
  setUp.setOption(idFnp, "--ids,--id", "Primers file", true, "IDs");

  setUp.setOption(lenCutOffsFnp, "--lenCutOffsFnp,--lenCutOffs", "Length cut offs per target, e.g. forSeekDeep/lenCutOffs.txt from genTargetInfoFromGenomes", true, "IDs");
  setUp.setOption(uniqueKmersPerTargetFnp, "--uniqueKmersPerTarget", "The per target k-mer table, e.g. forSeekDeep/uniqueKmers.tab.txt.gz or allKmers.tab.txt.gz from genTargetInfoFromGenomes", true, "IDs");
  bfs::path trimPrimersFnp;
  bfs::path trimKmersPerTargetFnp;
  bfs::path trimLenCutOffsFnp;
  setUp.setOption(trimPrimersFnp, "--trimPrimers",
      "Primers, in the same format as --id, that mark off regions within the amplified targets to trim the reads to. With these, --id, --lenCutOffs and --uniqueKmersPerTarget "
      "describe the amplified targets and are only used to work out which target a read is, split reads of several targets and drop foldbacks; a read doesn't need the amplified primers or the amplified "
      "target's length. Each read is cut to every region of its target it covers, and only the trimming primers, the region's length and the quality are checked. "
      "A target can have several regions, they can overlap, and each region is output as its own target. Regions are matched to the amplified targets by their k-mers", false, "Trim Primers");
  setUp.setOption(trimKmersPerTargetFnp, "--trimKmersPerTarget", "The k-mers per trimming region, from genTargetInfoFromGenomes run on --trimPrimers", !trimPrimersFnp.empty(), "Trim Primers");
  setUp.setOption(trimLenCutOffsFnp, "--trimLenCutOffs", "The length cut offs per trimming region, from genTargetInfoFromGenomes run on --trimPrimers", !trimPrimersFnp.empty(), "Trim Primers");
  const bool trimToRegions = !trimPrimersFnp.empty();
  uint32_t trimRegionMinKmers = 10;
  setUp.setOption(trimRegionMinKmers, "--trimRegionMinKmers", "The k-mer hits needed to find a trimming region in a read", false, "Trim Primers");
  uint32_t trimRegionSlop = 40;
  setUp.setOption(trimRegionSlop, "--trimRegionSlop", "Bases kept past each trimming primer for finding it, for reads whose first or last k-mer is further in because of sequencing errors", false, "Trim Primers");
  setUp.setOption(refSeqDir, "--refSeqDir", "reference Seq Dir, should have a .fasta file named by each target", false, "IDs");

  //in and out
  setUp.processDefaultReader(VecStr{"--fastq", "--fastqgz", "--fasta", "--fastagz", "--bam"});
  setUp.processDirectoryOutputName(true);
  if("sample" == sampleName){
    sampleName = bfs::basename(setUp.pars_.ioOptions_.firstName_);
    if(njh::endsWith(setUp.pars_.ioOptions_.firstName_.string(), ".gz")){
      sampleName = bfs::path(sampleName).replace_extension("").string();
    }
    //for unaligned bams (e.g. from nanopore dorado demux), default to the SM of the read group(s) in the header
    if(SeqIOOptions::inFormats::BAM == setUp.pars_.ioOptions_.inFormat_){
      setSampleNameFromBamHeader(setUp.pars_.ioOptions_.firstName_, sampleName, "--sampleName");
    }
  }
  setUp.setOption(sampleName, "--sampleName", "sample name, if input is a bam this defaults to the SM of the read group in the header otherwise the input file's basename");
  setUp.setOption(rename, "--rename", "rename input extracted reads");

  // filtering
  setUp.setOption(corePars.smallFragmentCutoff, "--minLenCutOff", "Hard cut off min length", false, "filtering");
  setUp.setOption(corePars.largeFragmentCutOff, "--maxLenCutOff", "Hard cut off max length", false, "filtering");

  uint32_t hardKmersPerTarget = 5;
  double kmerFracPerTarget = 0.10;
  setUp.setOption(hardKmersPerTarget, "--hardKmersPerTarget", "hard occurencet cut off Kmers Per Target", false, "filtering");
  setUp.setOption(kmerFracPerTarget, "--kmerFracPerTarget", "at minimum have to have this fraction of kmers in a Target", false, "filtering");

  //running
  setUp.setOption(numThreads, "--numThreads", "number of threads");
  bool noInvertedKmerIndex = false;
  setUp.setOption(noInvertedKmerIndex, "--noInvertedKmerIndex",
      "Probe every target's k-mer set for every read k-mer instead of using one k-mer -> target(s) index. The index gives identical counts, counting a k-mer shared between targets toward each of them exactly as the probe does; this exists to fall back if that is ever in doubt");
  const bool useInvertedKmerIndex = !noInvertedKmerIndex;

  MultihitSplitter::Pars splitPars;
  bool noMultihitSplit = false;
  setUp.setOption(noMultihitSplit, "--noMultihitSplit",
      "Don't split reads that hit several targets; by default each target's stretch of the read (found by k-mers unique to that target) is extracted as its own read, e.g. amplicons ligated end to end", false, "Multihit Split");
  splitPars.split_ = !noMultihitSplit;
  setUp.setOption(splitPars.minKmers_, "--multihitSplitMinKmers", "Unique k-mer hits a stretch of a multihit read needs to count as a target's amplicon", false, "Multihit Split");
  setUp.setOption(splitPars.maxGap_, "--multihitSplitMaxGap", "Largest gap between a target's unique k-mer hits within one stretch of a multihit read", false, "Multihit Split");
  setUp.setOption(splitPars.primerSlop_, "--multihitSplitPrimerSlop", "Sequence kept on each side of a stretch past the longest primer, so the primers are still found", false, "Multihit Split");
  setUp.setOption(splitPars.oneEnd_, "--multihitSplitOneEnd",
      "Multihit reads demultiplexed on one end's barcode (barcodeEnds=1 in the name): none leaves them as multihit, nearest keeps only the amplicon nearest that barcode (needs assignedOn and direction in the name), all keeps every piece. "
      "Pieces of one end reads are tagged barcodeEnds=1, so they only add to a sample-target that also has both ends reads (processClusters --minBothEndsBarcodeReads); some of these pieces come from other samples' molecules, and without that filter they put 268-444 reads on 5 targets in ligation Pool 1's no template controls. "
      "With the filter, all kept the controls clean and raised final reads by 38% (run 2), 67% (Pool 1) and 13% (Pool 2). Reads without barcodeEnds are always split", false, "Multihit Split");
  bool noKmerStretchTrim = false;
  setUp.setOption(noKmerStretchTrim, "--noKmerStretchTrim",
      "Don't trim reads of one target to that target's stretch of the read. By default a read assigned to a single target is cut back to the stretch "
      "covered by k-mers unique to that target (the same stretches the multihit split finds) plus the target's primer and --kmerStretchTrimSlop on each side, "
      "so sequence that belongs to no target, e.g. a human off-target product ligated onto the amplicon, no longer keeps the primers outside the primer search window. "
      "Only reads that wouldn't pass as they are (primers, length and quality) get trimmed, so a read that passes without trimming is never changed. "
      "On 200k reads, ligation Pool 2 passed 32.5% more reads with it (as many as --primerWithinStart 2000), mostly amplicons with human sequence ligated on, and run 3 (PCR barcoding) 0.5% more", false, "K-mer Stretch Trim");
  const bool kmerStretchTrim = !noKmerStretchTrim;
  uint32_t kmerStretchTrimSlop = 40;
  setUp.setOption(kmerStretchTrimSlop, "--kmerStretchTrimSlop",
      "Bases kept past each primer when trimming to a target's k-mer stretch, for reads whose first or last unique k-mer is further in because of sequencing errors. "
      "The primer ends up about this far into the trimmed read, reads that have been cut are always searched at least this far in", false, "K-mer Stretch Trim");
  AdapterTrimmer::Pars adapterPars;
  bool noAdapterTrim = false;
  setUp.setOption(noAdapterTrim, "--noAdapterTrim",
      "Don't cut adapters out of reads. By default every read is searched, in both directions, for the nanopore rapid, ligation and PCR barcoding adapters and the Illumina TruSeq and Nextera adapters "
      "(or the --adapters given), cut at each adapter found and the longest piece kept, before the read is matched to a target. "
      "Nanopore rapid kit reads whose transposase cut inside an amplicon start with the rapid adapter and part of the amplicon, without the cut these end up as false haplotypes that start with adapter", false, "Adapter Trim");
  adapterPars.trim_ = !noAdapterTrim;
  bfs::path adaptersFnp;
  setUp.setOption(adaptersFnp, "--adapters", "A fasta of the adapter sequences to cut out of reads instead of the defaults", false, "Adapter Trim");
  setUp.setOption(adapterPars.kmerLen_, "--adapterKmerLen", "How many bases in a row have to match an adapter (at most 32)", false, "Adapter Trim");
  setUp.setOption(adapterPars.minKmerHits_, "--adapterMinKmerHits", "How many of an adapter's k-mers a read needs to have that adapter", false, "Adapter Trim");
  if(!adaptersFnp.empty()){
    adapterPars.adapters_.clear();
    SeqInput adapterReader(SeqIOOptions::genFastaIn(adaptersFnp));
    for(const auto & adapter : adapterReader.readAllReads<seqInfo>()){
      adapterPars.adapters_.emplace_back(adapter.seq_);
    }
  }
  if(!njh::in(splitPars.oneEnd_, VecStr{"none", "nearest", "all"})){
    std::stringstream ss;
    ss << __PRETTY_FUNCTION__ << ", error " << "--multihitSplitOneEnd should be none, nearest or all, not " << splitPars.oneEnd_ << "\n";
    throw std::runtime_error { ss.str() };
  }


  //primers
  corePars.pDetPars.useMotif_ = false;
  corePars.pDetPars.allowable_.lqMismatches_ = 5;
  corePars.pDetPars.allowable_.hqMismatches_ = 5;
  corePars.pDetPars.allowable_.distances_.query_.coverage_ = .50;
  corePars.pDetPars.allowable_.largeBaseIndel_ = .99;
  corePars.pDetPars.allowable_.oneBaseIndel_ = 4;
  corePars.pDetPars.allowable_.twoBaseIndel_ = 2;
  corePars.pDetPars.primerWithin_ = 50;
  bool primerToUpperCase = false;
  setUp.setOption(primerToUpperCase, "--primerUpper",
                  "Leave primers in upper case", false, "Primer");
  corePars.pDetPars.primerToLowerCase_ = !primerToUpperCase;
  setUp.setOption(corePars.pDetPars.allowable_.distances_.query_.coverage_, "--primerCoverage",
                  "Amount of primers found", false, "Primer");
  setUp.setOption(corePars.pDetPars.allowable_.hqMismatches_, "--primerNumOfMismatches",
                  "Number of Mismatches to allow in primers", false, "Primer");
  setUp.setOption(corePars.pDetPars.allowable_.oneBaseIndel_, "--primerOneBaseIndels",
                  "Number Of One base indels to allow in primers", false, "Primer");
  setUp.setOption(corePars.pDetPars.allowable_.twoBaseIndel_, "--primerTwoBaseIndels",
                  "Number Of Two base indels to allow in primers", false, "Primer");
  setUp.setOption(corePars.pDetPars.allowable_.largeBaseIndel_, "--primerMoreThan2BaseIndels",
                  "Number Of greater than two base indels to allow in primers", false, "Primer");
  setUp.setOption(corePars.pDetPars.primerWithin_, "--primerWithinStart",
                  "How far into the read the primers are searched for. Reads that were cut to their target (split reads of several targets, --noKmerStretchTrim, --trimPrimers) are searched at least 50 bases in, or as far in as the cut left the primers, "
                  "so this doesn't need to be raised for sequence in front of the primers; kept low it keeps down false positives, and joined amplicons from passing whole",
                  false, "Primers");
  setUp.setOption(corePars.pDetPars.primerStart_, "--primerSearchStart",
                  "By default the primers are searched at the very beginning of seq, use this flag to start the search here",
                  false, "Primers");
  //copy over the regular determined primers to the back end primers
  corePars.backEndpDetPars.primerWithin_ = corePars.pDetPars.primerWithin_;
  corePars.backEndpDetPars.primerToLowerCase_ = corePars.pDetPars.primerToLowerCase_ ;

  corePars.backEndpDetPars.allowable_.hqMismatches_ = corePars.pDetPars.allowable_.hqMismatches_;
  corePars.backEndpDetPars.allowable_.lqMismatches_ = corePars.pDetPars.allowable_.lqMismatches_;
  corePars.backEndpDetPars.allowable_.distances_.query_.coverage_ = corePars.pDetPars.allowable_.distances_.query_.coverage_;
  corePars.backEndpDetPars.allowable_.largeBaseIndel_ = corePars.pDetPars.allowable_.largeBaseIndel_;
  corePars.backEndpDetPars.allowable_.oneBaseIndel_ = corePars.pDetPars.allowable_.oneBaseIndel_;
  corePars.backEndpDetPars.allowable_.twoBaseIndel_ = corePars.pDetPars.allowable_.twoBaseIndel_;
  corePars.backEndpDetPars.useMotif_ = corePars.pDetPars.useMotif_;

  corePars.qPars_.qualCheck_ = 30;
  corePars.qPars_.qualCheckCutOff_ = 0.15;
  setUp.setOption(corePars.qPars_.qualCheck_, "--qualCheckLevel",
                  "Bin qualities at this quality to do filtering on fraction above this", false, "Post Processing");
  setUp.setOption(corePars.qPars_.qualCheckCutOff_, "--qualCheckCutOff",
                  "The fractions of bases that have to be above the qualCheckLevel to be kept", false, "Post Processing");
  corePars.primIdsPars.mPars_.allowableErrors_ = 1;
  setUp.setOption(corePars.primIdsPars.mPars_.allowableErrors_, "--MIDAllowableErrors", "Number of errors to allow in MIDs, more errors allowed can lead to mis-sorting of reads by MIDs", false, "MID");
  setUp.finishSetUp(std::cout);
  // run log
  setUp.startARunLog(setUp.pars_.directoryName_);
  // parameter file
  setUp.writeParametersFile(setUp.pars_.directoryName_ + "parametersUsed.txt", false, false);

  // create Primers and MIDs
  PrimersAndMids ids(idFnp);

  if(ids.getTargets().empty()){
    std::stringstream ss;
    ss << __PRETTY_FUNCTION__ << ", error " << "no targets read in from " << idFnp << "\n";
    throw std::runtime_error { ss.str() };
  }

  //init determinators
  ids.initPrimerDeterminator();

  if(ids.containsMids()) {
    ids.initMidDeterminator(corePars.primIdsPars.mPars_);
  }
  //read in extra info
  ids.addLenCutOffs(lenCutOffsFnp);
  if(!refSeqDir.empty()){
    ids.addRefSeqs(refSeqDir);
  }
  uint32_t extractionKmer = ids.addKmerCounts(uniqueKmersPerTargetFnp, ignoreMissingKmerSets);
  //the table can hold all kmers per target (overlapping targets share kmers) or only unique kmers. All kmers pick the
  //winning target and kmers unique to a target count how many targets a read hits. A shared kmer counts toward every
  //target holding it.
  const njhseq::KmerCollectionReadMatcher kmerMatcher(ids.allKmersPerTarget_, extractionKmer, useInvertedKmerIndex);
  const auto & kmerTargetNames = kmerMatcher.names();

  const SimpleKmerHash splitHasher;
  //the trimming regions, output is by region instead of by amplified target
  std::unique_ptr<PrimersAndMids> trimIds;
  std::list<KmerSetInvertedIndex::SetMap> trimRegionKmers; //the indexes are kept by the address of these
  std::unordered_map<std::string, std::vector<TrimRegion>> trimRegionsPerTarget;
  if(trimToRegions){
    if(ids.containsMids()){
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << ", error " << "--trimPrimers can't be used with MIDs in --id, the MIDs are next to the amplified primers which aren't looked for when trimming to regions" << "\n";
      throw std::runtime_error { ss.str() };
    }
    trimIds = std::make_unique<PrimersAndMids>(trimPrimersFnp);
    if(trimIds->getTargets().empty()){
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << ", error " << "no regions read in from " << trimPrimersFnp << "\n";
      throw std::runtime_error { ss.str() };
    }
    if(trimIds->containsMids()){
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << ", error " << "--trimPrimers shouldn't have MIDs, " << trimPrimersFnp << "\n";
      throw std::runtime_error { ss.str() };
    }
    trimIds->initPrimerDeterminator();
    trimIds->addLenCutOffs(trimLenCutOffsFnp);
    const uint32_t trimKmer = trimIds->addKmerCounts(trimKmersPerTargetFnp, ignoreMissingKmerSets);
    if(trimKmer != extractionKmer){
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << ", error " << "the k-mer length of --trimKmersPerTarget, " << trimKmer << ", should be the same as --uniqueKmersPerTarget, " << extractionKmer << "\n";
      throw std::runtime_error { ss.str() };
    }
    //a region is within the amplified target that has most of its k-mers, or in more than one when the amplified targets overlap, not all of a
    //region's k-mers have to be found since a genome's region can be there without all of its amplified target (e.g. an amplified primer variant)
    OutputStream outRegions(njh::files::make_path(setUp.pars_.directoryName_, "trimRegions.tab.txt"));
    outRegions << "region\ttarget\tdirection\tregionKmers\tregionKmersInTarget\tregionKmersInTargetFrac" << std::endl;
    VecStr regionsWithoutTarget;
    auto regionNames = trimIds->getTargets();
    njh::sort(regionNames);
    for(const auto & regionName : regionNames){
      if(!njh::in(regionName, trimIds->allKmersPerTarget_)){
        continue;
      }
      const auto & regionKmers = trimIds->allKmersPerTarget_.at(regionName);
      //a region's primers can be in the other direction to the amplified primers
      std::unordered_set<uint64_t> regionKmersRevComp;
      for(const auto kmer : regionKmers){
        regionKmersRevComp.emplace(splitHasher.hash(splitHasher.revCompReverseHash(kmer)));
      }
      const auto countShared = [](const std::unordered_set<uint64_t> & kmers, const std::unordered_set<uint64_t> & targetKmers){
        uint32_t count = 0;
        for(const auto kmer : kmers){
          if(njh::in(kmer, targetKmers)){
            ++count;
          }
        }
        return count;
      };
      std::map<std::string, std::pair<uint32_t, bool>> shared; //count and if it's in the reverse direction
      uint32_t mostShared = 0;
      for(const auto & targetKmers : ids.allKmersPerTarget_){
        const uint32_t forwardCount = countShared(regionKmers, targetKmers.second);
        const uint32_t revCompCount = countShared(regionKmersRevComp, targetKmers.second);
        shared[targetKmers.first] = {std::max(forwardCount, revCompCount), revCompCount > forwardCount};
        mostShared = std::max({mostShared, forwardCount, revCompCount});
      }
      const double total = static_cast<double>(regionKmers.size());
      TrimRegion region;
      region.name_ = regionName;
      if(nullptr != trimIds->targets_.at(regionName).lenCuts_){
        region.maxSpan_ = trimIds->targets_.at(regionName).lenCuts_->maxLenChecker_.maxLen_;
      }
      for(const auto & fwd : trimIds->pDeterminator_->primers_.at(regionName).fwds_){
        region.forwardPrimerLen_ = std::max<uint32_t>(region.forwardPrimerLen_, len(fwd.info_));
      }
      for(const auto & rev : trimIds->pDeterminator_->primers_.at(regionName).revs_){
        region.reversePrimerLen_ = std::max<uint32_t>(region.reversePrimerLen_, len(rev.info_));
      }
      trimRegionKmers.emplace_back();
      trimRegionKmers.back()[regionName] = regionKmers;
      region.index_ = KmerSetInvertedIndex::get(trimRegionKmers.back());
      bool assigned = false;
      for(const auto & target : shared){
        const uint32_t count = target.second.first;
        const double frac = total == 0 ? 0 : count / total;
        if(count > 0 && ((count == mostShared && frac >= 0.5) || frac >= 0.9)){
          assigned = true;
          region.revComp_ = target.second.second;
          trimRegionsPerTarget[target.first].emplace_back(region);
          outRegions << regionName << "\t" << target.first << "\t" << (region.revComp_ ? "reverse" : "forward") << "\t" << regionKmers.size() << "\t" << count << "\t" << frac << std::endl;
        }
      }
      if(!assigned){
        regionsWithoutTarget.emplace_back(regionName);
      }
    }
    if(!regionsWithoutTarget.empty()){
      std::stringstream ss;
      ss << __PRETTY_FUNCTION__ << ", error " << "these --trimPrimers regions aren't within any of the amplified targets of --id (none of them has at least half of the region's k-mers): " << njh::conToStr(regionsWithoutTarget, ", ") << "\n";
      throw std::runtime_error { ss.str() };
    }
  }
  //what the output is for, the trimming regions or the amplified targets
  PrimersAndMids & outIds = trimToRegions ? *trimIds : ids;
  if(useInvertedKmerIndex){
    auto invIndex = njhseq::KmerSetInvertedIndex::get(ids.allKmersPerTarget_);
    std::cout << "kmerSetExclusivity\tsets\t" << invIndex->nSets_
        << "\ttotalKmers\t" << invIndex->totalKmers_
        << "\tsharedKmers\t" << invIndex->nSharedKmers_
        << "\texclusive\t" << njh::boolToStr(0 == invIndex->nSharedKmers_) << std::endl;
  }
  //for splitting multihit reads, positions of k-mers unique to one target are looked up in the inverted index
  std::shared_ptr<const KmerSetInvertedIndex> splitIndex;
  std::vector<uint32_t> matcherToSplitIdx(kmerTargetNames.size(), 0);
  std::vector<uint32_t> splitMaxSpan;
  if(splitPars.split_ || kmerStretchTrim || trimToRegions){
    splitIndex = KmerSetInvertedIndex::get(ids.allKmersPerTarget_);
    std::unordered_map<std::string, uint32_t> splitIdxByName;
    for(uint32_t idx = 0; idx < splitIndex->setNames_.size(); ++idx){
      splitIdxByName[splitIndex->setNames_[idx]] = idx;
    }
    for(uint32_t idx = 0; idx < kmerTargetNames.size(); ++idx){
      matcherToSplitIdx[idx] = splitIdxByName.at(kmerTargetNames[idx]);
    }
    //a stretch can't be longer than its target's longest allowed amplicon, which keeps two copies of a target in a row apart
    splitMaxSpan.assign(splitIndex->setNames_.size(), corePars.largeFragmentCutOff);
    for(uint32_t idx = 0; idx < splitIndex->setNames_.size(); ++idx){
      const auto & name = splitIndex->setNames_[idx];
      if(njh::in(name, ids.targets_) && nullptr != ids.targets_.at(name).lenCuts_){
        splitMaxSpan[idx] = ids.targets_.at(name).lenCuts_->maxLenChecker_.maxLen_;
      }
    }
  }


  //for trimming to a target's k-mer stretch, the longest forward and reverse primer of each target, the primers sit just outside the stretch
  std::unordered_map<std::string, std::pair<uint32_t, uint32_t>> primerLensPerTarget;
  for(const auto & primer : ids.pDeterminator_->primers_){
    auto & lens = primerLensPerTarget[primer.first];
    for(const auto & fwd : primer.second.fwds_){
      lens.first = std::max<uint32_t>(lens.first, len(fwd.info_));
    }
    for(const auto & rev : primer.second.revs_){
      lens.second = std::max<uint32_t>(lens.second, len(rev.info_));
    }
  }
  //the unique k-mer hits needed to find a target's stretch for the trim, by index order. Targets too short to have --multihitSplitMinKmers
  //k-mer positions between their primers need half their positions instead (at least 5), the insert length is estimated as the middle of the
  //target's length cutoffs minus its primers
  std::vector<MultihitSplitter::Pars> trimSegmentPars;
  if(nullptr != splitIndex){
    trimSegmentPars.assign(splitIndex->setNames_.size(), splitPars);
    for(uint32_t idx = 0; idx < splitIndex->setNames_.size(); ++idx){
      const auto & name = splitIndex->setNames_[idx];
      if(!njh::in(name, ids.targets_) || nullptr == ids.targets_.at(name).lenCuts_ || !njh::in(name, primerLensPerTarget)){
        continue;
      }
      const auto & lenCuts = *ids.targets_.at(name).lenCuts_;
      const auto & primerLens = primerLensPerTarget.at(name);
      const uint64_t ampliconLen = (static_cast<uint64_t>(lenCuts.minLenChecker_.minLen_) + lenCuts.maxLenChecker_.maxLen_) / 2;
      const uint64_t primersLen = primerLens.first + primerLens.second;
      const uint64_t positions = ampliconLen > primersLen + extractionKmer ? ampliconLen - primersLen - extractionKmer + 1 : 1;
      if(positions < splitPars.minKmers_){
        trimSegmentPars[idx].minKmers_ = std::max<uint32_t>(5, (positions + 1) / 2);
      }
    }
  }
  //how far in the primers are searched for in reads that were cut to their target. The cut leaves the primer about the slop into the read (plus the
  //difference between the longest primer, which split reads are cut with, and the shortest), a window smaller than that loses most of the split
  //reads, e.g. --primerWithinStart 24 passed 48 pieces per 100 reads of PRJNA1188139 instead of 146. Never less than the default window
  const uint32_t defaultPrimerWithin = 50;
  uint64_t minPrimerSize = std::numeric_limits<uint32_t>::max();
  for(const auto & lens : primerLensPerTarget){
    minPrimerSize = std::min<uint64_t>({minPrimerSize, lens.second.first, lens.second.second});
  }
  const uint32_t cutReadPrimerWithin = std::max<uint32_t>({corePars.pDetPars.primerWithin_, defaultPrimerWithin, kmerStretchTrimSlop,
    static_cast<uint32_t>(ids.pDeterminator_->getMaxPrimerSize() - minPrimerSize) + splitPars.primerSlop_,
    trimToRegions ? trimRegionSlop : 0});
  auto cutReadForwardPars = corePars.pDetPars;
  cutReadForwardPars.primerWithin_ = cutReadPrimerWithin;
  auto cutReadBackPars = corePars.backEndpDetPars;
  cutReadBackPars.primerWithin_ = cutReadPrimerWithin;

  // set up input
  seqInfo seq;
  SeqInput reader(setUp.pars_.ioOptions_);
  reader.openIn();



  //qual check off;
  std::unique_ptr<ReadChecker> qualChecker = std::make_unique<ReadCheckerQualCheck>(corePars.qPars_.qualCheck_,
                                                                                    corePars.qPars_.qualCheckCutOff_, true);

  //set up extraction counts outputs
  ExtractionStator masterCounts;
  //counts per target and MID, keyed by the name of the extracted file, <target><MID>, or <target>-MIDunknown/-MIDUnrecognized when the MID wasn't determined
  ExtractionStator perMidCounts;
  OutOptions outCountsOpts(njh::files::make_path(setUp.pars_.directoryName_, "extractionProfile.tab.txt"));
  OutputStream outCounts(outCountsOpts);

  OutOptions outStatsOpts(njh::files::make_path(setUp.pars_.directoryName_, "extractionStats.tab.txt"));
  OutputStream outStats(outStatsOpts);

  std::unordered_map<std::string, uint32_t> readsPerSet;
  std::unordered_map<std::string, uint32_t> readsPerSetRevComp;
  //pieces cut from multihit reads, per target; the reads they came from are counted in readsPerSet under "multihitSplit"
  std::unordered_map<std::string, uint32_t> piecesPerSet;
  std::unordered_map<std::string, uint32_t> piecesPerSetRevComp;
  uint64_t multihitSplitFoldback = 0;
  //reads of one target trimmed to the target's k-mer stretch, per target, and the bases cut off
  std::unordered_map<std::string, uint64_t> kmerTrimmedPerSet;
  uint64_t kmerTrimmedFrontBases = 0;
  uint64_t kmerTrimmedBackBases = 0;
  //reads with an adapter cut out, and the bases cut off
  const AdapterTrimmer adapterTrimmer(adapterPars);
  uint64_t adapterTrimmed = 0;
  uint64_t adapterTrimmedBases = 0;
  //when trimming to regions, reads of the region's target that don't cover the region, per region, and reads that passed for at least one
  //region, per amplified target
  std::unordered_map<std::string, uint64_t> notCoveredPerRegion;
  std::unordered_map<std::string, uint64_t> readsPassingAnyRegion;

  VecStr names = ids.getTargets();
  MultiSeqIO seqOut;
  auto initialExtractionDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"extractedReads"});
  auto failedExtractionDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"failedExtractedReads"});
  for(const auto & name : names){
    auto seqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(initialExtractionDir, name) );
    seqOut.addReader(name, seqOutOpts);
  }

  for(const auto & name : outIds.getTargets()){
    if(ids.containsMids()) {
      for(const auto & midname : getVectorOfMapKeys(ids.mids_)) {
        auto seqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(name, midname) ) );
        seqOut.addReader(njh::pasteAsStr(name, midname, "-passed"), seqOutOpts);

        auto failedSeqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(failedExtractionDir, name + midname) );
        seqOut.addReader(njh::pasteAsStr(name, midname, "-failed"), failedSeqOutOpts);
      }
      seqOut.addReader(njh::pasteAsStr(name, "-MIDUnrecognized", "-failed"),
                       SeqIOOptions::genFastqOutGz(
                         njh::files::make_path(failedExtractionDir, name+ "-MIDUnrecognized")));
      seqOut.addReader(njh::pasteAsStr(name, "-MIDunknown", "-failed"),
                 SeqIOOptions::genFastqOutGz(
                   njh::files::make_path(failedExtractionDir, name + "-MIDunknown")));
    } else {
      auto seqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, njh::pasteAsStr(name, sampleName) ) );
      seqOut.addReader(njh::pasteAsStr(name, "-passed"), seqOutOpts);

      auto failedSeqOutOpts = SeqIOOptions::genFastqOutGz(njh::files::make_path(failedExtractionDir, name) );
      seqOut.addReader(njh::pasteAsStr(name, "-failed"), failedSeqOutOpts);

    }
  }


  seqOut.addReader("undetermined", SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, "undetermined")));
  seqOut.addReader("multihit", SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, "multihit")));
  if(splitPars.split_){
    seqOut.addReader("multihitSplit", SeqIOOptions::genFastqOutGz(njh::files::make_path(setUp.pars_.directoryName_, "multihitSplit")));
  }

  //set up aligner
  // creating aligner
  // create aligner for primer identification
  auto scoreMatrix = substituteMatrix::createDegenScoreMatrixNoNInRef(
      setUp.pars_.generalMatch_, setUp.pars_.generalMismatch_);

  setUp.pars_.gapInfo_.gapLeftQueryOpen_ = 0;
  setUp.pars_.gapInfo_.gapLeftQueryExtend_ = 0;

  setUp.pars_.gapInfo_.gapRightQueryOpen_ = 0;
  setUp.pars_.gapInfo_.gapRightQueryExtend_ = 0;

  setUp.pars_.gapInfo_.gapLeftRefOpen_ = 0;
  setUp.pars_.gapInfo_.gapLeftRefExtend_ = 0;

  setUp.pars_.gapInfo_.gapRightRefOpen_ = 0;
  setUp.pars_.gapInfo_.gapRightRefExtend_ = 0;

  gapScoringParameters gapPars(setUp.pars_.gapInfo_);
  KmerMaps emptyMaps;
  bool countEndGaps = false;
  uint64_t maxPrimerSize = ids.pDeterminator_->getMaxPrimerSize();
  uint64_t maxReadSize =  std::max<uint64_t>(maxPrimerSize, outIds.pDeterminator_->getMaxPrimerSize()) * 4 + cutReadPrimerWithin;
  aligner alignObj(maxReadSize, gapPars, scoreMatrix, emptyMaps, setUp.pars_.qScorePars_, countEndGaps, false);
  alignObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);

  concurrent::AlignerPool alnPool(alignObj, numThreads);
  alnPool.initAligners();



  std::mutex mut;
  std::function<void()> readInComp = [&reader, &ids, &readsPerSet,&readsPerSetRevComp,&mut,
                                      &seqOut,&sampleName,&rename,
                                      &corePars,&alnPool,
                                      &qualChecker,
                                      &masterCounts,
                                      &perMidCounts,
                                      &hardKmersPerTarget,
                                      &kmerFracPerTarget,
                                      &kmerMatcher,
                                      &kmerTargetNames,
                                      &splitPars,
                                      &splitIndex,
                                      &splitHasher,
                                      &splitMaxSpan,
                                      &trimSegmentPars,
                                      &matcherToSplitIdx,
                                      &extractionKmer,
                                      &maxPrimerSize,
                                      &piecesPerSet,
                                      &piecesPerSetRevComp,
                                      &multihitSplitFoldback,
                                      &kmerStretchTrim,
                                      &kmerStretchTrimSlop,
                                      &primerLensPerTarget,
                                      &kmerTrimmedPerSet,
                                      &kmerTrimmedFrontBases,
                                      &kmerTrimmedBackBases,
                                      &adapterTrimmer,
                                      &adapterTrimmed,
                                      &adapterTrimmedBases,
                                      &outIds,
                                      &cutReadForwardPars,
                                      &cutReadBackPars,
                                      &trimToRegions,
                                      &trimRegionsPerTarget,
                                      &trimRegionMinKmers,
                                      &trimRegionSlop,
                                      &notCoveredPerRegion,
                                      &readsPassingAnyRegion
                                      ]() {
    auto maxMidSize = ids.getMaxMIDSize() + 2;
    seqInfo seq;
    std::unordered_map<std::string, uint32_t> readsPerSetCurrent;
    std::unordered_map<std::string, uint32_t> readsPerSetRevCompCurrent;
    std::unordered_map<std::string, uint32_t> piecesPerSetCurrent;
    std::unordered_map<std::string, uint32_t> piecesPerSetRevCompCurrent;
    uint64_t multihitSplitFoldbackCurrent = 0;
    std::unordered_map<std::string, uint64_t> kmerTrimmedPerSetCurrent;
    uint64_t kmerTrimmedFrontBasesCurrent = 0;
    uint64_t kmerTrimmedBackBasesCurrent = 0;
    uint64_t adapterTrimmedCurrent = 0;
    uint64_t adapterTrimmedBasesCurrent = 0;
    std::unordered_map<std::string, uint64_t> notCoveredPerRegionCurrent;
    std::unordered_map<std::string, uint64_t> readsPassingAnyRegionCurrent;
    auto trimRegionPars = splitPars;
    trimRegionPars.minKmers_ = trimRegionMinKmers;

    ExtractionStator masterCountsCurrent;
    ExtractionStator perMidCountsCurrent;

    auto kmerMatch = kmerMatcher.genMatch();

    auto currentAlignerObj = alnPool.popAligner();
    //primer, MID, length and quality checks for a read assigned to winnerSet, writing it to that target's passed or failed output
    auto extractTarget = [&](seqInfo & seq, const std::string & winnerSet, ExtractionStator::extractCase extractorCase, bool cutToTarget) -> bool {
        const auto & forwardPars = cutToTarget ? cutReadForwardPars : corePars.pDetPars;
        const auto & backPars = cutToTarget ? cutReadBackPars : corePars.backEndpDetPars;
        //only check min len if INVERSECHIMERA wasn't true
        if(seq.on_) {
          //min len
          outIds.targets_.at(winnerSet).lenCuts_->minLenChecker_.checkRead(seq);
          if(!seq.on_){
            extractorCase = ExtractionStator::extractCase::MINLENBAD;
          }
        }

        std::string MIDunknownName = ids.containsMids() ? "-MIDunknown" : "";
        std::string perMidName = njh::pasteAsStr(winnerSet, MIDunknownName);
        if(seq.on_){
          std::string frontPrimerName = "unrecognized";
          std::string backPrimerName = "unrecognized";
          std::string midName;
          bool passesMID = false;
          if(ids.containsMids()) {
            //if also contains mids check for those, the previous extraction by kmer step would have re-oriented everything into forward primer direction
            auto positionForwardResults = outIds.pDeterminator_->determineBestForwardPrimerPosFront(seq, forwardPars, *currentAlignerObj, VecStr{winnerSet});
            seq.reverseComplementRead(false, true);
            auto positionReverseResults = outIds.pDeterminator_->determineBestReversePrimerPosFront(seq, backPars, *currentAlignerObj, VecStr{winnerSet});
            seq.reverseComplementRead(false, true);
            frontPrimerName = positionForwardResults.primerName_;
            backPrimerName = positionReverseResults.primerName_;
            //check for mid only if the read passed for both primers
            if(winnerSet == positionForwardResults.primerName_ && winnerSet == positionReverseResults.primerName_) {
              MidDeterminator::MidDeterminePars midpars;
              midpars.searchStart_ = positionForwardResults.start_ > maxMidSize ? positionForwardResults.start_ - maxMidSize : 0;
              midpars.searchStop_ = midpars.searchStart_ + maxMidSize;
              midpars.allowableErrors_ = corePars.primIdsPars.mPars_.allowableErrors_;
              auto midResults = ids.mDeterminator_->searchRead(seq, midpars, midpars);
              auto midResultsProcessed = ids.mDeterminator_->processSearchRead(seq, midResults);
              if(midResultsProcessed.case_ == MidDeterminator::ProcessedRes::PROCESSED_CASE::MATCH) {
                //process the primers locations
                //this will both trim and determine the primers, redundant to do this twice for the forward primer with the trimming of the MID it's easier to just recheck
                //the modification below only makes sense if there's only 1 MID and it's attached to the forward primer which is the most common scernario
                auto modifiedSearchPars = forwardPars;
                modifiedSearchPars.primerStart_ = 0;
                modifiedSearchPars.primerWithin_ = outIds.pDeterminator_->getMaxPrimerSize();
                frontPrimerName = outIds.pDeterminator_->determineForwardPrimer(seq, modifiedSearchPars, *currentAlignerObj, VecStr{winnerSet});
                seq.reverseComplementRead(false, true);
                backPrimerName = outIds.pDeterminator_->determineWithReversePrimer(seq, backPars, *currentAlignerObj, VecStr{winnerSet});
                seq.reverseComplementRead(false, true);
                midName = midResultsProcessed.midName_;
                passesMID = true;
                perMidName = njh::pasteAsStr(winnerSet, midName);
              } else {
                perMidName = njh::pasteAsStr(winnerSet, "-MIDUnrecognized");
              }
            }
          } else {
            //this will both trim and determine the primers
            frontPrimerName = outIds.pDeterminator_->determineForwardPrimer(seq, forwardPars, *currentAlignerObj, VecStr{winnerSet});
            seq.reverseComplementRead(false, true);
            backPrimerName = outIds.pDeterminator_->determineWithReversePrimer(seq, backPars, *currentAlignerObj, VecStr{winnerSet});
            seq.reverseComplementRead(false, true);
            passesMID = true;
          }
          //passes if the front and end primer both match the expected primer based on the winder set from the kmer matching
          bool passesFront = frontPrimerName == winnerSet;
          bool passesBack = backPrimerName == winnerSet;
          if (!passesBack && !passesFront) {
            extractorCase = ExtractionStator::extractCase::FAILEDBOTHPRIMERS;
          } else if(!passesFront){
            extractorCase = ExtractionStator::extractCase::BADFORWARD;
          } else if(!passesBack){
            extractorCase = ExtractionStator::extractCase::BADREVERSE;
          }
          if (passesBack && passesFront) {
            //check for MID, if no mids the pass MIDs is always true
            if(!passesMID) {
              //failed first min length pass
              seqOut.openWrite(njh::pasteAsStr(winnerSet, "-MIDUnrecognized", "-failed"), seq);
              extractorCase = ExtractionStator::extractCase::BADMID;
            } else {
              //min len
              outIds.targets_.at(winnerSet).lenCuts_->minLenChecker_.checkRead(seq);
              if(!seq.on_){
                extractorCase = ExtractionStator::extractCase::MINLENBAD;
              }

              //max len
              if(seq.on_){
                outIds.targets_.at(winnerSet).lenCuts_->maxLenChecker_.checkRead(seq);
                if(!seq.on_){
                  extractorCase = ExtractionStator::extractCase::MAXLENBAD;
                }
              }

              //quality
              if(seq.on_){
                qualChecker->checkRead(seq);
                if(!seq.on_){
                  extractorCase = ExtractionStator::extractCase::QUALITYFAILED;
                }
              }
              if (seq.on_) {
                //passed all filters
                seqOut.openWrite(njh::pasteAsStr(njh::pasteAsStr(winnerSet, midName), "-passed"), seq);
              } else {
                //failed
                seqOut.openWrite(njh::pasteAsStr(njh::pasteAsStr(winnerSet, midName), "-failed"), seq);
              }
            }
          } else {
            seq.name_.append(njh::pasteAsStr("[", "frontPrimerName=", frontPrimerName, ";", "backPrimerName=", backPrimerName, ";","]"));
            //failed primer check
            seqOut.openWrite(njh::pasteAsStr(njh::pasteAsStr(winnerSet, MIDunknownName), "-failed"), seq);
          }
        } else {
          //failed first min length pass
          seqOut.openWrite(njh::pasteAsStr(njh::pasteAsStr(winnerSet, MIDunknownName), "-failed"), seq);
        }
        masterCountsCurrent.increaseCounts(winnerSet, seq.name_, extractorCase);
        if(ids.containsMids()){
          perMidCountsCurrent.increaseCounts(perMidName, seq.name_, extractorCase);
        }
        return ExtractionStator::extractCase::GOOD == extractorCase;
    };
    //cut a read, already in the target's forward direction, to each of the target's trimming regions it covers
    auto extractTrimRegions = [&](const seqInfo & seq, const std::string & target, ExtractionStator::extractCase extractorCase){
      const auto regions = trimRegionsPerTarget.find(target);
      if(regions == trimRegionsPerTarget.end()){
        return;
      }
      bool passedAny = false;
      seqInfo seqRevComp; //only made if a region is in the other direction
      for(const auto & region : regions->second){
        if(ExtractionStator::extractCase::GOOD != extractorCase){
          masterCountsCurrent.increaseCounts(region.name_, seq.name_, extractorCase);
          continue;
        }
        if(region.revComp_ && seqRevComp.seq_.empty()){
          seqRevComp = seq;
          seqRevComp.reverseComplementRead(false, true);
        }
        const auto & regionSeq = region.revComp_ ? seqRevComp : seq;
        const auto found = region.find(regionSeq.seq_, splitHasher, extractionKmer, trimRegionPars);
        if(0 == found.second){
          ++notCoveredPerRegionCurrent[region.name_];
          continue;
        }
        const uint32_t frontPad = region.forwardPrimerLen_ + trimRegionSlop;
        const uint32_t backPad = region.reversePrimerLen_ + trimRegionSlop;
        const uint32_t start = found.first > frontPad ? found.first - frontPad : 0;
        const uint32_t end = std::min<uint64_t>(len(regionSeq), static_cast<uint64_t>(found.second) + backPad);
        auto piece = regionSeq.getSubRead(start, end - start);
        if(extractTarget(piece, region.name_, ExtractionStator::extractCase::GOOD, true)){
          passedAny = true;
        }
      }
      if(passedAny){
        ++readsPassingAnyRegionCurrent[target];
      }
    };

    while(reader.readNextReadLock(seq)){
      ++masterCountsCurrent.totalReadCount_;
      if(len(seq) < corePars.smallFragmentCutoff){
        ++masterCountsCurrent.smallFrags_;
        continue;
      }
      if (len(seq) > corePars.largeFragmentCutOff) {
        ++masterCountsCurrent.largeFrags_;
        continue;
      }
      //the read length cut offs are on the read as sequenced, a short amplicon whose read is under --minLenCutOff once its adapters are cut is
      //still kept, the target's own length cut offs check what's left
      if(adapterTrimmer.pars_.trim_){
        const auto kept = adapterTrimmer.find(seq.seq_);
        if(kept.adapterFound_){
          ++adapterTrimmedCurrent;
          adapterTrimmedBasesCurrent += len(seq) - (kept.end_ - kept.start_);
          seq = seq.getSubRead(kept.start_, kept.end_ - kept.start_);
        }
      }
      //all kmers per target pick the winner, kmers unique to a target count the targets hit
      kmerMatch.reset();
      kmerMatcher.addKmers(seq.seq_, true, kmerMatch);
      kmerMatcher.countFound(kmerMatch);

      //filter
      const double nInputKmers = static_cast<double>(kmerMatch.kmers.size());
      for(auto * found : {&kmerMatch.found, &kmerMatch.foundRevComp, &kmerMatch.uniqueFound, &kmerMatch.uniqueFoundRevComp}){
        for (auto & count : *found) {
          if(count < hardKmersPerTarget || static_cast<double>(count)/nInputKmers <= kmerFracPerTarget) {
            count = 0;
          }
        }
      }
      kmerMatcher.determineWinner(kmerMatch);
      std::string winnerSet = kmerMatch.hasWinner() ? kmerTargetNames[kmerMatch.winnerIdx] : "undetermined";
      const bool winnerRevComp = kmerMatch.winnerRevComp;
      // doesn't take into account inverse hits, 1 target with large hits in forward and reverse
      uint32_t target_hit_count = njhseq::KmerCollectionReadMatcher::countCollectionsHit(kmerMatch, hardKmersPerTarget);
      bool multi_hit = target_hit_count > 1;
      if (multi_hit) {
        winnerSet = "multihit";
      }
      if(multi_hit && splitPars.split_){
        //stretches of the targets that made this a multihit
        std::vector<char> allowed(splitIndex->setNames_.size(), 0);
        for(uint32_t idx = 0; idx < kmerTargetNames.size(); ++idx){
          if(kmerMatch.uniqueFound[idx] >= hardKmersPerTarget || kmerMatch.uniqueFoundRevComp[idx] >= hardKmersPerTarget){
            allowed[matcherToSplitIdx[idx]] = 1;
          }
        }
        const auto segments = MultihitSplitter::findSegments(seq.seq_, *splitIndex, splitHasher, extractionKmer, allowed, splitMaxSpan, splitPars);
        std::set<uint32_t> segmentTargets;
        for(const auto & seg : segments){
          segmentTargets.emplace(seg.setIdx_);
        }
        MetaDataInName inMeta;
        const bool hasMeta = MetaDataInName::nameHasMetaData(seq.name_);
        if(hasMeta){
          inMeta = MetaDataInName(seq.name_);
        }
        const bool oneEndRead = hasMeta && inMeta.containsMeta("barcodeEnds") && "1" == inMeta.getMeta("barcodeEnds");
        if(segmentTargets.size() > 1 && !(oneEndRead && "none" == splitPars.oneEnd_)){
          //a target seen again in the other orientation means the molecule was read through and back (a foldback),
          //so only the stretches before it turns are separate molecules
          std::vector<MultihitSplitter::Segment> kept;
          std::unordered_map<uint32_t, bool> orientationSeen;
          bool foldback = false;
          for(const auto & seg : segments){
            auto seen = orientationSeen.find(seg.setIdx_);
            if(seen != orientationSeen.end() && seen->second != seg.revComp_){
              foldback = true;
              break;
            }
            orientationSeen[seg.setIdx_] = seg.revComp_;
            kept.emplace_back(seg);
          }
          if(foldback){
            ++multihitSplitFoldbackCurrent;
          }
          //demultiplexed on one end's barcode, the amplicon nearest that barcode is the one it labels best
          if("nearest" == splitPars.oneEnd_ && oneEndRead && inMeta.containsMeta("assignedOn") && inMeta.containsMeta("direction") && kept.size() > 1){
            const bool barcodeAtFront = foldback || (("back" != inMeta.getMeta("assignedOn")) == ("rev" != inMeta.getMeta("direction")));
            const auto nearest = barcodeAtFront ? kept.front() : kept.back();
            kept = {nearest};
          }
          //a piece has the read's barcode at one of its own ends at most, so it doesn't count as support from both ends
          //(processClusters --minBothEndsBarcodeReads), otherwise a few joined molecules from other samples let bleed through, e.g. into no template controls
          std::string keptMeta;
          if(hasMeta && inMeta.containsMeta("barcodeEnds")){
            MetaDataInName barcodeEndsMeta;
            barcodeEndsMeta.addMeta("barcodeEnds", 1);
            keptMeta = barcodeEndsMeta.createMetaName();
          }
          ++readsPerSetCurrent["multihitSplit"];
          seqOut.openWrite("multihitSplit", seq);
          const uint64_t pad = maxPrimerSize + splitPars.primerSlop_;
          const auto threadId = estd::to_string(std::this_thread::get_id());
          for(const auto & seg : kept){
            const uint32_t start = seg.start_ > pad ? seg.start_ - pad : 0;
            const uint32_t end = std::min<uint64_t>(len(seq), seg.end_ + pad);
            auto piece = seq.getSubRead(start, end - start);
            const std::string & pieceSet = splitIndex->setNames_[seg.setIdx_];
            if(seg.revComp_){
              ++piecesPerSetRevCompCurrent[pieceSet];
            } else {
              ++piecesPerSetCurrent[pieceSet];
            }
            if(rename){
              piece.name_ = njh::pasteAsStr(sampleName, ".", pieceSet, ".mhSplit", piecesPerSetCurrent[pieceSet] + piecesPerSetRevCompCurrent[pieceSet], ".", threadId, keptMeta);
              if(seg.revComp_){
                piece.reverseComplementRead(false, true);
                piece.name_.append("_Comp");
              }
            } else {
              piece.name_ = njh::pasteAsStr(seq.name_, "_mhSplit", start, "-", end);
              if(seg.revComp_){
                piece.reverseComplementRead(true, true);
              }
            }
            seqOut.openWrite(pieceSet, piece);
            if(trimToRegions){
              extractTrimRegions(piece, pieceSet, ExtractionStator::extractCase::GOOD);
            } else {
              extractTarget(piece, pieceSet, ExtractionStator::extractCase::GOOD, true);
            }
          }
          continue;
        }
      }
      if(winnerRevComp){
        ++readsPerSetRevCompCurrent[winnerSet];
        seq.reverseComplementRead(true, true);
      } else {
        ++readsPerSetCurrent[winnerSet];
      }

      if(rename){
        //keep how many of the read's ends had the sample's barcode (nanopore demultiplexing), for clustering to count
        std::string keptMeta;
        if(MetaDataInName::nameHasMetaData(seq.name_)){
          MetaDataInName inMeta(seq.name_);
          if(inMeta.containsMeta("barcodeEnds")){
            MetaDataInName barcodeEndsMeta;
            barcodeEndsMeta.addMeta("barcodeEnds", inMeta.getMeta("barcodeEnds"));
            keptMeta = barcodeEndsMeta.createMetaName();
          }
        }
        auto threadId = estd::to_string(std::this_thread::get_id());
        seq.name_ = njh::pasteAsStr(sampleName, ".",winnerSet, ".", readsPerSetCurrent[winnerSet] + readsPerSetRevCompCurrent[winnerSet], ".", threadId, keptMeta);
        if(winnerRevComp){
          seq.name_.append("_Comp");
        }
      }

      seqOut.openWrite(winnerSet, seq);
      if("undetermined" == winnerSet){
        ++masterCountsCurrent.readsUnrecBarcode_;
        continue;
      }
      if("multihit" == winnerSet){
        ++masterCountsCurrent.multihit_;
        continue;
      }

      auto extractorCase = ExtractionStator::extractCase::GOOD;

      //check for inverse hit (has hit of the same target in both the forward and reverse directions)
      if (kmerMatch.found[kmerMatch.winnerIdx] >= hardKmersPerTarget && kmerMatch.foundRevComp[kmerMatch.winnerIdx]  >= hardKmersPerTarget) {
        extractorCase = ExtractionStator::extractCase::INVERSECHIMERA;
        seq.on_ = false;
      }
      //sequence of no target on either side of the amplicon (e.g. a human off-target product ligated onto it) pushes the primers
      //out of the primer search window, so cut the read back to the winner's stretch of unique k-mers plus its primers and some slop;
      //the read is already in the winner's forward direction, a read with no stretch (too few unique k-mers) is left as is
      bool cutToTarget = false;
      //when trimming to regions the regions are found on their own, so here the stretch is only for whether the target reaches the read's ends
      if(seq.on_ && (kmerStretchTrim || trimToRegions)){
        std::vector<char> allowed(splitIndex->setNames_.size(), 0);
        allowed[matcherToSplitIdx[kmerMatch.winnerIdx]] = 1;
        const auto segments = MultihitSplitter::findSegments(seq.seq_, *splitIndex, splitHasher, extractionKmer, allowed, splitMaxSpan, trimSegmentPars[matcherToSplitIdx[kmerMatch.winnerIdx]]);
        uint32_t stretchStart = std::numeric_limits<uint32_t>::max();
        uint32_t stretchEnd = 0;
        //several stretches in the same direction (an error dense patch, or copies of the target) are kept together, only the outside is trimmed
        for(const auto & seg : segments){
          if(!seg.revComp_){
            stretchStart = std::min(stretchStart, seg.start_);
            stretchEnd = std::max(stretchEnd, seg.end_);
          }
        }
        if(stretchEnd > 0){
          const auto & primerLens = primerLensPerTarget.at(winnerSet);
          const uint32_t frontPad = primerLens.first + (trimToRegions ? std::max<uint32_t>(kmerStretchTrimSlop, corePars.pDetPars.primerWithin_) : kmerStretchTrimSlop);
          const uint32_t backPad = primerLens.second + (trimToRegions ? std::max<uint32_t>(kmerStretchTrimSlop, corePars.pDetPars.primerWithin_) : kmerStretchTrimSlop);
          const uint32_t trimStart = stretchStart > frontPad ? stretchStart - frontPad : 0;
          const uint32_t trimEnd = std::min<uint64_t>(len(seq), static_cast<uint64_t>(stretchEnd) + backPad);
          bool trim = trimStart > 0 || trimEnd < len(seq);
          if(trim && !trimToRegions){
            //a read that would pass as it is (both primers, then the length and quality checks) isn't touched, so trimming never changes or
            //loses a read that passes without it (on short amplicons the first unique k-mer can start well past the primer and the trim would clip it).
            //Finding both primers alone isn't enough, a reverse primer matched by chance in the extra sequence extracts a read that's too long
            auto untrimmed = seq;
            bool passesAsIs = winnerSet == ids.pDeterminator_->determineForwardPrimer(untrimmed, corePars.pDetPars, *currentAlignerObj, VecStr{winnerSet});
            if(passesAsIs){
              untrimmed.reverseComplementRead(false, true);
              passesAsIs = winnerSet == ids.pDeterminator_->determineWithReversePrimer(untrimmed, corePars.backEndpDetPars, *currentAlignerObj, VecStr{winnerSet});
            }
            if(passesAsIs){
              ids.targets_.at(winnerSet).lenCuts_->minLenChecker_.checkRead(untrimmed);
              if(untrimmed.on_){
                ids.targets_.at(winnerSet).lenCuts_->maxLenChecker_.checkRead(untrimmed);
              }
              if(untrimmed.on_){
                qualChecker->checkRead(untrimmed);
              }
              passesAsIs = untrimmed.on_;
            }
            trim = !passesAsIs;
          }
          if(trim){
            if(!trimToRegions){
              ++kmerTrimmedPerSetCurrent[winnerSet];
              kmerTrimmedFrontBasesCurrent += trimStart;
              kmerTrimmedBackBasesCurrent += len(seq) - trimEnd;
              seq = seq.getSubRead(trimStart, trimEnd - trimStart);
              cutToTarget = true;
            }
            //like a multihit piece, the amplicon no longer reaches the cut end, so that end's barcode doesn't vouch for it and the read
            //doesn't count as support from both ends (processClusters --minBothEndsBarcodeReads); kept at 2, rescued molecules joined to
            //another sample's molecule put DHFR and CRT calls (18-25 reads) in run 2 and run 3 diversity panel samples
            if(MetaDataInName::nameHasMetaData(seq.name_)){
              MetaDataInName trimmedMeta(seq.name_);
              if(trimmedMeta.containsMeta("barcodeEnds") && "2" == trimmedMeta.getMeta("barcodeEnds")){
                trimmedMeta.addMeta("barcodeEnds", 1, true);
                trimmedMeta.resetMetaInName(seq.name_);
              }
            }
          }
        }
      }
      if(trimToRegions){
        extractTrimRegions(seq, winnerSet, extractorCase);
      } else {
        extractTarget(seq, winnerSet, extractorCase, cutToTarget);
      }
    }
    {
      std::lock_guard<std::mutex> lockGuard(mut);

      masterCounts.addOtherExtractorCounts(masterCountsCurrent);
      perMidCounts.addOtherExtractorCounts(perMidCountsCurrent);

      for(const auto & readsPerSetCount : readsPerSetCurrent){
        readsPerSet[readsPerSetCount.first] += readsPerSetCount.second;
      }
      for(const auto & readsPerSetRevCompCount : readsPerSetRevCompCurrent){
        readsPerSetRevComp[readsPerSetRevCompCount.first] += readsPerSetRevCompCount.second;
      }
      for(const auto & count : piecesPerSetCurrent){
        piecesPerSet[count.first] += count.second;
      }
      for(const auto & count : piecesPerSetRevCompCurrent){
        piecesPerSetRevComp[count.first] += count.second;
      }
      multihitSplitFoldback += multihitSplitFoldbackCurrent;
      for(const auto & count : kmerTrimmedPerSetCurrent){
        kmerTrimmedPerSet[count.first] += count.second;
      }
      kmerTrimmedFrontBases += kmerTrimmedFrontBasesCurrent;
      kmerTrimmedBackBases += kmerTrimmedBackBasesCurrent;
      adapterTrimmed += adapterTrimmedCurrent;
      adapterTrimmedBases += adapterTrimmedBasesCurrent;
      for(const auto & count : notCoveredPerRegionCurrent){
        notCoveredPerRegion[count.first] += count.second;
      }
      for(const auto & count : readsPassingAnyRegionCurrent){
        readsPassingAnyRegion[count.first] += count.second;
      }
    }
  };

  njh::concurrent::runVoidFunctionThreaded(readInComp, numThreads);


  uint64_t totalReadsProcessed = 0;
  for(const auto & readsPerSetCount : readsPerSet){
    totalReadsProcessed += readsPerSetCount.second;
  }
  for(const auto & readsPerSetRevCompCount : readsPerSetRevComp){
    totalReadsProcessed += readsPerSetRevCompCount.second;
  }

  uint64_t totalReadsPassedAllFilters = 0;


  outCounts << "sample\ttotalReadsProcessed\ttarget\tcount\tfrac\tforwardCount\tfracForward";
  outCounts << "\tminLenFailed\tminLenFailedFrac";
  outCounts << "\tinverseChimeraFailed\tinverseChimeraFailedFrac";
  outCounts << "\tmaxLenFailed\tmaxLenFailedFrac";
  outCounts << "\tqualityFailed\tqualityFailedFrac";
  outCounts << "\tforwardPrimerFailed\tforwardPrimerFailedFrac";
  outCounts << "\treversePrimerFailed\treversePrimerFailedFrac";
  outCounts << "\tbothForRevPrimerFailed\tbothForRevPrimerFailedFrac";
  outCounts << "\tMIDFailed\tMIDFailedFrac";
  outCounts << "\ttotalFailed\ttotalFailedFrac";
  outCounts << "\tpassed\tpassedFrac";
  outCounts << "\tfromMultihitSplit";
  outCounts << "\tkmerStretchTrimmed";
  if(trimToRegions){
    outCounts << "\ttrimRegionNotCovered\ttrimRegionTargets";
  }
  outCounts << std::endl;


  uint64_t totalExtractedAllTargets = 0;
  uint64_t totalExtractedAllTargetsForward = 0;
  uint64_t totalExtractedUndetermined = readsPerSet["undetermined"] + readsPerSetRevComp["undetermined"];
  uint64_t totalExtractedMultihit = readsPerSet["multihit"] + readsPerSetRevComp["multihit"];


  //the amplified targets the output is for, themselves unless trimming to regions, then the reads of the targets a region is in
  std::unordered_map<std::string, VecStr> targetsPerOutput;
  for(const auto & setName : ids.getTargets()){
    if(!trimToRegions){
      targetsPerOutput[setName] = VecStr{setName};
    }
    totalExtractedAllTargets += readsPerSet[setName] + readsPerSetRevComp[setName] + piecesPerSet[setName] + piecesPerSetRevComp[setName];
    totalExtractedAllTargetsForward += readsPerSet[setName] + piecesPerSet[setName];
  }
  for(const auto & target : trimRegionsPerTarget){
    for(const auto & region : target.second){
      targetsPerOutput[region.name_].emplace_back(target.first);
    }
  }

  for(const auto & setName : outIds.getTargets()){
    uint64_t fromMultihitSplit = 0;
    uint64_t totalExtracted = 0;
    uint64_t totalExtractedForward = 0;
    for(const auto & target : targetsPerOutput[setName]){
      fromMultihitSplit += piecesPerSet[target] + piecesPerSetRevComp[target];
      totalExtracted += readsPerSet[target] + readsPerSetRevComp[target] + piecesPerSet[target] + piecesPerSetRevComp[target];
      totalExtractedForward += readsPerSet[target] + piecesPerSet[target];
    }

    uint64_t totalBad = masterCounts.counts_[setName][false].minLenBad_ + masterCounts.counts_[setName][true].minLenBad_ +
          masterCounts.counts_[setName][false].maxLenBad_ + masterCounts.counts_[setName][true].maxLenBad_+
          masterCounts.counts_[setName][false].inverse_chimera_ + masterCounts.counts_[setName][true].inverse_chimera_+
          masterCounts.counts_[setName][false].qualityFailed_ + masterCounts.counts_[setName][true].qualityFailed_ +
          masterCounts.counts_[setName][false].badForward_ + masterCounts.counts_[setName][true].badForward_ +
          masterCounts.counts_[setName][false].badReverse_ + masterCounts.counts_[setName][true].badReverse_ +
          masterCounts.counts_[setName][false].failedBothPrimers_ + masterCounts.counts_[setName][true].failedBothPrimers_+
          masterCounts.counts_[setName][false].badmid_ + masterCounts.counts_[setName][true].badmid_ +
          notCoveredPerRegion[setName];

    outCounts << sampleName
              << "\t" << totalReadsProcessed
              << "\t" << setName
              << "\t" << totalExtracted
              << "\t" << static_cast<double>(totalExtracted) / static_cast<double>(totalReadsProcessed)
              << "\t" << totalExtractedForward
              << "\t" << (totalExtracted == 0 ? 0 : static_cast<double>(totalExtractedForward) / static_cast<double>(totalExtracted));
    //minlen
    outCounts<< "\t" << masterCounts.counts_[setName][false].minLenBad_ + masterCounts.counts_[setName][true].minLenBad_
    << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].minLenBad_ + masterCounts.counts_[setName][true].minLenBad_)/static_cast<double>(totalBad));
    //inverse chimera
    outCounts<< "\t" << masterCounts.counts_[setName][false].inverse_chimera_ + masterCounts.counts_[setName][true].inverse_chimera_
    << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].inverse_chimera_ + masterCounts.counts_[setName][true].inverse_chimera_)/static_cast<double>(totalBad));
    //maxlen
    outCounts << "\t" << masterCounts.counts_[setName][false].maxLenBad_ + masterCounts.counts_[setName][true].maxLenBad_
    << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].maxLenBad_ + masterCounts.counts_[setName][true].maxLenBad_)/static_cast<double>(totalBad));
    //quality
    outCounts << "\t" << masterCounts.counts_[setName][false].qualityFailed_ + masterCounts.counts_[setName][true].qualityFailed_
              << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].qualityFailed_ + masterCounts.counts_[setName][true].qualityFailed_)/static_cast<double>(totalBad));
    //failed forward
    outCounts << "\t" << masterCounts.counts_[setName][false].badForward_ + masterCounts.counts_[setName][true].badForward_
              << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].badForward_ + masterCounts.counts_[setName][true].badForward_)/static_cast<double>(totalBad));
    //bad reverse
    outCounts << "\t" << masterCounts.counts_[setName][false].badReverse_ + masterCounts.counts_[setName][true].badReverse_
              << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].badReverse_ + masterCounts.counts_[setName][true].badReverse_)/static_cast<double>(totalBad));
    //both bad reverse and failed forward
    outCounts << "\t" << masterCounts.counts_[setName][false].failedBothPrimers_ + masterCounts.counts_[setName][true].failedBothPrimers_
              << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].failedBothPrimers_ + masterCounts.counts_[setName][true].failedBothPrimers_)/static_cast<double>(totalBad));
    //failed MID
    outCounts << "\t" << masterCounts.counts_[setName][false].badmid_ + masterCounts.counts_[setName][true].badmid_
              << "\t" << (totalBad == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].badmid_ + masterCounts.counts_[setName][true].badmid_)/static_cast<double>(totalBad));
    //bad
    outCounts << "\t" << totalBad
              << "\t" << (totalExtracted == 0 ? 0 : static_cast<double>(totalBad) / static_cast<double>(totalExtracted));
    //good
    outCounts << "\t" << masterCounts.counts_[setName][false].good_ + masterCounts.counts_[setName][true].good_
              << "\t" << (totalExtracted == 0 ? 0 : static_cast<double>(masterCounts.counts_[setName][false].good_ + masterCounts.counts_[setName][true].good_) / static_cast<double>(totalExtracted))
              << "\t" << fromMultihitSplit
              << "\t" << kmerTrimmedPerSet[setName];
    if(trimToRegions){
      outCounts << "\t" << notCoveredPerRegion[setName]
                << "\t" << njh::conToStr(targetsPerOutput[setName], ",");
    }
    outCounts << std::endl;
    totalReadsPassedAllFilters += masterCounts.counts_[setName][false].good_ + masterCounts.counts_[setName][true].good_;
  }
  {
    uint64_t totalExtracted = readsPerSet["undetermined"] + readsPerSetRevComp["undetermined"];
    outCounts << sampleName
              << "\t" << totalReadsProcessed
              << "\t" << "undetermined"
              << "\t" << totalExtracted
              << "\t" << (totalReadsProcessed == 0 ? 0 : static_cast<double>(totalExtracted) / static_cast<double>(totalReadsProcessed))
              << "\t" << readsPerSet["undetermined"]
              << "\t" << (totalExtracted == 0 ? 0 : static_cast<double>(readsPerSet["undetermined"]) / static_cast<double>(totalExtracted))
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0
        << "\t" << 0;
    if(trimToRegions){
      outCounts << "\t" << 0 << "\t" << "";
    }
    outCounts << std::endl;
  }

  outStats << "sampleName\ttotalReadsProcessed\tfailedMinLen_" << corePars.smallFragmentCutoff << "\tfailedMinLenFrac\tfailedMaxLen_" << corePars.largeFragmentCutOff << "\tfailedMaxLenFrac\tundetermined\tundeterminedFrac\tmultihit\tmultihitFrac\textracted\textractedFrac\textractedForward\textractedForwardFrac\tpassed\tpassedFrac\tmultihitSplit\tmultihitSplitFrac\tpiecesFromMultihitSplit\tmultihitSplitFoldback\tkmerStretchTrimmed\tkmerStretchTrimmedFrontBases\tkmerStretchTrimmedBackBases" << (trimToRegions ? "\treadsPassingAnyTrimRegion" : "") << "\tadapterTrimmed\tadapterTrimmedBases" << std::endl;
  //reads that hit several targets and were split into pieces, the pieces are counted in extracted and passed
  const uint64_t totalMultihitSplit = readsPerSet["multihitSplit"];
  uint64_t totalPiecesFromMultihitSplit = 0;
  for(const auto & count : piecesPerSet){
    totalPiecesFromMultihitSplit += count.second;
  }
  for(const auto & count : piecesPerSetRevComp){
    totalPiecesFromMultihitSplit += count.second;
  }
  uint64_t totalKmerTrimmed = 0;
  for(const auto & count : kmerTrimmedPerSet){
    totalKmerTrimmed += count.second;
  }
  uint32_t grand_total_reads = totalReadsProcessed + masterCounts.smallFrags_ + masterCounts.largeFrags_;
  outStats << sampleName
					 << "\t" << grand_total_reads
					 << "\t" << masterCounts.smallFrags_
					 << "\t" << static_cast<double>(masterCounts.smallFrags_) / static_cast<double>(grand_total_reads)
           << "\t" << masterCounts.largeFrags_
           << "\t" << static_cast<double>(masterCounts.largeFrags_) / static_cast<double>(grand_total_reads)
					 << "\t" << totalExtractedUndetermined
					 << "\t" << static_cast<double>(totalExtractedUndetermined) / static_cast<double>(totalReadsProcessed)
           << "\t" << totalExtractedMultihit
           << "\t" << static_cast<double>(totalExtractedMultihit) / static_cast<double>(totalReadsProcessed)
					 << "\t" << totalExtractedAllTargets
					 << "\t" << static_cast<double>(totalExtractedAllTargets) / static_cast<double>(grand_total_reads)
					 << "\t" << totalExtractedAllTargetsForward
					 << "\t" << static_cast<double>(totalExtractedAllTargetsForward) / static_cast<double>(totalExtractedAllTargets)
           << "\t" << totalReadsPassedAllFilters
           << "\t" << static_cast<double>(totalReadsPassedAllFilters) / static_cast<double>(grand_total_reads)
           << "\t" << totalMultihitSplit
           << "\t" << static_cast<double>(totalMultihitSplit) / static_cast<double>(totalReadsProcessed)
           << "\t" << totalPiecesFromMultihitSplit
           << "\t" << multihitSplitFoldback
           << "\t" << totalKmerTrimmed
           << "\t" << kmerTrimmedFrontBases
           << "\t" << kmerTrimmedBackBases;
  if(trimToRegions){
    uint64_t totalReadsPassingAnyRegion = 0;
    for(const auto & count : readsPassingAnyRegion){
      totalReadsPassingAnyRegion += count.second;
    }
    outStats << "\t" << totalReadsPassingAnyRegion;
  }
  outStats << "\t" << adapterTrimmed
           << "\t" << adapterTrimmedBases;
  outStats << std::endl;

  if(ids.containsMids()){
    //the counts per target per MID, name matches the extracted file <target><MID> so it can be linked to the sample the MID belongs to
    //reads failing before the MID is looked for (inverse chimera, min length, either primer) are counted under MIDunknown, reads whose MID wasn't found under MIDUnrecognized
    OutputStream outPerMid(njh::files::make_path(setUp.pars_.directoryName_, "extractionProfilePerMID.tab.txt"));
    outPerMid << "sample\ttarget\tmid\tname\ttotal"
              << "\tpassed\tpassedFrac"
              << "\tminLenFailed\tmaxLenFailed\tqualityFailed\tinverseChimeraFailed"
              << "\tforwardPrimerFailed\treversePrimerFailed\tbothForRevPrimerFailed\tMIDFailed"
              << "\ttotalFailed" << std::endl;
    auto midNames = getVectorOfMapKeys(ids.mids_);
    njh::sort(midNames);
    midNames.emplace_back("MIDUnrecognized");
    midNames.emplace_back("MIDunknown");
    for(const auto & setName : ids.getTargets()){
      for(const auto & midName : midNames){
        const bool undeterminedMid = "MIDUnrecognized" == midName || "MIDunknown" == midName;
        const std::string name = undeterminedMid ? njh::pasteAsStr(setName, "-", midName) : njh::pasteAsStr(setName, midName);
        //both read directions
        ExtractionStator::extractCounts counts;
        if(njh::in(name, perMidCounts.counts_)){
          for(const auto & dir : perMidCounts.counts_.at(name)){
            counts.good_ += dir.second.good_;
            counts.bad_ += dir.second.bad_;
            counts.minLenBad_ += dir.second.minLenBad_;
            counts.maxLenBad_ += dir.second.maxLenBad_;
            counts.qualityFailed_ += dir.second.qualityFailed_;
            counts.inverse_chimera_ += dir.second.inverse_chimera_;
            counts.badForward_ += dir.second.badForward_;
            counts.badReverse_ += dir.second.badReverse_;
            counts.failedBothPrimers_ += dir.second.failedBothPrimers_;
            counts.badmid_ += dir.second.badmid_;
          }
        }
        const uint32_t total = counts.good_ + counts.bad_;
        outPerMid << sampleName
                  << "\t" << setName
                  << "\t" << midName
                  << "\t" << name
                  << "\t" << total
                  << "\t" << counts.good_
                  << "\t" << (0 == total ? 0 : static_cast<double>(counts.good_) / static_cast<double>(total))
                  << "\t" << counts.minLenBad_
                  << "\t" << counts.maxLenBad_
                  << "\t" << counts.qualityFailed_
                  << "\t" << counts.inverse_chimera_
                  << "\t" << counts.badForward_
                  << "\t" << counts.badReverse_
                  << "\t" << counts.failedBothPrimers_
                  << "\t" << counts.badmid_
                  << "\t" << counts.bad_
                  << std::endl;
      }
    }
  }

  return 0;
}

} // namespace njhseq
