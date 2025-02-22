//
// Created by Nicholas Hathaway on 7/23/24.
//

#include "SeekDeepUtilsRunner.hpp"
#include "SeekDeep/utils.h"


#include <njhseq/objects/kmer/KmerGatherer.hpp>
#include <njhseq/objects/kmer/SimpleKmerHash.hpp>


namespace njhseq {

int SeekDeepUtilsRunner::genTargetInfoFromGenomes(const njh::progutils::CmdArgs & inputCommands) {
	extractBetweenSeqsPars pars;

  KmerGatherer::KmerGathererPars countPars;
  countPars.noRevComp_ = true;
  countPars.kmerLength_ = 19;
  countPars.entropyFilter_ = 1.20;
	uint32_t minOverlap  = 10;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pars.primersFile, "--primers", "A file that contains three columns, target,forwardPrimer,reversePrimer 5` to 3` directions, same file as the input to SeekDeep", true);
	pars.verbose_ = setUp.pars_.verbose_;
	pars.debug_ = setUp.pars_.debug_;
	pars.setUpCoreOptions(setUp, true);
  countPars.numThreads_ = pars.pars.numThreads_;
	setUp.setOption(minOverlap, "--minOverlap", "Minimum overlap for stitching");

  setUp.setOption(countPars.kmerLength_, "--uniqKmerLength", "kmer Length");
  setUp.setOption(countPars.allowableCharacters_, "--allowableCharactersForUniqKmer",
                  "Only count kmers with these allowable Characters");
  setUp.setOption(countPars.entropyFilter_, "--entropyFilterForUniqKmer", "entropy Filter cut off, exclusive, will only keep kmers abovet this entropy level");


  setUp.finishSetUp(std::cout);

	njh::sys::requireExternalProgramThrow("bowtie2");
	njh::sys::requireExternalProgramThrow("samtools");

	PrimersAndMids ids(pars.primersFile);

	if(ids.getTargets().empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error in reading in target primers file " << pars.primersFile << "\n";
		ss << "Make sure there is a line that starts with target in file" << "\n";
		throw std::runtime_error{ss.str()};
	}
	ids.initPrimerDeterminator();

	extractBetweenSeqs(ids, pars);

	setUp.startARunLog(pars.outputDirPars.dirName_.string());


	auto forSeekDeepDir = njh::files::make_path(pars.outputDirPars.dirName_, "forSeekDeep");

	njh::files::makeDir(njh::files::MkdirPar{forSeekDeepDir});
	auto refSeqsDir = njh::files::make_path(forSeekDeepDir, "refSeqs");
	njh::files::makeDir(njh::files::MkdirPar{refSeqsDir});
	OutOptions lenCutOffsOpts(njh::files::make_path(forSeekDeepDir, "lenCutOffs.txt"));
	OutputStream lenCutOffsOut(lenCutOffsOpts);
	lenCutOffsOut << "target\tminlen\tmaxlen" << "\n";

	OutOptions overlapStatusOpts(njh::files::make_path(forSeekDeepDir, "overlapStatuses.txt"));
	OutputStream overlapStatusOut(overlapStatusOpts);
	overlapStatusOut << "target\tstatus" << "\n";
	for(const auto & tar : ids.getTargets()){
		auto primersRemovedFnp = njh::files::make_path(pars.outputDirPars.dirName_, tar, tar + "_primersRemoved.fasta");
		auto extractedSeqsFnp = njh::files::make_path(pars.outputDirPars.dirName_, tar, tar + ".fasta");
		auto primersRemovedFinalFnp = njh::files::make_path(refSeqsDir, tar + ".fasta");
		if(bfs::exists(primersRemovedFnp)){
			{
				auto finalSeqOpts = SeqIOOptions::genFastaInOut(primersRemovedFnp, primersRemovedFinalFnp);
				SeqIO reader(finalSeqOpts);
				reader.openIn();
				reader.openOut();
				seqInfo seq;
				while(reader.readNextRead(seq)){
					reader.write(seq);
				}
			}
			{
				std::vector<uint32_t> readLengths;
				SeqInput reader(SeqIOOptions::genFastaIn(extractedSeqsFnp));
				reader.openIn();
				seqInfo seq;
				while(reader.readNextRead(seq)){
					readLengths.emplace_back(len(seq));
				}
				auto minlen = vectorMinimum(readLengths);
				auto maxlen = vectorMaximum(readLengths);
				lenCutOffsOut << tar
						<< "\t" << (minlen > pars.minLenCutOffSizeExpand ? minlen - pars.minLenCutOffSizeExpand : 0)
						<< "\t" << maxlen + pars.maxLenCutOffSizeExpand << std::endl;
				uint32_t finalMaxSize = maxlen + pars.barcodeSize;
				uint32_t finalMinSize = minlen + pars.barcodeSize;

				uint32_t maxInsertSize = 2 * pars.pairedEndLength - minOverlap;

				std::string status;
				std::set<std::string> statuses;
				if(finalMaxSize > maxInsertSize){
					statuses.emplace("NoOverLap");
				} else {
					if(finalMaxSize >= pars.pairedEndLength){
						statuses.emplace("R1EndsInR2");
					}
					if(finalMinSize < pars.pairedEndLength){
						statuses.emplace("R1BeginsInR2");
					}

					if((finalMaxSize >= pars.pairedEndLength && finalMinSize < pars.pairedEndLength) ||
						(uAbsdiff(finalMaxSize, pars.pairedEndLength) < 10 || uAbsdiff(finalMinSize, pars.pairedEndLength) < 10)) {
						//add PerfectOverlap as a possible overlap status if reference seqs are above and below the expected paired end reads
						//or if the min or max size is within 10 bases of the paired end
						statuses.emplace("PerfectOverlap");
					}
				}
				overlapStatusOut << tar << "\t" << njh::conToStr(statuses, ",") << std::endl;
			}
		} else {
			std::cerr << "Warning, no sequences extracted for " << tar << std::endl;
		}
	}

  {

    KmerGatherer kGather(countPars);

    std::unordered_map<std::string, std::set<std::string>> fastasForSet;
    for(const auto & tar : ids.getTargets()){
      auto fastaFnp = njh::files::make_path(refSeqsDir, tar + ".fasta");
      if(bfs::exists(fastaFnp)){
				fastasForSet[tar].emplace(fastaFnp.string());
      }
    }
    std::vector<bfs::path> fastaFiles;
    for(const auto & seqSet : fastasForSet){
      for(const auto & fnp : seqSet.second){
				fastaFiles.emplace_back(fnp);
      }
    }
    std::map<std::string, std::set<uint64_t>> kmersPerSet;

		std::function<bool(const std::string &)> seqCheck = [&countPars](const std::string &k) {
			return std::all_of(k.begin(), k.end(),
												 [&countPars](char base) { return njh::in(base, countPars.allowableCharacters_); });
		};


    {
      auto allKmers = kGather.getUniqueKmersSetHashWithFiltersFromFastas(fastaFiles);
      setUp.rLog_.logCurrentTime("condense");
      setUp.rLog_.runLogFile_.flush();
      njh::concurrent::LockableQueue<std::string> seqSetNamesQueue(getVectorOfMapKeys(fastasForSet));
      for(const auto & name : fastasForSet){
        kmersPerSet[name.first] = std::set<uint64_t>{};
      }
      std::function<void()> condenseKmers = [&seqSetNamesQueue,&allKmers,&fastasForSet,&kmersPerSet](){
        std::string name;
        while(seqSetNamesQueue.getVal(name)){
          SimpleKmerHash hasher;
          for(const auto & fasta : fastasForSet.at(name)){
            for(const auto & k : allKmers.at(fasta)){
              kmersPerSet[name].emplace(k);
            }
          }
        }
      };
      njh::concurrent::runVoidFunctionThreaded(condenseKmers, countPars.numThreads_);
    }
    std::map<std::string, std::set<uint64_t>> uniqueKmersFinal;
    setUp.rLog_.logCurrentTime("compare");
    setUp.rLog_.runLogFile_.flush();

		{
    	SimpleKmerHash hasher;
    	OutputStream out(njh::files::make_path(forSeekDeepDir, "allKmers.tab.txt.gz"));
    	OutputStream outInfo(njh::files::make_path(forSeekDeepDir, "allKmersCounts.tsv"));
    	outInfo << "target\tKmerCount" << std::endl;
    	for(const auto & kmersForSet : kmersPerSet){
    		outInfo << kmersForSet.first << "\t" << kmersForSet.second.size() << std::endl;
    		for(const auto & kmer : kmersForSet.second){
    			out << kmersForSet.first
							<< "\t" << hasher.reverseHash(kmer) << "\n";
    		}
    	}
		}
    for(const auto & kmersForSet : kmersPerSet){
      uniqueKmersFinal[kmersForSet.first] = std::set<uint64_t>{};
    }
    {

      auto namesFound = getVectorOfMapKeys(kmersPerSet);
      if(namesFound.size() > 1){
        njh::concurrent::LockableQueue<std::string> seqSetNamesQueue(namesFound);
        std::function<void()> compareKmers = [&seqSetNamesQueue,&kmersPerSet,&uniqueKmersFinal](){
          std::string name;
          while(seqSetNamesQueue.getVal(name)){
            std::set<uint64_t> uniqueKmers;
            uint32_t count = 0;
            for(const auto & otherSet : kmersPerSet){
              if(otherSet.first == name){
                continue;
              }
              if(0 == count){
                std::vector<uint64_t> notShared;
                std::set_difference(
                    kmersPerSet.at(name).begin(), kmersPerSet.at(name).end(),
                    otherSet.second.begin(), otherSet.second.end(),
                    std::back_inserter(notShared));
                uniqueKmers = njh::vecToSet(notShared);
              }else{
                std::vector<uint64_t> notShared;
                std::set_difference(
                    uniqueKmers.begin(), uniqueKmers.end(),
                    otherSet.second.begin(), otherSet.second.end(),
                    std::back_inserter(notShared));
                uniqueKmers = njh::vecToSet(notShared);
              }
              ++count;
            }
            uniqueKmersFinal[name] = uniqueKmers;
          }
        };
        njh::concurrent::runVoidFunctionThreaded(compareKmers, countPars.numThreads_);
      }else if( namesFound.size() == 1){
        uniqueKmersFinal[namesFound.front()] = kmersPerSet[namesFound.front()];
      }
    }

    {
    	SimpleKmerHash hasher;
    	OutputStream out(njh::files::make_path(forSeekDeepDir, "uniqueKmers.tab.txt.gz"));
    	OutputStream outInfo(njh::files::make_path(forSeekDeepDir, "uniqueKmersCounts.tsv"));
    	outInfo << "target\tuniqueKmerCount" << std::endl;
    	for(const auto & kmersForSet : uniqueKmersFinal){
    		outInfo << kmersForSet.first << "\t" << kmersForSet.second.size() << std::endl;
    		for(const auto & kmer : kmersForSet.second){
    			out << kmersForSet.first
							<< "\t" << hasher.reverseHash(kmer) << "\n";
    		}
    	}
    }
  }

	return 0;
}


} // namespace njhseq

