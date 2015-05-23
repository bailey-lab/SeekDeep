#!/usr/bin/env python

import sys,os,argparse,pysam,ucscgenome
from collections import defaultdict

def oneCharStr(character, number):
    return character * number

def rearrangeSeqs(gDna, cigarTuples):
    '''
        rearrange the genomeic dna by adding - for insertions (1) or skipped region (3) 
    '''
    ret = gDna
    pos = 0
    numBases = 0
    for cTup in cigarTuples:
        if(cTup[0]  == 0):
            numBases += cTup[1]
        elif (cTup[0] == 1 or cTup[0] == 3):
            ret = ret[0:pos] + oneCharStr("-",cTup[1]) + ret[pos:]
            pos += numBases
            numBases = 0
        elif (cTup[0] == 2):
            numBases += cTup[1]
    return ret
          
        
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str, required = True, help = "Bam file")
    parser.add_argument('--outFile', type=str, required = True, help = "Out filename")
    parser.add_argument('--libraryName', type=str, required = True, help = "Name of library in bamfile")
    parser.add_argument('--overWrite',dest = 'overWrite', action = 'store_true', default = False, help = "Over write outfile if it exists", )
    return parser.parse_args()
    

def main():
    args = parse_args()
    if(os.path.exists(args.outFile) and not args.overWrite):
        sys.stdout.write("Error, file "  + args.outFile + " already exist, use --overWrite to overwrite the file")
        return 1
    yeastG = ucscgenome.Genome("sacCer3")
    samfile = pysam.AlignmentFile(args.file, "rb")
    totalBases = 0
    baseRates = defaultdict(lambda: defaultdict(int))
    
    for read in samfile.fetch():
        gDna = yeastG[samfile.getrname(read.reference_id)]
        gDna = gDna[read.reference_start:read.reference_end]
        gDna = rearrangeSeqs(gDna, read.cigartuples).upper()
        mapped = read.query_alignment_sequence.upper()
        for pos in xrange(read.query_alignment_sequence):
            if(gDna[pos] != "-" and mapped[pos] != "-"):
                totalBases += 1
                baseRates[gDna[pos]][mapped[pos]] += 1
    samfile.close()
    yeastG.close()
    outFile = open(args.outFile,"w")
    outFile.write("\t".join(["library","refBase", "seqBase","freq", "rate"]))
    for refBase,refBaseCounts in baseRates.iteritems():
        for seqBase,seqBaseCount in refBaseCounts.iteritems():
            outFile.write("\t".join([args.libraryName,refBase, seqBase,seqBaseCount, str(seqBaseCount/float(totalBases))]) + "\n")

main() 
            
            

