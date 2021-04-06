#! /usr/bin/env python

from Bio import SeqIO
import sys
import random
from collections import Counter

gener =  SeqIO.parse(sys.argv[1], 'fasta')
genome = gener.__next__()
out = sys.argv[2]

verbose = True

dmut = dict()

GENERATION = 10
MUT_RAT = 1E-4
MUT_VAR = 1E-5
MULT_GEN = 2
PROBA_RES = 0.04

N_INDIV = 100
LD = 0.2

resistantMutList = [e.split('\t')[0:2] for e in open('mutationResistant.txt').readlines()]
resistantMutList = [(e[0], int(e[1])) for e in resistantMutList]


dmut[0] = dict()
dmut[0]['0'] = list()
#for resName, resPos in resistantMutList:
#    dmut[0][resName] = [resPos]

    

MUT_NUM = round(MUT_RAT * len(genome))
MUT_NUM_VAR = round(MUT_VAR * len(genome))

for i in range(1, GENERATION + 1):
    dmut[i] = dict()
    for name, mutList in dmut[i-1].items():
        for offspring in range(round(random.gauss(MULT_GEN, 1))):
            nMut = random.randint(MUT_NUM - MUT_NUM_VAR, MUT_NUM + MUT_NUM_VAR)
            mutPos = random.sample(range(0, len(genome)), nMut) # Generate a list of MUT numbers (positions) between start and end of genome
            offName = str(offspring)
            for resName, resPos in resistantMutList:
                if random.random() < PROBA_RES and resPos not in mutPos and resPos not in mutList:
                    offName += resName
                    mutPos += [resPos]
            dmut[i][name + offName] = list(set(mutList + mutPos))


totalPos = list()
for n, l in dmut[GENERATION].items():
    totalPos += l


fn = len(dmut[GENERATION])
if verbose:
    print("Individuals", fn)
    resNames = [e[0] for e in resistantMutList]
    for res in resNames:
        print("Individuals with resistance", res, len([e for e in  dmut[GENERATION].keys() if resName in e]))
    print("Individuals wildtype", len([e for e in  dmut[GENERATION].keys() if 'Res' not in e]))

    print("Individuals:")
    for k, v in dmut[GENERATION].items():
        print(k, len(v))



if verbose:
    dist = dict()
    c = Counter(totalPos)
    for k, v in c.items():
        try:
            dist[v/fn] += 1
        except KeyError:
            dist[v/fn] = 1

    print("Allele frequency distribution:")
    for k, v in sorted(dist.items()):
        print(k, v)


if verbose:
    print("Example indiv:")
    for e in dmut[GENERATION].values():
        print(e)
        break

CHROM = genome.name
IND_LIST = random.sample(list(dmut[GENERATION].keys()), k=N_INDIV)
f = "ACGTACG"
mutated = set()
HEADER = """##fileformat=VCFv4.2
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID={NAME},length={LEN}>""".format(NAME=CHROM, LEN=len(genome))
CHR_HEAD = '\t'.join("#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT".split())  + '\t'+ '\t'.join(IND_LIST)

if verbose:
    print("Indiv kept")
    for e in IND_LIST:
          print(e)
          

with open(out, 'w') as outFile:
    print(HEADER, file=outFile)
    print(CHR_HEAD, file=outFile)
    for indiv in IND_LIST:
        mutated.update(set(dmut[GENERATION][indiv]))

    for e in sorted(mutated):
        REF = genome[e]
        ALT = f[f.index(REF) + random.randint(1, 3)]
        geno = list()
        for indiv in IND_LIST:
            if e in dmut[GENERATION][indiv]:
                geno.append('1')
            else:
                geno.append('0')
        # e + 1 because the POSITION are 1-based: https://samtools.github.io/hts-specs/VCFv4.3.pdf
        print(CHROM, e + 1, '.', REF, ALT, '.', '.', '.', 'GT', '\t'.join(geno), sep='\t', file=outFile)
