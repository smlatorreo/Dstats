import vcfpytools
import pandas as pd
from numpy import sqrt, std, nan
from sys import argv, stderr, exit
import argparse

#### Parse options ####
parser = argparse.ArgumentParser()
parser.add_argument('--vcf', action='store', help='The path for the VCF file')
parser.add_argument('--popfile', action = 'store', help = 'Tab separated file with 4-taxa configurations:X <TAB> Y <TAB> Test <TAB> OUT. Incompatible with --pops')
parser.add_argument('--pops', action = 'store', help = 'Comma separated string file with 4-taxa configurations: X,Y,Test.OUT. Incompatible with --popfile')
parser.add_argument('--block', action = 'store', type = int, default = 5000000, help = 'Block size in basepairs. Default=5000000')
parser.add_argument('--sites', action = 'store_true', help = 'Prints out the location of ABBA and BABA sites')
parser.add_argument('--popstats_comp', action = 'store_true', help = 'Reverts the D sign so results are comparable with popstats')
args = parser.parse_args()

#### Checks ####
if args.vcf == None:
    exit('ERROR: No VCF file was provided (--vcf)')
if args.popfile != None and args.pops != None:
    exit('ERROR: Incompatible use of options --popfile and --pops')

#### Taxa configuration ####
if args.popfile != None:
    dstat_configs = [line.strip().split('\t') for line in open(args.popfile, 'r').readlines()]
if args.pops != None:
    dstat_configs = [args.pops.split(',')]
print(dstat_configs)
configs = {'1001':'ABBA','0110':'ABBA',
           '1010':'BABA','0101':'BABA'}

#### D function ####
def Dstat(ABBA, BABA):
    D = (ABBA - BABA) / (ABBA + BABA)
    if args.popstats_comp == True:
        D = (BABA - ABBA) / (ABBA + BABA)
    return(D)

#### Count ABBA and BABA sites in the whole genome keeping their positional information ####
all_samples = list(set([i for l in dstat_configs for i in l]))
counts = {'_'.join(config):{'ABBA':[], 'BABA':[], 'Total':0} for config in dstat_configs}
for i in vcfpytools.get_genotypes_hap(vcf_file=args.vcf, filtered = False, binary=True, samples=all_samples): # The VCF must be filtered
    genotypes = ''.join(i['Genotypes'])
    geno_sample = {all_samples[idx]:g for idx, g in enumerate(genotypes)}
    for config in dstat_configs:
        config_str = '_'.join(config)
        config_genotype = ''.join([geno_sample[ID] for ID in config])
        if '.' not in config_genotype:
            counts[config_str]['Total'] += 1
            try:
                counts[config_str][configs[config_genotype]].append(i['Position'])
            except KeyError:
                continue

# Convert lists into pandas DataFrames
counts_new = {config:{} for config in counts}
for config in counts:
    counts_new[config]['ABBA'] = pd.DataFrame(counts[config]['ABBA'])
    counts_new[config]['BABA'] = pd.DataFrame(counts[config]['BABA'])
counts = counts_new

# Function to compute D per block
def D_block(counts_AB, block_contig, block_start, block_end):
    ABBA_jack = sum(~((counts_AB['ABBA'].iloc[:,0] == block_contig)&(counts_AB['ABBA'].iloc[:,1] >= block_start)&(counts_AB['ABBA'].iloc[:,1] < block_end )))
    BABA_jack = sum(~((counts_AB['BABA'].iloc[:,0] == block_contig)&(counts_AB['BABA'].iloc[:,1] >= block_start)&(counts_AB['BABA'].iloc[:,1] < block_end )))
    D_block = Dstat(ABBA_jack, BABA_jack)
    return D_block

# Compute global D and D per block
D = {config:{'D_global':Dstat(len(counts[config]['ABBA']),len(counts[config]['BABA'])),'D_jacks':[]} for config in counts}

contig_lengths = vcfpytools.get_chr_lengths(args.vcf)
block_size = args.block
for contig in contig_lengths:
    for bl in range(0, contig_lengths[contig], block_size):
        for config in counts:
            D[config]['D_jacks'].append(D_block(counts[config], contig, bl, bl+block_size))

print('X\tY\tTest\tOUT\tD\tStd_err\tZ\tN_ABBA\tN_BABA\tBlock_size\tN_blocks')
#### Compute D per block and print result###
#Algorithm:
#1. Define a number of blocks (T) based on the block size
#2. Iterate T times
#        Removing systematically each block
#        Compute D for the remaining dataset
#        pseudo_D = (D_genomic * T) - (D_jackknifed * (T-1))
for config in D:
    N_blocks = len(D[config]['D_jacks'])
    pseudovalues = [(D[config]['D_global'] * N_blocks) - (i * (N_blocks-1)) for i in D[config]['D_jacks']]
    D_err = std(pseudovalues) / sqrt(len(pseudovalues))
    D_Z = D[config]['D_global'] / D_err
    N_ABBAs = len(counts[config]['ABBA'])
    N_BABAs = len(counts[config]['BABA'])
    print('\t'.join(config.split('_')), D[config]['D_global'], D_err, D_Z, N_ABBAs, N_BABAs, block_size, N_blocks, sep = '\t')
