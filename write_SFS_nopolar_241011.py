import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pdb
import pickle
import argparse
import sys

# python /home/tc557/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/1000Genomes_30X/phased_vcf/SFS_241003/write_SFS_nopolar_241011.py -chrom 2 -pop YRI -out_prefix /tmp/deleteme123

def get_B_sequence(B_file):
    B_data = pd.read_csv(B_file, header = None,sep='\t')
    B_stat = np.array(B_data.loc[:,1:3])
    b_length = int(B_stat[:,1][-1])
    B_sequence = np.zeros(b_length,dtype=float) - 1
    endcheck=False

    z = np.copy(B_stat[:,2])
    z2 = np.copy(B_stat[:,0:2])
    B_stat[:,0]=z
    B_stat[:,1:]=z2

    # B_vals = np.unique(B_stat[:,0])
    # B = [np.where(B_vals==i)[0][0] for i in B_stat[:,0]]
    # B_stat[:,0] = B
    prev=0
    for k in range(0,B_stat.shape[0]):
    # k=0
    # while B_sequence[-1]!=-1:
        if B_stat[k,1]!=prev:
            print(f'\tProblem. There is an unannotated gap in B_file={B_file} at line={k}; nearest_index={prev}.Aborting',flush=True)
            sys.exit()
        # zB_sequence[int(B_stat[k,1]):int(B_stat[k,2])] = B_stat[k,0]

        B_sequence[int(B_stat[k,1]):int(B_stat[k,2])] = B_stat[k,0]
        prev=B_stat[k,2]


    if len(np.where(B_sequence<0)[0])>0:
        print(f'\tProblem. There are negative B values in B_sequence.',flush=True)
        sys.exit()
    return B_sequence


parser = argparse.ArgumentParser(description="Set inputs and outputs")
parser.add_argument('-chrom','--chrom',help='Chromosome',required=True,type=int)
# parser.add_argument('-neanderthal','--neanderthal',help='neanderthal genome',required=True,type=str)
parser.add_argument('-pop','--pop',help='1KGP population',required=True,type=str)
parser.add_argument('-out_prefix','--out_prefix',help='Output prefix',required=True,type=str)

args = parser.parse_args()
zargs = dir(args)
zargs = [zarg for zarg in zargs if zarg[0]!='_']
for zarg in zargs:
    print(f'{zarg} is ',end='')
    exec(f'{zarg}=args.{zarg}')
    exec(f'print(args.{zarg})')

bad_chars = ['.','-','N','W','R','Y']

print(f'Loading human reference fasta')
GRCh38_fasta = '/home/tc557/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/HGDP/GRCh38_ref/Homo_sapiens_assembly38.fasta'
with open(GRCh38_fasta,'r') as f:
    lines = f.readlines()
# get start and stop line of chrom of interest
print(f'Getting chrom={chrom} of interest...')
if chrom==22:
    nextchrom='X'
else:
    nextchrom=chrom+1
for i,h in enumerate(lines):
    if f'chr{chrom} ' in h:
        start_line_index = i
    elif f'chr{nextchrom} ' in h:
        end_line_index = i
        break
# write string of fasta
print(f'Writing human reference fasta to string...')
human_fasta = ''
for i in range(start_line_index+1,end_line_index):
    human_fasta += lines[i].split('\n')[0]
# convert fasta to list for easy manipulation
humanref_list = [i for i in human_fasta]
badhuman_indices = np.array([i for i,j in enumerate(humanref_list) if j in bad_chars ])
del human_fasta

dict_letters = {'A':'A','C':'C','G':'G','T':'T','a':'A','c':'C','g':'G','t':'T','N':'N','-':'N','n':'N','W':'N','R':'N','Y':'N','M':'N','K':'N','S':'N'}
def capitalise_ref(i,dict_letters):
    return dict_letters[i]


# pop='YRI'
print(f'Loading human allele counts')
filename = f'/home/tc557/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/1000Genomes_30X/phased_vcf/SFS_241003/processed_241003/pop{pop}_chrom{chrom}_CHROM_POS_REF_ALT_AC.txt.gz'
columns = ['CHROM','POS','REF','ALT','AC']
stats = pd.read_csv(filename,names=columns,sep=' ')
stats['POS'] = stats['POS']-1 # make zero based

positions = np.array(stats['POS'])
reference = list(stats['REF'])
alt = list(stats['ALT'])
allele_count = np.array(stats['AC'])
zlength = positions[-1]+1
seq_ac = np.zeros(zlength)
seq_ac[positions] = allele_count
num_derived = allele_count.max()

print(f'Writing human reference sequence...')
seq_ref = ['N']*zlength
for i,j in enumerate(positions):
    seq_ref[j] = reference[i]
abba = np.array([i for i,j in enumerate(positions) if seq_ref[j]!=humanref_list[j]]) # check human ref matches ref fasta
if len(abba)>0:
    print(f'ERROR! The humanref does not always equal the ref from vcf, which is problematic. Aborting.')
    sys.exit()
seq_ref = humanref_list

print(f'Writing human alternate sequence...')
seq_alt = ['N']*zlength
for i,j in enumerate(positions):
    seq_alt[j] = alt[i]



# abba = [i for i,j in enumerate(positions) if seq_ref[j]!=chimpref_list[j] and chimpref_list[j] not in bad_chars and seq_ref not in bad_chars] # where does human ref not match chimp
# get B sequence
print(f'Loading B sequence...',flush=True)
B_sequence = get_B_sequence(f'/home/tc557/Murphy_Bstat_hg20/chr{chrom}_nomissing_hg38.bed')

print(f'Loading mappability mask (assuming given regions are callable)...')
# get callable sequence
print(f'Loading mappability mask...',flush=True)
uncallable = pd.read_csv(f'/home/tc557/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/HGDP/GRCh38_ref/strict_perchrom/chr{chrom}_negative.bed',header=None,delimiter='\t')
uncallable = np.array(uncallable.iloc[:,1:3])
seq_mask = np.ones(len(B_sequence))
for jj in uncallable:
    zjstart = int(jj[0])
    zjend = int(jj[1])
    seq_mask[zjstart:zjend] = -1


num_masked_bases = len(seq_mask[seq_mask==-1])
print(f'number of masked bases = {num_masked_bases} out of {len(seq_mask)} total bases')

ndxx = []
minlength = zlength
print(f'Writing site counts...')
# for i in np.arange(0,minlength,1):
cc=0
zlen = len(seq_ref)
for i in range(0,len(seq_ref)):
    if cc%1000000==0:
        print(f'i={cc} out of {zlen}')
    if seq_mask[i]==-1 or B_sequence[i]<0.8 or seq_ref[i]=='N':
        cc+=1
        continue
    else:
        ndxx.append(seq_ac[i])
    cc+=1

filename = f'{out_prefix}_allfiles.pickle'
# with open(filename,'wb') as f:
    # pickle.dump(files_dict,f)

xdom = np.arange(0,num_derived+2,1)
for sitetype in ['ndxx']: # ,'nd00','nd11','nd01','nd10','ndx1','ndx0','nd1x','nd0x']:
    exec(f'zcounts_array = {sitetype}')
    zcounts, x_ = np.histogram(zcounts_array,bins=xdom)
    filename = f'{out_prefix}_{sitetype}.txt.gz'
    np.savetxt(filename,zcounts,fmt='%i')
    print(f'Saved counts to {filename}')