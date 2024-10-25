# infer a demographic model with mushi

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mushi
import argparse
import pdb
import pickle
import sys

# usage 
# python /home/tc557/falsifying_bottleneck/infer_mushi_model_241011.py -pop YRI -out_prefix /home/tc557/testingmushidelete_ -trend_1 0 -trend_2 10 -ridge 750 -mu 1.25e-08 -pts 200 -most_ancient_gens 2e+5 -maxiter 300 -folded False

parser = argparse.ArgumentParser(description="Set inputs and outputs")
parser.add_argument('-pop','--pop',help='1KGP population',required=True,type=str)
parser.add_argument('-out_prefix','--out_prefix',help='Output prefix',required=True,type=str)
parser.add_argument('-ridge','--ridge',help='ridge penalty',required=True,type=float)
parser.add_argument('-trend_1','--trend_1',help='trend_1',required=int,type=float)
parser.add_argument('-trend_2','--trend_2',help='trend_2',required=True,type=float)
parser.add_argument('-mu','--mu',help='mutation rate per gen per bp',required=False,type=float,default=1.25e-08)
parser.add_argument('-pts','--pts',help='Pts used for inference',required=False,type=int,default=200)
parser.add_argument('-most_ancient_gens','--most_ancient_gens',help='Most ancient number of generations',required=False,type=float,default=2e+05)
parser.add_argument('-maxiter','--maxiter',help='Maximum number of iterations',required=False,type=int,default=300)
parser.add_argument('-folded','--folded',help='Use folded or unfoled SFS',required=False,type=str,default="False")


args = parser.parse_args()
zargs = dir(args)
zargs = [zarg for zarg in zargs if zarg[0]!='_']
for zarg in zargs:
    print(f'{zarg} is ',end='')
    exec(f'{zarg}=args.{zarg}')
    exec(f'print(args.{zarg})')

zpop=pop
if folded=="False":
    folded=False
elif folded=="True":
    folded=True
else:
    print(f'Invalid folded parameter. Must be True or False')
    sys.exit()

SFSfile = f'/home/tc557/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/1000Genomes_30X/phased_vcf/SFS_241003/processed_241011/nopolar_pop{zpop}_allchrs_ndxx.txt.gz'
observed_sfs = np.loadtxt(SFSfile)
observed_sfs_one_to_nminusone = observed_sfs[1:-1]
# plt.plot(sfs_one_to_nminusone)
ksfs = mushi.kSFS(observed_sfs_one_to_nminusone)
# print(f'sfs.sum() ={sfs.sum()}')
ksfs.infer_eta(mu0=observed_sfs.sum()*mu, trend_kwargs=(trend_1, trend_2), ridge_penalty = ridge, pts=pts,ta=most_ancient_gens,folded=folded, max_iter=maxiter, verbose=True)

gen=30
_observed_data, expected_sfs, expected_sfs_lower, expected_sfs_upper, ancestral_state_misidentification_rate = ksfs.return_statistics(folded=folded)

mushi_model = np.column_stack([ksfs.eta.change_points*gen, ksfs.eta.y[0:-1]/2])

parameters = {
'mushi_model':mushi_model,
'observed_SFS':observed_sfs_one_to_nminusone, # this should be the same as below if folded=False, otherwise will be collapsed
'observed_SFS_mushi':_observed_data,
'expected_sfs':expected_sfs,
'expected_sfs_lower':expected_sfs_lower,
'expected_sfs_upper':expected_sfs_upper,
'ancestral_state_misidentification_rate':ancestral_state_misidentification_rate # this will be None if folded=True
}
filename = f'{out_prefix}model_parameters.pickle'
with open(filename,'wb') as f:
    pickle.dump(parameters,f)
print(f'Saved parameters to {filename}')

plt.figure()
plt.plot(range(1,len(_observed_data)+1),_observed_data,label="Observed",linewidth=4,marker='o',markersize=10)
plt.plot(range(1,len(_observed_data)+1),expected_sfs,label="Expected",linewidth=4,marker='o',markersize=10,alpha=0.4)
plt.xlabel('Number derived alleles')
plt.ylabel('Frequency')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.tick_params(which='major',length=20)
plt.tick_params(which='minor',length=10)
ztitle = f'{zpop}; {trend_1}; {trend_2}; {ridge}; {folded}'
plt.title(ztitle)
filename = f'{out_prefix}SFS_expected_observed.png'
plt.savefig(filename)
print(f'Saved SFS figure to {filename}')

# mushi model plot
plt.figure()
plt.plot(mushi_model[:,0],mushi_model[:,1],label=zpop,linewidth=4)


# psmc model plots, see http://localhost:1231/notebooks/twin_peaks/twinpeaks_240813/plots_for_paper_240813.ipynb
mu=1.25e-08
final_params_file_pan = '/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813/D_64/b_100/spread1_0.01/spread2_50/muoverr_1.5/iterations30/thresh_1/thetafixed_0.0008/popsample_YRI_NA18488/final_parameters.txt'
final_params_pan = np.loadtxt(final_params_file_pan)
time_array = list(final_params_pan[:,1])
time_array.insert(0,0)
time_array = np.array(time_array)
plt.stairs(edges=(time_array/mu)*gen,values=(1/final_params_pan[:,2])/mu,label=zpop,linewidth=4,linestyle="solid",baseline=None,color="red")

plt.xscale('log')
plt.xlabel('Time')
plt.ylabel('$N(t)$')
plt.legend()
plt.tick_params(which='major',length=20)
plt.tick_params(which='minor',length=10)
plt.title(ztitle)
filename = f'{out_prefix}inferred_model.png'
plt.savefig(filename)
print(f'Saved N(t) figure to {filename}')