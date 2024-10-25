# take the processed allele counts and write the SFS
for pop in YRI;
        do
        for chrom in {1..22}
        # for chrom in 2;
                do
                # human_ref=chimp
                outfile=/home/tc557/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/1000Genomes_30X/phased_vcf/SFS_241003/processed_241011/nopolar_pop${pop}_chr${chrom}_ndxx.txt.gz
                commandline="python /home/tc557/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/1000Genomes_30X/phased_vcf/SFS_241003/write_SFS_nopolar_241011.py -chrom ${chrom} -pop ${pop} -out_prefix /home/tc557/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/1000Genomes_30X/phased_vcf/SFS_241003/processed_241011/nopolar_pop${pop}_chr${chrom}"
        echo "no polarisation"
                if [ -f $outfile ]; then
                        echo $outfile exists
                else
                        echo $outfile DOES NOT exist
                        echo $commandline
                        # sbatch --wrap="${commandline}" -A DURBIN-SL2-CPU -p icelake-himem --mem=50G -t 01:00:00 -c 1
                fi
                echo
        done
done