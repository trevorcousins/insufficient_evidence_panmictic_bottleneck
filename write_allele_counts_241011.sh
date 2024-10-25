# process allele counts 
for pop in YRI;
        do
        for chrom in {1..22}
                do
                commandline="bcftools view /home/tc557/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/1000Genomes_30X/phased_vcf/1kGP_high_coverage_Illumina.chr${chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz -s `cat /home/tc557/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/1000Genomes_30X/metadata/sample_pop.txt | grep \${pop} |  awk '{print $1}' | paste -sd ',' -` | bcftools query -f \"%CHROM %POS %REF %ALT %AC\n\" | awk '(length(\$3)<2)' | awk '(length(\$4)<2)' | gzip -c > /home/tc557/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/1000Genomes_30X/phased_vcf/SFS_241003/processed_241003/241004_pop${pop}_chrom${chrom}_CHROM_POS_REF_ALT_AC.txt.gz"
                echo $commandline
                sbatch --wrap="${commandline}" -A DURBIN-SL2-CPU -p cclake --mem=50G -t 03:00:00 -c 1
                echo
        done
done