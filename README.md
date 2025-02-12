# insufficient_evidence_panmictic_bottleneck

A recent paper in Science proposed that humans went through a very strong population bottleneck around 1 million years ago. We do not believe the data is in support of this. This repo describes some analysis we did, demonstrating that simpler panmictic models without the bottleneck better fit the data. 

## Real data processing

In the paper, there are three different data sets that we looked at. All of these are based on data fom the 1000 Genomes Project (1KGP), but they have different processing. We call these three different processings "GRCh38, neutral, unpolarized", "GRCh37, neutral, unpolarized", and "GRCh37, noncoding, polarized”. 

### GRCh38, neutral, unpolarized

Download 1KGP VCFs from https://www.internationalgenome.org/category/vcf/ . We used a B-map from Murphy et al (eLife, 2023), which was lifted over to GRCh38, and only kept regions with b > 0.8. We filtered out regions according to a strict mappability mask, as downloaded from ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38 .<br>
Write allele counts at each positions with write_allele_counts_241011.sh<br>
(we did not polarise the alleles using a chimp or gorilla sequence, which is clearly the wrong thing to do, however mushi can model the rate of mididentification very well. Note that standard annotations of the ancestral allele based on multiple-sequence alignment leave a lot to be desired, see Sup Fig 2)<br>
Write the SFS with write_SFS_submission_nopolarisation_241011.sh which calls write_SFS_nopolar_241011.py (this includes a strict mappabiltiy mask, and a B-map from Murphy et al., 2023)<br>

### GRCh37, neutral, unpolarized

This processing is analogous to that in "GRCh38, neutral, unpolarized", except we used the GRCh37 VCFs as downloaded from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ . The mappability mask was downloaded from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks .

### GRCh37, neutral, unpolarized

Annotations (for exonic regions etc) were downloaded from https://www.gencodegenes.org/human/release_35lift37.html .<br>
Coding regions (as described by HAVANA or ENSEMBL) were obtained with :
```
zcat < gencode.v35lift37.annotation.gff3.gz | awk '$3 == "CDS" {print}' | grep -v "pseudogene" | grep "HAVANA" | bedtools sort -i stdin -g human.hg19.genome | gzip > HAVANA.gff.gz
zcat < gencode.v35lift37.annotation.gff3.gz | awk '$3 == "CDS" {print}' | grep -v "pseudogene" | grep "ENSEMBL" | bedtools sort -i stdin -g human.hg19.genome | gzip > ENSEMBL.gff.gz

bedtools complement -i ENSEMBL.gff.gz -g human.hg19.genome > hg19.ENSEMBL.complement.bed
bedtools complement -i HAVANA.gff.gz -g human.hg19.genome > hg19.HAVANA.complement.bed
```
and subsequently masked from analysis. We used the ancestral allele annotation from the 2015 1KGP paper, which was obtained from the VCF with : 
```
bcftools query -f "%CHROM %POS %REF %ALT %AA\n" ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | awk '(length($3)<2)' | awk '(length($4)<2)' | awk '{$5 = substr($5, 1, length($5)-3); print}'| gzip -c > 241004_chrom${chrom}_CHROM_POS_REF_ALT_AA.txt.gz
```

## Running mushi

Using each of the 3 data sets, we inferred a model with mushi using searching_mushi_model_241011.sh, which calls infer_mushi_model_241011.py<br>

## Running FitCoal

Using each of the 3 data sets, we inferred a model with FitCoal using infer_FitCoal_241011.sh<br>

## Analysis 

Various other bits of analysis are done in plots_for_paper_241025_upload.ipynb<br>

## Preprint

Our preprint is available to read at https://www.biorxiv.org/content/10.1101/2024.10.21.619456v1
