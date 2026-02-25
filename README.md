# Pipeline for SRKW analysis of MHC gene diversity
[will add more text here later] 


## Pulling MHC Genes
Extract all CDS regions from the genome so that I can annotate them.
- Run the `cds_maker.py` python script that will pull the gene regions from the genome based on [Marty’s KW.gene.gff file](https://www.nature.com/articles/s41559-023-01995-0)

Translate the extracted CDS regions into proteins for BLAST. You could blast this whole `KW_pep_out.fa` file for complete annotations if you want!
```
seqkit translate KW_cds_out.fa > KW_pep_out.fa 
```

I am just going to use the annotations on NCBI of [MHC genes from another annotation](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9733/), extract them from my `KW_pep_out.fa`, and then use their positions to find the MHC genes out of `vcfs_JKL.vcf`
- There were 17 genes containing “Histocompatibility” in gene details, so just pulled them all as proteins
```
module load blast+/2.12.0
makeblastdb -in mhc_prot_ncbi.faa -dbtype prot -out mhc_db
blastp -query KW_pep_out.fa -db mhc_db -out mhc_hits -outfmt 6

#Filter mhc_hits for >40% identity match and e-value >1e-5
awk '$3 >= 40 && $11 < 1e-5' mhc_hits > mhc_hits_filtered
```
Okay now I have this `mhc_hits_filtered file`, which contains positions of MHC genes in `KW.gene.gff`. I can use this to find the MHC genes in the vcf files!

### 1. Extract MHC coordinates from KW.gene.gff
We want the mRNA line for each hit so that we have the full gene start/end ranges.

```
#First finding unique gene IDs from the `mhc_hits_filtered` file.
cut -f1 mhc_hits_filtered | sort -u > mhc_gene_ids

#Find the CDS coordinates from these mhc_gene_ids from the gff file.
grep -f mhc_gene_ids KW.gene.gff | awk '$3=="CDS"' > mhc_cds_features.gff
```

### 2. Format `mhc_coords.temp` into a bed file format
VCFtools needs (car, start, end) format
Bed uses 0-based starts, so this will also subtract 1 from GFF start coordinates
```
awk -v OFS="\t" '{print $1, $4-1, $5, $9}' mhc_cds_features.gff > mhc_exons.bed
```

### 3. Finally, I will use bcftools to extract SNPs that fall within these MHC regions

Pull SNPs in the `kw_151.snp.final.vcf.gz` from `mhc_targets.bed` for MHC gene diversity analysis

```
bcftools index kw_151.snp.final.vcf.bgzf.gz #index vcf
bcftools view -R mhc_exons.bed kw_151.snp.final.vcf.bgzf.gz -Oz -o mhc_exon_variants.vcf.gz
bcftools index mhc_exon_variants.vcf.gz #index output
```
### 4. Filter the SNPs 

I only want to include biallelic sites with high quality scores. This decreases noise downstream and seems to be appropriate according to findings in other cetaceans (Heimeier et al.). 

```
view -m2 -M2 -v snps mhc_exon_variants.vcf.gz -Oz -o mhc_exon_biallelic_snps.vcf.gz # get biallelic

bcftools filter -e 'QUAL<80 || MQ<40 || MQRankSum>12.5 || MQRankSum<-12.5 || ReadPosRankSum>8 || ReadPosRankSum<-8' \
  mhc_exon_biallelic_snps.vcf.gz -Oz -o mhc_exon_filtered.vcf.gz
bcftools view -f PASS mhc_exon_filtered.vcf.gz -Oz -o mhc_exon_filtered_PASS.vcf.gz #keep only passing sites
bcftools index mhc_exon_filtered_PASS.vcf.gz

# Sanity check: get number of sites
bcftools view -H mhc_exon_filtered_PASS.vcf.gz | wc -l 
#257 - only lost 6 sites from mhc_exon_biallelic_snps.vcf
```

I also want to calculate depth for each genotype and remove sites with low sequencing depth (<5x). 

```
bcftools filter -e 'MEAN(FORMAT/DP) < 5' mhc_exon_filtered_PASS.vcf.gz -Oz -o mhc_depth_filtered.vcf.gz
bcftools index mhc_depth_filtered.vcf.gz

bcftools view -H mhc_depth_filtered.vcf.gz | wc -l
#254 - only lost 3 with lower than 5x coverage
```

For heterozygous sites, I will filter for allelic imbalance where one allele represents <10% or >90% of reads.

First, export per-site AD (ref/alt counts) for each sample with vcftools.

```
vcftools --gzvcf mhc_depth_filtered.vcf.gz --get-INFO AD --out mhc_AD
```

Now vcftools can calculate allelic balance per site or per genotype. Per site is a lot more strict (one bad samples ruins the site for everyone) and per genotype is looser (only removes that site for the one sample that is imbalanced). I will try both strict and loose. 

Starting with filtering per site:


## Genetic differentiation between samples
<img width="761" height="436" alt="image" src="https://github.com/user-attachments/assets/e9849c97-8088-4331-b44a-6f79fe316d9f" />


## Nucleotide diversity

Get nucleotide diversity using vcftools across each population to be analyzed separately in R
```
bcftools view   -S <(bcftools query -l mhc_depth_filtered.vcf.gz | grep -E '^[J|K|L|calf|neo]')   mhc_depth_filtered.vcf.gz   -o mhc_srkw.vcf #SRKW - N = 
bcftools view   -S <(bcftools query -l mhc_depth_filtered.vcf.gz | grep -E '^[T|MBH]')   mhc_depth_filtered.vcf.gz   -o mhc_TKW.vcf #Transients - N = 
bcftools view   -S <(bcftools query -l mhc_depth_filtered.vcf.gz | grep -E '^[P|AKW13]')   mhc_depth_filtered.vcf.gz   -o mhc_AKW.vcf #ARKW - N = 
bcftools view   -S <(bcftools query -l mhc_depth_filtered.vcf.gz | grep -E '^[o|AWK04|NG]')   mhc_depth_filtered.vcf.gz   -o mhc_off.vcf #offshore - N = 


vcftools --vcf mhc_srkw.vcf --site-pi --out srkw_pi
vcftools --vcf mhc_TKW.vcf --site-pi --out TKW_pi
vcftools --vcf mhc_AKW.vcf --site-pi --out AKW_pi
vcftools --vcf mhc_off.vcf --site-pi --out off_pi
```

I then created `pos2gene.py` to take the output from `#_pi` and figure out which annotated MHC genes the site positions are in. I ignored multiple isoforms because we cant actually tell which isoform is correct for each individual, so they are just duplicated scores.

# Move to data analysis in `mhc_diversity.R`

### Plots for mean π across all genes by population:
<img width="900" height="900" alt="image" src="https://github.com/user-attachments/assets/7f0dd1e9-dd2f-40d4-8823-1eaee9912f94" />


### Plots for each population's π by gene:
<img width="1520" height="1512" alt="image" src="https://github.com/user-attachments/assets/6d1ca385-dc58-4384-9601-2f978aeb62d1" />


### Plot of rolling average π: 
<img width="2706" height="1512" alt="image" src="https://github.com/user-attachments/assets/687a6674-dfb8-4b5d-88d1-e93d864a615e" />

## Focusing only on MHC class II
### Plot of rolling average π:
<img width="2764" height="1010" alt="image" src="https://github.com/user-attachments/assets/3840e938-44f3-4206-aa84-19631185f34c" />

### Plot of mean π for each population: 
<img width="1252" height="1010" alt="image" src="https://github.com/user-attachments/assets/3f99ae1d-a27f-4226-ada3-487fa76bd774" />






