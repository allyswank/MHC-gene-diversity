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

For heterozygous sites, I will filter for allelic imbalance where one allele represents <10% or >90% of reads. I am filtering AB per genotype rather than per site so that I don't lose good data across all samples just because a site is bad in one sample.

```
# Calculate AB with gatk (AF is the same as AB here)
gatk VariantAnnotator -R genome.fa -V mhc_depth_filtered.vcf.gz -O mhc_depth_filtered_AB.vcf.gz -A AlleleFraction

# Use bcftools to mask heterozygotes with allelic imbalance
bcftools view -i 'GT="0/1" && (FMT/AF>0.1 || FMT/AF<0.9)' mhc_depth_filtered_AB.vcf.gz > mhc_ab_filtered.vcf

# Look at per site missingness
bcftools +fill-tags mhc_ab_filtered.vcf.gz -- -t F_MISSING | bcftools view -O z -o mhc_ab_filtered_tags.vcf.gz
bcftools query -f '%CHROM\t%POS\t%F_MISSING\n' mhc_ab_filtered_tags.vcf.gz
    # Only 2 sites with 12% missingness - all other sites <2% missingness.
```



## Genetic differentiation between samples
<img width="1401" height="737" alt="image" src="https://github.com/user-attachments/assets/c3880387-e894-460a-b39a-c368f0d6ff64" />
This is using the mhc_ab_filtered.vcf.gz file.


## Haplotype inference

First I will need to take my multi-sample vcf and create one diploid consensus fasta per gene per individual. I will then merge all exons within class II loci into one sequence per sample and make a fasta file for each gene where each fasta sequence is from one individual. At this point, the data are still unphased and heterozygous sites are still unresolved. 

```
# make a gatk genome dict file and index the vcf
gatk CreateSequenceDictionary -R genome.fa -O genome.dict
gatk IndexFeatureFile --input mhc_ab_filtered.vcf.gz

samples=$(bcftools query -l mhc_ab_filtered.vcf.gz) #this will loop through each sample

for sample in $samples; do
    gatk SelectVariants -V mhc_ab_filtered.vcf.gz -sn $sample -O ${sample}.vcf.gz
    gatk FastaAlternateReferenceMaker -R genome.fa -V ${sample}.vcf.gz -O ${sample}.fa -L mhc_exons.bed --use-iupac-sample $sample
done

python3 make_loci_fastas

# Fix headers for PHASE
sed -i '/^>/s/ /_/g' DOB.fa 
sed -i '/^>/s/.fa/_/g' DOB.fa
sed -i '/^[^>]/s/\./N/g' DOB.fa # change "." to "N" in sequences
```

Do a quick check to make sure the fasta's look right for each loci (Ex.: `seqkit stats DRA.fa`, all loci should be same number of bp's), and then run PHASE to resolve heterozygote allele frequencies. 

```
# convert these fasta alignments to PHASE format using SeqPHASE
for gene in `ls *.fa`; do
    python3 ../SeqPHASE/seqphase1.py -1 ${gene} -o ${gene}
done

# run PHASE
for gene in `ls *.fa`; do
    ../phase/PHASE ${gene}.inp ${gene}.out 100 1 100 -X10 -MR
done

# convert phased output back to fasta format
for gene in `ls *.fa`; do
    python3 ../SeqPHASE/seqphase2.py -i ${gene}.out -c ${gene}.const -o ${gene}_phased.fa
done
```

Now that phasing is complete, we can do final filtering and begin data analysis. I will start by removing unresolved alleles and those with a frequency less than 2% of the dataset.

```
# find any unresolved alleles
for gene in `ls *_phased.fa`; do
  seqkit grep -s -r -p "[NKMSWBDHV.]" ${gene} > ${gene}_unresolved.fa
done

# keep only alleles that are fully resolved for the rest of analysis (keeping R [Purine (A or G)] and Y [Pyrimidine (C or T)])
for gene in `ls *_phased.fa`; do
  seqkit grep -vs -r -p "[NKMSWBDHV.]" ${gene} > ${gene}_clean.fa
done

# get rid of alleles if the individuals had one but not both of their alleles filtered
for gene in *_clean.fa; do
  grep "^>" "$gene" | sed 's/^>_//' | cut -d'_' -f1 | sort | uniq -c | awk '$1==2 {print $2}' > keep_ids
  awk '{print "^_"$1"_"}' keep_ids > keep_patterns #make sure id will match fasta header
  seqkit grep -n -r -f keep_patterns "$gene" > "${gene}_paired.fa"
done

#I also renamed these files so that they are just (gene-name)_paired.fa 
for f in *_paired.fa; do
    base="${f%%.*}"
    mv "$f" "${base}_paired.fa"
done

# identify unique alleles and count how many times we see each one
for gene in `ls *_paired.fa`; do
  seqkit rmdup -s -D ${gene}_dup_counts ${gene} > ${gene}_unique.fa
done
#unique.fa = unique alleles
#dup_counts = counts per allele

# sum total alleles
for gene in *_paired.fa; do
  total=$(awk '{sum+=$1} END {print sum}' "${gene}_dup_counts")
  echo "$gene total alleles: $total"
done

# calculate frequency per allele
for gene in *_paired.fa; do
  awk '{freq=$1/'"$total"'; print $1, $2, freq}' "${gene}_dup_counts" > "${gene}_freqs"
done

# check remaining number of alleles for each gene
for f in *_freqs; do
  echo "$f alleles: $(wc -l < "$f")"
done

# check that all individuals only have a/b alleles (diploid)
for f in *_paired.fa; do
  echo $f
  grep "^>" $f | sed 's/^>_//' | cut -d'_' -f1 | sort | uniq -c | awk '$1!=2'
done
```
Haplotype network analysis in `mhc_haplotypes.R`!
<img width="970" height="649" alt="Screenshot 2026-04-06 at 11 30 13 PM" src="https://github.com/user-attachments/assets/f86676ef-495b-4a5a-b851-8f9492823e1b" />


## Identify non-synonomous mutations
I want to do this after phasing. Can then evaluate how specific regions might carry functional changes in a population.


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
<img width="760" height="588" alt="image" src="https://github.com/user-attachments/assets/aed8eaac-bf47-4ca1-8299-5314c048fdb1" />

### Plots for each population's π by gene:
<img width="760" height="588" alt="image" src="https://github.com/user-attachments/assets/658e6adf-8445-4c19-91ea-f4c1a46a8949" />

## Focusing only on MHC class II
### Plot of mean π for each population: 
<img width="760" height="588" alt="image" src="https://github.com/user-attachments/assets/81084078-c8fd-4ead-94e5-ad10707d3463" />
<img width="760" height="588" alt="image" src="https://github.com/user-attachments/assets/9b0b60e5-106b-4837-bfd6-a0e005bbb21b" />






