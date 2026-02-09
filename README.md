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

#Find the mRNA coordinates from these mhc_gene_ids from the gff file.
grep -f mhc_gene_ids KW.gene.gff | awk '$3=="mRNA"' | cut -f1,4,5,9 > mhc_coords.tmp
```

### 2. Format `mhc_coords.temp` into a bed file format
VCFtools needs (car, start, end) format
Bed uses 0-based starts, so this will also subtract 1 from GFF start coordinates
```
awk -v OFS="\t" '{print $1, $2-1, $3, $4}' mhc_coords.tmp > mhc_targets.bed
```

### 3. Finally, I will use bcftools to extract SNPs that fall within these MHC regions

Isolate only vcf’s from SRKWs, and then pull SNPs from `mhc_targets.bed`.
```
bcftools view   -S <(bcftools query -l vcf_files | grep -E '^[JKL]')   vcf_files   -o vcfs_JKL.vcf
```
```
bcftools index vcfs_JKL.vcf.gz #index vcf
bcftools view -R mhc_targets.bed vcfs_JKL.vcf.gz -o mhc_variants.vcf
```

## Nucleotide diversity

Create separate files for J, K, and L matrilines
```
bcftools view   -S <(bcftools query -l mhc_variants.vcf | grep -E '^[J]')   mhc_variants.vcf   -o mhc_J.vcf
bcftools view   -S <(bcftools query -l mhc_variants.vcf | grep -E '^[K]')   mhc_variants.vcf   -o mhc_K.vcf
bcftools view   -S <(bcftools query -l mhc_variants.vcf | grep -E '^[L]')   mhc_variants.vcf   -o mhc_L.vcf
```

Get nucleotide diversity for each matrilineal family
```
vcftools --vcf mhc_J.vcf --site-pi --out J_pi
vcftools --vcf mhc_K.vcf --site-pi --out K_pi
vcftools --vcf mhc_L.vcf --site-pi --out L_pi
```

I then created `pos2gene.py` to take the outputs from `J_pi`, `K_pi`, and `L_pi` and figure out which annotated MHC genes the site positions are in. I ignored multiple isoforms because we cant actually tell which isoform is correct for each individual, so they are just duplicated scores.

# Move to data analysis in `mhc_diversity.R`
