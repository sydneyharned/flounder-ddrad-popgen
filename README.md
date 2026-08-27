# flounder-ddrad-popgen
ddRAD-seq population genomics of southern flounder (Paralichthys lethostigma) across NC and TX estuaries. Stacks · PLINK · ADMIXTURE · hierfstat · R

# Southern Flounder ddRAD Population Structure: North Carolina & Texas

ddRAD-seq based population genomics of juvenile southern flounder (*Paralichthys lethostigma*) across estuarine nursery habitats in North Carolina and Texas, assessing spatial and temporal genetic structure relevant to fisheries management.

> Harned, S., Burford Reiskind, M. et al. (*in review*). An assessment of spatial and temporal population structure in southern flounder (*Paralichthys lethostigma*) in North Carolina. Preprint at https://doi.org/10.64898/2026.04.24.720543

---

## Background

Southern flounder are a commercially important and declining estuarine species managed under active fisheries review in the southeastern U.S. Juveniles spend 1–2 years in shallow estuarine nurseries before moving offshore. This study uses ddRAD-seq to characterize fine-scale population structure across coastal nursery habitats in North Carolina (sampled across multiple years: 2014–2016, 2022, 2023) and Texas (2014–2016), and assesses temporal stability of genetic structure within sites sampled across years.

---

## Sampling Design

### North Carolina

Juvenile flounder were collected by NC Division of Marine Fisheries (NCDMF) via the P120 estuarine trawl survey (3.2 m headrope, 6.4 mm bar mesh body, 3.2 mm cod end; 1 min tow at 1.1 m/s). In 2016, additional sampling used an otter trawl or 2-meter beam trawl. Sampling years: 2014, 2015, 2016, 2022, 2023.

Sampling locations were grouped into four geographic regions:

| Region | Code | Environmental character |
|--------|------|------------------------|
| Albemarle Sound | AL | Low temp, high freshwater variation, single narrow inlet (Oregon Inlet) |
| Pamlico Sound | PS | Low, stable salinity/temperature; enclosed; limited tidal exchange |
| Neuse River | NR | Intermediate; higher flushing rate; more variable than Pamlico |
| New River (south) | NR | Higher salinity/temperature; greater tidal input |

Sites within each region varied by year. Full sampling schema by year in supplemental Table S1.

### Texas

Juvenile flounder collected by Texas Parks and Wildlife (2014–2016) from three wild locations and one hatchery:

| Site | Method |
|------|--------|
| East Bay | 1- or 2-meter beam trawl (1–2 min at 1–2 knots); hand-pulled in restricted areas |
| Bastrop Bayou | Same as above |
| Chocolate Bayou | Same as above |
| Sea Center Texas Hatchery (Lake Jackson) | 2016 only — hatchery vs. wild comparison |

---

## Pipeline Overview

```
Fin clips → DNA extraction (Qiagen DNeasy)
        │
        ▼
ddRAD library preparation (SphI + MluCI)
        │
        ▼
Illumina NovaSeq (single-end, 90 bp)
        │
        ▼
1. Quality control                    [FastQC]
        │
        ▼
2. Demultiplexing & SNP calling       [Stacks v1.24: process_radtags, denovo.pl]
        │
        ▼
3. Species verification               [Treemix v1.13]
        │
        ▼
4. Filtering                          [Stacks populations, PLINK v1.19, DartR (R)]
        │
        ├── TX vs. NC structure       [hierfstat: pairwise FST]
        │
        ├── Texas analysis            [hierfstat, PCA: adegenet + ggplot2]
        │
        └── North Carolina analysis   [hierfstat, PCA: adegenet, ADMIXTURE v1.3]
```

---

## Methods

### 1. Quality control

```bash
fastqc raw_reads/*.fastq.gz -o fastqc_output/
# Criterion: phred score > 33
```

### 2. Demultiplexing and SNP calling

```bash
# Demultiplex variable-length barcodes
process_radtags -p ./raw/ -o ./demultiplexed/ \
    -b barcodes.txt \
    --inline_null \
    -r -c -q

# Trim reads to uniform length
# (trimmed to 90 bp prior to de novo assembly)

# De novo SNP calling
denovo_map.pl -T 8 -o ./stacks_output/ \
    --samples ./demultiplexed/ \
    --popmap popmap.txt \
    -m 3 -M 2 -n 2
# m=3: minimum stack depth
# M=2: mismatches allowed between loci within individual
# n=2: mismatches allowed between loci in catalog
```

Total loci identified across NC + TX: 3,643,301; mean per-sample coverage: 22.2×

### 3. Species verification

```bash
# Treemix v1.13 — maximum-joining trees per sampling year
# Bootstrap: 1000 replicates
# Purpose: confirm all samples are P. lethostigma
# (4 P. dentatus individuals identified in 2023 and removed)
treemix -i input.treemix.gz -o treemix_output -bootstrap -k 500 -noss
```

### 4. Dataset filtering

Two Stacks populations runs were used depending on analysis:

```bash
# Broad dataset (TX vs. NC, Texas analyses): loci in ≥1 population
populations -P ./stacks_output/ -M popmap.txt \
    -p 1 -r 0.6 \
    --write-random-snp \
    --max-obs-het 0.8 \
    --plink --vcf

# NC-only dataset: stricter locus sharing (loci in ≥2 populations)
populations -P ./stacks_output/ -M popmap_NC.txt \
    -p 2 -r 0.6 \
    --write-random-snp \
    --max-obs-het 0.8 \
    --plink --vcf
# Final NC dataset: 25,085 loci
```

**Key filtering parameter decisions:**

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `-r` (min within-group genotyping rate) | 0.6 | Balances locus retention with completeness |
| `--write-random-snp` | 1 SNP per locus | Avoids pseudoreplication within RAD loci |
| `--max-obs-het` | 0.8 | Removes likely paralogous loci |
| `-p` (populations sharing) | 1 (broad) / 2 (NC) | Stricter for NC to remove rare, location-specific variants |

```bash
# PLINK filtering
plink --file populations.plink \
    --maf 0.01 \
    --geno 0.5 \
    --mind 0.3 \       # adjusted to 0.5 depending on dataset average missingness
    --make-bed \
    --out filtered_snps

# Format conversion for downstream analyses
# PGDSpider v2.1.1.0: PLINK → genetix format

# HWE filtering in R
library(dartR)
gl <- gl.filter.hwe(gl, method = "out-all")
# "out-all": removes loci deviating from HWE across all populations
```

Multiple SNP datasets were generated optimized for specific analyses (see sections 2.5–2.7 in manuscript).

### 5. TX vs. NC structure

```r
library(hierfstat)

# Pairwise FST per sampling year (2014, 2015, 2016)
# Individuals grouped by state and year (e.g., NC_2014, TX_2015)
fst <- pairwise.WCfst(dat)   # Weir & Cockerham 1984

# 95% confidence intervals via locus bootstrapping
ci <- boot.ppfst(dat, nboot = 1000, quant = c(0.025, 0.975))
```

### 6. Texas analysis

```r

library(hierfstat); library(adegenet); library(ggplot2)

# Summary statistics per year
stats <- basic.stats(dat)   # HO, HE, FIS

# PCA (no a priori group assignment required)
# Retain first 2 PCs (greatest variance explained; no further structure in subsequent PCs)
pca <- dudi.pca(tab(genind_obj, NA.method = "mean"), scannf = FALSE, nf = 2)

# Temporal analysis: Bastrop Bayou (sampled 2014, 2015, 2016)
# Pairwise FST + PCA as above
```

### 7. North Carolina analysis

```r

library(hierfstat); library(adegenet); library(ggplot2)

# Summary statistics grouped by region × year
stats <- basic.stats(dat)

# Pairwise FST: by specific location × year, and by pooled region × year
fst <- pairwise.WCfst(dat)

# PCA
pca <- dudi.pca(tab(genind_obj, NA.method = "mean"), scannf = FALSE, nf = 2)

# ADMIXTURE: K=1 to n (n = total sampling locations in dataset)
# Run outside R via command line, then plot Q matrices in ggplot2
# Optimal K selected by lowest cross-validation error
```

```bash
# ADMIXTURE (run per K value)
for K in $(seq 1 17); do
    admixture --cv filtered_snps.bed $K | tee log_K${K}.out
done

grep "CV" log_K*.out | sort -k3 -n   # identify best K
```

### 8. Temporal structure

Sites sampled in ≥3 years were analyzed for temporal genetic stability:

| Region | Site | Years sampled |
|--------|------|--------------|
| Pamlico Sound | Swanquarter | 2014, 2015, 2016 |
| Neuse River | Hancock Creek | 2014, 2016, 2022, 2023 |
| Neuse River | Slocum Creek | 2014, 2016, 2023 |
| Neuse River | Clubfoot Creek | 2014, 2016, 2022, 2023 |
| New River | Mill Creek | 2014, 2015, 2016 |
| New River | Virginia Creek | 2014, 2016, 2023 |

For each temporal dataset: pairwise FST (hierfstat) + PCA (adegenet).


## Dependencies

### Command-line tools
| Tool | Version | Purpose |
|------|---------|---------|
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | — | Read quality assessment |
| [Stacks](https://catchenlab.life.illinois.edu/stacks/) | 1.24 | Demultiplexing and de novo SNP calling |
| [PLINK](https://www.cog-genomics.org/plink/) | 1.19 | MAF, missingness, and individual filtering |
| [Treemix](https://bitbucket.org/nygcresearch/treemix/) | 1.13 | Species verification |
| [ADMIXTURE](https://dalexander.github.io/admixture/) | 1.3 | Ancestry estimation |
| [PGDSpider](http://www.cmpg.unibe.ch/software/PGDSpider/) | 2.1.1.0 | File format conversion |

### R packages
| Package | Version | Purpose |
|---------|---------|---------|
| [hierfstat](https://cran.r-project.org/package=hierfstat) | 0.5-11 | FST, summary statistics |
| [adegenet](https://github.com/thibautjombart/adegenet) | 1.3-1 | PCA |
| [ade4](https://cran.r-project.org/package=ade4) | 1.7-23 | Multivariate analysis |
| [dartR](https://cran.r-project.org/package=dartR) | — | HWE filtering |
| [ggplot2](https://ggplot2.tidyverse.org/) | — | Visualization |

---

## Citation

Harned, S., Burford Reiskind, M.O. et al. (*in review*). An assessment of spatial and temporal population structure in southern flounder (*Paralichthys lethostigma*) in North Carolina. Preprint: https://doi.org/10.64898/2026.04.24.720543

**Data collection supported by:** NC Division of Marine Fisheries (P120 survey); Texas Parks and Wildlife

---

## Contact

Sydney Harned · spharned@ncsu.edu · [sydneyharned.com](https://sydneyharned.com)
