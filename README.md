# Human Multi-omics Pipeline: Рабочий TODO List

## **Phase 1: Environment Setup** [Week 1]

### **1.1 System Requirements Check**
- [ ] Verify available storage space (minimum 4TB)
- [ ] Check RAM availability (minimum 64GB recommended)
- [ ] Ensure stable internet connection (for large downloads)
- [ ] Confirm Linux/Unix environment or WSL on Windows

###  **1.2 Install Core Tools**
```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install -y wget curl git build-essential

# Install conda/mamba
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh

# Install bioinformatics tools
conda install -c bioconda -c conda-forge \
  bedtools=2.31.0 \
  samtools=1.17 \
  deeptools=3.5.4 \
  homer=4.11 \
  macs2=2.2.9.1 \
  ucsc-bigwigtobedgraph=377 \
  ucsc-bigwigaverageoverbed=377
```

### **1.3 Create Project Structure**
```bash
mkdir -p human_multiomics_pipeline/{data,scripts,results,logs}
cd human_multiomics_pipeline

# Data subdirectories
mkdir -p data/{reference,encode,roadmap,conservation,validation,processed}
mkdir -p data/encode/{atac_seq,metadata}
mkdir -p data/roadmap/{histone_peaks,signal_tracks,chromhmm_states}
mkdir -p data/conservation/{phylop,phastcons,gerp}
```

---

## **Phase 2: Reference Data** [Week 1]

### **2.1 Download Human Genome (hg38)**
```bash
cd data/reference/

# Human genome sequence
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa

# Chromosome sizes
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
```
**Expected size:** ~3GB compressed, ~3.2GB uncompressed

### **2.2 Download Gene Annotations**
```bash
# GENCODE annotations (current version)
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz
gunzip gencode.v44.basic.annotation.gtf.gz
```
**Expected size:** ~50MB

### **2.3 Process Annotations**
```bash
# Extract TSS positions
awk '$3=="gene"' gencode.v44.basic.annotation.gtf | \
awk 'BEGIN{OFS="\t"} {
  gsub(/[";]/, "", $10); gsub(/[";]/, "", $14);
  if($7=="+") print $1, $4-1, $4, $10, $6, $7;
  else print $1, $5-1, $5, $10, $6, $7
}' > gene_tss.bed

# Create promoter regions (±2kb)
awk 'BEGIN{OFS="\t"} {
  start = $2-2000; if(start < 0) start = 0;
  end = $3+2000;
  print $1, start, end, $4, $5, $6
}' gene_tss.bed > promoter_regions.bed
```

### **2.4 Download Blacklisted Regions**
```bash
wget -c https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz
```

---

## **Phase 3: Validation Datasets** [Week 1]

### **3.1 VISTA Enhancer Database**
```bash
cd data/validation/

# VISTA enhancers (preprocessed hg38 coordinates)
wget -c https://enhancer.lbl.gov/cgi-bin/imagedb3.pl?form=presentation_form&show=1&organism=Human&sequence=1 \
  -O vista_raw.html

# Parse VISTA coordinates (manual processing required)
# Alternative: Use community-curated version
wget -c https://raw.githubusercontent.com/ENCODE-DCC/atac-seq-pipeline/master/genome/hg38/hg38_vista_enhancers.bed.gz
gunzip hg38_vista_enhancers.bed.gz
```
**Note:** VISTA requires manual curation - check coordinates manually

### **3.2 FANTOM5 Enhancers**
```bash
# FANTOM5 human enhancers
wget -c http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/F5.hg38.enhancers.bed.gz
gunzip F5.hg38.enhancers.bed.gz
mv F5.hg38.enhancers.bed fantom5_enhancers.bed
```

### **3.3 ENCODE cCREs**
```bash
# ENCODE candidate cis-Regulatory Elements
wget -c https://www.encodeproject.org/files/ENCFF503GCK/@@download/ENCFF503GCK.bed.gz \
  -O encode_ccres_hg38.bed.gz
gunzip encode_ccres_hg38.bed.gz

# Filter for enhancer-like elements
awk '$4 ~ /[Ee]nhancer/ || $4 ~ /[Dd]istal/' encode_ccres_hg38.bed > encode_enhancers.bed
```

---

## **Phase 4: Conservation Data** [Week 1-2]

### **4.1 PhyloP Conservation Scores**
```bash
cd data/conservation/phylop/

# PhyloP 100-way scores (bigWig format)
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw
```
**Expected size:** ~8GB

### **4.2 PhastCons Elements**
```bash
cd ../phastcons/

# PhastCons conserved elements
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons100way/hg38.phastCons100way.bed.gz
gunzip hg38.phastCons100way.bed.gz

# PhastCons scores
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons100way/hg38.phastCons100way.bw
```

### **4.3 GERP++ Data**
```bash
cd ../gerp/

# GERP++ constrained elements (alternative source if main fails)
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/gerpppElements.txt.gz \
  || echo "GERP++ data unavailable, proceeding without"
```

---

## **Phase 5: ENCODE ATAC-seq Data** [Week 2-3]

### **5.1 Priority Cell Lines Setup**
```bash
cd data/encode/atac_seq/

# Define priority cell lines
declare -a CELL_LINES=("GM12878" "K562" "HepG2" "A549" "HUVEC" "H1-hESC")
```

### **5.2 ENCODE File Downloads**
**Working ENCODE file accessions (verified 2024):**

```bash
# GM12878 ATAC-seq
wget -c https://www.encodeproject.org/files/ENCFF356YES/@@download/ENCFF356YES.bed.gz \
  -O GM12878_ENCFF356YES_peaks.bed.gz

wget -c https://www.encodeproject.org/files/ENCFF833POA/@@download/ENCFF833POA.bigWig \
  -O GM12878_ENCFF833POA_signal.bigWig

# K562 ATAC-seq  
wget -c https://www.encodeproject.org/files/ENCFF038DDS/@@download/ENCFF038DDS.bed.gz \
  -O K562_ENCFF038DDS_peaks.bed.gz

wget -c https://www.encodeproject.org/files/ENCFF045OAB/@@download/ENCFF045OAB.bigWig \
  -O K562_ENCFF045OAB_signal.bigWig

# HepG2 ATAC-seq
wget -c https://www.encodeproject.org/files/ENCFF031FSF/@@download/ENCFF031FSF.bed.gz \
  -O HepG2_ENCFF031FSF_peaks.bed.gz
```

### **5.3 Process ATAC-seq Data**
```bash
# Uncompress all bed files
find . -name "*.bed.gz" -exec gunzip {} \;

# Merge all ATAC-seq peaks
cat *_peaks.bed | sort -k1,1 -k2,2n > ../processed/all_atac_peaks.bed
bedtools merge -i ../processed/all_atac_peaks.bed -d 200 > ../processed/merged_atac_peaks.bed
```

---

## **Phase 6: Roadmap Epigenomics Data** [Week 2-4]

### **6.1 Histone Modification Peaks**
```bash
cd data/roadmap/histone_peaks/

# Core histone marks
declare -a MARKS=("H3K4me1" "H3K4me3" "H3K27ac" "H3K27me3" "H3K9me3" "H3K36me3")

for mark in "${MARKS[@]}"; do
  echo "Downloading $mark..."
  mkdir -p $mark
  
  # Download consolidated peaks
  wget -r -np -nH --cut-dirs=6 -A "*.broadPeak.gz" \
    "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/$mark/" \
    -P $mark/ 2>/dev/null
  
  # Uncompress
  find $mark/ -name "*.gz" -exec gunzip {} \;
done
```
**Expected size:** ~500GB total

### **6.2 ChromHMM States**
```bash
cd ../chromhmm_states/

# Representative epigenomes for major cell types
declare -a EPIGENOMES=("E003" "E116" "E118" "E119" "E122" "E123" "E127")

for epi in "${EPIGENOMES[@]}"; do
  wget -c "http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/${epi}_15_coreMarks_segments.bed.gz"
done

gunzip *.gz
```

---

## **Phase 7: Data Processing & Integration** [Week 4]

### **7.1 Create Analysis Intervals**
```bash
cd data/processed/

# Create 200bp windows across genome
bedtools makewindows -g ../reference/hg38.chrom.sizes -w 200 > genome_200bp_windows.bed

# Remove blacklisted regions
bedtools subtract -a genome_200bp_windows.bed \
  -b ../reference/hg38-blacklist.v2.bed > clean_200bp_windows.bed

# Intersect with accessible regions
bedtools intersect -a clean_200bp_windows.bed \
  -b merged_atac_peaks.bed -u > accessible_intervals.bed
```

### **7.2 Process Histone Modifications**
```bash
# Merge peaks for each histone mark
for mark in H3K4me1 H3K4me3 H3K27ac H3K27me3; do
  echo "Processing $mark..."
  
  find ../roadmap/histone_peaks/$mark -name "*.broadPeak" -exec cat {} \; | \
  cut -f1-3 | sort -k1,1 -k2,2n > ${mark}_all_peaks.bed
  
  bedtools merge -i ${mark}_all_peaks.bed -d 200 > ${mark}_merged_peaks.bed
done
```

### **7.3 Feature Extraction**
```bash
# Extract ATAC-seq signal for each interval
bigWigAverageOverBed ../encode/atac_seq/GM12878_*_signal.bigWig \
  accessible_intervals.bed atac_signal_values.txt

# Extract conservation scores
bigWigAverageOverBed ../conservation/phylop/hg38.phyloP100way.bw \
  accessible_intervals.bed phylop_values.txt

# Calculate distances to TSS
bedtools closest -a accessible_intervals.bed -b ../reference/gene_tss.bed -d > intervals_tss_distance.txt
```

---

## **Phase 8: Scoring & Prediction** [Week 5]

### **8.1 Create Feature Matrix**
```bash
# Combine all features into matrix using awk
awk 'BEGIN{OFS="\t"; print "chr","start","end","atac_signal","h3k4me1","h3k27ac","h3k4me3","conservation","dist_tss"} 
FNR==NR{atac[NR]=$2; next}
FNR==NR{cons[NR]=$2; next}
FNR==NR{dist[NR]=$NF; next}
{
  # Binary overlap with histone marks
  h3k4me1 = (system("bedtools intersect -a <(echo \""$1"\t"$2"\t"$3"\") -b H3K4me1_merged_peaks.bed -u | wc -l") > 0) ? 1 : 0
  h3k27ac = (system("bedtools intersect -a <(echo \""$1"\t"$2"\t"$3"\") -b H3K27ac_merged_peaks.bed -u | wc -l") > 0) ? 1 : 0
  h3k4me3 = (system("bedtools intersect -a <(echo \""$1"\t"$2"\t"$3"\") -b H3K4me3_merged_peaks.bed -u | wc -l") > 0) ? 1 : 0
  
  print $1,$2,$3,atac[FNR],h3k4me1,h3k27ac,h3k4me3,cons[FNR],dist[FNR]
}' atac_signal_values.txt phylop_values.txt intervals_tss_distance.txt accessible_intervals.bed > feature_matrix.txt
```

### **8.2 Apply Scoring Algorithm**
```bash
# Calculate composite scores
awk 'BEGIN{OFS="\t"} NR>1 {
  # Normalize features (simple min-max scaling)
  atac_norm = ($4 > 50) ? 1 : $4/50
  cons_norm = ($8 > 2) ? 1 : $8/2
  
  # Distance penalty
  dist_penalty = ($9 < 2000) ? 0.5 : 1.0
  
  # Composite score
  score = (atac_norm * 0.30) + ($5 * 0.20) + ($6 * 0.25) + (cons_norm * 0.15) + (dist_penalty * 0.10)
  
  print $1,$2,$3,score,$4,$5,$6,$7,$8,$9
}' feature_matrix.txt | sort -k4,4nr > scored_intervals.txt
```

### **8.3 Select Top Candidates**
```bash
# Select high-confidence predictions (score > 0.7)
awk '$4 > 0.7' scored_intervals.txt > high_confidence_enhancers.bed

# Get top 5000 predictions
head -n 5000 scored_intervals.txt > top_5000_predictions.bed
```

---

## **Phase 9: Validation & Analysis** [Week 5-6]

### **9.1 VISTA Validation**
```bash
cd data/processed/

# Calculate overlap with VISTA enhancers
bedtools intersect -a top_5000_predictions.bed -b ../validation/hg38_vista_enhancers.bed -u > vista_validated.bed

# Calculate metrics
total_predicted=$(wc -l < top_5000_predictions.bed)
vista_overlap=$(wc -l < vista_validated.bed)
total_vista=$(wc -l < ../validation/hg38_vista_enhancers.bed)

echo "Precision: $(echo "scale=3; $vista_overlap/$total_predicted" | bc)"
echo "Recall: $(echo "scale=3; $vista_overlap/$total_vista" | bc)"
```

### **9.2 Cross-validation with Other Datasets**
```bash
# FANTOM5 validation
bedtools intersect -a top_5000_predictions.bed -b ../validation/fantom5_enhancers.bed -u > fantom5_validated.bed

# ENCODE cCREs validation  
bedtools intersect -a top_5000_predictions.bed -b ../validation/encode_enhancers.bed -u > encode_validated.bed

# Calculate consensus (present in multiple databases)
bedtools intersect -a vista_validated.bed -b fantom5_validated.bed -u | \
bedtools intersect -a - -b encode_validated.bed -u > consensus_validated.bed
```

### **9.3 Motif Analysis**
```bash
# Extract sequences for top predictions
bedtools getfasta -fi ../reference/hg38.fa -bed top_5000_predictions.bed > top_predictions.fa

# HOMER motif analysis
findMotifsGenome.pl top_5000_predictions.bed hg38 motif_results/ \
  -size 200 -mask -len 8,10,12 -p 8
```

---

## **Phase 10: Results & Reporting** [Week 6]

### **10.1 Generate Summary Statistics**
```bash
cd results/

# Create final results summary
cat > pipeline_summary.txt << EOF
Human Multi-omics Enhancer Prediction Results
=============================================

Input Data:
- ATAC-seq peaks: $(wc -l < ../data/processed/merged_atac_peaks.bed)
- Analysis intervals: $(wc -l < ../data/processed/accessible_intervals.bed)
- VISTA enhancers: $(wc -l < ../data/validation/hg38_vista_enhancers.bed)

Predictions:
- Total predictions: $(wc -l < ../data/processed/top_5000_predictions.bed)
- High-confidence (>0.7): $(wc -l < ../data/processed/high_confidence_enhancers.bed)

Validation:
- VISTA overlap: $(wc -l < ../data/processed/vista_validated.bed)
- FANTOM5 overlap: $(wc -l < ../data/processed/fantom5_validated.bed)
- Consensus validated: $(wc -l < ../data/processed/consensus_validated.bed)
EOF
```

### **10.2 Create Final Output Files**
```bash
# Copy key results to results directory
cp ../data/processed/top_5000_predictions.bed ./final_enhancer_predictions.bed
cp ../data/processed/vista_validated.bed ./vista_validated_enhancers.bed
cp ../data/processed/consensus_validated.bed ./high_confidence_validated.bed
cp ../data/processed/feature_matrix.txt ./
```

---

## **Critical Checkpoints**

### **Data Integrity Checks**
- [ ] Verify hg38.fa file size (~3.2GB)
- [ ] Check ATAC-seq files downloaded correctly
- [ ] Validate BED file formats (3+ columns)
- [ ] Confirm conservation files accessible

### **Quality Control**
- [ ] ATAC-seq peak counts reasonable (>10K per cell line)
- [ ] Feature matrix has expected number of rows
- [ ] Scoring produces reasonable distribution (0-1 range)
- [ ] VISTA overlap >500 enhancers (expected minimum)

### **Expected Performance Targets**
- [ ] Precision vs VISTA: >60%
- [ ] Recall vs VISTA: >50%
- [ ] F1-score: >55%
- [ ] Total runtime: <1 week

---

## **Troubleshooting Common Issues**

### **Large File Downloads**
```bash
# Resume interrupted downloads
wget -c [URL]

# Use aria2 for faster parallel downloads
aria2c -x 4 -s 4 [URL]
```

### **Memory Issues**
```bash
# Sort large files externally
sort -T /tmp -k1,1 -k2,2n large_file.bed > sorted_file.bed

# Process chromosomes individually if needed
for chr in {1..22} X Y; do
  grep "^chr${chr}" input.bed > chr${chr}_subset.bed
done
```

### **Missing Dependencies**
```bash
# Install missing tools
conda install -c bioconda [tool_name]

# Check tool versions
bedtools --version
samtools --version
```

---

## **Final Verification Checklist**

- [ ] All downloads completed successfully
- [ ] Feature matrix contains expected number of features
- [ ] Predictions generated and ranked
- [ ] Validation performed against multiple datasets  
- [ ] Results files created in results/ directory
- [ ] Pipeline performance meets target metrics
- [ ] All intermediate files preserved for reproducibility

**Estimated Total Time:** 4-6 weeks
**Storage Used:** ~3-4TB
**Key Output:** final_enhancer_predictions.bed with scored regulatory elements
