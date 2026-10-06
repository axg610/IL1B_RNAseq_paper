```powershell
ssh alex.gao1@arc.ucalgary.ca
```

## setup hpc

```bash SETUP TOOLS
salloc -c 8 --mem=32GB -t 05:00:00
export PATH="/home/alex.gao1/anaconda3/bin":$PATH
export PATH="/home/alex.gao1/tools/samtools-1.21":$PATH
conda init
conda activate genomics-tools

# cd /home/alex.gao1/tools
# apptainer pull docker://memesuite/memesuite

# cd /work/newton_lab/ag_analysis/ref_seq
# wget https://hocomoco13.autosome.org/final_bundle/hocomoco13/H13CORE/formatted_motifs/H13CORE_meme_format.meme
```

```bash SETUP DIRECTORIES
mkdir -p /work/newton_lab/ag_analysis/a549-il1b-bud-timecourse/xstreme
cd /work/newton_lab/ag_analysis/a549-il1b-bud-timecourse/xstreme

# move in /genelists from local...
```

## define and pull promoters for each kinetic cluster
```bash

GENE_DIR=/work/newton_lab/ag_analysis/a549-il1b-bud-timecourse/xstreme/genelists
PROMOTER_DIR=/work/newton_lab/ag_analysis/a549-il1b-bud-timecourse/xstreme/promoters
TEMP_DIR=/work/newton_lab/ag_analysis/a549-il1b-bud-timecourse/xstreme/_temp

GTF=/work/newton_lab/ag_analysis/ref_seq/hg38.gtf
GENOME=/work/newton_lab/ag_analysis/ref_seq/hg38.fa
# SIZES=/work/newton_lab/ag_analysis/ref_seq/genome.sizes

UP=1000
DOWN=100
# -1000, +100 seems to be common?
# PMC7778956  

mkdir -p "$PROMOTER_DIR" "$TEMP_DIR"

for list in "$GENE_DIR"/*.txt; do
    base=$(basename "$list" .txt)
    echo "Processing cluster: $base"

    # filter GTF for transcripts matching genes in this cluster
    awk 'NR==FNR {gsub("\r",""); genes[$1]; next} 
         $3=="transcript" {
             match($0, /gene_name "([^"]+)"/, arr)
             if (arr[1] in genes) print
         }' "$list" "$GTF" > "$TEMP_DIR/${base}_transcripts.gtf"

    # extract one TSS per gene
    awk -v UP=$UP -v DOWN=$DOWN 'BEGIN{OFS="\t"} 
         {
           match($0, /gene_name "([^"]+)"/, arr); gene=arr[1]
           if ($7 == "+") {
               tss=$4
               start=(tss-UP<0?0:tss-UP)
               end=tss+DOWN
           } else {
               tss=$5
               start=(tss-DOWN<0?0:tss-DOWN)
               end=tss+UP
           }
           print $1, start, end, gene, ".", $7
         }' "$TEMP_DIR/${base}_transcripts.gtf" \
         | sort -k4,4 -u > "$TEMP_DIR/${base}_tss_window.bed"

    # filter to canonical chromosomes
    awk '$1 ~ /^chr([0-9]+|X|Y|M)$/' "$TEMP_DIR/${base}_tss_window.bed" \
        > "$TEMP_DIR/${base}_tss_window_canonical.bed"

    # extract FASTA sequences
    bedtools getfasta -fi "$GENOME" \
                      -bed "$TEMP_DIR/${base}_tss_window_canonical.bed" \
                      -fo "$PROMOTER_DIR/${base}_promoters.fa" \
                      -s -name

done
```

## script for running xstreme
```bash 
cat > run_xstreme.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=xstreme
#SBATCH --output=/work/newton_lab/ag_analysis/a549-il1b-bud-timecourse/xstreme/logs/xstreme_%j.out
#SBATCH --error=/work/newton_lab/ag_analysis/a549-il1b-bud-timecourse/xstreme/logs/xstreme_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

FASTA="$1"
BACKGROUND="$2"
OUTDIR="$3"

DB="/work/newton_lab/ag_analysis/ref_seq/H13CORE_meme_format.meme"
MEME="/home/alex.gao1/tools/memesuite_latest.sif"

mkdir -p "$OUTDIR"

BASENAME=$(basename "$FASTA" .fa)
RESULTDIR="$OUTDIR/${BASENAME}"

apptainer exec --bind /work:/work "$MEME" xstreme \
        --p "$FASTA" \
        --n "$BACKGROUND" \
        --oc "$RESULTDIR" \
        --m "$DB" \
        --dna --mea-only

# Extract main result files
if [ -f "$RESULTDIR/xstreme.html" ]; then
    cp "$RESULTDIR/xstreme.html" "$OUTDIR/${BASENAME}.html"
    cp "$RESULTDIR/xstreme.tsv" "$OUTDIR/${BASENAME}.txt"
    echo "EXTRACTED: $BASENAME"
else
    echo "WARNING: $RESULTDIR/xstreme.html not found"
fi
EOF
```

## start xstreme runs
```bash

# db is HOCOLOCO V13
# background is promoters from all clusters combined

BACKGROUND="./promoters/all_promoters.fa"
OUTPUT="./xstreme_promoters"

mkdir -p "$OUTPUT"

for f in ./promoters/*.fa; do

    FASTA="$f"

    if [[ "$FASTA" == "./promoters/all_promoters.fa" ]]; then
        echo "Skipping $FASTA"
        continue
    fi

    sbatch run_xstreme.sh "$FASTA" "$BACKGROUND" "$OUTPUT"

done
```

