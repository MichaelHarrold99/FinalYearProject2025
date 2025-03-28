Description:
#   This script is part of my final-year bioinformatics project investigating
#   potential hypermutator genes across multiple bacterial genomes. It uses the
#   NCBI Datasets CLI to download genomic data from a specified BioProject, merges
#   all the genome FASTA files, creates a BLAST database, and runs TBLASTN queries
#   for selected protein sequences.
#
# Usage:
#   1. Make sure both the NCBI Datasets CLI and BLAST+ tools (makeblastdb, tblastn) 
#      are installed and available in your PATH.
#   2. Place your query protein files in the directory specified by QUERY_DIR.
#   3. Update the variables (e.g., BIOPROJECT, QUERY_PROTEINS) if needed.
#   4. Run this script in a Unix-based shell (e.g., ./hypermutator_tblastn_pipeline.sh).
#!/bin/bash
BIOPROJECT="PRJNA325248"
COMBINED_FASTA="combined_genomes.fna"
BLAST_DB_NAME="my_combined_db"
WORK_DIR="$HOME/${BIOPROJECT}_workflow"
QUERY_PROTEINS=("pa01_mutL.faa"  "pa01_mutS_protein.faa" "pa01_uvrD.faa")
QUERY_DIR="$HOME/queries"
OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq"

mkdir -p "$WORK_DIR"
cd "$WORK_DIR" || { echo "Cannot change to directory $WORK_DIR. Exiting."; exit 1; }

echo "Downloading BioProject $BIOPROJECT data..."
datasets download genome accession  --assembly-source "GenBank" "$BIOPROJECT" --dehydrated --filename "${BIOPROJECT}_dehydrated.zip"
unzip -o "${BIOPROJECT}_dehydrated.zip"

echo "Rehydrating dataset..."
datasets rehydrate --directory "$WORK_DIR"

echo "Collecting genome FASTA files..."
find ncbi_dataset/data/ -type f -name "*.fna" -exec mv {} . \;

echo "Combining genome FASTA files into ${COMBINED_FASTA}..."
cat GCA*.fna > "$COMBINED_FASTA"

echo "Creating BLAST database: $BLAST_DB_NAME..."
makeblastdb -in "$COMBINED_FASTA" -dbtype nucl -parse_seqids -out "$BLAST_DB_NAME" -title "Combined Genome Database for $BIOPROJECT"

for QUERY_FILE in "${QUERY_PROTEINS[@]}"; do
  if [[ -f "$QUERY_DIR/$QUERY_FILE" ]]; then
    BASE_NAME="${QUERY_FILE%.*}"
    OUTPUT="${BASE_NAME}_tblastn_results.txt"
    echo "Running TBLASTN for $QUERY_FILE..."
    tblastn -query "$QUERY_DIR/$QUERY_FILE" -db "$BLAST_DB_NAME" -outfmt "$OUTFMT" -out "$OUTPUT" -evalue 1e-5 -num_threads 4
  else
    echo "Warning: Query file '$QUERY_DIR/$QUERY_FILE' not found. Skipping."
  fi
done

echo "Workflow complete."
