#!/bin/bash

# =========================================================
# === get RefSeq genes with position for Manhattan plot ===
# =========================================================

# add gene biotype
cd /.../03_data_visualization/01_fMRI/Manhattan_QQ_Plots/

# download refseq identifiers
wget -O README.txt https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/README.txt
wget -O GRCh37_latest_genomic.gff.gz https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz

# get chromomsome translation and remove ^M
wget -O GRCh37_NCBI2UCSC.txt https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh37_NCBI2UCSC.txt
cat -v GRCh37_NCBI2UCSC.txt | sed "s/\^M//g" > GRCh37_NCBI2UCSC.tmp.txt; \mv GRCh37_NCBI2UCSC.tmp.txt GRCh37_NCBI2UCSC.txt

# rename chr1..22 X Y MT
awk -F'\t' 'NR==FNR { gsub(/chr/,""); print $0 } $2=="M" { $2=="MT"}' OFS='\t' GRCh37_NCBI2UCSC.txt > GRCh37_NCBI2UCSC_chr122XYMT.txt

# get chromsome (translate), gene, biotype, description, synonyms from resources/GRCh37_latest_genomic.gff.gz
awk -F'\t' 'NR==FNR { chr[$1]=$2; next }
!($1 in chr) { next } $9!~/Name/ || $9!~/biotype/ { next }
{ gene=$9; sub(/.*Name=/,"",gene); sub(/[;].*/,"",gene)}
{ biotype=$9; sub(/.*gene_biotype=/,"",biotype); sub(/[;].*/,"",biotype) } $9!~/biotype/
{ description=$9; sub(/.*description=/,"",description); sub(/[;].*/,"",description) } $9!~/description/ { description="" }
{ synonym=$9; sub(/.*gene_synonym=/,"",synonym); sub(/[;].*/,"",synonym) } $9!~/synonym/ { synonym="" }
{ chrom=chr[$1] } !($1 in chr) { chrom="" } { print $1, $2, $3, $4, $5, $6, $7, $8, chrom, gene, biotype, description, synonym, $9 }' OFS='\t' GRCh37_NCBI2UCSC_chr122XYMT.txt <(gzip -dc GRCh37_latest_genomic.gff.gz) > GRCh37_latest_genomic.edit.gff

# sort by name, start, and stop coordinate and remove duplicate genes
cat GRCh37_latest_genomic.edit.gff | sort -k10,10 -k4,4n -k5,5n > GRCh37_latest_genomic.edit.tmp.gff
\mv GRCh37_latest_genomic.edit.tmp.gff GRCh37_latest_genomic.edit.gff
awk -F'\t' '!($10 in gene) { print; gene[$10]; next}' OFS='\t' GRCh37_latest_genomic.edit.gff > GRCh37_latest_genomic.edit.tmp.gff
\mv GRCh37_latest_genomic.edit.tmp.gff GRCh37_latest_genomic.edit.gff

# extract synonyms
awk -F'\t' '{ print $13 }' GRCh37_latest_genomic.edit.gff > GRCh37_latest_genomic.edit.synonyms.gff

# gzip
chmod 770 *
gzip -f GRCh37_latest_genomic.edit.gff
