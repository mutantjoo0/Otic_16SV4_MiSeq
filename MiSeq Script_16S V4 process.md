# ONR 16S V4 Illumina MiSeq amplicon sequencing analysis

Sequencing data were processed with USEARCH (ver.10.0.240 x64) with UPARSE OTU picking method and subsequent analyses were performed in QIIME (ver.1.8).

# Process with USEARCH

## Merge paired end reads
```
#decompress the reads
gunzip *.gz

#make output dir for merged fastq file
mkdir fastqmerged

#join the forward and reverse reads for each sample, write a single file (merged.fq) will take all samples in it
./usearch10.0.240_i86linux64 -fastq_mergepairs *R1*.fastq -relabel @ -fastq_maxdiffs 10 -fastqout fastqmerged/merged.fq -fastq_merge_maxee 1.0 -fastq_minmergelen 200 -fastq_maxmergelen 300
```

## Dereplicate
```
#remove duplicate sequences (to reduce computation time downstream)
./usearch10.0.240_i86linux64 -fastx_uniques fastqmerged/merged.fq -fastqout fastqmerged/uniques_combined_merged.fastq -sizeout
```

## Remove singletons
```
#remove singletons and reorder data by the length of the sequence
./usearch10.0.240_i86linux64 -sortbysize fastqmerged/uniques_combined_merged.fastq -fastqout fastqmerged/nosigs_uniques_combined_merged.fastq -minsize 2
```

## Precluster 
```
#preclust the sequences
./usearch10.0.240_i86linux64 -cluster_fast fastqmerged/nosigs_uniques_combined_merged.fastq -centroids_fastq fastqmerged/denoised_nosigs_uniques_combined_merged.fastq -id 0.9 -maxdiffs 5 -abskew 10 -sizein -sizeout -sort size
```

## Reference-based OTU picking 
```
#download Silva database (silva 128 qiime release)
refer-> http://www.metagenomics.wiki/tools/16s/qiime/otu-clustering/silva

#perform closed-refernce OTU pick against Silva database. Keep all sequences that could not hit the Silva database. 
./usearch10.0.240_i86linux64 -usearch_global fastqmerged/denoised_nosigs_uniques_combined_merged.fastq -id 0.97 -db /mnt/home/leejooy5/sequence/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta -strand plus -uc fastqmerged/ref_seqs.uc -dbmatched fastqmerged/closed_reference.fasta -notmatchedfq fastqmerged/failed_closed.fq
```

## OTU clustering with chimera filtering
```
#reorder seqeuences that faile to hit the Silva database
./usearch10.0.240_i86linux64 -sortbysize fastqmerged/failed_closed.fq -fastaout fastqmerged/sorted_failed_closed.fq

#de novo pick and chimera check
./usearch10.0.240_i86linux64 -cluster_otus fastqmerged/sorted_failed_closed.fq -minsize 2 -otus fastqmerged/denovo_otus.fasta -relabel OTU_dn_ -uparseout fastqmerged/denovo_out.up
```

## Combine rep sets between de novo and reference-based OTU picking
```
#join the two files that have the sequences for each OTU
 cat fastqmerged/denovo_otus.fasta fastqmerged/closed_reference.fasta > full_rep_set.fna
```

## Map rep set to pre-dereplicated sequences 
```
#map our OTU sequences back to the original dataset to determine the abundance of each OTU in each sample. Write a OTU table (species by sample table)
./usearch10.0.240_i86linux64 -usearch_global fastqmerged/merged.fq -db fastqmerged/full_rep_set.fna  -strand plus -id 0.97 -uc OTU_map.uc -otutabout OTU_table.txt -biomout OTU_jsn.bio
```

# Analysis with QIIME
```
#assign taxonomy to OTUs
assign_taxonomy.py -i fastqmerged/full_rep_set.fna -o taxonomy -r /mnt/home/leejooy5/sequence/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta -t /mnt/home/leejooy5/sequence/SILVA_128_QIIME_release/taxonomy/16S_only/97/consensus_taxonomy_7_levels.txt

#add taxonomy to OTU table
biom add-metadata -i otu_jsn.biom -o otu_table_tax.biom --observation-metadata-fp=taxonomy/full_rep_set_tax_assignments.txt --sc-separated=taxonomy --observation-header=OTUID,taxonomy

#summarize the taxonomy at different levels 
summarize_taxa.py -i otu_table_tax.biom -o taxa_summary

#align sequences to a database
align_seqs.py -i fastqmerged/full_rep_set.fna -o alignment -t /mnt/home/leejooy5/sequence/SILVA_128_QIIME_release/rep_set_aligned/97/97_otus_aligned.fasta

#filter alignment, removes excess gaps
filter_alignment.py -i alignment/full_rep_set_aligned.fasta -o filtered_alignment

#make a tree
make_phylogeny.py -i filtered_alignment/full_rep_set_aligned_pfiltered.fasta -o rep_set.tre

#summarize the OTU table to look at sequenceing depth (min depth=2313)
biom summarize-table -i otu_table_tax.biom -o otu_tablesum.txt

#filter out low sample
filter_samples_from_otu_table.py -i otu_table_tax.biom -o otu_table_tax_filt.biom -n 1000

#rarefy the OTU table (normalized the table to the lowest depth)
single_rarefaction.py -i otu_table_tax_filt.biom -o single_rare.biom -d 2313

#calculate alpha diversiy
alpha_diversity.py -m shannon,observed_species,simpson_e,PD_whole_tree -i single_rare.biom -o alpha_diversty.txt -t rep_set.tre


#calculate beta diversity and make coordinates
beta_diversity.py -i single_rare.biom -o beta_div -m bray_curtis,weighted_unifrac,unweighted_unifrac -t rep_set.tre
principal_coordinates.py -i beta_div -o coords
