# Mapping the reads

In case you get reads from multiple sequencing lanes you have to consolidate them in a single file per individual

```
# Create a directory to store concatenated individual reads.
mkdir /data/scratch/projects/punim1528/cat_ind_reads

# Iterate through unique identifiers extracted from a list of filenames.
for ind in $(cat /data/gpfs/projects/punim1528/a_minax/reads/list_names.txt | grep -oE '^\w{2,3}\_\w{1}\_\w{2}' | sort | uniq); do 

    # Retrieve and concatenate the first set of filtered reads for the current identifier.
    reads1=$(ls /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/clean_ready_to_assemble/*1.gz | grep $ind)
    zcat $reads1 > /data/scratch/projects/punim1528/cat_ind_reads/"$ind".1.gz

    # Retrieve and concatenate the second set of filtered reads for the current identifier.
    reads2=$(ls /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/clean_ready_to_assemble/*2.gz | grep $ind)
    zcat $reads2 > /data/scratch/projects/punim1528/cat_ind_reads/"$ind".2.gz
done
```

Because our experimental design consisted in multiple biological replicates, we need to create a file with this information. In our case

Mapeamos las lecturas a la referencia:

```
module load conda
source activate /home/fabianc.salgado/data/POCtrinity
module load perl
module load rsem

for i in $(seq 1 9)
do
align_and_estimate_abundance.pl --seqType fq \
--left /home/fabianc.salgado/shared/paula_torres/gasteracantha/filter_reads/remove_overrep/rmoverrep_blacklist_paired_unaligned_PTUR00$i.left.fa \
--right /home/fabianc.salgado/shared/paula_torres/gasteracantha/filter_reads/remove_overrep/rmoverrep_blacklist_paired_unaligned_PTUR00$i.right.fa \
--transcripts /home/fabianc.salgado/shared/paula_torres/gasteracantha/Repeat/run_repeatmasker/rmblast/run_rmblast_2000/trinity_cdhit_2000.fasta.masked \
--est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --output_dir /home/fabianc.salgado/shared/paula_torres/gasteracantha/mapeo/mapeo_2000/PTUR00"$i".RSEM
done
```
