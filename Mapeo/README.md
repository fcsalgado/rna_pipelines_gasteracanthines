# Mapeo

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
