# Identifying coding regions using ```TransDecoder``` and ```hmmer```

First, we will use [TransDecoder](https://github.com/TransDecoder/TransDecoder) for two porpuses, identify dentify long open reading frames (ORFs) using ```TransDecoder.LongOrfs``` make final predictions about which ORFs in our transcripts are real with ```TransDecoder.Predict```.
[hmmer](http://hmmer.org/) will be use for searching sequence databases for sequence homologs.


```bash
#!/bin/bash
#SBATCH -N 1 # Número de nodos
#SBATCH -n 16 # Número de núcleos
#SBATCH -t 4-23:00 # Límite de tiempo (D-HH:MM)       
#SBATCH -o transdeco.out # Salida STDOUT
#SBATCH -e transdeco.err # Salida STDERR
# mail alert at start, end and abortion of execution  
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=fsalgadoroa@student.unimelb.edu.au

module load TransDecoder/5.5.0
module load HMMER/3.3.2

TransDecoder.LongOrfs -t /data/scratch/projects/punim1528/assembly_M3_A/coding_regions/trinity_cdhit_aminax_3M_A.fasta

hmmscan --cpu 16 --domtblout /data/scratch/projects/punim1528/assembly_M3_A/coding_regions/pfam.domtblout /data/scratch/projects/punim1528/assembly_M3_A/coding_regions/db/Pfam-A.hmm trinity_combine.fasta.transdecoder_dir/longest_orfs.pep

TransDecoder.Predict -t /data/scratch/projects/punim1528/assembly_M3_A/coding_regions/trinity_cdhit_aminax_3M_A.fasta --retain_pfam_hits /data/scratch/projects/punim1528/assembly_M3_A/coding_regions/pfam.domtblout --cpu 16

```
