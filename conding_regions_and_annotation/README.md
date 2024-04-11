## Identifying coding regions using ```TransDecoder``` and ```hmmer```

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

## Remove transcriptome redundancy

De novo transcriptome assemblies are intricate, encompassing both biological diversity through alternately spliced transcripts and nucleotide sequence variability, alongside technical challenges such as fragmented transcripts. Because we want to qunatify gene expression at the gene level, we aim to simplify much of this complexity by distilling it down to a single transcriptome, where each gene is represented by one transcript. This will enable us to accurately quantify gene expression against a standardized backdrop. During this step, we will group transcripts based on their amino acid sequences and choose a singular representative from each cluster. 

We will use [vsearch](https://github.com/torognes/vsearch) to cluster transcripts with similar sequences and choose one _centroid_ transcript the represent them all. We will use a identity thershold of 90%: 

```bash
#!/bin/bash
#SBATCH -N 1 # Número de nodos
#SBATCH -n 16 # Número de núcleos
#SBATCH -t 4-23:00 # Límite de tiempo (D-HH:MM)       
#SBATCH -o vsearch.out # Salida STDOUT
#SBATCH -e vsearch.err # Salida STDERR
# mail alert at start, end and abortion of execution  
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=fsalgadoroa@student.unimelb.edu.au

module load VSEARCH/2.22.1

vsearch --threads 16 --log LOGFile --cluster_fast /data/scratch/projects/punim1528/assembly_M3_A/coding_regions/trinity_cdhit_aminax_3M_A.fasta.transdecoder.cds --id 0.90 --centroids centroids.fasta --uc clusters.uc

```

