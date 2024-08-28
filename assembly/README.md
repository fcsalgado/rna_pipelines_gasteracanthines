# Transcriptome Assembly

To assemble a de novo transcriptome for the Christmas spider, we used two different de novo transcriptome assembly algorithms with all the sequences and default parameters: Trinity v2.15.1 (Grabherr et al., 2011), which uses a k-mer size of 25, and rnaSPAdes v3.15.5 (Bushmanova et al., 2019) with selected k-mer sizes of 49 and 73. 

## Trinity Assembly
```bash
#!/bin/bash
#SBATCH -N 1 # Número de nodos
#SBATCH -n 20 # Número de núcleos
#SBATCH -t 4-23:00 # Límite de tiempo (D-HH:MM)
#SBATCH -o trinity.out # Salida STDOUT
#SBATCH -e trinity.err # Salida STDERR
# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=fsalgadoroa@student.unimelb.edu.au

module load Trinity/2.15.1 #trinityis super picky with the dependencies. Keep that in mind
module load Salmon/1.5.2

# --min_contig_length is critical parameter, since it will filter away many reads. For Austracantha I left --min_contig_length 300

Trinity --seqType fq --SS_lib_type RF --max_memory 50G --CPU 20 --min_contig_length 300 --left /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/clean_ready_to_assemble/*.1.gz  --right /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/clean_ready_to_assemble/*.2.gz --output /data/scratch/projects/punim1528/trinity_output --full_cleanup
```

## rnaSPAdes
```bash
#!/bin/bash
#SBATCH -N 1 # Número de nodos
#SBATCH -n 20 # Número de núcleos
#SBATCH -t 4-23:00 # Límite de tiempo (D-HH:MM)       
#SBATCH -o spades.out # Salida STDOUT
#SBATCH -e spades.err # Salida STDERR
#SBATCH --mem=200G
# mail alert at start, end and abortion of execution  
#SBATCH --mail-type=ALL
# send mail to this address
#SBATCH --mail-user=fsalgadoroa@student.unimelb.edu.au

module load SPAdes/3.15.5

rnaspades.py -t 20 --dataset  /data/scratch/projects/punim1528/config.yaml -o /data/scratch/projects/punim1528/spades_aligment
```

## Evigene
```bash
/home/fcsalgado/software/evigene/scripts/rnaseq/trformat.pl -out /scratch4/a_minax/transcriptomes_evigene/trinity.transcripts.fasta -log -in /scratch4/a_minax/assembly_M1_M3_A/trinity_M1_M3_A.Trinity.fasta

/home/fcsalgado/software/evigene/scripts/rnaseq/trformat.pl -out /scratch4/a_minax/transcriptomes_evigene/spades.transcripts.fasta -log -in /scratch4/a_minax/spades/transcripts.fasta

cat /scratch4/a_minax/transcriptomes_evigene/*.fasta > /scratch4/a_minax/transcriptomes_evigene/transcript_joint.fa

cd /scratch4/a_minax/transcriptomes_evigene/
```


# Transcriptome Assembly Quality Assessment

This is really important, especially because we do not have a reference. We will follow all the steps recommended [here](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Assembly-Quality-Assessment). 

## The 'Gene' Contig Nx Statistic

This is the classic statistic for assemblies. Because this is de novo assembly and most of the transcripts are small, we want N50 values not too high. Check [here](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats) for details

```bash
module load conda
source activate /home/fabianc.salgado/data/POCtrinity

/datacnmat01/biologia/biologia.evolutiva/fabianc.salgado/POCtrinity/opt/trinity-2.9.1/util/TrinityStats.pl /home/fabianc.salgado/shared/paula_torres/gasteracantha/trinity/trinity_without_2000/Trinity_2000.fasta
```

## Run BUSCO

This command is evaluating the completeness of a transcriptome assembly (in FASTA format) using BUSCO. It utilizes the arachnid lineage dataset (arachnida_odb10) and runs in transcriptome mode with parallel processing on 10 CPU cores. The results will be stored in the specified output directory
```
module load python3.6/3.6.6
python3 /home/fabianc.salgado/shared/paula_torres/gasteracantha/busco/busco/bin/busco \
--config /home/fabianc.salgado/shared/paula_torres/gasteracantha/busco/busco/config/config.ini \
--in /home/fabianc.salgado/shared/paula_torres/gasteracantha/cdhit/trinity_cdhit_2000.fasta -o gasteracantha_2000_busco_out \
-l arachnida_odb10 -m tran -c 10 -f

#!/bin/sh
#SBATCH -N 1 # Número de nodos
#SBATCH -n 10 # Número de núcleos
#SBATCH -t 5-23:00 # Límite de tiempo (D-HH:MM)       
#SBATCH -o busco.out # Salida STDOUT
#SBATCH -e busco.err # Salida STDERR
# mail alert at start, end and abortion of execution  
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=fsalgadoroa@student.unimelb.edu.au

module load BUSCO/5.4.5
busco --in /data/scratch/projects/punim1528/assembly_300/trinity_output.Trinity.fasta -o /data/scratch/projects/punim1528/busco/austracantha_busco_out -l arachnida_odb10 -m tran -c 10 -f

```

Most of your assembly (>80%) needs to match the database. There is a high probability that you ended up with a high number of duplicates. This may be because, for instance, a high level of duplication may be explained by a recent whole duplication event (biological) or a chimeric assembly of haplotypes. Let's discard the technical part and remove the redundant transcripts.

## CDHIT

Remove redundant transcripts with [CD-HIT](https://sites.google.com/view/cd-hit)

```
module load CD-HIT/4.8.1

cd-hit-est -i /data/scratch/projects/punim1528/assembly_300/trinity_output.Trinity.fasta -o /data/scratch/projects/punim1528/trinity_cdhit_aminax -c 0.95 -M 32000 -T 16
```

**After this repeat the BUSCO analyses to check if the matching improves**

## Reads representation in the assembly

Ideally, your transcriptome assembly should encompass approximately 80% or more of the input RNA-Seq reads. The unassembled reads that remain are likely associated with transcripts that are expressed at low levels, lacking adequate coverage for assembly, or are of poor quality and abnormal

```bash
#!/bin/sh
#SBATCH -N 1 # Número de nodos
#SBATCH -n 2 # Número de núcleos
#SBATCH -t 5-23:00 # Límite de tiempo (D-HH:MM)
#SBATCH -o abundance.out # Salida STDOUT
#SBATCH -e abundance.err # Salida STDERR
# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=fsalgadoroa@student.unimelb.edu.au

module load Bowtie2/2.4.5
module load SAMtools/1.16.1

#build bowtie database
bowtie2-build /data/scratch/projects/punim1528/assembly_M3_A/trinity_output.Trinity.fasta /data/scratch/projects/punim1528/assembly_M3_A/trans_M3_A.fna

#map reads to database adults
IFS=$'\n'
cat /data/scratch/projects/punim1528/cat_ind_reads/samples_A.txt | while read ind; do

f1=$(echo ${ind} |cut -f 3)
f2=$(echo ${ind} |cut -f 4)
code=$(echo ${ind} |cut -f 4 | grep -oE "\w+\_\w+\_\w+")

bowtie2 -p 10 -q --no-unal -k 20 -x /data/scratch/projects/punim1528/assembly_M3_A/trans_M3_A.fna -1 /data/scratch/projects/punim1528/cat_ind_reads/"$f1" -2 /data/scratch/projects/punim1528/cat_ind_reads/"$f2" 2> /data/scratch/projects/punim1528/assembly_M3_A/align_stats_"$code".txt | samtools view -@10 -Sb -o bowtie2.bam; done


#map reads to database m3
IFS=$'\n'
cat /data/scratch/projects/punim1528/cat_ind_reads/samples_m3.txt | while read ind; do

f1=$(echo ${ind} |cut -f 3)
f2=$(echo ${ind} |cut -f 4)
code=$(echo ${ind} |cut -f 4 | grep -oE "\w+\_\w+\_\w+")

bowtie2 -p 10 -q --no-unal -k 20 -x /data/scratch/projects/punim1528/assembly_M3_A/trans_M3_A.fna -1 /data/scratch/projects/punim1528/cat_ind_reads/"$f1" -2 /data/scratch/projects/punim1528/cat_ind_reads/"$f2" 2> /data/scratch/projects/punim1528/assembly_M3_A/align_stats_"$code".txt | samtools view -@10 -Sb -o bowtie2.bam; done
```

Check the output of each pair of reads saved in _align_stats_"$ind".txt_

## TransRate

```bash
#run script on cluster: 

#!/bin/sh
#SBATCH -N 1 # Número de nodos
#SBATCH -n 32 # Número de núcleos
#SBATCH -t 5-23:00 # Límite de tiempo (D-HH:MM)
#SBATCH -o abundance.out # Salida STDOUT
#SBATCH -e abundance.err # Salida STDERR
# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=fsalgadoroa@student.unimelb.edu.au

#concatenate right and left reads in a single file
zcat *1.gz > left_reads.gz
zcat *2.gz > right_reads.gz

module load Transrate/1.0.3

transrate --assembly /data/scratch/projects/punim1528/assembly_M3_A/trinity_cdhit_aminax_3M_A.fasta --left left_reads.gz --right right_reads.gz --threads 32

rm -f left_reads.gz right_reads.gz 


```
