
# Transcriptome Assembly

We will use trinity to assemble the transcriptome

```
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

# Transcriptome Assembly Quality Assessment

This is really important, especially because we do not have a reference. We will follow all the steps recommended [here](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Assembly-Quality-Assessment). 

## The 'Gene' Contig Nx Statistic

This is the classic statistic for assemblies. Because this is de novo assembly and most of the transcripts are small, we want N50 values not too high. Check [here](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats) for details

```
module load conda
source activate /home/fabianc.salgado/data/POCtrinity

/datacnmat01/biologia/biologia.evolutiva/fabianc.salgado/POCtrinity/opt/trinity-2.9.1/util/TrinityStats.pl /home/fabianc.salgado/shared/paula_torres/gasteracantha/trinity/trinity_without_2000/Trinity_2000.fasta
```




## Corremos BUSCO
```
module load python3.6/3.6.6
python3 /home/fabianc.salgado/shared/paula_torres/gasteracantha/busco/busco/bin/busco \
--config /home/fabianc.salgado/shared/paula_torres/gasteracantha/busco/busco/config/config.ini \
--in /home/fabianc.salgado/shared/paula_torres/gasteracantha/cdhit/trinity_cdhit_2000.fasta -o gasteracantha_2000_busco_out \
-l arachnida_odb10 -m tran -c 10 -f
```

## CDHIT

Corremos Cdhit para remover transcritos redundantes:

```
/datacnmat01/biologia/biologia.evolutiva/shared/cdhit/cd-hit-est  \
-i /home/fabianc.salgado/shared/paula_torres/gasteracantha/trinity/trinity_without_2000/Trinity_2000.fasta \
-o trinity_cdhit_2000 -c 0.95 -M 32000 -T 16
```

## Quitamos elementos repetitivos

RepeatModeler:

```
module load perl/5.30.3 

perl /home/fabianc.salgado/shared/paula_torres/gasteracantha/Repeat/RepeatModeler/RepeatModeler/RepeatModeler \
-database gasteracantha -engine rmblast -pa 16
```
RepeatMasker:

```
perl /home/fabianc.salgado/shared/paula_torres/gasteracantha/Repeat/RepeatMasker/RepeatMasker \
-pa 16 -gff -lib consensi.fa.classified trinity_cdhit_2000.fasta
```
