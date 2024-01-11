# Sequence Filtering

Initially, it is necessary to filter the reads as they may contain adapters and other artifacts. The filtering protocol used in this case is based on one published by Harvard, which can be found [aqui](https://github.com/harvardinformatics/TranscriptomeAssemblyTools).

## Examining Sequence Quality (FastQC)
The first step with raw reads is to run FastQC to assess the quality of the sequences. On the cluster, it is executed as follows:

```

module load fastqc/0.11.7

fastqc /home/reads/*.gz
```

Luego de analizar la calidad de las lecturas procedemos a filtrar las secuencias.

## Remove Erroneous k-mers

Since state-of-the-art transcriptome assemblers build assemblies using a graph-based approach from k-mers, erroneous k-mers can negatively impact the assemblies. It is recommended to remove them. For this purpose, we use [rCorrector](https://github.com/mourisl/Rcorrector), a tool specifically designed for RNA-seq data. In addition to demonstrating superior performance in parallel comparisons with other kmer-based read error correction tools, it includes labels in the fastq output indicating whether the read has been corrected or if an error has been detected but cannot be corrected.

```
# in some cases companies run samples in m ultiple lanes, with long IDs, in this case you can create a list of files with following command replacing the regular expression depending on the situation

ls *.gz | sed -E "s/\_\R[1,2]\.\w+\.\w+$//g" | sort | uniq > list_names.txt

#!/bin/bash
#SBATCH -N 1 # Número de nodos
#SBATCH -n 10 # Número de núcleos
#SBATCH -t 4-23:00 # Límite de tiempo (D-HH:MM)
#SBATCH -o rcorrector.out # Salida STDOUT
#SBATCH -e rcorrector.err # Salida STDERR
# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=fsalgadoroa@student.unimelb.edu.au

module load Anaconda3/2023.07-2
source activate /home/fsalgadoroa/.conda/envs/rcorrector
#mkdir /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads

cat /data/gpfs/projects/punim1528/a_minax/reads/list_names.txt | while read ind;
do
/data/gpfs/projects/punim1528/a_minax/scripts/rcorrector/run_rcorrector.pl -t 10 -1 /data/gpfs/projects/punim1528/a_minax/reads/"$ind"_R1.fastq.gz -2 /data/gpfs/projects/punim1528/a_minax/reads/"$ind"_R2.fastq.gz -od /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads; done
```

The output fastq files will include "cor" in their names. The corrected reads will have a "cor" suffix in their labels.

```
@ERR1101637.57624042/1 l:165 m:203 h:218 cor
GTACAACCCTTCCAACCTCCACCGTCTTATATACGAAGCGCCTTGAGTGTGTGTGTGCATGAGCCAAAGGGAATACCG
+
<@@D=D2ACDAHB?:<8EGE;B;FE<?;??DBDD))0)8@GDF@C)))5-;5CAE=..7;@DEEC@C;A9?BB=@>@B
```
The values l and h in the read headers indicate the counts of k-mers for the lowest, median, and highest frequency k-mers that occur in the read. There will also be some reads that cannot be corrected.

These uncorrectable reads often contain Ns or represent other low-complexity sequences. There is a Python script to perform this task. It also removes the "cor" label from the headers of corrected sequences, as it may cause issues for downstream tools:

```
#!/bin/bash
#SBATCH -N 1 # Número de nodos
#SBATCH -n 1 # Número de núcleos
#SBATCH -t 4-23:00 # Límite de tiempo (D-HH:MM)
#SBATCH -o rcorrector.out # Salida STDOUT
#SBATCH -e rcorrector.err # Salida STDERR
# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=fsalgadoroa@student.unimelb.edu.au

module load Anaconda3/2023.07-2
source activate /home/fsalgadoroa/.conda/envs/rcorrector
#mkdir /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads

cat /data/gpfs/projects/punim1528/a_minax/reads/list_names.txt | while read ind; do
python /data/gpfs/projects/punim1528/a_minax/scripts/TranscriptomeAssemblyTools/utilities/FilterUncorrectabledPEfastq.py -1 /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/"$ind"_R1.cor.fq.gz -2 /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/"$ind"_R2.cor.fq.gz -s "$ind"; done

```

## Remove sequence adaptors

```
#!/bin/bash
#SBATCH -N 1 # Número de nodos
#SBATCH -n 2 # Número de núcleos
#SBATCH -t 4-23:00 # Límite de tiempo (D-HH:MM)
#SBATCH -o trim.out # Salida STDOUT
#SBATCH -e trim.err # Salida STDERR
# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=fsalgadoroa@student.unimelb.edu.au

module load cutadapt/4.2
module load FastQC/0.12.1-Java-11

cat /data/gpfs/projects/punim1528/a_minax/reads/list_names.txt | while read ind; do
/data/gpfs/projects/punim1528/a_minax/scripts/TrimGalore-0.6.10/trim_galore --cores 2 --paired --retain_unpaired --phred33 --output_dir /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/trimmed_reads --length 36 -q 5 --stringency 1 -e 0.1 --fastqc /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/filtered_uncorrectable/unfixrm_"$ind"_R1.cor.fq.gz /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/filtered_uncorrectable/unfixrm_"$ind"_R2.cor.fq.gz; done
```
## Mapping against SILVA database to remove contamination 

```
# Download SILVA (https://www.arb-silva.de) sequences 

wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_LSUParc_tax_silva.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSUParc_tax_silva.fasta.gz

# unzip and concatenate 

gzip -d *.gz

cat *.fasta > silva.db

rm *.fasta

#remplace U for T because we are working with cDNA sequences

sed -i 's/U/T/g' silva.db

# index the fasta file to be used in bowtie, this can take ages to run and it can be really heavy. Use multiple processors

bowtie2-build /data/gpfs/projects/punim1528/a_minax/silva_db/silva.db /data/gpfs/projects/punim1528/a_minax/silva_db/SILVA

#use bowtie to map agains silva and keep those sequences that are not in the database (unpaired) 

#!/bin/bash
#SBATCH -N 1 # Número de nodos
#SBATCH -n 2 # Número de núcleos
#SBATCH -t 4-23:00 # Límite de tiempo (D-HH:MM)
#SBATCH -o trim.out # Salida STDOUT
#SBATCH -e trim.err # Salida STDERR
# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=fsalgadoroa@student.unimelb.edu.au

module load Bowtie2/2.4.5

cat /data/gpfs/projects/punim1528/a_minax/reads/list_names.txt | while read ind; do
bowtie2 --quiet --very-sensitive-local \
--phred33  -x /data/gpfs/projects/punim1528/a_minax/silva_db/SILVA -1 /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/trimmed_reads/unfixrm_"$ind"_R1.cor_val_1.fq.gz -2 /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/trimmed_reads/unfixrm_"$ind"_R2.cor_val_2.fq.gz --threads 6 \
--met-file "$ind"_bowtie2_metrics.txt \
--al-conc-gz /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/silva_pair/"$ind"_blacklist_paired_aligned.fq.gz \
--un-conc-gz /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/silva_pair/"$ind"_blacklist_paired_unaligned.fq.gz  \
--al-gz /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/silva_pair/"$ind"_blacklist_unpaired_aligned.fq.gz \
--un-gz /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/silva_pair/"$ind"_blacklist_unpaired_unaligned.fq.gz;
done
```

We are interested in those sequences with the suffix 'blacklist_paired_unaligned.fq.gz' because they did not align with the database, indicating that these reads are free of contamination

## Remove over-represented reads

1. Run Fastqc with the sequence with the suffix 'blacklist_paired_unaligned.fq.gz'

2. Let's extract the data we are intersted in from this Fastqc run

```
mkdir /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/silva_pair/fastqc

for file in $(ls -1 /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/silva_pair/*.gz | xargs -n 1 basename); do

pp=$(echo "$file" | sed "s/.gz//g") 

unzip -p /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/silva_pair/"$pp"_fastqc.zip "$pp"_fastqc/fastqc_data.txt > /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/silva_pair/fastqc/"$pp"_fastqc.txt; done
```

3. Now let's run the script _RemoveFastqcOverrepSequenceReads.py_ to remove the over-represented reads

```
#!/bin/bash
#SBATCH -N 1 # Número de nodos
#SBATCH -n 2 # Número de núcleos
#SBATCH -t 4-23:00 # Límite de tiempo (D-HH:MM)
#SBATCH -o trim.out # Salida STDOUT
#SBATCH -e trim.err # Salida STDERR
# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=fsalgadoroa@student.unimelb.edu.au

cat /data/gpfs/projects/punim1528/a_minax/reads/list_names.txt | while read ind; do
python /data/gpfs/projects/punim1528/a_minax/scripts/TranscriptomeAssemblyTools/utilities/RemoveFastqcOverrepSequenceReads.py \
-1 /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/silva_pair/"$ind"_blacklist_paired_unaligned.fq.1.gz \
-2 /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/silva_pair/"$ind"_blacklist_paired_unaligned.fq.2.gz \
-fql /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/silva_pair/fastqc/"$ind"_blacklist_paired_unaligned.fq.1_fastqc.txt \
-fqr /data/gpfs/projects/punim1528/a_minax/reads/filtered_reads/silva_pair/fastqc/"$ind"_blacklist_paired_unaligned.fq.2_fastqc.txt;
done
```


