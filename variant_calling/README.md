# Variant calling

## Supertranscripts

[following this](https://github.com/trinityrnaseq/trinityrnaseq/wiki/SuperTranscripts)

```bash
module load Trinity/2.15.1
module load Salmon/1.5.2

/data/gpfs/projects/punim1528/a_minax/scripts/trinityrnaseq/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py --trinity_fasta centroids.fasta
```

## obtain multiple BAMs

```bash
#!/bin/bash
#SBATCH -N 1 # Número de nodos
#SBATCH -n 30 # Número de núcleos
#SBATCH -t 4-23:00 # Límite de tiempo (D-HH:MM)
#SBATCH --array=1-18%9
#SBATCH -o abundance.out # Salida STDOUT
#SBATCH -e abundance.err # Salida STDERR
#SBATCH --mem=200G
# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=fsalgadoroa@student.unimelb.edu.au

module load BCFtools/1.15.1
module load SAMtools/1.16.1
module load STAR/2.7.10b

readarray -t FILES1 < <(cut -f 3 samples_nospace.txt)
readarray -t FILES2 < <(cut -f 4 samples_nospace.txt)

# Get the file corresponding to the current task ID
FILE1=${FILES1[$SLURM_ARRAY_TASK_ID - 1]}
FILE2=${FILES2[$SLURM_ARRAY_TASK_ID - 1]}
prefix=$(echo $FILE1 | grep -Eo "\w+\_\w+\_\w+")
# Uncompress the file

STAR --runThreadN 15 --genomeDir /data/scratch/projects/punim1528/variant_calling/star_genome_idx  --runMode alignReads  --twopassMode Basic  --alignSJDBoverhangMin 12  --outSAMtype BAM SortedByCoordinate  --limitBAMsortRAM 199117913013  --readFilesIn /data/scratch/projects/punim1528/cat_ind_reads/ungzip/$FILE1 /data/scratch/projects/punim1528/cat_ind_reads/ungzip/$FILE2 --outFileNamePrefix /data/scratch/projects/punim1528/variant_calling/BAMs/$prefix

samtools index /data/scratch/projects/punim1528/variant_calling/BAMs/${prefix}Aligned.sortedByCoord.out.bam

bcftools mpileup -Ou -f /data/scratch/projects/punim1528/assembly_M3_A/trinity_genes.fasta /data/scratch/projects/punim1528/variant_calling/BAMs/${prefix}Aligned.sortedByCoord.out.bam | bcftools call -mv -Ov -o /data/scratch/projects/punim1528/variant_calling/vcfs/$prefix.vcf
```

