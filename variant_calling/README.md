# Variant calling

## Supertranscripts

[following this](https://github.com/trinityrnaseq/trinityrnaseq/wiki/SuperTranscripts)

```bash
module load Trinity/2.15.1
module load Salmon/1.5.2

/data/gpfs/projects/punim1528/a_minax/scripts/trinityrnaseq/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py --trinity_fasta centroids.fasta
```
## call SNPs using GATK following pipeline written by trinity authors
```bash
#!/bin/bash
#SBATCH -N 1 # Número de nodos
#SBATCH -n 25 # Número de núcleos
#SBATCH -t 6-23:00 # Límite de tiempo (D-HH:MM)
#SBATCH -o call_gatk.out # Salida STDOUT       
#SBATCH -e call_gatk.err # Salida STDERR       
#SBATCH --array=1-10%5
#SBATCH --mem=200G
# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=fsalgadoroa@student.unimelb.edu.au

module load Trinity/2.15.1
module load BCFtools/1.15.1
module load SAMtools/1.16.1
module load STAR/2.7.10b
module load GATK/4.3.0.0-Java-11
module load Trinity/2.15.1
module load picard/3.0.0-Java-17

export PICARD_HOME=/apps/easybuild-2022/easybuild/software/Core/picard/3.0.0-Java-17
export GATK_HOME=/apps/easybuild-2022/easybuild/software/Compiler/GCCcore/11.3.0/GATK/4.3.0.0-Java-11/

readarray -t FILES1 < <(cut -f 3 /data/scratch/projects/punim1528/variant_calling_m1_m3_a/adults.txt)
readarray -t FILES2 < <(cut -f 4 /data/scratch/projects/punim1528/variant_calling_m1_m3_a/adults.txt)

# Get the file corresponding to the current task ID
FILE1=${FILES1[$SLURM_ARRAY_TASK_ID - 1]}
FILE2=${FILES2[$SLURM_ARRAY_TASK_ID - 1]}
prefix=$(echo $FILE1 | cut -f 8 -d "/" | grep -Eo "\w+\_\w+\_\w+")

mkdir /data/scratch/projects/punim1528/variant_calling_m1_m3_a/gatk_vcfs/$prefix

/data/gpfs/projects/punim1528/a_minax/scripts/trinityrnaseq/Analysis/SuperTranscripts/AllelicVariants/run_variant_calling.py \
           --st_fa /data/scratch/projects/punim1528/variant_calling_m1_m3_a/gatk_vcfs/trinity_genes.fasta \
           --st_gtf /data/scratch/projects/punim1528/variant_calling_m1_m3_a/gatk_vcfs/trinity_genes.gtf \
           -p $FILE1 $FILE2 \
           -t 5 \
           -m 181559536565 \
           -o /data/scratch/projects/punim1528/variant_calling_m1_m3_a/gatk_vcfs/$prefix
```

## Merge VCF files

Because the VCF files were created separately, we need to merge them. First, let's put them all in a single folder with the correct names for the files and samples.

```bash
for folder in $(ls -d G*/ | sed "s/\///g"); do
IFS=$'\n'
cp $folder/filtered_output.vcf vcfs/$folder.vcf
sed -i "s/sample/$folder/g" vcfs/$folder.vcf;done
```

Now let's merge them

```bash
#!/bin/bash
#SBATCH -N 1 # Número de nodos
#SBATCH -n 2 # Número de núcleos
#SBATCH -t 4-23:00 # Límite de tiempo (D-HH:MM)
#SBATCH -o abundance.out # Salida STDOUT
#SBATCH -e abundance.err # Salida STDERR
#SBATCH --mem=50G
# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=fsalgadoroa@student.unimelb.edu.au

module load BCFtools/1.15.1

# Directory containing the VCF files
vcf_dir="/data/scratch/projects/punim1528/variant_calling_m1_m3_a/vcfs"

# Output merged VCF file
merged_vcf="merged_samples.vcf"

# Create an array to hold the file paths
vcf_files=()

# Loop over all VCF files in the directory
for vcf in "$vcf_dir"/*.vcf; do
  #compress 
  bgzip "$vcf"
  vcf="${vcf}.gz"
  # Index each VCF file
  bcftools index "$vcf"
  # Add the file path to the array
  vcf_files+=("$vcf")
done

# Merge the VCF files
bcftools merge -o "$merged_vcf" -O z "${vcf_files[@]}"

```
## obtain multiple BAMs and VCF files

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

# create index
STAR --runThreadN --runThreadN 15 --runMode genomeGenerate --genomeDir /data/scratch/projects/punim1528/variant_calling_m1_m3_a/star_genome_idx --genomeFastaFiles /data/scratch/projects/punim1528/variant_calling_m1_m3_a/trinity_genes.fasta --genomeSAindexNbases 10 --sjdbGTFfile /data/scratch/projects/punim1528/variant_calling_m1_m3_a/trinity_genes.gtf --sjdbOverhang 150  --limitGenomeGenerateRAM 199117913013

readarray -t FILES1 < <(cut -f 3 cat_ind_reads/salmon_m1_m3_a/all_reads.txt)
readarray -t FILES2 < <(cut -f 4 cat_ind_reads/salmon_m1_m3_a/all_reads.txt)

# Get the file corresponding to the current task ID
FILE1=${FILES1[$SLURM_ARRAY_TASK_ID - 1]}
FILE2=${FILES2[$SLURM_ARRAY_TASK_ID - 1]}
prefix=$(echo $FILE1 | cut -f 8 -d "/" | grep -Eo "\w+\_\w+\_\w+")

STAR --runThreadN 15 --genomeDir /data/scratch/projects/punim1528/variant_calling/star_genome_idx  --runMode alignReads  --twopassMode Basic  --alignSJDBoverhangMin 12  --outSAMtype BAM SortedByCoordinate  --limitBAMsortRAM 199117913013  --readFilesIn $FILE1 $FILE2 --outFileNamePrefix /data/scratch/projects/punim1528/variant_calling_m1_m3_a/BAMs/$prefix

samtools index /data/scratch/projects/punim1528/variant_calling_m1_m3_a/BAMs/${prefix}Aligned.sortedByCoord.out.bam

bcftools mpileup -Ou -f /data/scratch/projects/punim1528/variant_calling_m1_m3_a/trinity_genes.fasta /data/scratch/projects/punim1528/variant_calling_m1_m3_a/BAMs/${prefix}Aligned.sortedByCoord.out.bam | bcftools call -mv -Ov -o /data/scratch/projects/punim1528/variant_calling_m1_m3_a/vcfs/$prefix.vcf
```

