# Differential expression analyses

Create a file with the experimental design a names of the files: 

```bash
linesw=$(ls /data/scratch/projects/punim1528/cat_ind_reads | sed -E "s/\.\w{1}.gz//g" | sort | uniq | grep -E "B\w{1}$" | wc -l)
for s in $(seq 1 $linesw); do  ind=$(ls /data/scratch/projects/punim1528/cat_ind_reads | sed -E "s/\.\w{1}.gz//g" | sort | uniq | grep -E "B\w{1}$" | awk -v s="$s" 'NR==s'); morph="black"; echo -e "$morph\t${morph}_$s\t$ind.1.gz\t$ind.2.gz\n" >> samples.txt;done
```
**If by any chance the reads are in different files, merged them in a single file with the following line:**
```bash
#!/bin/bash
#SBATCH -N 1 # Número de nodos
#SBATCH -n 7 # Número de núcleos
#SBATCH -t 1-23:00 # Límite de tiempo (D-HH:MM)
#SBATCH -o joint.out # Salida STDOUT
#SBATCH -e joint.err # Salida STDERR
#SBATCH --array=1-14%7
#SBATCH --mem=50G
# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=fsalgadoroa@student.unimelb.edu.au

readarray -t ind < <(cut -f 3 /data/scratch/projects/punim1528/processed_reads/molt3_v2/codes_unique.txt)

# Get the file corresponding to the current task ID
ind=${ind[$SLURM_ARRAY_TASK_ID - 1]}

seq1=$(ls /data/scratch/projects/punim1528/processed_reads/molt3_v2/clean/*.1 | grep $ind)
seq2=$(ls /data/scratch/projects/punim1528/processed_reads/molt3_v2/clean/*.2 | grep $ind)
cat $seq1 > /data/scratch/projects/punim1528/processed_reads/molt3_v2/joint/$ind.fq.1
gzip /data/scratch/projects/punim1528/processed_reads/molt3_v2/joint/$ind.fq.1
cat $seq2 > /data/scratch/projects/punim1528/processed_reads/molt3_v2/joint/$ind.fq.2
gzip /data/scratch/projects/punim1528/processed_reads/molt3_v2/joint/$ind.fq.2
```


## Map and count transcripts with Salmon

```bash
#!/bin/bash
#SBATCH -N 1 # Número de nodos
#SBATCH -n 2 # Número de núcleos
#SBATCH -t 4-23:00 # Límite de tiempo (D-HH:MM)       
#SBATCH -o abundance.out # Salida STDOUT
#SBATCH -e abundance.err # Salida STDERR
# mail alert at start, end and abortion of execution  
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=fsalgadoroa@student.unimelb.edu.au

module load Perl/5.34.1
module load Trinity/2.15.1
module load Bowtie2/2.4.5
#module load bioinfo/express/1.5.1
module load kallisto/0.48.0
module load RSEM/1.3.3
module load Salmon/1.5.2
module load SAMtools/1.16.1


#sheck format samples.txt
# salmon
/data/gpfs/projects/punim1528/a_minax/scripts/align_and_estimate_abundance.pl --transcripts /data/scratch/projects/punim1528/assembly_M1_M3_A/nonredundant_M1_M3_A.fasta --seqType fq --samples_file /data/scratch/projects/punim1528/cat_ind_reads/salmon_m1_m3_a/molt3v2.txt --est_method salmon --coordsort_bam --trinity_mode --prep_reference --output_dir /data/scratch/projects/punim1528/cat_ind_reads/salmon_m1_m3_a
```



