# Mapping the reads

In case you get reads from multiple sequencing lanes you have to consolidate them in a single file per individual

```bash
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

Because our experimental design consisted in multiple biological replicates, we need to create a file with this information. In our case the last part of the name of the individuals `"B\w{1}$"` has information about the morph. let's create the file samples.txt followwing the format `cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq` for the black morph

```bash
# Count the number of lines matching the pattern 'B\w{1}$' in the directory listing.
linesw=$(ls /data/scratch/projects/punim1528/cat_ind_reads | sed -E "s/\.\w{1}.gz//g" | sort | uniq | grep -E "B\w{1}$" | wc -l)

# Iterate over a sequence of numbers from 1 to the number of lines counted.
for s in $(seq 1 $linesw); do 

    # Extract a specific line that matches the pattern 'B\w{1}$' from the directory listing.
    ind=$(ls /data/scratch/projects/punim1528/cat_ind_reads | sed -E "s/\.\w{1}.gz//g" | sort | uniq | grep -E "B\w{1}$" | awk -v s="$s" 'NR==s')

    # Define the variable 'morph' as "black".
    morph="black"

    # Append formatted output to the file 'samples.txt', including values of 'morph', 's', 'ind.1.gz', and 'ind.2.gz'.
    echo -e "$morph\t${morph}_$s\t$ind.1.gz\t$ind.2.gz\n" >> samples.txt;

done
```

Let's map the reads to the the tracscritome using Salmon

```
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
/data/gpfs/projects/punim1528/a_minax/scripts/align_and_estimate_abundance.pl --transcripts /data/scratch/projects/punim1528/assembly_300/trinity_cdhit_aminax.fasta --seqType fq --samples_file samples.txt --est_method salmon  --coordsort_bam --trinity_mode --prep_reference --output_dir /data/scratch/projects/punim1528
```
