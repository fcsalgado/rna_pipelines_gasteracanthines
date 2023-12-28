# Filtrado de secuencias 

Inicialmente es necesario filtrar los reads porque pueden tener adaptadores y demás, el protocolo 
de filtrado utilizado en este caso se basa en uno publicado por Harvard que se encuentra [aqui](https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html).

## Examinar calidad de las secuencias (FastQC)
Lo primero que se hace con los reads crudos es correr FastQC para verificar la calidad de las secuencias, en el cluster se corre así:

```

module load fastqc/0.11.7

fastqc /home/reads/*.gz
```

Luego de analizar la calidad de las lecturas procedemos a filtrar las secuencias.

## Remover los k-mers erróneos

Debido a que los ensambladores de transcriptomas de última generación utilizan construyen el ensamblaje con un enfoque de grafos a partir de k-mers, los k-mers erróneos pueden afectar negativamente a los ensamblajes, se recomienda quitarlos. 
Para quitarlos utilizamos [rCorrector](https://github.com/mourisl/Rcorrector), una herramienta diseñada específicamente para datos de RNA-seq. 
Además de tener un rendimiento superior en comparaciones en paralelo de herramientas de corrección de errores de lectura basadas en kmer, 
incluye etiquetas en la salida fastq que indican si la lectura se ha corregido o si se ha detectado que contiene un error, pero no se puede corregir.

Creamos un directorio para correr rCorrector, en este caso como las secuencias tienen un código de números, lo más eficiente es correrlo 
con un ciclo for de 1 a 9 donde los nombres de las secuencias son van del PTUR001 al PTUR009.

```
#!/bin/bash

#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -t 20-12:30:30
#SBATCH -o corrector.out
#SBATCH -e error_corrector.err

module load perl

for i in $(seq 1 9)
do
perl /opt/ohpc/pub/apps/rcorrector/run_rcorrector.pl -t 12 -1 \
/home/fabianc.salgado/shared/paula_torres/gasteracantha/secuencias/PTUR00$i.left.fq \
-2 /home/fabianc.salgado/shared/paula_torres/gasteracantha/secuencias/PTUR00$i.right.fq

###########alternative based on other names###############

# in some cases companies run samples in m ultiple lanes, with long IDs, in this case you can create a list of files with following command replacing the regular expression depending oin the situation

ls *.gz | sed -E "s/\_\R[1,2]\.\w+\.\w+$//g" | sort | uniq > list_names.txt

# then the above loop would change for

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

Los archivos fastq de salida incluirán "cor" en sus nombres. Las lecturas corregidas tendrán un sufijo "cor" en sus etiquetas:

```
@ERR1101637.57624042/1 l:165 m:203 h:218 cor
GTACAACCCTTCCAACCTCCACCGTCTTATATACGAAGCGCCTTGAGTGTGTGTGTGCATGAGCCAAAGGGAATACCG
+
<@@D=D2ACDAHB?:<8EGE;B;FE<?;??DBDD))0)8@GDF@C)))5-;5CAE=..7;@DEEC@C;A9?BB=@>@B
```
Los valores l y h en los encabezados de lectura indican los recuentos de kmer , para los kmers de frecuencia más baja, 
mediana y más alta que ocurren en la lectura. También van a aparacer unas lecturas que no se pueden corregir.

Estas lecturas que no se pueden corregir a menudo están plagadas de N o representan otras secuencias de baja complejidad. 
Hay un script en Python para realizar esta tarea. También elimina la etiqueta "cor" de los encabezados de las secuencias corregidas, ya que pueden causar problemas a las herramientas posteriores:

```
module load python

for i in $(seq 1 9)
do
python /opt/ohpc/pub/apps/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py \ 
-1 /home/fabianc.salgado/shared/paula_torres/gasteracantha/filter_reads/rcorrector/PTUR00$i.left.cor.fq \
-2 /home/fabianc.salgado/shared/paula_torres/gasteracantha/filter_reads/rcorrector/PTUR00$i.right.cor.fq -s PTUR00$i
done
```

## Quitamos adaptadores

```
module load python3.6

for i in $(seq 1 9)
do
/home/fabianc.salgado/shared/paula_torres/gasteracantha/TrimGalore-0.6.0/trim_galore \
 --paired --retain_unpaired --phred33 --output_dir trimmed_reads_PTUR00$i --length 36 -q 5 \
 --stringency 1 -e 0.1 /home/fabianc.salgado/shared/paula_torres/gasteracantha/filter_reads/rcorrector/PTUR00$i.left.cor.fq \
   /home/fabianc.salgado/shared/paula_torres/gasteracantha/filter_reads/rcorrector/PTUR00$i.right.cor.fq 
done
```

## Removemos secuencias sobre-representadas 

```
module load python

for i in $(seq 1 9)
do
python /opt/ohpc/pub/apps/TranscriptomeAssemblyTools/RemoveFastqcOverrepSequenceReads.py \
-1 /home/fabianc.salgado/shared/paula_torres/gasteracantha/filter_reads/silva/reads_assembly/blacklist_paired_unaligned_PTUR00$i.fq.1.gz \
-2 /home/fabianc.salgado/shared/paula_torres/gasteracantha/filter_reads/silva/reads_assembly/blacklist_paired_unaligned_PTUR00$i.fq.2.gz \
-fql /home/fabianc.salgado/shared/paula_torres/gasteracantha/filter_reads/fastqc/fastqc_files/PTUR00"$i"_fastqc_data_1.txt \
-fqr /home/fabianc.salgado/shared/paula_torres/gasteracantha/filter_reads/fastqc/fastqc_files/PTUR00"$i"_fastqc_data_2.txt
done
```

Volvemos a correr fastQc
