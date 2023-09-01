
# Ensamblaje 

El ensamblaje se corre con Trinity:

```
module load conda
source activate /home/fabianc.salgado/data/POCtrinity
module load perl/5.32.0

/datacnmat01/biologia/biologia.evolutiva/fabianc.salgado/POCtrinity/opt/trinity-2.9.1/Trinity \
 --seqType fq --SS_lib_type RF --max_memory 225G --CPU 32 --min_contig_length 2000 \ 
 --left /home/fabianc.salgado/shared/paula_torres/gasteracantha/filter_reads/remove_overrep/*.left.fa \ 
 --right /home/fabianc.salgado/shared/paula_torres/gasteracantha/filter_reads/remove_overrep/*.right.fa \
 --output /home/fabianc.salgado/shared/paula_torres/gasteracantha/trinity/trinity_without_2000
```

## Calculo de estadísticas del ensamblaje

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
