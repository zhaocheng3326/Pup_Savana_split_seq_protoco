#!/bin/bash
DIR=/home/chenzh/My_project/Savana_split_seq_protocol
TMPD=$DIR/tmp_data
SRC=$DIR/src
DATA=$DIR/data
RE=$DIR/results
DOC=$DIR/doc
BDOC=$DIR/big_doc
BIN=$DIR/bin

cd $DIR



#' smt2 seq mapping
STAR_REF=/home/chenzh/Genome_new/Mouse/RefSeq/clean_chr_refdata/star
GTF=/home/chenzh/Genome_new/Mouse/RefSeq/clean_chr_refdata/genes/genes.gtf
RSEM_REF=/home/chenzh/Genome_new/Mouse/RefSeq/clean_chr_refdata/rsem/rsem
SMT_DIR=$TMPD/smt_data




#SAMPLE=Sample_101C 




for SAMPLE in `cat $BIN/smt.sample.ID`
#for SAMPLE in `cat $BIN/temp.ID`
do
	SAMPLE_FASTQ=$DATA/smt_data/${SAMPLE}.fastq.gz
	SAMPLE_OUT_DIR=${SMT_DIR}/${SAMPLE}
	echo ${SAMPLE}
	mkdir -p ${SAMPLE_OUT_DIR}
	cd ${SAMPLE_OUT_DIR}
	STAR --runThreadN 4 --outFilterMultimapNmax 1 --genomeDir ${STAR_REF} --readFilesIn ${SAMPLE_FASTQ} --quantMode TranscriptomeSAM --outSAMstrandField intronMotif --readFilesCommand zcat --outSAMtype BAM Unsorted  --outFileNamePrefix ${SAMPLE}.
	rsem-calculate-expression --no-bam-output --single-cell-prior --alignments -p 4 ${SAMPLE}.Aligned.toTranscriptome.out.bam ${RSEM_REF} ${SAMPLE}.rsem -q
done


cd ${SMT_DIR}	
FILES=`ls ${SMT_DIR}/Sample_*/*.rsem.genes.results|tr '\n' ' '`
python $SRC/merge_rsem_data.py -i $FILES -o ${SMT_DIR}/merge -f .rsem.genes.results


cd ${SMT_DIR}	
\ls ${SMT_DIR}/*/*Log.final.out > multi.qc.file.list
\ls -d ${SMT_DIR}/*/*.rsem.stat >> multi.qc.file.list
#conda activate Vel
multiqc -f  -dd 1 --flat  -o multiQC -l multi.qc.file.list


#SAMPLE=M_12_82


SMALL_DIR=$TMPD/smallseq_data
N=86 ## 101-8-2-1 -4(adapter length 5) 
cutada=$DOC/cutadapt_3prime.fa ### From paper
RGB=$BDOC/sc_smallRNA_annotation/Mouse/RefSeq/clean_chr_refdata/bowtie/bowtie
REF=$BDOC/sc_smallRNA_annotation/Mouse/RefSeq/clean_chr_refdata/fasta/genome.fa
ANDB=$DOC/small.anno.bed
ANOD=$DOC/annotation.report.order


for SAMPLE in `cat $BIN/smallseq.sample.ID`
do
	SAMPLE_FASTQ=$DATA/smallseq_data/${SAMPLE}.fastq.gz
	SAMPLE_OUT_DIR=${SMALL_DIR}/${SAMPLE}
	echo ${SAMPLE}
	mkdir -p ${SAMPLE_OUT_DIR}; cd ${SAMPLE_OUT_DIR}
	zcat ${SAMPLE_FASTQ} > ${SAMPLE}.fastq




	umi_tools extract --bc-pattern=NNNNNNNN -I ${SAMPLE}.fastq -S ${SAMPLE}.extract.fastq -L ${SAMPLE}.extract.log
	cutadapt -a file:${cutada}   -e 0.1 -o 5 -m 18 -M ${N} -u 2 -o ${SAMPLE}.extract.cutadp.fastq  ${SAMPLE}.extract.fastq  >${SAMPLE}.extract.cutadp.log  ## trim 2bp from begining
	\rm ${SAMPLE}.fastq ${SAMPLE}.extract.fastq


#' mapping
	bowtie  -a --best --strata -v 2 -m 50  -S -q -p 2 $RGB ${SAMPLE}.extract.cutadp.fastq 2>bowie.log |python $SRC/filter_sam_mismatch.py|samtools view -Sb -F 4 > ${SAMPLE}.bowtie.bam 
	samtools sort ${SAMPLE}.bowtie.bam  -o ${SAMPLE}.bowtie.sort.bam
	\rm ${SAMPLE}.bowtie.bam   ${SAMPLE}.extract.cutadp.fastq 


	#' deduplication
	samtools index ${SAMPLE}.bowtie.sort.bam
	umi_tools dedup --method directional --output-stats dedup.log1  --read-length -I ${SAMPLE}.bowtie.sort.bam -S ${SAMPLE}.bowtie.sort.dedup.temp.bam --random-seed 123 > /dev/null
	samtools sort ${SAMPLE}.bowtie.sort.dedup.temp.bam -o ${SAMPLE}.bowtie.sort.dedup.bam
	picard SamToFastq INPUT=${SAMPLE}.bowtie.sort.dedup.bam FASTQ=${SAMPLE}.bowtie.sort.dedup.bam.fastq.temp  2>> log.txt
	seqkit fx2tab ${SAMPLE}.bowtie.sort.dedup.bam.fastq.temp|sort |uniq |seqkit tab2fx -w 0 > ${SAMPLE}.bowtie.sort.dedup.bam.fastq

#' remapping
	bowtie  -a --best --strata -v 2 -m 50  -S -q -p 2 $RGB ${SAMPLE}.bowtie.sort.dedup.bam.fastq 2>bowie.dedup.log |python $SRC/filter_sam_mismatch.py|samtools view -Sb -F 4 > ${SAMPLE}.bowtie.dedup.reMap.bam  ### recover the multiple mapping results
	samtools sort ${SAMPLE}.bowtie.dedup.reMap.bam -o ${SAMPLE}.bowtie.dedup.reMap.sort.bam
	\rm ${SAMPLE}.bowtie.sort.bam ${SAMPLE}.bowtie.sort.dedup.temp.bam ${SAMPLE}.bowtie.sort.dedup.bam.fastq.temp dedup.log1_edit_distance.tsv dedup.log1_per_umi_per_position.tsv dedup.log1_per_umi.tsv ${SAMPLE}.bowtie.dedup.reMap.bam


#' featureCounting
	intersectBed -a <(bam2bed <${SAMPLE}.bowtie.dedup.reMap.sort.bam|awk -F "\t" '{OFS="\t";if (($3-$2) <  200) print $1,int(($2+$3)/2)-1,int(($2+$3)/2),$4}') -b ${ANDB} -wa -wb |cut -f 4,11,12 |sort |uniq > ${SAMPLE}_tmpCount.detail.out
	python $SRC/count_smallRNA.ordered.exp.py -a ${ANOD} -c ${SAMPLE}_tmpCount.detail.out -o ${SAMPLE}


done




