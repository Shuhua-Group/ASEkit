sample=ERR188035
samtools view -hb /picb/humpopg-bigdata5/huangke/ASE_pipeline_res/GBR68/WASP.bam/${sample}/${sample}Aligned.sortedByCoord.out.WASP.bam -L region.bed | samtools sort -n >  ${sample}.Chr22.bam
bedtools bamtofastq  -i ${sample}.Chr22.bam -fq ${sample}.Chr22.R.fastq -fq2 ${sample}.Chr22.L.fastq
bgzip ${sample}.Chr22.R.fastq
bgzip  ${sample}.Chr22.L.fastq

