sample=HG00155
#samtools view -hb /picb/humpopg-bigdata5/huangke/ASE_pipeline_res/GBR68/WASP.bam/${sample}/${sample}Aligned.sortedByCoord.out.WASP.bam -L region.bed > ${sample}.Chr22.bam
#bamToFastq -i ${sample}.Chr22.bam -fq ${sample}.Chr22.fastq
bcftools view -R ./region.bed /picb/humpopg-bigdata5/huangke/ASE_pipeline_res/GBR68/GBR68.chr.vcf/${sample}/${sample}.het.vcf.gz | bgzip > ${sample}.Chr22.vcf.gz
