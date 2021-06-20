ASEkit Calling \
   --sample ./sample.info.txt \
        --rnaseq ./RNA.test.data \
        --vcf ./vcf.test.data \
        --index ~/hg19 \
		--STAR ~/STAR2.7/source/STAR \
        --process 4 \
		--outdir ./out_res

#example
# --index /data/hg19  ( reference genome index produced by STAR)
# --STAR /software/STAR2.7/source/STAR  (STAR software filepath if you can not install STAR)
