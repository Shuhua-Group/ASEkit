# -*- coding: utf-8 -*-
"""

Created on October 2020

@email: huangke@shangahitech.edu.cn

"""

import pandas as pd
import os
import re
import sys
#from multiprocessing import Pool
import threading
import numpy as np
import argparse
import queue


description= \
        "Description:\n\n " + \
        "Input file:           \
        1. VCF files directory  \
        2. RNA fastq files directory\n" + \
        "The function of this script is to call RNA reads of every heterozygosity SNVs\n" + \
        "The method of this script is that: \
        1. STAR mapping using WASP model \
        2. Remove WASP filed reads from bam  \
        3. Phaser caculate allelic reads from bam"




parser = argparse.ArgumentParser()



parser.add_argument('Calling',
        help='call ASE')

parser.add_argument('--sample',
        dest='sample',
        type=str,
        required=True,
        help='samples file that first column is DNAid and second column is RNAid')

parser.add_argument('--rnaseq',
        dest='rnaseq',
        type=str,
        required=True,
        help='RNA fastq directory')

parser.add_argument('--vcf',
        dest='vcf',
        type=str,
        required=True,
        help='VCF file directory')

parser.add_argument('--process',
        dest='process',
        type=int,
        required=True,
        help='Number of parallel processes')

parser.add_argument('--index',
        dest='index',
        type=str,
        required=True,
        help='Reference genome index produced by STAR')
parser.add_argument('--outdir',
        dest='outdir',
        type=str,
        required=True,
        help='output directory')

## optional
parser.add_argument('--STAR',
         dest='STAR',
         type=str,
         default='Null',
         help='The filepath of STAR software; if you can not install STAR successfully')


args = parser.parse_args()

## get phaser.py file path
def phaser_path():
    filepath=os.path.join(os.path.split(os.path.realpath(__file__))[0],'src/phaser.py')
    return filepath
def star_WASP_filter(sample_int):
    data=pd.read_table(args.sample)
    data.columns=['DNA_name','RNA_name']
    RNAseq_file_list=list()
    vcf_file_list=list()
    RNA_samplename=data['RNA_name'][sample_int]
    DNA_samplename=data['DNA_name'][sample_int]
    print(RNA_samplename)
    print(DNA_samplename)
    RNAseq_file=os.listdir(args.rnaseq)
    vcf_file=os.listdir(args.vcf)
    for filename in RNAseq_file:
        RNAseq_file_list.append(os.path.join(args.rnaseq,filename))
    for filename in vcf_file:
        if re.findall('tbi',filename):
            pass
        else:
            vcf_file_list.append(os.path.join(args.vcf,filename))
    RNA_matching = [RNA_filepath for RNA_filepath in RNAseq_file_list if RNA_samplename  in RNA_filepath]
    vcf_matching = [vcf_filepath for vcf_filepath in vcf_file_list if DNA_samplename  in vcf_filepath]


    if  re.findall('.gz',vcf_matching[0]):
        os.system('gunzip -d '+vcf_matching[0])
        vcf_matching[0]=vcf_matching[0][:-3]
    else:
        vcf_matching[0]=vcf_matching[0]
    print(RNA_matching)
    print(vcf_matching)
    os.system('mkdir -p '+os.path.join(args.outdir,RNA_samplename))
    out_put=os.path.join(args.outdir,RNA_samplename,RNA_samplename)
    #index:/picb/humpopg-bigdata5/ningzhilin/RNA_genotype/WASP/hg19
    #STAR2.7:/picb/humpopg-bigdata5/huangke/software/STAR2.7/source/STAR
    if args.STAR=='Null':
        order1= 'STAR  --runThreadN 4 --genomeDir '+args.index +'  \
		         --readFilesIn '+RNA_matching[0] +' '+RNA_matching[1] +'\
                         --twopassMode Basic  --readFilesCommand zcat  --waspOutputMode SAMtag --varVCFfile '+vcf_matching[0]+'  \
                          --outSAMtype BAM SortedByCoordinate  \
			--outFileNamePrefix '+out_put
    else:
       order1= args.STAR+'  --runThreadN 4 --genomeDir '+args.index +'  \
                         --readFilesIn '+RNA_matching[0] +' '+RNA_matching[1] +'\
                         --twopassMode Basic  --readFilesCommand zcat  --waspOutputMode SAMtag --varVCFfile '+vcf_matching[0]+'  \
                          --outSAMtype BAM SortedByCoordinate  \
                        --outFileNamePrefix '+out_put
    order2='samtools view -h '+out_put+'Aligned.sortedByCoord.out.bam | \
            grep -v "vW:i:[2-7]"  | samtools sort > ' + out_put+'Aligned.sortedByCoord.out.WASP.bam'
    order3='samtools index '+out_put+'Aligned.sortedByCoord.out.WASP.bam'
    
    order4='bgzip '+vcf_matching[0]
    order5=' tabix -f '+vcf_matching[0]+'.gz'
    order6=' python2 '+phaser_path() +'  \
            --vcf '+vcf_matching[0]+'.gz' +' \
            --bam '+out_put+'Aligned.sortedByCoord.out.WASP.bam \
         --paired_end 1 --mapq 255 --baseq 10 --threads 4 --sample '+DNA_samplename +' --o '+out_put 
    print(order1)
    os.system(order1)
    os.system(order2)
    os.system(order3)
    os.system(order4)
    os.system(order5)
    os.system(order6)
    sys.stdout.flush()
def main():
    import concurrent.futures
    data_sample=pd.read_table(args.sample)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.process) as executor:
        executor.map(star_WASP_filter, range(len(data_sample)))
if __name__ == "__main__":
    main()
