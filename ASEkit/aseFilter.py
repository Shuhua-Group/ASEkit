# -*- coding: utf-8 -*-
"""

Created on October 2020

@email: huangke@shangahitech.edu.cn

"""

import pandas as pd
import os
from scipy import stats
import statsmodels.stats.multitest as FDR
import argparse
pd.options.mode.chained_assignment=None
parser = argparse.ArgumentParser()



description= \
        "Description:\n\n " + \
        "Input :           \
        1.  total reads cutoff \
        2. FDR cutff    \
        3. AI (allelic imblance) value cutoff\n" + \
        "The function of this script is to filter ASE loci meeting conditional"


parser.add_argument('Filter',
        help='ASE filter')

parser.add_argument('--rawdir',
        type=str,required=True,
        help='result directory produced by aseCalling')

parser.add_argument('--totalreads',
        dest='totalreads',
        type=int,
        required=True,
        help='exonic variants total reads cutoff,default 10',default=10)

parser.add_argument('--fdr',
        dest='fdr',
        type=float,
        required=True,
        help='Multiple testing binomial pvalue by BH,default 0.05',default=0.05)

parser.add_argument('--AIvalue',
        dest='AIvalue',
        type=float,
        required=True,
        help='allelic imbalance value',default=0.2)
parser.add_argument('--sample',
        dest='sample',
        type=str,
        required=True,
        help='sample info file, two columns with DNA and RNA id')

parser.add_argument('--outdir',
        dest='outdir',
        type=str,
        required=True,
        help='output directory'
        )

args = parser.parse_args()

'''Function to filter exonic SNP,detect ASE loci'''
def aseFilter():
    data=pd.read_table(args.sample,sep='\s+')
    data.columns=['DNAid','RNAid']
    data_ase_list=list()
    total_info_list=list()
    data_AI_list=list()
    os.system('mkdir '+args.outdir)
    
    ## for every sample
    for i in range(0,len(data)):
        samplename=data['RNAid'].iloc[i]
        print(samplename)
        filepath=os.path.join('./',args.rawdir,samplename,samplename+'.allelic_counts.txt')
        data_sample=pd.read_table(filepath)

        ## generate all coding site file including reads info
        data_sample=data_sample.loc[(data_sample['totalCount']>=args.totalreads)]
        data_sample=data_sample.reset_index(drop=True)
        data_total=data_sample
        data_total[samplename+'_ref|alt']=data_total['refCount'].astype(str)+'|'+data_total['altCount'].astype(str)
        total_info_list.append(data_total[['variantID',samplename+'_ref|alt']])
        
        ## generate all ASE site reads info file
        binom_pvalue_list=list()
        data_ase=data_sample
        for i in range(len(data_ase)):
            pvalue=stats.binom_test(data_ase['refCount'][i],n=data_ase['totalCount'][i],p=0.5)
            binom_pvalue_list.append(pvalue)
        data_ase.insert(loc=8,column='binom_pvalue',value=binom_pvalue_list)
        fdr_list=list(FDR.multipletests(binom_pvalue_list,method='fdr_bh')[1])
        data_ase.insert(loc=9,column='FDR',value=fdr_list)
        AI_list=((data_ase['refCount']/data_ase['totalCount'])-0.5).abs()
        data_ase.insert(loc=10,column='AI',value=AI_list)
        data_ase=data_ase.loc[(data_ase['FDR']<args.fdr)& (data_ase['AI']>args.AIvalue)]
        data_ase.to_csv(os.path.join('./',args.rawdir,samplename,samplename+'.filter.txt'),index=None,sep='\t')
        data_ase[samplename+'_ref|alt']=data_ase['refCount'].astype(str)+'|'+data_ase['altCount'].astype(str)
        data_ase_list.append(data_ase[['variantID',samplename+'_ref|alt']])
        data_ase[samplename+'_AI']=AI_list
        data_AI_list.append(data_ase[['variantID',samplename+'_AI']])


    ## merge two type result file
    data_info=total_info_list[0]
    for j in range(1,len(total_info_list)):
        data_info=pd.merge(data_info,total_info_list[j],on=['variantID'],how='outer')
    data_ase=data_ase_list[0]
    for j in range(1,len(data_ase_list)):
        data_ase=pd.merge(data_ase,data_ase_list[j],on=['variantID'],how='outer')
    data_ase_info=pd.merge(data_info,data_ase['variantID'],how='inner')
    
    ## output only ASE site megre file
    data_ase.to_csv(os.path.join(args.outdir,'ase.site.merge.txt'),index=None,sep='\t')
    ## output ASE site and allelic imblance of population level
    data_ase_info.to_csv(os.path.join(args.outdir,'population.ase.txt'),index=None,sep='\t')

    ##output AI value file
    data_AI=data_AI_list[0]
    for n in range(1,len(data_AI_list)):
        data_AI=pd.merge(data_AI,data_AI_list[n],on=['variantID'],how='outer')
    data_AI.to_csv(os.path.join(args.outdir,'population.ase.AI.txt'),index=None,sep='\t')
def main():
    ##run aseFilter
    aseFilter()

if __name__ == '__main__':
    main()


