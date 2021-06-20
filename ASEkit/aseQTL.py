import argparse
import pandas as pd
import sys
import gzip as gz
import os
import scipy.stats as stats
import statsmodels.stats.multitest
from multiprocessing import Pool
parser = argparse.ArgumentParser()

parser.add_argument('aseQTL',
                help='Compare candidate regulatory effects of homozygous and heterozygous SNP')

parser.add_argument("--ase",
        dest='ase',
        type=str,
        required=True,
        help='population information level ase file produced by aseFilter')
parser.add_argument("--vcf",
        dest='vcf',
        type=str,
        required=True,
        help='vcf filepath file according by chromosome')

parser.add_argument("--hetNumber",
        dest='hetNumber',
        type=int,
        default=8,
        required=True,
        help='minimum heterozygous exonic locus number to run aseQTL')

parser.add_argument("--cisRegion",
        dest='cisRegion',
        type=int,
        required=True,
        default='100000',
        help=' SNP within this distance are considered candidate regultory SNP')


parser.add_argument('--process',
        type=int,
        required=False,
        default=3,
       help='multiprocessing according chromosome')

parser.add_argument("--outdir",
        type=str,
        required=False,
        help="aseQTL result file output")
args = parser.parse_args()


''' Function to detect vcf filepath file format'''
def vcf_file():
    data=pd.read_table(args.vcf,sep='\s+',dtype=str)
    data.columns=['chr','filepath']
    ## whether include 'chr'
    if len(data['chr'][0])>1:
        data['chr']=data['chr'].str[3:]
    else:
        pass
    return data

'''Function to count the number of '#' line '''
def vcf_skip_line(geno_file):
    i=0
    file_data=gz.open(geno_file,'r')
    for line in file_data:
         line=line.decode()
         i=i+1
         if line[0:2]=='##':
             i=i+1
         else:
             break
    file_data.close()

    return i
'''Function to produce aseQTL format file by genes'''
def aseQTLformat(rank):
    
    vcf_filepath=vcf_file()
    geno_file=vcf_filepath['filepath'][rank]
    Chr=vcf_filepath['chr'][rank]
    #print('ase filepath',args.ase)
    print('caculate vcf filepath'+'\t',geno_file)

    ## remain at least (hetNumber) individuals heterozygous in ASE locus
    ## add chromosome and position information
    data_ase=pd.read_table(args.ase,dtype=str)
    chr_list=list(data_ase['variantID'].str.split('_').str[0])
    pos_list=list(data_ase['variantID'].str.split('_').str[1])
    data_ase.insert(loc=1,column='Chr',value=chr_list)
    data_ase.insert(loc=2,column='POS',value=pos_list)
    data_ase=data_ase.loc[data_ase['Chr']==str(Chr)]
    data_ase['sum']=data_ase.iloc[:,3:].notnull().sum(axis=1)
    data_ase=data_ase.loc[data_ase['sum']>=args.hetNumber]
    data_ase=data_ase.drop(['sum'],axis=1)
    
    ## replace genotype by 0,1,2 format
    data_geno=pd.read_table(geno_file,dtype=str,skiprows=int((vcf_skip_line(geno_file)-1)/2))
    if len(data_geno.iloc[0][10])==3:
        if data_geno.iloc[0][10][1]=='/':
            data_geno=data_geno.replace({'0/0':'0','0/1':'1','1/0':'1','1/1':'2'})
        else:
            data_geno=data_geno.replace({'0|0':'0','0|1':'1','1|0':'1','1|1':'2'})
    data_geno['POS']=data_geno['POS'].astype(int)
    os.system('mkdir -p '+os.path.join('./aseQTL.format.file/',str(Chr)))
    
    for i in range(len(data_ase)):
        POS=int(data_ase['POS'].iloc[i])
        ## upstream and downstream cis region
        ## output format:
        ## AI   SNP1    SNP2    SNP3

        ## cis region SNP 
        data_ase_geno=data_geno.loc[(data_geno['POS']>=POS-args.cisRegion)&(data_geno['POS']<=POS+args.cisRegion)]
        data_ase_format=data_ase_geno.T.iloc[9:]
        data_ase_format.columns=list(data_ase_geno.T.iloc[1])
        data_ase_format.insert(loc=0,column='allele imblance',value=list(data_ase.iloc[i][3:]))
        data_ase_format=data_ase_format[data_ase_format['allele imblance'].notnull()]
        
        ##caulcate AI value:|(ref/ref+alt)-0.5|
        AI=(data_ase_format['allele imblance'].str.split('|').str[0].astype(int)/(data_ase_format['allele imblance'].str.split('|').str[1].astype(int)+data_ase_format['allele imblance'].str.split('|').str[0].astype(int))-0.5).abs()
        data_ase_format.insert(loc=1,column='AI_value',value=list(AI))
        Chr=str(Chr)
        ## produce result file by genes
        data_ase_format.to_csv(os.path.join('./aseQTL.format.file/',Chr,Chr+'_'+str(POS)+'.format.txt'),index=None,sep='\t')



'''Function to get produce chromosome number by aseFormat function'''
def chr_number():
    j=0
    for dir_number in os.walk('./aseQTL.format.file/'):
        j=j+1
        if j==1:
            return dir_number[1]
        else:
            break


''' Function to caculate local FDR'''
def aseQTL_local_fdr(df):
    data_ase=df
    fdr_list=multitest.multipletests(data_ase['pvalue'],method='fdr_bh')[1]
    data_ase['fdr']=fdr_list
    return data_ase
    #data_aseQTL.to_csv('ASEqtl.all.Chr.res.txt',index=None,sep='\t')

''' Function to run aseQTL association'''
def aseQTL(Chr):
    Chr=str(Chr)
    df_list=list()
    i=0
    print('run aseQTL Chr'+Chr)
    
    for filename in os.listdir(os.path.join('./aseQTL.format.file/',Chr)):
        filepath=os.path.join('./aseQTL.format.file/',Chr,filename)
        ASE_pos=filepath.split('/')[-1].split('.')[0]
        SNP_list=list()
        ASE_site_list=list()
        pvalue_list=list()
        data=pd.read_table(filepath)

        ## skip all AI value is 0.5
        if len(data.drop_duplicates(['AI_value']))==1:
            continue
        
        ##  skip no SNP 
        if len(data.columns)==2:
            continue
        
        ## split cadidate SNPs into two groups: 'heterozygous and homozygous'
        for SNP in list(data.columns[2:]):
            corr_data=data[['AI_value',SNP]]
            het_corr_data=corr_data.loc[corr_data[SNP]==1]    
            
            ## drop homo SNP
            if len(het_corr_data)==0:
                continue
            homo_corr_data=corr_data.loc[corr_data[SNP]!=1]
            if len(homo_corr_data)==0:
                continue

            ## mannwhiteyu statistic single side test
            pvalue=stats.mannwhitneyu(het_corr_data['AI_value'],homo_corr_data['AI_value'],alternative='greater')[1]
            SNP_list.append(Chr+':'+SNP)
            ASE_site_list.append(ASE_pos)
            pvalue_list.append(pvalue)
        df=pd.DataFrame({'ASE_site':ASE_site_list,'SNP':SNP_list,'pvalue':pvalue_list})
        fdr_list=statsmodels.stats.multitest.multipletests(df['pvalue'],method='fdr_bh')[1]
        df['fdr']=fdr_list
        df_list.append(df)
    
    ## concat result and produce result file
    df_total=pd.concat(df_list)
    df_total['Chr']=Chr
    
    ##local FDR
    col=['ASE_site','Chr','SNP','pvalue','fdr']
    df_total.to_csv(os.path.join(args.outdir,'ASEqtl.Chr'+Chr+'.res.txt'),index=None,sep='\t',columns=col)

'''Function to run aseQTLformat in parallel'''
def pool_aseQTLfromat():
    vcf_filepath=pd.read_table(args.vcf,sep='\s+')
    p=Pool(args.process)
    for rank in range(len(vcf_filepath)):
        p.apply_async(aseQTLformat,(rank,))
    p.close()
    p.join()

'''Function to run aseQTL  in parallel'''
def pool_aseQTL():
    p=Pool(args.process)
    chr_list=chr_number()
    for Chr in chr_list:
        p.apply_async(aseQTL,(Chr,))
    p.close()
    p.join()

def main():
    ## mkdir 
    if  os.path.exists(args.outdir)==False:
        try:
            os.system('mkdir '+args.outdir)
        except:
            pass
    else:
        pass
    
    ##run aseQTLformat in parallel
    pool_aseQTLfromat()
    ##run aseQTL in parallel
    pool_aseQTL()

if __name__ == "__main__":
    main()

