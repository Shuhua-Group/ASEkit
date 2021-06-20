# _*_ coding: utf-8 _*_
########## PACKAGES REQUIRED##########
if(require("devtools")){
    print("devtools is loaded correctly")
} else {
    print("trying to install devtools")
    install.packages("devtools",repos = "https://cloud.r-project.org")
    if(require('devtools')){
        print("devtools installed and loaded")
    } else {
        stop("could not install devtools")
    }
}

if(require("getopt")){
    print("getopt is loaded correctly")
} else {
    print("trying to install getopt")
    install.packages("getopt",repos = "https://cloud.r-project.org")
    if(require('getopt')){
        print("getopt installed and loaded")
    } else {
        stop("could not install getopt")
    }
}

if(require("eagle")){
    print("eagle is loaded correctly")
} else {
    print("trying to install eagle")
    devtools::install_github("davidaknowles/eagle",upgrade="never",repos = "https://cloud.r-project.org")
    if(require('eagle')){
        print("eagle installed and loaded")
    } else {
        stop("could not install eagle")
    }
}

########## PACKAGES ##########
library('getopt')
library('eagle')
library('parallel')
########## INPUT ##########
command=matrix(c( 
    'help', 'h', 0,'loical', '显示此帮助信息',
    'ase_result', 'a', 1,'character', 'population.ase file',
    'phenotype_file', 'p', 1, 'character', 'phenotype file',
    'output_prefix', 'o', 1, 'character', 'output prefix'),
    byrow=T, ncol=5
)
args=getopt(command)

if (!is.null(args$help) || is.null(args$ase_result) || is.null(args$phenotype_file) || is.null(args$output_prefix)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}

## file format transform
ase_res_format_trans <- function(ase_res){
    mat_size <- dim(ase_res)
    ref_mat <- matrix(NA,mat_size[1],mat_size[2])
    alt_mat <- matrix(NA,mat_size[1],mat_size[2])
    geno_mat <- matrix(0,mat_size[1],mat_size[2])
    for(i in 1:mat_size[1]){
        for(j in 1:mat_size[2]){
            if(ase_res[i,j] != ''){
                res_temp <- as.numeric(strsplit(ase_res[i,j],'|',fixed = TRUE)[[1]])
                ref_mat[i,j] <- res_temp[1]
                alt_mat[i,j] <- res_temp[2]
                geno_mat[i,j] <- 1
            }
        }
    }
    return(list(ref_mat,alt_mat,geno_mat))
}
########## DATA READ IN ##########

## ase result
#args <- list(ase_result='ase_res.txt',phenotype_file='pheno.txt')
ase_res <- read.table(args$ase_result,header=T,row.names=1,stringsAsFactors=F,sep='\t')
ase_res_temp <- ase_res_format_trans(ase_res)
ref_count <- ase_res_temp[[1]]
alt_count <- ase_res_temp[[2]]
genotype <- ase_res_temp[[3]]
rownames(genotype) <- rownames(ase_res)
heterozygosity <- t(apply(genotype,1,function(x)ifelse(x==1,TRUE, FALSE)))
total_count <- ref_count+alt_count
## phenotype
phenotype <- read.table(args$phenotype_file,header=T,row.names=1,stringsAsFactors=F)


########## MAIN FUNCTION ##########




phenotype_cal <- function(line){
    pheno_now <- line
    pos.del= which(pheno_now > (mean(pheno_now)+3*sd(pheno_now)) | pheno_now < (mean(pheno_now)-3*sd(pheno_now))) ## 离群点
    if(length(pos.del)>0){ 
        environmentVar <- pheno_now[-pos.del]
        alt <- t(alt_count[,-pos.del])
        ref <- t(ref_count[,-pos.del])
        het <- t(heterozygosity[,-pos.del])
        eqtlGenotypes <- t(genotype[,-pos.del])
        totalreads <- t(total_count[,-pos.del])
    }else{
        environmentVar= pheno_now 
        alt <- t(alt_count)
        ref <- t(ref_count)
        het <- t(heterozygosity)
        eqtlGenotypes <- t(genotype)
        totalreads <- t(total_count)
    }
    ########## Default parameters ##########
    count.cutoff=5 # monoallelic if fewer reads than this...
    prop.cutoff=0.01 # ...or lower proportion than this mapping to one allele
    prop.mono.cutoff=0.5 # maximum proportion of hets who are allowed to how monoallelic expression
    minSampleSize=20
    minGroupSize=10 # minimum testable individuals in an environmental group (e.g. smokers)
    minModel=T     
    alt.list <- list()
    n.list <- list()
    x.null <- list()
    x.full <- list()
    original.index <- list()
    # iterate over exonic SNPs
    non.problematic.counter=1
    for (snp.index in 1:dim(eqtlGenotypes)[2]) {
        if (snp.index %% 1000 == 0){
              print(snp.index)}
        valid=het[,snp.index] & totalreads[,snp.index]>5 ## 对每个样本的布尔值 杂合&total reads > 5 
        valid[which(is.na(valid))]=FALSE ## NA设为False
        # check we have at least 20 valid samples 若该位点符合条件的个体数量小于minSampleSize(20)，舍去，进入下一个点
        if (sum(na.omit(valid))<minSampleSize)
            next
        ## 只保留符合条件的个体
        a=alt[valid,snp.index]
        r=ref[valid,snp.index]
        heteq=eqtlGenotypes[valid,snp.index]==1  
        x=environmentVar[valid]
        # check not too many are in one group (e.g. female, non-smokers)
        if ((length(x)-max(table(x)))<minGroupSize) 
            next
        n=a+r
        # check there isn't too much mono-allelic expression 是否太偏向一个allel
        min.a.r=apply(cbind(a,r),1,min)
        is.mono=(min.a.r<count.cutoff)|((as.double(min.a.r)/n)<prop.cutoff)  
        if (mean(is.mono)>prop.mono.cutoff)
            next
        ## 将较少的allel设为alt
        alt.list[[non.problematic.counter]]=if (minModel) pmin(a,r) else a
        n.list[[non.problematic.counter]]=n
        original.index[[non.problematic.counter]]=snp.index
        num.samples=length(x)
        #print(num.samples)  
        x.full[[non.problematic.counter]]=cbind(env=x,intercept=numeric(num.samples)+1.0)
        x.null[[non.problematic.counter]]=cbind(intercept=numeric(num.samples)+1.0)
        # add eQTL heterozygosity
        if (!any(is.na(heteq))) if ((sd(heteq)>0) & (min(table(heteq))>5)){ 
            x.full[[non.problematic.counter]]=cbind(x.full[[non.problematic.counter]],eqtl=heteq)
            x.null[[non.problematic.counter]]=cbind(x.null[[non.problematic.counter]],eqtl=heteq)
            # lots of conditions to decide whether to add an interaction terms
            if (min(table(heteq))>=10) { # only if at least 10 hets and 10 non-hets
                e1=x[heteq==1.0]
                e0=x[heteq==0.0]    
                if ((length(e1)-max(table(e1)))>=10 & (length(e0)-max(table(e0)))>=10){ 
                    if (length(unique(x))>2 | min(table(x,heteq))>=5){ 
                        temp=cbind(x.full[[non.problematic.counter]],interaction=x*heteq)
                        if (det(t(temp) %*% temp) > 1e-8){
                            x.full[[non.problematic.counter]]=temp
                        }
                    }
                }
            }
        }
    non.problematic.counter=non.problematic.counter+1
    }
    original.index=unlist(original.index)
    #---------------- run the model --------------------------
    library('eagle')
    s=eagle.settings()
    s$debug=F
    s$rev.model=2 # local regression
    s$normalised.depth=scale(log10(unlist(lapply(n,mean)))) # only required for rev.model=3
    s$max.iterations=10000
    s$convergence.tolerance=.001
    s$coeff.regulariser=0.1
    s$learn.rev=T
    # learnt parameters from the DGN dataset
    s$rep.global.shape=1.0
    s$rep.global.rate=0.0033
    s$traceEvery=1
    system.time( res <- eagle.helper(alt.list,n.list,x.full,x.null,s) ) 
    res <- t(rbind(colnames(eqtlGenotypes)[original.index],res$p.values))
    return(res)
}

cl <- makeCluster(getOption("cl.cores", length(rownames(phenotype))))
clusterExport(cl,varlist=c('genotype','phenotype','heterozygosity','ref_count','alt_count','total_count'),envir = environment())
result <- parApply(cl, phenotype, 1, phenotype_cal)
for(i in 1:length(rownames(phenotype))){
    res <- data.frame(result[[i]])
    colnames(res) <- c('SNP','p_value')
    res$p_value <- as.numeric(as.character(res$p_value))
    res$FDR <-p.adjust(res$p_value,method = 'BH')
    res <- res[order(res$FDR),]
    output <- paste(args$output_prefix,'/',rownames(phenotype)[i],'.ASE.association.txt',sep='')
    write.table(res,output,quote=F,sep='\t',row.names=F,col.names=T)
}
