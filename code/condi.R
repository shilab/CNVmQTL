setwd("~/condiCNVmQTL");
library(cpgen);
library(data.table);

cnvmQTL <- read.table("1KGgenotype-Conradgenotype/mQTL_CIS_TRANS_P100_twoDatasets.nominalPvalue.overlapRemoved.txt",header = T);
snpmQTL <- read.table("LDAnalysis/SNPmqtls-Bell-Zhang-Banovich.txt",header =F, sep="\t")
SNP_init <- as.data.frame(read.table("LDAnalysis/new.genotype.GT.update.txt", header=T, sep="\t"));
SNP <- SNP_init[complete.cases(SNP_init),];
snpmQTL_sort = snpmQTL[order( snpmQTL[,2], snpmQTL[,3] ),];

CNV_geno_KG <- read.table("GenotypeCNV/CNV-matrix-diffGenotypeKG-Bell.txt",header = T);
CNV_geno_Conrad <- read.table("GenotypeCNVConrad/CNV-matrix-GenotypeConrad-Bell.txt",header = T);

methy_KG <- read.table("GenotypeCNV/Methylation-matrix-genotype-Bell.txt",header = T);
methy_Conrad <- read.table("GenotypeCNVConrad/Methylation-matrix-GenotypeConrad-Bell.txt",header = T)

res.s <- NULL;
window=1000000;
for(i in 1:nrow(cnvmQTL)){
    cat(i,"\t");
    if (cnvmQTL$Source[i] == "1000GPgenotype"){
        tmp_geno = CNV_geno_KG[CNV_geno_KG$Clone %in% paste0(cnvmQTL[i,1],":",cnvmQTL[i,2],"-",cnvmQTL[i,3]),5:ncol(CNV_geno_KG)];
        tmp_start = as.numeric(as.matrix(cnvmQTL[i,2])) - 1000000;
        tmp_end = as.numeric(as.matrix(cnvmQTL[i,3])) + 1000000;

        tmp_snp = snpmQTL_sort[which(snpmQTL_sort$V2 == cnvmQTL[i,1]),];
        tmp_index = tmp_snp[intersect(which(tmp_snp$V3 <= tmp_end),which(tmp_snp$V3 >= tmp_start)),];
        tmp_SNP_mat = SNP[which(SNP$ID %in% tmp_index$V1),4:ncol(SNP)];

        tmp_methy = methy_KG[which(methy_KG$TargetID %in% cnvmQTL$CpG[i]),5:ncol(methy_KG)];

        comm_sample = intersect(intersect(colnames(tmp_geno),colnames(tmp_SNP_mat)),colnames(tmp_methy));

        X = tmp_SNP_mat[,which(colnames(tmp_SNP_mat) %in% comm_sample)];
        M = tmp_geno[,colnames(tmp_geno) %in% comm_sample];
        y = tmp_methy[,colnames(tmp_methy) %in% comm_sample];       

        sum(colnames(X) == colnames(M));
        sum(colnames(M) == colnames(y));
        
        X <- t(as.matrix(X));
        X <- matrix(as.numeric(X),nrow = length(comm_sample),byrow=F);
        M <- t(as.matrix(M));
        M <- as.matrix(as.numeric(M));
        y = as.numeric(matrix(y,ncol=1));

        if (length(unique(sort(M))) > 1){
          # estimates = cGWAS.emmax(y, M=M, X= X, A = NULL,verbose = TRUE)
          estimates = cGWAS(y, M=M, X= X,verbose = TRUE)
          result = cbind(as.matrix(cnvmQTL$Source[i]), estimates$beta, estimates$se, estimates$p_value)
          colnames(result) = c("CNVsource", "beta", "se", "p_value")
          res.s = rbind(res.s,result)
          rm("result");
        }
        if (length(unique(sort(M))) == 1){
            result = cbind(as.matrix(cnvmQTL$Source[i]), "NA", "NA", "NA")
            colnames(result) = c("CNVsource", "beta", "se", "p_value")
            res.s = rbind(res.s,result)
            rm("result")
        }
    }
    
    if (cnvmQTL$Source[i] == "ConradGenotype"){
        tmp_geno = CNV_geno_Conrad[CNV_geno_Conrad$Clone %in% paste0(cnvmQTL[i,1],":",cnvmQTL[i,2],"-",cnvmQTL[i,3]),5:ncol(CNV_geno_Conrad)];
        tmp_start = as.numeric(as.matrix(cnvmQTL[i,2])) - 1000000;
        tmp_end = as.numeric(as.matrix(cnvmQTL[i,3])) + 1000000;

        tmp_snp = snpmQTL_sort[which(snpmQTL_sort$V2 == cnvmQTL[i,1]),];
        tmp_index = tmp_snp[intersect(which(tmp_snp$V3 <= tmp_end),which(tmp_snp$V3 >= tmp_start)),];
        tmp_SNP_mat = SNP[which(SNP$ID %in% tmp_index$V1),4:ncol(SNP)];

        tmp_methy = methy_Conrad[which(methy_Conrad$TargetID %in% cnvmQTL$CpG[i]),5:ncol(methy_Conrad)];

        comm_sample = intersect(intersect(colnames(tmp_geno),colnames(tmp_SNP_mat)),colnames(tmp_methy));

        X = tmp_SNP_mat[,colnames(tmp_SNP_mat) %in% comm_sample];
        M = tmp_geno[,colnames(tmp_geno) %in% comm_sample];       
        y = tmp_methy[,colnames(tmp_methy) %in% comm_sample];        

        sum(colnames(X) == colnames(M));
        sum(colnames(M) == colnames(y));

        X <- t(as.matrix(X));
        X <- matrix(as.numeric(X),nrow = length(comm_sample),byrow=F);
        M <- t(as.matrix(M));
        M <- as.matrix(as.numeric(M));
        y = as.numeric(matrix(y,ncol=1));

        if (length(unique(sort(M))) > 1){
            # estimates = cGWAS.emmax(y, M=M, X= X, A = NULL, verbose = TRUE);
            estimates = cGWAS(y, M=M, X= X,verbose = TRUE)
            result = cbind(as.matrix(cnvmQTL$Source[i]), estimates$beta, estimates$se, estimates$p_value);
            colnames(result) = c("CNVsource", "beta", "se", "p_value");
            res.s = rbind(res.s,result);
            rm("result");
        }
        if (length(unique(sort(M))) == 1){
            result = cbind(as.matrix(cnvmQTL$Source[i]), "NA", "NA", "NA");
            colnames(result) = c("CNVsource", "beta", "se", "p_value");
            res.s = rbind(res.s,result);
            rm("result");
        }
        
    }
}

write.table(res.s,"conditional_res_cGWAS",quote = F,col.names = T,row.names = F,sep="\t");
