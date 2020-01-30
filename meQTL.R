setwd("~/Desktop/methy/")
source("./mxeqtl.R")
library(MatrixEQTL)
methy <- read.table("./GenotypeCNVConrad/Methylation-matrix-GenotypeConrad-Bell.txt",header = T);
methy_expre = methy[,5:ncol(methy)];

gene_mat <-NULL;
for( sl in 1:nrow(methy_expre)) {
  cat(sl,"\t");
  mat = methy_expre[sl,];
  mat = t(apply(mat, 1, rank, ties.method = "average"));
  gene_mat = rbind(gene_mat,mat);
  rm(mat);
}

colnames(gene_mat) = colnames(methy)[5:ncol(methy)]
rownames(gene_mat) = methy[,1]
write.table(gene_mat,"./GenotypeCNVConrad/Methylation-matrix-GenotypeConrad-Bell_mat",quote = F,col.names = T,row.names = T,sep="\t")


me<-mxeqtl('./GenotypeCNVConrad/CNV_mat','./GenotypeCNVConrad/cnvloc.txt','./GenotypeCNVConrad/Methylation-matrix-GenotypeConrad-Bell_mat','./cpgloc.txt',covariates="",cis_output_file='./GenotypeCNVConrad/methy_cis',cis_pval=0.05,
             model="linear", MAF=0, cis_dist=1000000, missing="NA",trans_output_file='./GenotypeCNVConrad/methy_trans',
             trans_pval=1e-5);
me<-mxeqtl('./GenotypeCNV/CNV_mat','./GenotypeCNV/cnvloc.txt','./GenotypeCNV/Methylation-matrix-genotype-Bell_mat','./cpgloc.txt',covariates="",cis_output_file='./GenotypeCNV/methy_cis',cis_pval=0.05,
            model="linear", MAF=0, cis_dist=1000000, missing="NA",trans_output_file='./GenotypeCNV/methy_trans',trans_pval=1e-5);



                                         trans_pval=1e-5);