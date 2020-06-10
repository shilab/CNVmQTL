setwd("~/HiC/");
#######Step 1: Convert the Raw data to normalized data:
chr1_raw <- read.table("chr1_10kb.RAWobserved",header=F);
chr1_raw <- as.matrix(chr1_raw);

chr1_KRnorm <- read.table("chr1_10kb.KRnorm",header=F);
chr1_KRnorm <- as.matrix(chr1_KRnorm);
resolution = 10000;

chr1_intensity <- matrix(".",nrow=nrow(chr1_raw),ncol = 6);
for(k in 1:nrow(chr1_raw)){
  print(k);
  i = chr1_raw[k,1]/resolution + 1;
  j = chr1_raw[k,2]/resolution + 1;
  if (chr1_KRnorm[i] != "NaN" && chr1_KRnorm[j] != "NaN"){
    left = c(chr1_raw[k,1],(chr1_raw[k,1] + resolution));
    right = c(chr1_raw[k,2],(chr1_raw[k,2] + resolution));
    intensity = chr1_raw[k,3]/(chr1_KRnorm[i]*chr1_KRnorm[j]);
    subline = c("chr1",left,right,intensity);
    chr1_intensity[k,] = subline;
  } 
  else{
    left = c(chr1_raw[k,1],(chr1_raw[k,1] + resolution));
    right = c(chr1_raw[k,2],(chr1_raw[k,2] + resolution));
    subline = c("chr1",left,right,"NaN");
    chr1_intensity[k,] = subline;
  }
}
write.table(chr1_intensity,"chr1_intensity.txt",quote = F,col.names = F,row.names = F,sep="\t");\


###Step 2: Extract the CpG NCBI36 coordinate of CNV-meQTL:
cpg_meqtl_id <- read.table("cpg_meQTL_id",header = F);
cpg_meqtl_id <- as.matrix(cpg_meqtl_id);

cpg_posi <- read.table("cpg_posi",header = F);
cpg_posi <- as.matrix(cpg_posi);

cpg_meqtl_posi = matrix(".",nrow = nrow(cpg_meqtl_id),ncol = ncol(cpg_posi));
for(i in 1:nrow(cpg_meqtl_id)){
  cpg_meqtl_posi[i,] = cpg_posi[which(cpg_meqtl_id[i,1] == cpg_posi[,1]),];
}
write.table(cpg_meqtl_posi,"cpg_meQTL_posi",quote = F,col.names = F,row.names = F,sep="\t");

##################Step 3: Project the b37 to b36 reference pair:
chr_b37 <- read.table("chr1_noNA.value.bed",header = F);
chr_37 <- as.matrix(chr_b37);
integer_I <- read.table("chr1_noNA.value.bed.b36.index",header = F);
integer_J <- read.table("chr1_noNA.value.bed.interval_j.b36.index",header = F);
integer_I <- as.matrix(integer_I);
integer_J <- as.matrix(integer_J);

chr_37_index_I <- matrix(paste(chr_37[,1],"-",chr_37[,2],"-",chr_37[,3],sep=""),ncol=1);
chr_37_index_J <- matrix(paste(chr_37[,1],"-",chr_37[,4],"-",chr_37[,5],sep=""),ncol=1);
integer_I_index <- unique(matrix(paste(integer_I[,1],"-",integer_I[,2],"-",integer_I[,3],sep=""),ncol=1));
integer_J_index <- unique(matrix(paste(integer_J[,1],"-",integer_J[,2],"-",integer_J[,3],sep=""),ncol=1));

orga_index_I <- matrix(",",nrow = nrow(integer_I_index),ncol(integer_I));
sub = matrix(paste(integer_I[,1],"-",integer_I[,2],"-",integer_I[,3],sep=""),ncol=1);
for(i in 1:nrow(integer_I_index)){
  if(length(which(sub[,1] == integer_I_index[i,1])) == 1){
    orga_index_I[i,] = integer_I[which(sub[,1] == integer_I_index[i,1]),];
  } else{
    subregion = integer_I[which(sub[,1] == integer_I_index[i,1]),];
    orga_index_I[i,1:3] = subregion[1,1:3];
    orga_index_I[i,4] = "->";
    orga_index_I[i,5] = subregion[1,1];
    orga_index_I[i,6] = min(subregion[,6:7]);
    orga_index_I[i,7] = max(subregion[,6:7]);
  }
}

orga_index_J <- matrix(",",nrow = nrow(integer_J_index),ncol(integer_J));
sub = matrix(paste(integer_J[,1],"-",integer_J[,2],"-",integer_J[,3],sep=""),ncol=1);
for(i in 1:nrow(integer_J_index)){
  if(length(which(sub[,1] == integer_J_index[i,1])) == 1){
    orga_index_J[i,] = integer_J[which(sub[,1] == integer_J_index[i,1]),];
  } else{
    subregion = integer_J[which(sub[,1] == integer_J_index[i,1]),];
    orga_index_J[i,1:3] = subregion[1,1:3];
    orga_index_J[i,4] = "->";
    orga_index_J[i,5] = subregion[1,1];
    orga_index_J[i,6] = min(subregion[,6:7]);
    orga_index_J[i,7] = max(subregion[,6:7]);
  }
}

res <- matrix(".",nrow(chr_37),7);
nul = NULL;
for(i in 1:nrow(chr_37)){
  #print(i);
  subindex_I = which(chr_37_index_I[i,1] == integer_I_index[,1]);
  subindex_J = which(chr_37_index_J[i,1] == integer_J_index[,1]);
  if(length(subindex_I) == 1 && length(subindex_J) == 1){
    res[i,1:3] = orga_index_I[subindex_I,5:7];
    res[i,4:6] = orga_index_J[subindex_J,5:7];
    res[i,7] = chr_37[i,6];
  }
  if(length(subindex_I) == 1 && length(subindex_J) == 0){
    res[i,1:3] = orga_index_I[subindex_I,5:7];
    res[i,7] = chr_37[i,6];
  }
  if(length(subindex_I) == 0 && length(subindex_J) == 1){
    res[i,4:6] = orga_index_J[subindex_J,5:7];
    res[i,7] = chr_37[i,6];
  }
  # if(length(subindex_I) > 1 || length(subindex_J) > 1){
  #   print(i);
  # }
}
write.table(res,"res_chr1",quote = F,sep="\t",row.names = F,col.names = F);

res_sub <- read.table("res_chr1_sub",header = F);
res_sub <- as.matrix(res_sub);
diff_I = as.numeric(res_sub[,3])-as.numeric(res_sub[,2])
diff_J = as.numeric(res_sub[,6])-as.numeric(res_sub[,5])


############Step 5: To overlap CNV_meQTL with Hi-C region:
HiC_bychr <- read.table("res_chr1.noNA.b36",header = F);
HiC_bychr = as.matrix(HiC_bychr);
HiC_bychr_bin_I = matrix(paste(HiC_bychr[,1],"-",HiC_bychr[,2],"-",HiC_bychr[,3],sep=""),ncol=1);
HiC_bychr_bin_I = format(HiC_bychr_bin_I, scientific=FALSE)
HiC_bychr_bin_J = matrix(paste(HiC_bychr[,4],"-",HiC_bychr[,5],"-",HiC_bychr[,6],sep=""),ncol=1);
HiC_bychr_bin_J = format(HiC_bychr_bin_J, scientific=FALSE)
new_HiC_bychr = cbind(HiC_bychr_bin_I,HiC_bychr_bin_J,HiC_bychr[,7]);
new_HiC_bychr = format(new_HiC_bychr, scientific=FALSE)
new_HiC_bychr = gsub(" ","",new_HiC_bychr);

cnveQTL <- read.table("CNVeQTL_2kb",header = F);
cnveQTL <- as.matrix(cnveQTL);
cnv_index = matrix(paste(cnveQTL[,1],"-",cnveQTL[,2],"-",cnveQTL[,3],sep=""),ncol=1);
cnv_index = gsub(" ","",cnv_index);
#cnv_index = format(cnv_index, scientific=FALSE)
cpg_index = matrix(paste(cnveQTL[,4],"-",cnveQTL[,5],"-",cnveQTL[,6],sep=""),ncol=1);
cpg_index = gsub(" ","",cpg_index);
#cpg_index = format(cpg_index, scientific=FALSE)

cnv_HiC_by_chr_I <- read.table("res_chr1.noNA.b36.cnv.inter.I",header = F);
cnv_HiC_by_chr_J <- read.table("res_chr1.noNA.b36.cnv.inter.J",header = F);
cnv_HiC_by_chr_I = as.matrix(cnv_HiC_by_chr_I); 
# cnv_HiC_by_chr_I_index = matrix(paste(cnv_HiC_by_chr_I[,1],"-",cnv_HiC_by_chr_I[,2],"-",cnv_HiC_by_chr_I[,3],sep=""),ncol=1);
# cnv_HiC_by_chr_I_index = format(cnv_HiC_by_chr_I_index, scientific=FALSE)
cnv_HiC_by_chr_J = as.matrix(cnv_HiC_by_chr_J);
# cnv_HiC_by_chr_J_index = matrix(paste(cnv_HiC_by_chr_J[,1],"-",cnv_HiC_by_chr_J[,2],"-",cnv_HiC_by_chr_J[,3],sep=""),ncol=1)
# cnv_HiC_by_chr_J_index = format(cnv_HiC_by_chr_J_index, scientific=FALSE)
cnv_HiC_by_chr = unique(rbind(cbind(cnv_HiC_by_chr_I[,1:6],cnv_HiC_by_chr_I[,10]),cnv_HiC_by_chr_J));
cnv_HiC_by_chr_qtlindex = matrix(paste(cnv_HiC_by_chr[,1],"-",cnv_HiC_by_chr[,2],"-",cnv_HiC_by_chr[,3],sep=""),ncol=1);
cnv_HiC_by_chr_qtlindex = format(cnv_HiC_by_chr_qtlindex, scientific=FALSE);
cnv_HiC_by_chr_qtlindex = gsub(" ","",cnv_HiC_by_chr_qtlindex);
cnv_HiC_by_chr_hicindex = matrix(paste(cnv_HiC_by_chr[,4],"-",cnv_HiC_by_chr[,5],"-",cnv_HiC_by_chr[,6],sep=""),ncol=1);
cnv_HiC_by_chr_hicindex = format(cnv_HiC_by_chr_hicindex, scientific=FALSE);
cnv_HiC_by_chr_hicindex = gsub(" ","",cnv_HiC_by_chr_hicindex);

cpg_HiC_by_chr_I <- read.table("res_chr1.noNA.b36.cpg.2kb.inter.I",header = F);
cpg_HiC_by_chr_J <- read.table("res_chr1.noNA.b36.cpg.2kb.inter.J",header = F);
cpg_HiC_by_chr_I = as.matrix(cpg_HiC_by_chr_I);
# cpg_HiC_by_chr_I_index = matrix(paste(cpg_HiC_by_chr_I[,1],"-",cpg_HiC_by_chr_I[,2],"-",cpg_HiC_by_chr_I[,3],sep=""),ncol=1);
# cpg_HiC_by_chr_I_index = format(cpg_HiC_by_chr_I_index, scientific=FALSE)
cpg_HiC_by_chr_J = as.matrix(cpg_HiC_by_chr_J);
# cpg_HiC_by_chr_J_index = matrix(paste(cpg_HiC_by_chr_J[,1],"-",cpg_HiC_by_chr_J[,2],"-",cpg_HiC_by_chr_J[,3],sep=""),ncol=1);
# cpg_HiC_by_chr_J_index = format(cpg_HiC_by_chr_J_index, scientific=FALSE);
cpg_HiC_by_chr = unique(rbind(cbind(cpg_HiC_by_chr_I[,1:6],cpg_HiC_by_chr_I[,10]),cpg_HiC_by_chr_J))
cpg_HiC_by_chr_qtlindex = matrix(paste(cpg_HiC_by_chr[,1],"-",cpg_HiC_by_chr[,2],"-",cpg_HiC_by_chr[,3],sep=""),ncol=1);
cpg_HiC_by_chr_qtlindex = format(cpg_HiC_by_chr_qtlindex, scientific=FALSE);
cpg_HiC_by_chr_qtlindex = gsub(" ","",cpg_HiC_by_chr_qtlindex);
cpg_HiC_by_chr_hicindex = matrix(paste(cpg_HiC_by_chr[,4],"-",cpg_HiC_by_chr[,5],"-",cpg_HiC_by_chr[,6],sep=""),ncol=1);
cpg_HiC_by_chr_hicindex = format(cpg_HiC_by_chr_hicindex, scientific=FALSE);
cpg_HiC_by_chr_hicindex = gsub(" ","",cpg_HiC_by_chr_hicindex);

res <- NULL;
for(i in 1:nrow(cnveQTL)){
  print(i);
  subcnv_HiC = unique(cnv_HiC_by_chr_hicindex[which(cnv_index[i,1] == cnv_HiC_by_chr_qtlindex[,1]),]);
  subcpg_HiC = unique(cpg_HiC_by_chr_hicindex[which(cpg_index[i,1] == cpg_HiC_by_chr_qtlindex[,1]),]);
  
  if(length(subcnv_HiC) >0 && length(subcpg_HiC) > 0){
    sub_res = NULL;
    for(j in 1:length(subcnv_HiC)){
      sub_tmp1 = new_HiC_bychr[which(subcnv_HiC[j] == new_HiC_bychr[,1]),];
      
      sub_res = NULL;
      for(k in 1:length(subcpg_HiC)){
        sub_tmp2 = sub_tmp1[which(subcpg_HiC[k] == sub_tmp1[,2]),];
        if(length(sub_tmp2) > 0){
          sub_res = rbind(sub_res,c(cnveQTL[i,],sub_tmp2));
        }else{
          sub_res = rbind(sub_res,c(cnveQTL[i,],matrix("NA",nrow = 1,ncol = 3)));
        }
      }
      res = rbind(res,sub_res);
    }
    
    sub_res = NULL;
    for(j in 1:length(subcnv_HiC)){
      sub_tmp1 = new_HiC_bychr[which(subcnv_HiC[j] == new_HiC_bychr[,2]),];
      
      sub_res = NULL;
      for(k in 1:length(subcpg_HiC)){
        sub_tmp2 = sub_tmp1[which(subcpg_HiC[k] == sub_tmp1[,1]),];
        if(length(sub_tmp2) > 0){
          sub_res = rbind(sub_res,c(cnveQTL[i,],sub_tmp2));
        }else{
          sub_res = rbind(sub_res,c(cnveQTL[i,],matrix("NA",nrow = 1,ncol = 3)));
        }
      }
      res = rbind(res,sub_res);
    }
    # subcnv_HiC = NULL;
    # subcpg_HiC = NULL;
  }
}
res = unique(res);
write.table(res,"res_chr1_veri_2kb",quote = F,sep="\t",col.names = F,row.names = F);


###Step 6: add CpG ID into final result:
cpgid <- read.table("cpgid_2kb",header = F);
cpgid = unique(as.matrix(cpgid));
# cpg_2kb = matrix(paste(cpgid[,2],"-",cpgid[,3],"-",cpgid[,4],sep=""),ncol=1);
# cpg_2kb = gsub(" ","",cpg_2kb);
# 
cpg_5kb_index = matrix(paste(cpgid[,1],"-",cpgid[,2],"-",cpgid[,3],sep=""),ncol=1);
cpg_5kb_index = gsub(" ","",cpg_5kb_index);

cpg_5kb_res <- read.table("cpg_2kb",header=F);
#cpg_2kb_res <- read.table("cpg_2kb",header = F);
cpg_5kb_res = as.matrix(cpg_5kb_res);
cpg_5kb_res_index = matrix(paste(cpg_5kb_res[,1],"-",cpg_5kb_res[,2],"-",cpg_5kb_res[,3],sep=""),ncol=1);
cpg_5kb_res_index = gsub(" ","",cpg_5kb_res_index);
#cpg_2kb_res = as.matrix(cpg_2kb_res);

idres = NULL;
for(i in 1:nrow(cpg_5kb_res)){
  index = which(cpg_5kb_res_index[i,1] == cpg_5kb_index[,1]);
  if(length(index) > 1) {print(i)}
  idres <- rbind(idres,cpgid[index,]);
}

######identify the gene of CpG:
cpggene <- read.table("cpggene",header=F);
cpgres <- read.table("cpgres",header=F);


#####Identify cis and trans:
index <- read.table("index",header = F);
index <- as.matrix(index);
index_combine <- gsub(" ","",matrix(paste(index[,1],"-",index[,2],"-",index[,3],"-",index[,4],sep=""),ncol=1));

alleqtl <- read.table("alleqtl",header = F);
alleqtl <- as.matrix(alleqtl);
alleqtl <- gsub(" ","",matrix(paste(alleqtl[,1],"-",alleqtl[,2],"-",alleqtl[,3],"-",alleqtl[,4],sep=""),ncol=1));

res <- matrix(".",nrow(alleqtl),1);
for(i in 1:nrow(alleqtl)){
  res[i,1]=index[which(alleqtl[i,1] == index_combine[,1]),5];
}
write.table(res,"res",col.names = F,row.names = F,sep="\t",quote = F);