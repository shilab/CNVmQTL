
#####MSE and MAE
setwd("~/Desktop");
regions <- read.table("regions",header = F);
regions <- as.matrix(regions);
dis <- matrix(".",nrow = nrow(regions),1);

for(i in 1:nrow(regions)){
  print(i);
  if (regions[i,1] > regions[i,4]) { dis[i,1]=regions[i,1]-regions[i,4];}
  if (regions[i,3] > regions[i,2]) { dis[i,1]=regions[i,3]-regions[i,2];}
}
write.table(dis,"dis",quote=F,sep="\t",col.names = F,row.names = F);

mQTL <- read.table("~/Desktop/mQTL",header = F);
dis <- read.table("~/Desktop/dis",header = F);
mQTL <- as.matrix(mQTL);
dis <- as.matrix(dis);
mQTL_index = matrix(gsub(" ","",paste(mQTL[,1],"-",mQTL[,2],"-",mQTL[,3],"-",mQTL[,4],"-",mQTL[,5],"-",mQTL[,11],"-",mQTL[,12])),ncol=1);
dis_index = matrix(gsub(" ","",paste(dis[,1],"-",dis[,2],"-",dis[,3],"-",dis[,7],"-",dis[,4],"-",dis[,5])),ncol=1);

res <- matrix(",",851,8);
for(i in 1:nrow(mQTL)){
  index = which(mQTL_index[i,1] == dis_index[,1]);
  if(length(index) ==1){
    res[i,] = dis[index,];
  }
}
write.table(res,"res",quote=F,sep="\t",col.names = F,row.names = F);

########Generate the CpG shared regions:
setwd("~/HiC/enrich/cnv_per/");
#files = list.files(path=".",pattern = "cnv_per_");
per_ids <- read.table("~/HiC/enrich/per_dis");
per_ids <- as.matrix(per_ids);

for(f in 1:1000){
  data <- read.table(paste("cnv_per_",f,sep=""),header = F);
  data <- as.matrix(data);
  newdata <- cbind(data,per_ids[,7]);
  
  cpg_per <- NULL;
  for(chr in 1:22){
    
    subchr=newdata[which(newdata[,1] == chr),];
    sub_per=per_ids[which(as.numeric(per_ids[,1]) == chr),];
    sub_per_cis = sub_per[which(sub_per[,7] == "cis"),];
    sub_per_trans = sub_per[which(sub_per[,7] == "trans"),,drop=FALSE];
    
    subdistance_cis = sub_per_cis[sample(1:nrow(sub_per_cis),replace = FALSE),6];
    sub_cpg_cis = sub_per_cis[sample(1:nrow(sub_per_cis)),8];
    index_cis = which(subchr[,4] == "cis");
    
    subdistance_trans = sub_per_trans[sample(1:nrow(sub_per_trans)),6];
    sub_cpg_trans = sub_per_trans[sample(1:nrow(sub_per_trans)),8];
    index_trans = which(subchr[,4] == "trans");
    
    sub_cpg_dis = matrix(".",nrow(subchr),2);
    sub_cpg_dis[index_cis,1]=subdistance_cis;
    sub_cpg_dis[index_cis,2]=sub_cpg_cis;
    sub_cpg_dis[index_trans,1]=subdistance_trans;
    sub_cpg_dis[index_trans,2]=sub_cpg_trans;
    
    per_newdata = cbind(subchr,sub_cpg_dis);
    subcpg_per = matrix(".",nrow(per_newdata),2);
    
    for(j in 1:nrow(subchr)){
      if(as.numeric(sub_cpg_dis[j,1]) != 0 && as.numeric(sub_per[j,6])!=0){
        if (sub_per[j,5] < sub_per[j,2] && sub_per[j,7] == "cis") {
          subcpg_per[j,2] = as.numeric(subchr[j,2]) - as.numeric(sub_cpg_dis[j,1]);
          subcpg_per[j,1] = as.numeric(subcpg_per[j,2]) - as.numeric(sub_cpg_dis[j,1]);
          if(as.numeric(subcpg_per[j,1]) <0 | as.numeric(subcpg_per[j,2])){
            subcpg_per[j,1] = as.numeric(subchr[j,3]) + as.numeric(sub_cpg_dis[j,1]);
            subcpg_per[j,2] = as.numeric(subcpg_per[j,1]) + as.numeric(sub_cpg_dis[j,2]);
          }
        }
        if (sub_per[j,3] < sub_per[j,4] && sub_per[j,7] == "cis") {
          subcpg_per[j,1] = as.numeric(subchr[j,3]) + as.numeric(sub_cpg_dis[j,1]);
          subcpg_per[j,2] = as.numeric(subcpg_per[j,1]) + as.numeric(sub_cpg_dis[j,2]);
        }
      }
      if(as.numeric(sub_per[j,6])==0 && as.numeric(sub_cpg_dis[j,1]) != 0){
        subcpg_per[j,1] = as.numeric(subchr[j,3]) + as.numeric(sub_cpg_dis[j,1]);
        subcpg_per[j,2] = as.numeric(subcpg_per[j,1]) + as.numeric(sub_cpg_dis[j,2]);
      }
      
      if(as.numeric(sub_cpg_dis[j,1]) == 0){
        integer = sample(as.numeric(subchr[j,2]):as.numeric(subchr[j,3]),1,replace = FALSE)
        subcpg_per[j,1] = integer;
        subcpg_per[j,2] = integer + as.numeric(sub_cpg_dis[j,2]);
      }
        
      if (sub_per[j,5] < sub_per[j,2] && sub_per[j,7] == "trans") {
        subcpg_per[j,2] = as.numeric(subchr[j,2]) - as.numeric(sub_cpg_dis[j,1]);
        subcpg_per[j,1] = as.numeric(subcpg_per[j,2]) - as.numeric(sub_cpg_dis[j,1]);
        if(as.numeric(subcpg_per[j,1] <0) | as.numeric(subcpg_per[j,2])){
          subcpg_per[j,1] = as.numeric(subchr[j,3]) + as.numeric(sub_cpg_dis[j,1]);
          subcpg_per[j,2] = as.numeric(subcpg_per[j,1]) + as.numeric(sub_cpg_dis[j,2]);
        }
      }
      
      if (sub_per[j,3] < sub_per[j,4] && sub_per[j,7] == "trans") {
        subcpg_per[j,1] = as.numeric(subchr[j,3]) + as.numeric(sub_cpg_dis[j,1]);
        subcpg_per[j,2] = as.numeric(subcpg_per[j,1]) + as.numeric(sub_cpg_dis[j,2]);
      }
    }
    subres = cbind(matrix(subchr[,1],nrow = nrow(subcpg_per),1),subcpg_per);
    cpg_per <- rbind(cpg_per,subres);
  }
  write.table(cpg_per,paste("cpg_per_5kb_cnv_per_",f,sep=""),quote=F,sep="\t",col.names = F,row.names = F);
}


