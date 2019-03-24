library(gplots)
library(viridis)

data <- read.csv("F:/Thesis/AlzheimersDisease/UtshaDas/AlzheimersDisease1/Dataset/GSE5281/HIP_T_UpDown.csv", row.names=1)


names_up <- vector('character')
fc_up <- vector('numeric')
lfc_up <- vector('numeric')
names_down <- vector('character')
fc_down <- vector('numeric')
lfc_down <- vector('numeric')

data_ = as.matrix(data)

for (i in 1:nrow(data)){
  vec1 = as.numeric(data[i:i,1:13])
  vec2 = as.numeric(data[i:i,14:23])
  
  m1 = mean(vec1)
  m2 = mean(vec2)
  
  data_[i,] = (data_[i,]-mean(data_[i,]))/(max(data_[i,]) - min(data_[i,]))
  FC = m2/m1
  
  logFC = log2(FC)
  
  if (logFC > 0)
  {
    names_up <- c(names_up, row.names(data)[i])
    fc_up <- c(fc_up, FC)
    lfc_up <- c(lfc_up, logFC)    
  }
  else
  {
    names_down <- c(names_down, row.names(data)[i])
    fc_down <- c(fc_down, FC)
    lfc_down <- c(lfc_down, logFC)    
  }
}


resUp = data.frame("Genes"=names_up, "FC Value"=fc_up, "log2 of FC"=lfc_up)
resDown = data.frame("Genes"=names_down, "FC Value"=fc_down, "log2 of FC"=lfc_down)

names = vector('character')
lfc = vector('numeric')

for (i in 1:length(names_up))
{
  names = c(names, names_up[i])
  lfc = c(lfc, lfc_up[i])
}
for (i in 1:length(names_down))
{
  names = c(names, names_down[i])
  lfc = c(lfc, lfc_down[i])
}
res_abs = data.frame(row.names=names, "log2 FC Value"=lfc)

save(res_abs, file="absfc.rdata")

concol = vector('character')

for (j in 1:ncol(data_))
{
  if (j <= 13)
    concol = c(concol, '#f92c32')
  else 
    concol = c(concol, '#127cd4')  
}  

heatmap.2(data_, main="DEGs", trace="none", margins=c(4,14), density="none", 
          cexRow=1, cexCol=0.2, col=viridis(150), ColSideColors=concol)
legend(1,3.50,legend=c("Group1","Group2"),fill=c('#f92c32','#127cd4'),cex=0.4)

write.csv(resUp, file="F:/Thesis/AlzheimersDisease/UtshaDas/AlzheimersDisease1/Dataset/GSE5281/HIP_T_upregulated.csv")
write.csv(resDown, file="F:/Thesis/AlzheimersDisease/UtshaDas/AlzheimersDisease1/Dataset/GSE5281/HIP_T_downregulated.csv")
