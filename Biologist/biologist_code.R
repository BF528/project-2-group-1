#7.1
dir1 <- "/../project/bf528/project_2/data/samples"
dir2 <- "/../projectnb/bf528/users/group1/project2"
# load the data file
P0_1 <- read.csv(file.path(dir2,"programmer_divya/P0_1_cufflinks/genes.fpkm_tracking"), sep = "\t", 
                 stringsAsFactors = F)
#P0_1 <- read.csv("C:/Users/41186/Desktop/gene_exp.diff",
#                 sep = "\t",
#                 stringsAsFactors = F)
P0_2 <- read.csv(file.path(dir1,"P0_2/genes.fpkm_tracking"), sep = "\t", 
                 stringsAsFactors = F)
AD_1 <- read.csv(file.path(dir1,"Ad_1/genes.fpkm_tracking"), sep = "\t", 
                 stringsAsFactors = F)
AD_2 <- read.csv(file.path(dir1,"Ad_2/genes.fpkm_tracking"), sep = "\t", 
                 stringsAsFactors = F)

# build the funtion: take the FPKM value from a gene list
subset.gene <- function(x1,x2,x3,x4,target.gene){
  result1 = rep(0,length(target.gene))
  result2 = rep(0,length(target.gene))
  result3 = rep(0,length(target.gene))
  result4 = rep(0,length(target.gene))
  for (i in c(1:length(target.gene))){
    if (target.gene[i] %in% x1$gene_short_name) 
      result1[i] = subset(x1,gene_short_name==
                            target.gene[i])$FPKM
    if (target.gene[i] %in% x2$gene_short_name)
      result2[i] = subset(x2,gene_short_name==
                            target.gene[i])$FPKM
    if (target.gene[i] %in% x3$gene_short_name)
      result3[i] = subset(x3,gene_short_name==
                            target.gene[i])$FPKM
    if (target.gene[i] %in% x4$gene_short_name)
      result4[i] = subset(x4,gene_short_name==
                            target.gene[i])$FPKM
  }
  result_P0 = (result1+result2)/2
  result_Ad = (result3+result4)/2
  names(result_P0) <- target.gene
  names(result_Ad) <- target.gene
  return(list("P0"=result_P0,"AD"=result_Ad))
}
sar_list <- c("Pdlim5","Pygm","Myoz2","Des","Csrp3",
              "Tcap","Cryab")
sar_gene.fpkm <- subset.gene(P0_1,P0_2,AD_1,AD_2,sar_list) 
sar_gene.fpkm

mito_list <- c("Mpc1","Prdx3","Acat1","Echs1","Slc25a11","Phyh")
mito_gene.fpkm <- subset.gene(P0_1,P0_2,AD_1,AD_2,mito_list)
mito_gene.fpkm

cellc_list <- c("Cdc7","E2f8","Cdk7","Cdc26","Cdc6","Cdc27",
                "E2f1","Bora","Cdc45","Rad51","Aurkb","Cdc23")
cellc_gene.fpkm <- subset.gene(P0_1,P0_2,AD_1,AD_2,cellc_list)
cellc_gene.fpkm

#plot function for genes highlighted in Fig 1D
gene_plot <- function(gene_list,title){
  # Gene_list: highlighted gene with fpkm
  # Remove 0 fpkm
  gene_list[[1]] = gene_list[[1]][!gene_list[[1]] %in% 0]
  gene_list[[2]] = gene_list[[2]][!gene_list[[2]] %in% 0]
  # x = 100 is the x axis position for P0, x = 200 is for Ad
  xdat <- c(100,200)
  # Create a random vector representing the point characters
  pch_vec <- sample(c(1:6,8,15:18),
                    length(gene_list[[1]]),
                    replace = F)
  # Create a random vector representing the color of lines
  color_vec <- sample(c("yellow","red","green",
                        "orange","lightseagreen",
                        "purple","orange4","moccasin",
                        "lightsteelblue1","darkgrey",
                        "deeppink"),
                      length(gene_list[[1]]),
                      replace = F)
  # Plot the first line
  plot(xdat,
       c(gene_list[[1]][1],gene_list[[2]][1]),
       type = "o",
       lty = 1, 
       col = "blue", 
       pch = 0, 
       xaxt = "n", 
       xlab = "", 
       ylab = "FPKM",
       main = title,
       xlim = c(95,280),
       ylim = c(min(min(gene_list[[1]]),min(gene_list[[2]])) - 5, 
                max(max(gene_list[[1]]),max(gene_list[[2]])) + 5),
       lwd = 3 )
  # Change the x ticks to P0 & Ad
  axis(side = 1, 
       at = c(100,200),
       labels = c("P0","Ad"))
  # Build a for loop to add the subsequent lines
  # The cols and pchs are taken from 
  # the random vectors we created above 
  for (i in c(2:length(gene_list[[1]]))){
    points(xdat, 
           c(gene_list[[1]][i], gene_list[[2]][i]),
           col = color_vec[i], 
           pch = pch_vec[i])
    lines(xdat,
          c(gene_list[[1]][i], gene_list[[2]][i]),
          col = color_vec[i], 
          lty = 1, 
          lwd = 3)
  }
  # add the legend 
  legend("topright",
         legend = names(gene_list[[1]]), 
         col = c("blue",color_vec[-1]),
         pch = c(0,pch_vec[-1]))
}
pdf(file.path(dir2,"biologist/figure.pdf"))
sar_plot <- gene_plot(sar_gene.fpkm,"Sarcomere")
mito_plot <- gene_plot(mito_gene.fpkm,"Mitochondria")
cellc_plot <- gene_plot(cellc_gene.fpkm,"Cell cycle")


#7.2
up <- read.csv(file.path(dir2,"analyst/up_cluster.csv"),skip = 1,
               stringsAsFactors = F,na.strings = "")
# GO terms for upregulated genes 
up <- up[1:2413,]  
# subset the clusters with enrichment score > -log10(0.05)
up <- up[complete.cases(up),]
up <- subset(up,Category!="Category")
# remove the redundant colnames
library(readxl)
up_paper_common <- read_xlsx(file.path(dir2,"biologist/304269r2_online_table_i.xlsx"),
                             sheet = 4, skip = 1)
up_paper_unique <- read_xlsx(file.path(dir2,"biologist/304269r2_online_table_i.xlsx"),
                             sheet = 8, skip = 1)
# Go terms for upregulated genes in vivo from paper
up_paper_common <- up_paper_common[complete.cases(up_paper_common),]
up_paper_unique <- up_paper_unique[complete.cases(up_paper_unique),]
up_paper <- rbind(up_paper_common,up_paper_unique)
up_paper <- subset(up_paper,Category %in% 
                     c("GOTERM_MF_FAT","GOTERM_BP_FAT","GOTERM_CC_FAT"))
up_select <- up$Term %in% up_paper$Term
up$Overlap <- ifelse(up_select,'Yes','No')
# create a new column with overlap information
write.csv(up,file = file.path(dir2,"biologist/TableS1.csv"))

down <- read.csv(file.path(dir2,"analyst/down_cluster.csv"),skip=1,
                 stringsAsFactors = F,na.strings = "")
# GO terms for downregulated genes
down <- down[1:2293,]
# subset the clusters with enrichment score bigger than -log10(0.05)
down <- down[complete.cases(down),]
down <- subset(down,Category!="Category")
# remove the redundant colnames
down_paper_common <- read_xlsx(file.path(dir2,"biologist/304269r2_online_table_i.xlsx"),
                               sheet = 5, skip = 1)
down_paper_unique <- read_xlsx(file.path(dir2,"biologist/304269r2_online_table_i.xlsx"),
                               sheet = 9, skip = 1)
# Go terms for upregulated genes from paper
down_paper_common <- down_paper_common[complete.cases(down_paper_common),]
down_paper_unique <- down_paper_unique[complete.cases(down_paper_unique),]
down_paper <- rbind(down_paper_common,down_paper_unique)
down_paper <- subset(down_paper,Category %in% 
                       c("GOTERM_MF_FAT","GOTERM_BP_FAT","GOTERM_CC_FAT"))
down_select <- down$Term %in% down_paper$Term
down$Overlap <- ifelse(down_select,'Yes','No')
# create a new column with overlap information
write.csv(down,file = file.path(dir2,"biologist/TableS2.csv"))
#7.3
P4_1 <- read.csv(file.path(dir1,"P4_1/genes.fpkm_tracking"), sep = "\t", header = T,
                 stringsAsFactors = F)
P4_2 <- read.csv(file.path(dir1,"P4_2/genes.fpkm_tracking"), sep = "\t", header = T,
                 stringsAsFactors = F)
P7_1 <- read.csv(file.path(dir1,"P7_1/genes.fpkm_tracking"), sep = "\t", header = T,
                 stringsAsFactors = F)
P7_2 <- read.csv(file.path(dir1,"P7_2/genes.fpkm_tracking"), sep = "\t", header = T,
                 stringsAsFactors = F)
exp_diff <- read.csv(file.path(dir2,"analyst/gene_exp.diff"), sep = "\t", header = T,
                     stringsAsFactors = F)
# Take the average of FPKM for duplicate gene
P0_1 <- P0_1[c("gene_short_name","FPKM")]
P0_2 <- P0_2[c("gene_short_name","FPKM")]
P4_1 <- P4_1[c("gene_short_name","FPKM")]
P4_2 <- P4_2[c("gene_short_name","FPKM")]
P7_1 <- P7_1[c("gene_short_name","FPKM")]
P7_2 <- P7_2[c("gene_short_name","FPKM")]
AD_1 <- AD_1[c("gene_short_name","FPKM")]
AD_2 <- AD_2[c("gene_short_name","FPKM")]
library(dplyr)
P0_1 <- P0_1 %>% group_by(gene_short_name) %>% summarise_all(sum)
P0_2 <- P0_2 %>% group_by(gene_short_name) %>% summarise_all(sum)
P4_1 <- P4_1 %>% group_by(gene_short_name) %>% summarise_all(sum)
P4_2 <- P4_2 %>% group_by(gene_short_name) %>% summarise_all(sum)
P7_1 <- P7_1 %>% group_by(gene_short_name) %>% summarise_all(sum)
P7_2 <- P7_2 %>% group_by(gene_short_name) %>% summarise_all(sum)
AD_1 <- AD_1 %>% group_by(gene_short_name) %>% summarise_all(sum)
AD_2 <- AD_2 %>% group_by(gene_short_name) %>% summarise_all(sum)

dat_merge <- merge(P0_1,P0_2,by="gene_short_name")
names(dat_merge)[c(2,3)] <- c("P0_1","P0_2")
dat_merge <- merge(dat_merge,P4_1,by="gene_short_name")
names(dat_merge)[4] <- "P4_1"
dat_merge <- merge(dat_merge,P4_2,by="gene_short_name")
names(dat_merge)[5] <- "P4_2"
dat_merge <- merge(dat_merge,P7_1,by="gene_short_name")
names(dat_merge)[6] <- "P7_1"
dat_merge <- merge(dat_merge,P7_2,by="gene_short_name")
names(dat_merge)[7] <- "P7_2"
dat_merge <- merge(dat_merge,AD_1,by="gene_short_name")
names(dat_merge)[8] <- "AD_1"
dat_merge <- merge(dat_merge,AD_2,by="gene_short_name")
names(dat_merge)[9] <- "AD_2"
names(dat_merge)
dim(dat_merge)
exp_diff <- exp_diff[order(exp_diff$q_value)[1:200],]
dat_merge <- dat_merge[dat_merge$gene_short_name %in% exp_diff$gene,]
dat_merge <- dat_merge[complete.cases(dat_merge),]
rownames(dat_merge)<-dat_merge$gene_short_name
heatmap(as.matrix(dat_merge[2:9]),
        cexCol=1.2,cexRow=0.1,labRow=F,
        ColSideColors=rep(c('cyan','deepskyblue',
                            'blue','blue4'),each=2)
        )
legend('topright',legend=c('P0','P4','P7','Ad'),
       fill=c('cyan','deepskyblue','blue','blue4'),
       cex=0.5
       )
dev.off()
