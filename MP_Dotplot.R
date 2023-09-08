library(Seurat)
# Set working directory
wd <- "H:/Manuscripts/2023_DS10_Cj/DESEQ2/Re_Analysis_MP" # MODIFY
setwd(wd)
filename <- "H:/Manuscripts/2023_DS10_Cj/DESEQ2/Re_Analysis_MP/MF.rds" # MODIFY
x <-readRDS(filename)
DimPlot(x)
######################
#############################################################Dot-Plot#############################################################
DotPlot(x, features = c("Hmox1","Ccl7","Apoe","Ccl2","Ifi207","Ednrb","Pf4","Mafb","Cd163","Mrc1",
                       "Lyz1","Fn1","Lpl","Retnla","Ear2","Tnip3","Flrt3","Fcrls","Il1rn","Ccr2",
                       "Plac8","Chil3","Thbs1","Gngt2","Clec4e","Lst1","Gsr","Ace","Msrb1","Prdx5",
                      "Aqp1","Ecm1", "Tppp3","Lyve1","Emp1","Fxyd2","C1qc","Vsig4","Tns1", "Timp2",
                      "Selenop","Ccl8","Maf", "C1qa","Gas6","Igfbp4", "F13a1","C1qb","Rnase4","Itm2b",
                      "Mmp12","Mmp19","Mmp14","Spp1","Axl","AA467197","Ctss","Lgals3","Anxa4","Psap","Ccl12"),
cols = c("blue", "red"), dot.scale = 6) + theme(axis.text.x = element_text(angle = 90, hjust=1))





