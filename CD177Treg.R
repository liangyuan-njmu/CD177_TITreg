##########################################################
#Figure1
##########################################################
#HCC CD4
library(Seurat)
library(ggplot2)
library(harmony)

data1=readRDS('../GSE98638_data.RDS')

data1<-NormalizeData(data1)
data1<-ScaleData(data1,features=rownames(data1))
data1<-FindVariableFeatures(data1,selection.method="vst",nfeatures=4000)
data1<-RunPCA(data1,features=VariableFeatures(data1),npcs=50)
data1<-RunHarmony(data1, group.by.vars="Patient",plot_convergence = F,reduction.save = "harmony",max.iter.harmony = 30) #迭代5次结束
data1<-FindNeighbors(data1,reduction = "harmony",dims=1:50)
data1<-FindClusters(data1,resolution=1)
data1<-RunUMAP(data1,dims=1:50,reduction = "harmony")

TsnePlot<-DimPlot(data1,raster=T, group.by="RNA_snn_res.0.5",reduction = "umap",pt.size = 4,label = T,label.size=4,cols=cols,raster.dpi=c(2048,2048))
ggsave("./data1_res0.5_umap.pdf",TsnePlot,width=11.5,height=10,units='cm')

TsnePlot<-DimPlot(data1,raster=T, group.by="RNA_snn_res.1",reduction = "umap",pt.size = 4,label = T,label.size=4,cols=cols,raster.dpi=c(2048,2048))
ggsave("./data1_res1_umap.pdf",TsnePlot,width=11.5,height=10,units='cm')

TsnePlot<-DimPlot(data1,raster=T, group.by="majorCluster",reduction = "umap",pt.size = 4,label = T,label.size=4,cols=cols,raster.dpi=c(2048,2048))
ggsave("./data1_majorCluster.pdf",TsnePlot,width=13,height=10,units='cm')

fea=FeaturePlot(data1,raster=T,features=c('CD4','CD8A','CD8B','FOXP3','CTLA4'),col=c("lightgray","red"),pt.size=2,reduction="umap",ncol=1,raster.dpi=c(1024,1024))
ggsave(paste0("./data1_","CD4_markers_fea.pdf"),fea,width=11,height=10*5,units='cm',limitsize = FALSE)

tmp=subset(data1,subset=RNA_snn_res.0.5 %in% c(0,6,3,13,12,1,7,9,2))
cells=colnames(tmp)[which(tmp@assays$RNA@data['CD8A',]==0 & tmp@assays$RNA@data['CD8B',]==0)]
CD4T=subset(tmp,cells=cells)
CD4T<-NormalizeData(CD4T)
CD4T<-ScaleData(CD4T,features=rownames(CD4T))
CD4T<-FindVariableFeatures(CD4T,selection.method="vst",nfeatures=4000)
CD4T<-RunPCA(CD4T,features=VariableFeatures(CD4T),npcs=50)
CD4T<-RunHarmony(CD4T, group.by.vars="Patient",plot_convergence = F,reduction.save = "harmony",max.iter.harmony = 30) #迭代5次结束
CD4T<-FindNeighbors(CD4T,reduction = "harmony",dims=1:50)
CD4T<-FindClusters(CD4T,resolution=2)
CD4T<-RunUMAP(CD4T,dims=1:50,reduction = "harmony")

TsnePlot<-DimPlot(CD4T,raster=T, group.by="RNA_snn_res.0.5",reduction = "umap",pt.size = 4,label = T,label.size=4,cols=cols,raster.dpi=c(2048,2048))
ggsave("./CD4T_res0.5_umap.pdf",TsnePlot,width=11.5,height=10,units='cm')
TsnePlot<-DimPlot(CD4T,raster=T, group.by="RNA_snn_res.2",reduction = "umap",pt.size = 8,label = T,label.size=4,cols=cols,raster.dpi=c(2048,2048))
ggsave("./CD4T_res2_umap.pdf",TsnePlot,width=11.5,height=10,units='cm')

fea=FeaturePlot(CD4T,raster=T,features=c('CD4','FOXP3','CTLA4'),col=c("lightgray","red"),pt.size=4,reduction="umap",ncol=1,raster.dpi=c(1024,1024))
ggsave(paste0("./CD4T_","CD4_markers_fea.pdf"),fea,width=11,height=10*3,units='cm',limitsize = FALSE)
VLN=VlnPlot(CD4T,pt.size=0,features =c('CD4','CD8A','CD8B','FOXP3','CTLA4'),group.by="RNA_snn_res.2",cols=cols)
ggsave(paste0("./CD4T_","_markers_vln.pdf"),VLN,width=50,height=30,units='cm',limitsize = FALSE)

HCCTreg=subset(CD4T,subset=RNA_snn_res.2 %in% c(0,4,5,10))
HCCTreg=subset(HCCTreg,subset=sampleType %in% c('NTR','PTR','TTR'))

badterm <- c(
  "^RPS","^RPL","^MT-","^MTMR","^MTND","^MTRN","^MTCO","^MRPL","^MRPS","^HBA","^HBB","^MTATP",
  "^MALAT1$", "^XIST$", "^XIST_intron$",
  "^AC[0-9]+\\.", "^AP[0-9]+\\.", "^AL[0-9]+\\.", "-AS[0-9]*$", "^MIR[0-9]", "^LINC[0-9]", "^SNOR", "^sno", "^SNR",
  "^TMSB", "^HIST", "^HSP", "^IFI", "^HLA-", "^ATP", "-", ":", "\\.", '^KIAA',
  "^DNAJ", '^FTH', '^FTL', '^LGALS'
)
badgenelt<-unique(grep(pattern=paste(c(badterm),collapse="|"),x=rownames(HCCTreg),perl=T,value=T))

diff_rs1=FindMarkers(HCCTreg,features=rownames(HCCTreg)[-which(rownames(HCCTreg) %in% badgenelt)],ident.1="TTR",ident.2='PTR',group.by="sampleType",logfc.threshold=0,min.pct=0.15,only.pos=F,pseudocount.use=F)
diff_rs1=diff_rs1[order(diff_rs1$avg_log2FC,decreasing=T),]
diff_rs1$gene=rownames(diff_rs1)

CD4T$type=ifelse(CD4T$RNA_snn_res.2 %in% c(0,4,5,10),'Treg','Others')
tmp=subset(CD4T,subset=sampleType %in% c('TTH','TTR'))
diff_rs2=FindMarkers(tmp,features=rownames(tmp)[-which(rownames(tmp) %in% badgenelt)],ident.1="Treg",ident.2='Others',group.by="type",logfc.threshold=0,min.pct=0.15,only.pos=F,pseudocount.use=F)
diff_rs2=diff_rs2[order(diff_rs2$avg_log2FC,decreasing=T),]
diff_rs2$gene=rownames(diff_rs2)

inte_genes=intersect(diff_rs1$gene,diff_rs2$gene)

diff_rs1=diff_rs1[inte_genes,]
diff_rs2=diff_rs2[inte_genes,]

df=data.frame(
    avg_log2FC_TP=diff_rs1$avg_log2FC,
    p_adj_TP=diff_rs1$p_val_adj,
    avg_log2FC_TO=diff_rs2$avg_log2FC,
    p_adj_TO=diff_rs2$p_val_adj,
    gene=inte_genes
)
df$sig=ifelse(df$p_adj_TP<0.001 & df$p_adj_TO<0.001 & df$avg_log2FC_TP>=2 & df$avg_log2FC_TO>=2,"Sig1",
        ifelse(df$p_adj_TP<0.001 & df$p_adj_TO<0.001 & df$avg_log2FC_TP<=(-2) & df$avg_log2FC_TO<=(-2),"Sig2", "Nosig"))
labels=c(
    #intersect(diff_rs1$gene[diff_rs1$p_val_adj<0.001][1:20],diff_rs2$gene[diff_rs2$p_val_adj<0.001][1:20]),
    #intersect(tail(diff_rs1[diff_rs1$p_val_adj<0.001,],20)$gene,tail(diff_rs2[diff_rs2$p_val_adj<0.001,],20)$gene)
    'CD177','IL1R2','TNFRSF13B','WLS','TNFRSF9','LAYN','CLEC7A','CCR8','TSPAN13','IL1R1','TNS3',
    'ANXA1','CCR7','KLF2','KLF3','S1PR1','TC2N','RUNX2'
)
df$label=ifelse(df$gene %in% labels,df$gene,"")
#df$label=ifelse(df$sig %in% c('Sig1','Sig2'),df$gene,"")
p<-ggplot(df,aes(x=avg_log2FC_TP,y=avg_log2FC_TO,color=sig,label=label))+
    ggrastr::rasterize(geom_point(size=2),dpi=150)+
    ggrepel::geom_text_repel(size=6,max.overlaps=Inf,nudge_x=0.25,nudge_y=0.25)+
    geom_hline(yintercept = c(2,-2), color="grey",linetype = 2,size = 0.5)+
    geom_vline(xintercept = c(2,-2), color="grey",linetype = 2,size = 0.5)+
    scale_color_manual(values=c("Sig1"="#EC478B",'Sig2'='#3E94D6',"NoSig"="lightgrey"))+
    xlab("avg_log2FC_TP")+ylab("avg_log2FC_TO")+
    theme_bw()+theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0.5),panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
ggsave("HCCTreg_TP_TO.pdf",p,width=5.5,height=4.5)

df1=df
HCC_TITreg_genes=df1$gene[df1$sig=='Sig1']

#CRC CD4 GSE146771
library(Seurat)
library(ggplot2)
library(harmony)

data2=readRDS('../../data/GSE146771_CRC/GSE146771_Smartseq2.RDS')

data2<-NormalizeData(data2)
data2<-ScaleData(data2,features=rownames(data2))
data2<-FindVariableFeatures(data2,selection.method="vst",nfeatures=4000)
data2<-RunPCA(data2,features=VariableFeatures(data2),npcs=50)
data2<-RunHarmony(data2, group.by.vars="Sample",plot_convergence = F,reduction.save = "harmony",max.iter.harmony = 30) #迭代5次结束
data2<-FindNeighbors(data2,reduction = "harmony",dims=1:50)
data2<-FindClusters(data2,resolution=0.5)
data2<-RunUMAP(data2,dims=1:50,reduction = "harmony")

TsnePlot<-DimPlot(data2,raster=T, group.by="RNA_snn_res.0.5",reduction = "umap",pt.size = 4,label = T,label.size=4,cols=cols,raster.dpi=c(2048,2048))
ggsave("./data2_res0.5_umap.pdf",TsnePlot,width=11.5,height=10,units='cm')

TsnePlot<-DimPlot(data2,raster=T, group.by="RNA_snn_res.1",reduction = "umap",pt.size = 4,label = T,label.size=4,cols=cols,raster.dpi=c(2048,2048))
ggsave("./data2_res1_umap.pdf",TsnePlot,width=11.5,height=10,units='cm')

fea=FeaturePlot(data2,raster=T,features=c('CD3D','CD3E','CD4','CD8A','CD8B','FOXP3','CTLA4'),col=c("lightgray","red"),pt.size=2,reduction="umap",ncol=1,raster.dpi=c(1024,1024))
ggsave(paste0("./data2_","CD4_markers_fea.pdf"),fea,width=11,height=10*7,units='cm',limitsize = FALSE)

tmp=subset(data2,subset=RNA_snn_res.0.5 %in% c(0,4,6))
cells=colnames(tmp)[which(tmp@assays$RNA@data['CD8A',]==0 & tmp@assays$RNA@data['CD8B',]==0)]
CD4T=subset(tmp,cells=cells)
CD4T<-NormalizeData(CD4T)
CD4T<-ScaleData(CD4T,features=rownames(CD4T))
CD4T<-FindVariableFeatures(CD4T,selection.method="vst",nfeatures=4000)
CD4T<-RunPCA(CD4T,features=VariableFeatures(CD4T),npcs=50)
CD4T<-RunHarmony(CD4T, group.by.vars="Sample",plot_convergence = F,reduction.save = "harmony",max.iter.harmony = 30) #迭代5次结束
CD4T<-FindNeighbors(CD4T,reduction = "harmony",dims=1:50)
CD4T<-FindClusters(CD4T,resolution=2)
CD4T<-RunUMAP(CD4T,dims=1:50,reduction = "harmony")

TsnePlot<-DimPlot(CD4T,raster=T, group.by="RNA_snn_res.0.5",reduction = "umap",pt.size = 4,label = T,label.size=4,cols=cols,raster.dpi=c(2048,2048))
ggsave("./CD4T_res0.5_umap.pdf",TsnePlot,width=11.5,height=10,units='cm')
TsnePlot<-DimPlot(CD4T,raster=T, group.by="RNA_snn_res.2",reduction = "umap",pt.size = 10,label = T,label.size=4,cols=cols,raster.dpi=c(2048,2048))
ggsave("./CD4T_res2_umap.pdf",TsnePlot,width=11.5,height=10,units='cm')

fea=FeaturePlot(CD4T,raster=T,features=c('CD4','FOXP3','CTLA4'),col=c("lightgray","red"),pt.size=10,reduction="umap",ncol=1,raster.dpi=c(2048,2048))
ggsave(paste0("./CD4T_","CD4_markers_fea.pdf"),fea,width=11,height=10*3,units='cm',limitsize = FALSE)
VLN=VlnPlot(CD4T,pt.size=0,features =c('CD4','CD8A','CD8B','FOXP3','CTLA4'),group.by="RNA_snn_res.2",cols=cols)
ggsave(paste0("./CD4T_","_markers_vln.pdf"),VLN,width=50,height=30,units='cm',limitsize = FALSE)

HCCTreg=subset(CD4T,subset=RNA_snn_res.2 %in% c(0,4,6))

badterm <- c(
  "^RPS","^RPL","^MT-","^MTMR","^MTND","^MTRN","^MTCO","^MRPL","^MRPS","^HBA","^HBB","^MTATP",
  "^MALAT1$", "^XIST$", "^XIST_intron$",
  "^AC[0-9]+\\.", "^AP[0-9]+\\.", "^AL[0-9]+\\.", "-AS[0-9]*$", "^MIR[0-9]", "^LINC[0-9]", "^SNOR", "^sno", "^SNR",
  "^TMSB", "^HIST", "^HSP", "^IFI", "^HLA-", "^ATP", "-", ":", "\\.", '^KIAA',
  "^DNAJ", '^FTH', '^FTL', '^LGALS'
)
badgenelt<-unique(grep(pattern=paste(c(badterm),collapse="|"),x=rownames(HCCTreg),perl=T,value=T))

diff_rs1=FindMarkers(HCCTreg,features=rownames(HCCTreg)[-which(rownames(HCCTreg) %in% badgenelt)],ident.1="T",ident.2='P',group.by="Tissue",logfc.threshold=0,min.pct=0.15,only.pos=F,pseudocount.use=F)
diff_rs1=diff_rs1[order(diff_rs1$avg_log2FC,decreasing=T),]
diff_rs1$gene=rownames(diff_rs1)
diff_rs1$avg_log2FC[1:2]=c(7.3,7.3)

CD4T$type=ifelse(CD4T$RNA_snn_res.2 %in% c(0,4,6),'Treg','Others')
tmp=subset(CD4T,subset=Tissue == 'T')
diff_rs2=FindMarkers(tmp,features=rownames(tmp)[-which(rownames(tmp) %in% badgenelt)],ident.1="Treg",ident.2='Others',group.by="type",logfc.threshold=0,min.pct=0.15,only.pos=F,pseudocount.use=F)
diff_rs2=diff_rs2[order(diff_rs2$avg_log2FC,decreasing=T),]
diff_rs2$gene=rownames(diff_rs2)

inte_genes=intersect(diff_rs1$gene,diff_rs2$gene)

diff_rs1=diff_rs1[inte_genes,]
diff_rs2=diff_rs2[inte_genes,]

df=data.frame(
    avg_log2FC_TP=diff_rs1$avg_log2FC,
    p_adj_TP=diff_rs1$p_val_adj,
    avg_log2FC_TO=diff_rs2$avg_log2FC,
    p_adj_TO=diff_rs2$p_val_adj,
    gene=inte_genes
)
df$sig=ifelse(df$p_adj_TP<0.001 & df$p_adj_TO<0.001 & df$avg_log2FC_TP>=2 & df$avg_log2FC_TO>=2,"Sig1",
        ifelse(df$p_adj_TP<0.001 & df$p_adj_TO<0.001 & df$avg_log2FC_TP<=(-2) & df$avg_log2FC_TO<=(-2),"Sig2", "Nosig"))
labels=c(
    #intersect(diff_rs1$gene[diff_rs1$p_val_adj<0.001][1:20],diff_rs2$gene[diff_rs2$p_val_adj<0.001][1:20]),
    #intersect(tail(diff_rs1[diff_rs1$p_val_adj<0.001,],20)$gene,tail(diff_rs2[diff_rs2$p_val_adj<0.001,],20)$gene)
    'CD177','IL1R2','TNFRSF13B','WLS','TNFRSF9','LAYN','CLEC7A','CCR8','TSPAN13','IL1R1','TNS3',
    'ANXA1','CCR7','KLF2','KLF3','S1PR1','TC2N','RUNX2'
)
df$label=ifelse(df$gene %in% labels,df$gene,"")
#df$label=ifelse(df$sig %in% c('Sig1','Sig2'),df$gene,"")
p<-ggplot(df,aes(x=avg_log2FC_TP,y=avg_log2FC_TO,color=sig,label=label))+
    ggrastr::rasterize(geom_point(size=2),dpi=150)+
    ggrepel::geom_text_repel(size=6,max.overlaps=Inf,nudge_x=0.25,nudge_y=0.25)+
    geom_hline(yintercept = c(2,-2), color="grey",linetype = 2,size = 0.5)+
    geom_vline(xintercept = c(2,-2), color="grey",linetype = 2,size = 0.5)+
    scale_color_manual(values=c("Sig1"="#EC478B",'Sig2'='#3E94D6',"NoSig"="lightgrey"))+
    xlab("avg_log2FC_TP")+ylab("avg_log2FC_TO")+
    theme_bw()+theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0.5),panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
ggsave("CRCTreg_TP_TO.pdf",p,width=5.5,height=4.5)

CRC_TITreg_genes=df2$gene[df2$sig=='Sig1']

#ccRCC CD4 GSE121636
library(Seurat)
library(ggplot2)
library(harmony)

N1=Read10X('~/Subject/CD177TITreg/data/GSE121636_ccRCC/ccRCC_N1/')
N1=CreateSeuratObject(N1,min.cells=10,min.features=200)
N1$SampleResource=rep("N1",ncol(N1))

N2=Read10X('~/Subject/CD177TITreg/data/GSE121636_ccRCC/ccRCC_N2/')
N2=CreateSeuratObject(N2,min.cells=10,min.features=200)
N2$SampleResource=rep("N2",ncol(N2))

N3=Read10X('~/Subject/CD177TITreg/data/GSE121636_ccRCC/ccRCC_N3/')
N3=CreateSeuratObject(N3,min.cells=10,min.features=200)
N3$SampleResource=rep("N3",ncol(N3))

P1=Read10X('~/Subject/CD177TITreg/data/GSE121636_ccRCC/ccRCC_P1/')
P1=CreateSeuratObject(P1,min.cells=10,min.features=200)
P1$SampleResource=rep("P1",ncol(P1))

P2=Read10X('~/Subject/CD177TITreg/data/GSE121636_ccRCC/ccRCC_P2/')
P2=CreateSeuratObject(P2,min.cells=10,min.features=200)
P2$SampleResource=rep("P2",ncol(P2))

P3=Read10X('~/Subject/CD177TITreg/data/GSE121636_ccRCC/ccRCC_P3/')
P3=CreateSeuratObject(P3,min.cells=10,min.features=200)
P3$SampleResource=rep("P3",ncol(P3))

T1=Read10X('~/Subject/CD177TITreg/data/GSE121636_ccRCC/ccRCC_T1/')
T1=CreateSeuratObject(T1,min.cells=10,min.features=200)
T1$SampleResource=rep("T1",ncol(T1))

T2=Read10X('~/Subject/CD177TITreg/data/GSE121636_ccRCC/ccRCC_T2/')
T2=CreateSeuratObject(T2,min.cells=10,min.features=200)
T2$SampleResource=rep("T2",ncol(T2))

T3=Read10X('~/Subject/CD177TITreg/data/GSE121636_ccRCC/ccRCC_T3/')
T3=CreateSeuratObject(T3,min.cells=10,min.features=200)
T3$SampleResource=rep("T3",ncol(T3))

data3=merge(N1,y=c(N2,N3,P1,P2,P3,T1,T2,T3))
data3[["percent.mt"]] <- PercentageFeatureSet(data3, pattern = "^MT-")

pdf('data3_QC.pdf',width=15,height=6)
VlnPlot(data3,pt.size=0,group.by="SampleResource", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

data3 <- subset(data3, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)

data3<-NormalizeData(data3)
data3<-ScaleData(data3,features=rownames(data3))
data3<-FindVariableFeatures(data3,selection.method="vst",nfeatures=4000)
data3<-RunPCA(data3,features=VariableFeatures(data3),npcs=50)
data3<-RunHarmony(data3, group.by.vars="SampleResource",plot_convergence = F,reduction.save = "harmony",max.iter.harmony = 30) #迭代5次结束
data3<-FindNeighbors(data3,reduction = "harmony",dims=1:50)
data3<-FindClusters(data3,resolution=0.5)
data3<-RunUMAP(data3,dims=1:50,reduction = "harmony")

TsnePlot<-DimPlot(data3,raster=T, group.by="RNA_snn_res.0.5",reduction = "umap",pt.size = 4,label = T,label.size=4,cols=cols,raster.dpi=c(2048,2048))
ggsave("./data3_res0.5_umap.pdf",TsnePlot,width=11.5,height=10,units='cm')

TsnePlot<-DimPlot(data3,raster=T, group.by="RNA_snn_res.1",reduction = "umap",pt.size = 4,label = T,label.size=4,cols=cols,raster.dpi=c(2048,2048))
ggsave("./data3_res1_umap.pdf",TsnePlot,width=11.5,height=10,units='cm')

fea=FeaturePlot(data3,raster=T,features=c('KLRG1','KLRF1','CD3D','CD3E','CD4','CD8A','CD8B','FOXP3','CTLA4'),col=c("lightgray","red"),pt.size=2,reduction="umap",ncol=1,raster.dpi=c(1024,1024))
ggsave(paste0("./data3_","CD4_markers_fea.pdf"),fea,width=11,height=10*9,units='cm',limitsize = FALSE)

tmp=subset(data3,subset=RNA_snn_res.0.5 %in% c(13,8,4,12,2,11))
cells=colnames(tmp)[which(tmp@assays$RNA@data['CD8A',]==0 & tmp@assays$RNA@data['CD8B',]==0)]
CD4T=subset(tmp,cells=cells)
CD4T=subset(CD4T,subset= RNA_snn_res.1 %in% c(0:7,10:11))
CD4T<-NormalizeData(CD4T)
CD4T<-ScaleData(CD4T,features=rownames(CD4T))
CD4T<-FindVariableFeatures(CD4T,selection.method="vst",nfeatures=4000)
CD4T<-RunPCA(CD4T,features=VariableFeatures(CD4T),npcs=50)
CD4T<-RunHarmony(CD4T, group.by.vars="SampleResource",plot_convergence = F,reduction.save = "harmony",max.iter.harmony = 30) #迭代5次结束
CD4T<-FindNeighbors(CD4T,reduction = "harmony",dims=1:50)
CD4T<-FindClusters(CD4T,resolution=1)
CD4T<-RunUMAP(CD4T,dims=1:50,reduction = "harmony")

TsnePlot<-DimPlot(CD4T,raster=T,order=T, group.by="RNA_snn_res.1",reduction = "umap",pt.size = 10,label = T,label.size=4,cols=cols,raster.dpi=c(2048,2048))
ggsave("./CD4T_res1_umap.pdf",TsnePlot,width=11.5,height=10,units='cm')

fea=FeaturePlot(CD4T,raster=T,order=T,features=c('CD4','FOXP3','CTLA4'),col=c("lightgray","red"),pt.size=10,reduction="umap",ncol=1,raster.dpi=c(2048,2048))
ggsave(paste0("./CD4T_","CD4_markers_fea.pdf"),fea,width=11,height=10*3,units='cm',limitsize = FALSE)
VLN=VlnPlot(CD4T,pt.size=0,features =c('KLRD1','XCR2','KLRG1','CD4','CD8A','CD8B','FOXP3','CTLA4'),group.by="RNA_snn_res.1",cols=cols)
ggsave(paste0("./CD4T_","_markers_vln.pdf"),VLN,width=50,height=30,units='cm',limitsize = FALSE)

HCCTreg=subset(CD4T,subset=RNA_snn_res.1 %in% c(7))

badterm <- c(
  "^RPS","^RPL","^MT-","^MTMR","^MTND","^MTRN","^MTCO","^MRPL","^MRPS","^HBA","^HBB","^MTATP",
  "^MALAT1$", "^XIST$", "^XIST_intron$",
  "^AC[0-9]+\\.", "^AP[0-9]+\\.", "^AL[0-9]+\\.", "-AS[0-9]*$", "^MIR[0-9]", "^LINC[0-9]", "^SNOR", "^sno", "^SNR",
  "^TMSB", "^HIST", "^HSP", "^IFI", "^HLA-", "^ATP", "-", ":", "\\.", '^KIAA',
  "^DNAJ", '^FTH', '^FTL', '^LGALS'
)
badgenelt<-unique(grep(pattern=paste(c(badterm),collapse="|"),x=rownames(HCCTreg),perl=T,value=T))

HCCTreg$SampleType=ifelse(HCCTreg$SampleResource %in% paste0('T',1:3),'Tumor',
                    ifelse(HCCTreg$SampleResource %in% paste0('P',1:3),'PBMC','Normal'))

diff_rs1=FindMarkers(HCCTreg,features=rownames(HCCTreg)[-which(rownames(HCCTreg) %in% badgenelt)],ident.1="Tumor",ident.2='PBMC',group.by="SampleType",logfc.threshold=0,min.pct=0.15,only.pos=F,pseudocount.use=F)
diff_rs1=diff_rs1[order(diff_rs1$avg_log2FC,decreasing=T),]
diff_rs1$gene=rownames(diff_rs1)
diff_rs1$avg_log2FC[1:12]=c(7.5,7.5)

CD4T$type=ifelse(CD4T$RNA_snn_res.1 %in% c(7),'Treg','Others')
CD4T$SampleType=ifelse(CD4T$SampleResource %in% paste0('T',1:3),'Tumor',
                    ifelse(CD4T$SampleResource %in% paste0('P',1:3),'PBMC','Normal'))
tmp=subset(CD4T,subset=SampleType == 'Tumor')
diff_rs2=FindMarkers(tmp,features=rownames(tmp)[-which(rownames(tmp) %in% badgenelt)],ident.1="Treg",ident.2='Others',group.by="type",logfc.threshold=0,min.pct=0.15,only.pos=F,pseudocount.use=F)
diff_rs2=diff_rs2[order(diff_rs2$avg_log2FC,decreasing=T),]
diff_rs2$gene=rownames(diff_rs2)

inte_genes=intersect(diff_rs1$gene,diff_rs2$gene)

diff_rs1=diff_rs1[inte_genes,]
diff_rs2=diff_rs2[inte_genes,]

df=data.frame(
    avg_log2FC_TP=diff_rs1$avg_log2FC,
    p_adj_TP=diff_rs1$p_val_adj,
    avg_log2FC_TO=diff_rs2$avg_log2FC,
    p_adj_TO=diff_rs2$p_val_adj,
    gene=inte_genes
)
df$sig=ifelse(df$p_adj_TP<0.001 & df$p_adj_TO<0.001 & df$avg_log2FC_TP>=2 & df$avg_log2FC_TO>=2,"Sig1",
        ifelse(df$p_adj_TP<0.001 & df$p_adj_TO<0.001 & df$avg_log2FC_TP<=(-2) & df$avg_log2FC_TO<=(-2),"Sig2", "Nosig"))
labels=c(
    #intersect(diff_rs1$gene[diff_rs1$p_val_adj<0.001][1:20],diff_rs2$gene[diff_rs2$p_val_adj<0.001][1:20]),
    #intersect(tail(diff_rs1[diff_rs1$p_val_adj<0.001,],20)$gene,tail(diff_rs2[diff_rs2$p_val_adj<0.001,],20)$gene)
    'CD177','IL1R2','TNFRSF13B','WLS','TNFRSF9','LAYN','CLEC7A','CCR8','TSPAN13','IL1R1','TNS3',
    'ANXA1','CCR7','KLF2','KLF3','S1PR1','TC2N','RUNX2'
)
df$label=ifelse(df$gene %in% labels,df$gene,"")
#df$label=ifelse(df$sig %in% c('Sig1','Sig2'),df$gene,"")
p<-ggplot(df,aes(x=avg_log2FC_TP,y=avg_log2FC_TO,color=sig,label=label))+
    ggrastr::rasterize(geom_point(size=2),dpi=150)+
    ggrepel::geom_text_repel(size=6,max.overlaps=Inf,nudge_x=0.25,nudge_y=0.25)+
    geom_hline(yintercept = c(2,-2), color="grey",linetype = 2,size = 0.5)+
    geom_vline(xintercept = c(2,-2), color="grey",linetype = 2,size = 0.5)+
    scale_color_manual(values=c("Sig1"="#EC478B",'Sig2'='#3E94D6',"NoSig"="lightgrey"))+
    xlab("avg_log2FC_TP")+ylab("avg_log2FC_TO")+
    theme_bw()+theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0.5),panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
ggsave("ccRCCTreg_TP_TO.pdf",p,width=5.5,height=4.5)

ccRCC_TITreg_genes=df3$gene[df3$sig=='Sig1']

#inte TITreg_genes
HCCTreg=readRDS('../HCC/HCC_Treg.RDS')
CRCTreg=readRDS('../CRC/CRC_Treg.RDS')
ccRCCTreg=readRDS('../ccRCC/ccRCC_Treg.RDS')

VLN=VlnPlot(HCCTreg,pt.size=0,features =c('CD177','LAYN','CCR8'),group.by="sampleType",cols=ggsci::pal_futurama("planetexpress")(5)[c(1,3,4)])
ggsave(paste0("./HCCTreg_integenes","_markers_vln.pdf"),VLN,width=50,height=30,units='cm',limitsize = FALSE)

VLN=VlnPlot(CRCTreg,pt.size=0,features =c('CD177','LAYN','CCR8'),group.by="Tissue",cols=ggsci::pal_futurama("planetexpress")(5)[c(1,3,4)])
ggsave(paste0("./CRCTreg_integenes","_markers_vln.pdf"),VLN,width=50,height=30,units='cm',limitsize = FALSE)

VLN=VlnPlot(ccRCCTreg,pt.size=0,features =c('CD177','LAYN','CCR8'),group.by="SampleType",cols=ggsci::pal_futurama("planetexpress")(5)[c(1,3,4)])
ggsave(paste0("./ccRCCTreg_integenes","_markers_vln.pdf"),VLN,width=50,height=30,units='cm',limitsize = FALSE)

plot_list=lapply(c('CD177','LAYN','CCR8'),function(x){
    df=data.frame(
        value=as.numeric(as.vector(HCCTreg@assays$RNA@data[x,])),
        group=factor(HCCTreg$sampleType,levels=c('PTR','NTR','TTR'))
    )
    ggplot(df,aes(x=group,y=value,color=group))+
        ggrastr::rasterise(geom_point(position = position_jitter(width = 0.3, height = 0),size=1,alpha = 0.4, shape = 21),dpi=150)+
        geom_boxplot(outlier.shape = NA,width=0.4,alpha=0)+
        scale_color_manual(values=ggsci::pal_futurama("planetexpress")(4)[c(1,3,4)])+
        ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test",comparisons = list(c("TTR","NTR"),c('TTR','PTR'),c('PTR','NTR')),size=4)+
        theme_bw()+theme(text = element_text(size = 8),axis.text.x = element_text(angle = 315, hjust = 0, vjust = 0),panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"), legend.position="bottom") +
        labs(x="", y = 'Exp.',title=x)
})
pdf('HCCTreg_integenes_boxplot.pdf',width=7,height=4)
cowplot::plot_grid(plotlist=plot_list,nrow=1)
dev.off()

plot_list=lapply(c('CD177','LAYN','CCR8'),function(x){
    df=data.frame(
        value=as.numeric(as.vector(CRCTreg@assays$RNA@data[x,])),
        group=factor(CRCTreg$Tissue,levels=c('P','N','T'))
    )
    ggplot(df,aes(x=group,y=value,color=group))+
        ggrastr::rasterise(geom_point(position = position_jitter(width = 0.3, height = 0),size=1,alpha = 0.4, shape = 21),dpi=150)+
        geom_boxplot(outlier.shape = NA,width=0.4,alpha=0)+
        scale_color_manual(values=ggsci::pal_futurama("planetexpress")(4)[c(1,3,4)])+
        ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test",comparisons = list(c("T","N"),c('T','P'),c('P','N')),size=4)+
        theme_bw()+theme(text = element_text(size = 8),axis.text.x = element_text(angle = 315, hjust = 0, vjust = 0),panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"), legend.position="bottom") +
        labs(x="", y = 'Exp.',title=x)
})
pdf('CRCTreg_integenes_boxplot.pdf',width=7,height=4)
cowplot::plot_grid(plotlist=plot_list,nrow=1)
dev.off()

plot_list=lapply(c('CD177','LAYN','CCR8'),function(x){
    df=data.frame(
        value=as.numeric(as.vector(ccRCCTreg@assays$RNA@data[x,])),
        group=factor(ccRCCTreg$SampleType,levels=c('PBMC','Tumor'))
    )
    ggplot(df,aes(x=group,y=value,color=group))+
        ggrastr::rasterise(geom_point(position = position_jitter(width = 0.3, height = 0),size=1,alpha = 0.4, shape = 21),dpi=150)+
        geom_boxplot(outlier.shape = NA,width=0.4,alpha=0)+
        scale_color_manual(values=ggsci::pal_futurama("planetexpress")(4)[c(1,4)])+
        ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test",comparisons = list(c("Tumor","PBMC")),size=4)+
        theme_bw()+theme(text = element_text(size = 8),axis.text.x = element_text(angle = 315, hjust = 0, vjust = 0),panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"), legend.position="bottom") +
        labs(x="", y = 'Exp.',title=x)
})
pdf('ccRCCTreg_integenes_boxplot.pdf',width=4,height=4)
cowplot::plot_grid(plotlist=plot_list,nrow=1)
dev.off()

##########################################################################################
#Figure2 Pancancer CD177+Treg https://www.science.org/doi/10.1126/science.abe6474
data2=readRDS('../data/PanCancer_TIT/data2.RDS')

Treg=subset(data2,subset=RNA_snn_res.0.5 %in% c('0','7','11','15'))

library(harmony)
Treg<-NormalizeData(Treg)
Treg<-FindVariableFeatures(Treg,selection.method="vst",nfeatures=4000)
all.genes<-rownames(Treg)
Treg<-ScaleData(Treg,features=all.genes)
Treg<-RunPCA(Treg,features=VariableFeatures(Treg),npcs=50)
Treg<-RunHarmony(Treg, group.by.vars="patient",plot_convergence = F,reduction.save = "harmony",max.iter.harmony = 30) #迭代5次结束
Treg<-FindNeighbors(Treg,reduction = "harmony",dims=1:50)
Treg<-FindClusters(Treg,resolution=0.5)
Treg<-RunUMAP(Treg,dims=1:50,reduction = "harmony")

Treg$RNA_snn_res.0.5=factor(Treg$RNA_snn_res.0.5,levels=0:7)
Treg$Cluster=paste0('C',Treg$RNA_snn_res.0.5)
clucol=c('#387EB8','#4DAF4A','#984EA2','#FF7F00','#FEFF33',"#FF1393",'#FFFFBF','#F781C0')#,'#BF5B17','#8DD3C7','#01FFFF','#FF1393')
names(clucol)=paste0('C',0:7)

TsnePlot<-DimPlot(Treg,raster=T, group.by="Cluster",reduction = "umap",pt.size = 2,label = T,label.size=6,cols=clucol,raster.dpi=c(1024,1024))
ggsave("./Treg_Cluster_umap.pdf",TsnePlot,width=12,height=10,units='cm')

TsnePlot<-DimPlot(Treg,raster=T, group.by="loc",reduction = "umap",pt.size = 2,label = T,label.size=6,cols=cols,raster.dpi=c(1024,1024))
ggsave("./Treg_loc_umap.pdf",TsnePlot,width=12,height=10,units='cm')

TsnePlot<-DimPlot(Treg,raster=T, group.by="accessory",reduction = "umap",pt.size = 2,label = T,label.size=6,cols=cols,raster.dpi=c(1024,1024))
ggsave("./Treg_accessory_umap.pdf",TsnePlot,width=12,height=10,units='cm')

TsnePlot<-DimPlot(Treg,raster=T, group.by="stype",reduction = "umap",pt.size = 2,label = T,label.size=6,cols=cols,raster.dpi=c(1024,1024))
ggsave("./Treg_stype_umap.pdf",TsnePlot,width=12,height=10,units='cm')

library("scales")
Treg$CD3D_scale <- rescale(Treg$CD3D, to = c(0, 5))
Treg$IL2RA_scale <- rescale(Treg$IL2RA, to = c(0, 5))
Treg$FOXP3_scale <- rescale(Treg$FOXP3, to = c(0, 5))
Treg$CD3E_scale <- rescale(Treg$CD3E, to = c(0, 5))
Treg$CD3G_scale <- rescale(Treg$CD3G, to = c(0, 5))
Treg$CD4_scale <- rescale(Treg$CD4, to = c(0, 5))
Treg$CD8A_scale <- rescale(Treg$CD8A, to = c(0, 5))
Treg$CD8B_scale <- rescale(Treg$CD8B, to = c(0, 5))

fea=FeaturePlot(Treg,raster=T,order=T,features=c('CD177','BATF'),col=c("lightgray","red"),pt.size=2,reduction="umap",ncol=1,raster.dpi=c(1024,1024))
ggsave(paste0("./Treg","_markers_fea.pdf"),fea,width=11,height=10*2,units='cm',limitsize = FALSE)

fea=FeaturePlot(Treg,raster=T,order=T,features=c('rna_CD8A','rna_CD8B','rna_CD4','rna_FOXP3','rna_CTLA4'),col=c("lightgray","red"),pt.size=2,reduction="umap",ncol=1,raster.dpi=c(1024,1024))
ggsave(paste0("./Treg","_markers2_fea.pdf"),fea,width=11,height=10*5,units='cm',limitsize = FALSE)

VLN=VlnPlot(Treg,pt.size=0,sort=T,features =c('CD177','BATF'),group.by="RNA_snn_res.0.5",cols=cols)
ggsave(paste0("./Treg_",'CD177',"_markers_vln.pdf"),VLN,width=20,height=10,units='cm',limitsize = FALSE)

VLN=VlnPlot(Treg,pt.size=0,sort=T,features =c('CD177','BATF'),group.by="loc",cols=cols)
ggsave(paste0("./Treg_",'CD177',"_loc_markers_vln.pdf"),VLN,width=20,height=10,units='cm',limitsize = FALSE)

VLN=VlnPlot(Treg,pt.size=0,sort=T,features =c('CD177','BATF'),group.by="cancerType",cols=cols)
ggsave(paste0("./Treg_",'CD177',"_cancerType_markers_vln.pdf"),VLN,width=20,height=10,units='cm',limitsize = FALSE)

#Markers
badterm <- c(
  "^RPS","^RPL","^MT-","^MTMR","^MTND","^MTRN","^MTCO","^MRPL","^MRPS","^HBA","^HBB","^MTATP",
  "^MALAT1$", "^XIST$", "^XIST_intron$",
  "^AC[0-9]+\\.", "^AP[0-9]+\\.", "^AL[0-9]+\\.", "-AS[0-9]*$", "^MIR[0-9]", "^LINC[0-9]", "^SNOR", "^sno", "^SNR",
  "^TMSB", "^HIST", "^HSP", "^IFI", "^HLA-", "^ATP", "-", ":", "\\.", '^KIAA',
  "^DNAJ", '^FTH', '^FTL', '^LGALS'
)
badgenelt<-unique(grep(pattern=paste(c(badterm),collapse="|"),x=rownames(Treg),perl=T,value=T))

future::plan("multicore", workers = 5)
library(dplyr)
Treg@active.ident<-factor(Treg$RNA_snn_res.0.5)
dt.markers <- FindAllMarkers(Treg, features=rownames(Treg)[-which(rownames(Treg) %in% badgenelt)],only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dt.markers$pct_diff<-dt.markers$pct.1-dt.markers$pct.2
dt.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) ->top50

write.table(top50,'Treg_Cluster_diffgenes_top50.txt',quote=F,row.names=F,col.names=T,sep='\t')

selected_genes=c(
    'TNFRSF9','TNFRSF4','TNFRSF18','CCR8','ICOS','CD177','LAYN',
    'LEF1','TXNIP','SELL','CCR7',
    'DUSP1','FOS','RGS1','CXCR4',
    'CXCL13','STMN1','TOP2A','MKI67','TYMS',
    'RSAD2','OAS3','MX2','MX1','OAS1',
    'NDFIP2','CISH','IL23R','CTSH','CCR5',
    'COCH','THADA','CTSL','CADM1','IKZF4',
    'TCF7','TC2N','CD40LG','GIMAP4','ANXA1'
)
Treg$type<-paste0(Treg$Cluster,"_",Treg$loc)
order_type<-c()
for(i in levels(factor(Treg$Cluster))){
    order_type<-c(order_type,paste0(i,"_",c('P','N','T')))
}

mat<-Treg@assays$RNA@data
a<-sapply(selected_genes,function(x){
    sapply(order_type,function(y){
        mean(as.numeric(as.vector(mat[x,which(Treg$type==y)])))
    })
})
t(a)->a

annotation_col = data.frame(
  Sample = factor(rep(c('P','N','T'),8),levels=c('P','N','T')),
  CellType=factor(rep(levels(factor(Treg$Cluster)),rep(3,8)),levels=levels(factor(Treg$Cluster)))
)
rownames(annotation_col) = colnames(a)

ann_colors = list(
  Sample = samplecol,
  CellType=clucol
)

library(pheatmap)
pdf("Treg_Cluster_slectedgenes_pheatmap.pdf",width=10,height=10)
pheatmap(a,scale="row",annotation_colors=ann_colors,annotation_col=annotation_col,cellwidth=15,cellheight=15,border_color=NA,angle_col=315,cluster_rows=F,cluster_cols=F,show_rownames=T,show_colnames=F,color=rev(heatmaply::RdBu(100)))
dev.off()


#CellType number
df=data.frame(
  value=as.numeric(as.vector(table(Treg$Cluster))),
  celltype=names(table(Treg$Cluster))
)
df2=df %>% arrange(desc(value))
df$celltype=factor(df$celltype,levels=rev(df2$celltype))
p=ggplot(df,aes(x=celltype,y=value))+
    geom_bar(stat = 'identity',fill='#80E6EC',width=0.5)+
    coord_flip() + 
    theme_bw()+theme(text = element_text(size = 8),axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0),panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"), legend.position="bottom") +
    labs(x="", y = "Cell Nums")
ggsave("CellType_CellNumbers.pdf",p,width=3,height=3)

#SampleType
sapply(levels(factor(Treg$loc)),function(x){
    sapply(levels(factor(Treg$Cluster)),function(y){
        round(length(which(Treg$loc==x & Treg$Cluster==y))/length(which(Treg$Cluster==y)),4)
    })
})->tmp
df<-data.frame(value=as.numeric(as.vector(unlist(tmp))),
			  sampletype=factor(rep(colnames(tmp),rep(nrow(tmp),ncol(tmp))),levels=c("P","N",'T')),
			  celltype=factor(rep(rownames(tmp),ncol(tmp)),levels=levels(factor(unique(Treg$Cluster)))))

df$celltype=factor(df$celltype,levels=rev(df2$celltype))
samplecol=ggsci::pal_futurama("planetexpress")(3)[c(3,1,2)]
names(samplecol)=levels(df$sampletype)

p=ggplot(df,aes(x=celltype,y=value,fill=sampletype))+
    geom_bar(stat = 'identity',width=0.5)+
    coord_flip() + 
    scale_fill_manual(values=samplecol)+
    theme_bw()+theme(text = element_text(size = 10),axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0),panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"), legend.position="bottom") +
    labs(x="", y = "Fraction")
ggsave("CellType_sampletype_fraction.pdf",p,width=3,height=4)

#cancertype
sapply(levels(factor(Treg$cancerType)),function(x){
    sapply(levels(factor(Treg$Cluster)),function(y){
        round(length(which(Treg$cancerType==x & Treg$Cluster==y))/length(which(Treg$Cluster==y)),4)
    })
})->tmp
df<-data.frame(value=as.numeric(as.vector(unlist(tmp))),
			  sampletype=factor(rep(colnames(tmp),rep(nrow(tmp),ncol(tmp))),levels=colnames(tmp)),
			  celltype=factor(rep(rownames(tmp),ncol(tmp)),levels=levels(factor(unique(Treg$Cluster)))))

df$celltype=factor(df$celltype,levels=rev(df2$celltype))
tumorcol=ggsci::pal_futurama("planetexpress")(7)[-c(1:3)]
names(tumorcol)=levels(df$sampletype)

p=ggplot(df,aes(x=celltype,y=value,fill=sampletype))+
    geom_bar(stat = 'identity',width=0.5)+
    coord_flip() + 
    scale_fill_manual(values=tumorcol)+
    theme_bw()+theme(text = element_text(size = 10),axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0),panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"), legend.position="bottom") +
    labs(x="", y = "Fraction")
ggsave("CellType_tumortype_fraction.pdf",p,width=3,height=4)

#
colours = colorRampPalette(c("grey89", "#FEF8E5","#FFEABC","#FDBC52","#F68523","#DD3226","#A31B3A","#5E2C81","#382F85","#28316C"))(1000)
fea=FeaturePlot(Treg,raster=T,order=T,features=c('CD177'),col=c("lightgray","red"),pt.size=2,reduction="umap",ncol=1,raster.dpi=c(1024,1024))+scale_color_gradientn(colors=colours)
ggsave(paste0("./Treg","_CD177_fea.pdf"),fea,width=11,height=10*1,units='cm',limitsize = FALSE)

VLN=VlnPlot(Treg,pt.size=0,sort=T,features =c('CD177'),group.by="Cluster",cols=clucol)
ggsave(paste0("./Treg_",'CD177',"_markers_vln.pdf"),VLN,width=12,height=8,units='cm',limitsize = FALSE)

C0=subset(Treg,subset=Cluster=='C0')
VLN=VlnPlot(C0,pt.size=0,sort=T,features =c('CD177'),group.by="loc",cols=samplecol)
ggsave(paste0("./C0_",'CD177',"_markers_vln.pdf"),VLN,width=12,height=8,units='cm',limitsize = FALSE)

x='CD177'
tmp_df=data.frame(
    value=as.numeric(as.vector(Treg@assays$RNA@data[x,])),
    sampletype=Treg$loc
)
p=ggplot(tmp_df,aes(x=sampletype,y=value,color=sampletype))+
        ggrastr::rasterise(geom_point(position = position_jitter(width = 0.3, height = 0),size=1,alpha = 0.5, shape = 21),dpi=150)+
        geom_boxplot(outlier.shape = NA,width=0.4,alpha=0)+
        scale_color_manual(values=samplecol)+
        ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test",comparisons = list(c('N','T'),c('N','P'),c('P','T')),size=3)+
        theme_bw()+theme(text = element_text(size = 8),axis.text.x = element_text(angle = 315, hjust = 0, vjust = 0),panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"), legend.position="bottom") +
        labs(x="", y = 'Expression',title=x)
ggsave('Treg_CD177_sampletype_boxplot.pdf',p,width=5,height=6)

#
library(ggpubr)

EffectorMolecules<-c('GZMB','ITGAE')
TFs<-c('BATF','PRDM1','BACH2','SATB1','rna_FOXP3')
Inhibitory_molecules<-c('CTLA4','GZMA','TNFSF14','GZMM')
Co_Stimulatory_Molecules<-c('CD27','CD70','CD177','ICOS','TNFRSF18','TNFRSF25','TNFRSF14','SEMA4D')
Activited<-c('IL7R','rna_IL2RA','CD69')
cytokines=c('IL10','CXCL13','TGFB1')

Naive=c('TCF7','CCR7','SELL','LEF1')
Activited<-c('IL7R','rna_IL2RA','CD69')
TFs<-c('BATF','PRDM1','BACH2','SATB1','rna_FOXP3')
Resident<-c('ITGAE','ITGA1','CXCR6',"CX3CR1",'CRTAM','ZNF683')
Co_Stimulatory_Molecules<-c('CD27','CD70','ICOS','TNFRSF18','TNFRSF25','TNFRSF14')

key_molecules=c(Naive,Activited,TFs,Resident,Co_Stimulatory_Molecules)
Treg$Cluster=factor(Treg$Cluster,levels=paste0('C',c(0,5,2,4,6,3,7,1)))
p<-DotPlot(Treg,cluster.idents=F,group.by="Cluster",features=key_molecules,cols=c("darkblue","yellow"))+scale_color_gradientn(colours=rev(heatmaply::RdBu(100)))+
    xlab("")+ylab("")+ggpubr::rotate_x_text(angle=45)
ggsave('Treg_function_genes_dotplot2.pdf',p,width=8,height=4)

library(ComplexHeatmap)
df=p$data
exp_mat<-df %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 
    
row.names(exp_mat) <- exp_mat$features.plot  
exp_mat <- exp_mat[,-1] %>% as.matrix()
percent_mat<-df %>% 
  select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() 
    
row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()
col_fun = circlize::colorRamp2(c(-1, 0, 2), heatmaply::RdBu(100)[c(1,30, 100)])

cell_fun = function(j, i, x, y, w, h, fill){
          grid.rect(x = x, y = y, width = w, height = h, 
                    gp = gpar(col = NA, fill = NA))
          grid.circle(x=x,y=y,r= percent_mat[i, j]/100 * min(unit.c(w, h)),
                      gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}
pdf('Treg_function_genes_ComplexHeatmap')
ComplexHeatmap::Heatmap(exp_mat,
        heatmap_legend_param=list(title="expression"),
        column_title = "clustered dotplot", 
        col=col_fun,
        cluster_rows = FALSE,
        clustering_method_columns='complete',
        cluster_columns = TRUE,
        rect_gp = gpar(type = "none"),
        cell_fun = cell_fun,
        border = "black")
dev.off()

#
Inhibitory=c('ENTPD1','TIGIT','IL10','LAG3','NT5E','HAVCR2','TGFB1','CTLA4','GZMA','TNFSF14','GZMM')
Resident<-c('ITGAE','ITGA1','CXCR6',"CX3CR1",'CRTAM','ZNF683')
Naive=c('TCF7','CCR7','SELL','LEF1')
Activited<-c('IL7R','rna_IL2RA','CD69')
Proliferation=c('STMN1','TOP2A','MKI67','CCNB1','TYMS','HMGB2','PCNA','DEK')

genelt=list(
    Inhibitory=Inhibitory,
    Resident=Resident,
    Naive=Naive,
    Activited=Activited,
    Proliferation=Proliferation
)
Treg=pochi::RunAUCell(Treg,assay="RNA",slot="counts",genesets=genelt, verbose=T)

Treg$Inhibitory=Treg@assays$AUC@data['Inhibitory',]
Treg$Resident=Treg@assays$AUC@data['Resident',]
Treg$Naive=Treg@assays$AUC@data['Naive',]
Treg$Activited=Treg@assays$AUC@data['Activited',]
Treg$Proliferation=Treg@assays$AUC@data['Proliferation',]

library(ggradar)
mean_rs=data.frame(
    CellType=Treg$Cluster,
    Inhibitory=Treg$Inhibitory,
    Resident=Treg$Resident,
    Naive=Treg$Naive,
    Activited=Treg$Activited,
    Proliferation=Treg$Proliferation
)
df=mean_rs %>% group_by(CellType) %>% summarise(Naive=mean(Naive),Inhibitory=mean(Inhibitory),Resident=mean(Resident),Activited=mean(Activited),Proliferation=mean(Proliferation))
scale_df=as.data.frame(apply(df[,-1],2,function(x) (x-min(x))/(max(x)-min(x))))
scale_df=scale_df %>% as_tibble(rownames = "group")
scale_df$group=df$CellType

pdf('Treg_cluster_5subtype_ggradarplot.pdf')
ggradar(
  scale_df,
  values.radar=c('0','0.5','1'),grid.min=0,grid.mid=0.5,grid.max=1,
  group.line.width=0.5,group.point.size=2,group.colours=clucol,
  background.circle.colour='white',gridline.mid.colour='grey',legend.position='bottom'
)
dev.off()

#function KEGG
badterm <- c(
  "^RPS","^RPL","^MT-","^MTMR","^MTND","^MTRN","^MTCO","^MRPL","^MRPS","^HBA","^HBB","^MTATP",
  "^MALAT1$", "^XIST$", "^XIST_intron$",
  "^AC[0-9]+\\.", "^AP[0-9]+\\.", "^AL[0-9]+\\.", "-AS[0-9]*$", "^MIR[0-9]", "^LINC[0-9]", "^SNOR", "^sno", "^SNR",
  "^TMSB", "^HIST", "^HSP", "^IFI", "^HLA-", "^ATP", "-", ":", "\\.", '^KIAA',
  "^DNAJ", '^FTH', '^FTL', '^LGALS'
)
badgenelt<-unique(grep(pattern=paste(c(badterm),collapse="|"),x=rownames(Treg),perl=T,value=T))

future::plan("multicore", workers = 5)
library(dplyr)
Treg@active.ident<-factor(Treg$Cluster)
dt.markers <- FindAllMarkers(Treg, features=rownames(Treg)[-which(rownames(Treg) %in% badgenelt)],only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dt.markers$pct_diff<-dt.markers$pct.1-dt.markers$pct.2
dt.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) ->top50

dt.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) ->top100

Enrich_analyse<-function(genelt,NumPathways=5,db=c('human','mouse'),enrichFunc=c('enrich','gse','enrichr'),DataBase=c('KEGG','GO','GMT'),GMT.dir=NULL){
    require(clusterProfiler)
    require(org.Hs.eg.db)
    require(fmsb)
    require(ggplot2)
    require(org.Mm.eg.db)

    if(db=="human"){
        OrgDb="org.Hs.eg.db"
        organism="hsa"
    }else if(db=="mouse"){
        OrgDb="org.Mm.eg.db"
        organism="mmu"
    }else{
        stop("DataBase ERROR!!!")
    }

    genelt=as.character(as.vector(unlist(genelt)))
    if(enrichFunc=="gse" & DataBase=="GO"){
        br<-bitr(genelt,fromType="SYMBOL",toType="ENTREZID",OrgDb=OrgDb)
        ge<-rev(1:length(br$ENTREZID))
        names(ge)<-br$ENTREZID
        ego<-gseGO(gene= ge,OrgDb= OrgDb,
            ont= "BP",pvalueCutoff  = 1,minGSSize=10,maxGSSize=500,verbose=T)
    }else if(enrichFunc=="gse" & DataBase=="KEGG"){
        br<-bitr(genelt,fromType="SYMBOL",toType="ENTREZID",OrgDb=OrgDb)
        ge<-rev(1:length(br$ENTREZID))
        names(ge)<-br$ENTREZID
        ego<-gseKEGG(gene= ge,organism= organism,
            pvalueCutoff  = 1,minGSSize=10,verbose=T)       
    }else if(enrichFunc=="enrich" & DataBase=="GO"){
        br <- bitr(genelt, fromType = "SYMBOL",toType = c("ENSEMBL", "ENTREZID"),OrgDb = OrgDb)
        ego <- enrichGO(gene= br$ENTREZID,OrgDb= OrgDb,ont= "BP",pAdjustMethod = "fdr",
            pvalueCutoff  = 1,qvalueCutoff  = 1,readable=TRUE)
    }else if(enrichFunc=="enrich" & DataBase=="KEGG"){
        br <- bitr(genelt, fromType = "SYMBOL",toType = c("ENSEMBL", "ENTREZID"),OrgDb = OrgDb)
        ego <- enrichKEGG(gene= br$ENTREZID,organism=organism,keyType= "kegg",pAdjustMethod = "fdr",
            pvalueCutoff  = 1,qvalueCutoff  = 1)
    }else if(enrichFunc=="enrichr" & DataBase=="GMT"){
        gmt<-clusterProfiler::read.gmt(GMT.dir)
        ego<-clusterProfiler::enricher(gene=genelt,pvalueCutoff=1,pAdjustMethod="fdr",minGSSize=10,qvalueCutoff =1,TERM2GENE=gmt)
    }else{
        print("Input wrong parameters......")
        break;
    }

    ego@result<-ego@result[!grepl("Melanogenesis|Viral|Bacterial|Influenza|Measles|sclerosis|diabetic|cardiomyocytes|Malaria|cardiomyopathy|myocarditis|lupus|Leishmaniasis|Asthma|arthritis|disease|cancer|leukemia|carcinoma|glioma|melanoma|infection|immunodeficiency|Amoebiasis|Hematopoietic|end|Hepatitis|carcinogenesis|diabetes|syndrome",ignore.case=T,ego@result$Description),]
    ego@result<-ego@result[which(ego@result$pvalue<0.05),]
    res1<-ego@result
    if(dim(res1)[1]<=NumPathways){
        res1<-res1
    }else{
        res1<-res1[1:NumPathways,]
    }
    df<-data.frame(value= -log10(res1$pvalue),pathway=factor(res1$Description,levels=rev(res1$Description)))
    return(df)
}

source("~/Biotools/Enrich_Radarchart.R")
tmp=lapply(levels(top50$cluster),function(x){
  genes=top50$gene[top50$cluster==x & top50$p_val_adj<0.001]
  Enrich_analyse(genelt=genes,NumPathways=5,enrichFunc=c('enrich'),db="human",DataBase=c('KEGG'))
})
df=reshape2::melt(tmp)
df$L1=rep(levels(top50$cluster),rep(5,length(levels(top50$cluster))))

tmp=sapply(unique(df$pathway),function(x){
  sapply(unique(df$L1),function(y){
    if(length(which(df$pathway==x & df$L1==y))!=0){
      df$value[df$pathway==x & df$L1==y]
    }else{
      0
    }
  })
})
colnames(tmp)=unique(df$pathway)
tmp=t(tmp)
library(pheatmap)
pdf("Treg_Cluster_top80_EnrichKEGG_pheatmap.pdf",width=8,height=8)
pheatmap(tmp,scale="row",cellwidth=15,cellheight=15,border_color=NA,angle_col=315,cluster_rows=T,cluster_cols=T,show_rownames=T,show_colnames=T,color=rev(heatmaply::RdBu(100)))
dev.off()

#Monocle2
devtools::load_all("~/R/library/4.2.0/monocle")

dat <- Seurat::as.CellDataSet(Treg)
dat <- estimateSizeFactors(dat)
dat <- detectGenes(dat, min_expr = 10)
fData(dat)$use_for_ordering <-fData(dat)$num_cells_expressed > 0.05 * ncol(dat)
expressed_genes <- row.names(subset(fData(dat),num_cells_expressed >= 10))
dat <- reduceDimension(dat,
                          max_components = 2,
                          norm_method = 'log',
                          num_dim = 20,
                          reduction_method = 'tSNE',
                          verbose = T,
                          check_duplicates=F)
dat <- clusterCells(dat,verbose = F)
clustering_DEG_genes <- differentialGeneTest(dat[expressed_genes,],
           fullModelFormulaStr = '~Cluster',
           cores = 10)

#dat <- setOrderingFilter(dat, ordering_genes = clustering_DEG_genes$gene_short_name[clustering_DEG_genes$use_for_ordering=="TRUE"])
dat <- setOrderingFilter(dat, ordering_genes = unique(top50$gene))
#dat <- setOrderingFilter(dat, ordering_genes = row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000])
dat <- reduceDimension(dat, method = 'DDRTree')
dat <- orderCells(dat)
saveRDS(dat,"Treg_monocle_top50.RDS")

dat$Cluster=Treg$Cluster
p<-plot_cell_trajectory(dat, cell_size=1, show_backbone=FALSE,cell_link_size = 1,color_by="Cluster",theta=0.8)+
	facet_wrap(~Cluster, nrow = 3)+scale_color_manual(values=clucol)
ggsave('./Treg_Cluster_wrap_top50.pdf',p)
p<-plot_cell_trajectory(dat, cell_size=0.5,show_tree=T, show_backbone=FALSE,cell_link_size = 1,color_by="Cluster")+
	scale_color_manual(values=clucol)
ggsave('./Treg_Cluster_top50.pdf',p,width=6,height=6)

samplecol=ggsci::pal_futurama("planetexpress")(3)[c(3,1,2)]
names(samplecol)=levels(factor(Treg$loc))
p<-plot_cell_trajectory(dat, cell_size=1, show_backbone=FALSE,cell_link_size = 1,color_by="loc")+
    scale_color_manual(values=samplecol)
ggsave('./Treg_SampleType_top50.pdf',p,useDingbat=F)

p<-plot_cell_trajectory(dat, cell_size=1, show_backbone=FALSE,cell_link_size = 1,color_by="Pseudotime")+viridis::scale_color_viridis()
ggsave('./Treg_Pseudotime_top50.pdf',p,useDingbat=F)

p<-plot_cell_trajectory(dat, cell_size=1,show_tree=T, show_backbone=FALSE,cell_link_size = 1,color_by="State",theta=0.8)+
	scale_color_manual(values=jjAnno::useMyCol('calm',n=6))
ggsave('./Treg_State_top50.pdf',p,width=6,height=6)

tmp=sapply(1:5,function(x){
  sapply(levels(factor(dat$Cluster)),function(y){
    length(which(dat$Cluster==y & dat$State==x))
  })
})
colnames(tmp)=paste0('State',1:5)

plot_lists=lapply(colnames(tmp),function(x){
  df=data.frame(
    value=as.numeric(as.vector(tmp[,x])),
    group=factor(rownames(tmp),levels=levels(factor(dat$Cluster)))
  )
  ggplot(df, aes(x="",y=value,fill=group))+
  geom_bar(stat="identity",width = 1)+
  scale_fill_manual(values=clucol)+labs(title=x)+
  coord_polar("y")+theme_minimal()+theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"))+Seurat::NoLegend()

})
pdf('Monocle2_State_CellType_pie.pdf')
cowplot::plot_grid(plotlist=plot_lists,nrow=1)
dev.off()

#
badterm <- c(
  "^RPS","^RPL","^MT-","^MTMR","^MTND","^MTRN","^MTCO","^MRPL","^MRPS","^HBA","^HBB","^MTATP",
  "^MALAT1$", "^XIST$", "^XIST_intron$",
  "^AC[0-9]+\\.", "^AP[0-9]+\\.", "^AL[0-9]+\\.", "-AS[0-9]*$", "^MIR[0-9]", "^LINC[0-9]", "^SNOR", "^sno", "^SNR",
  "^TMSB", "^HIST", "^HSP", "^IFI", "^HLA-", "^ATP", "-", ":", "\\.", '^KIAA',
  "^DNAJ", '^FTH', '^FTL', '^LGALS'
)
badgenelt<-unique(grep(pattern=paste(c(badterm),collapse="|"),x=rownames(Treg),perl=T,value=T))

ct=as.matrix(Treg@assays$RNA@counts)[-which(rownames(Treg) %in% badgenelt),]
pct=apply(ct,1,function(x){
    length(which(x>0))/ncol(ct)
})
genes=rownames(ct)[which(pct>=0.15)]

p_clu5=plot_genes_branched_heatmap(dat[genes,],
                branch_point = 2,
                num_clusters = 5,
                cores = 1,
                return_heatmap=T,
                show_rownames = F)
pdf('Treg_pct0.15_Branch2_clu5_pheatmap.pdf')
p_clu5$ph_res
dev.off()

p_clu4=plot_genes_branched_heatmap(dat[genes,],
                branch_point = 2,
                num_clusters = 4,
                cores = 1,
                return_heatmap=T,
                show_rownames = F)
pdf('Treg_pct0.15_Branch2_clu4_pheatmap.pdf')
p_clu4$ph_res
dev.off()

dt.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) ->top100
p_clu4=plot_genes_branched_heatmap(dat[unique(c(unique(top100$gene),VariableFeatures(Treg))),],
                branch_point = 2,
                num_clusters = 4,
                cores = 1,
                return_heatmap=T,
                show_rownames = F)
pdf('Treg_Branch2_clu4_pheatmap.pdf')
p_clu4$ph_res
dev.off()

genes=unique(c(unique(top100$gene),VariableFeatures(Treg)))

map_df=data.frame(
  Cluster=as.numeric(as.vector(sort(cutree(p_clu4$ph_res$tree_row, k=4)))),
  Gene=names(sort(cutree(p_clu4$ph_res$tree_row, k=4)))
)
map_df[map_df$Cluster==2,]

source("~/Biotools/Enrich_Radarchart.R")
for(i in 1:4){
  genes=map_df$Gene[map_df$Cluster==i]
  assign(paste0('df',i), OneGroup_Radarchart(genelt=genes,enrichFunc=c('enrich'),db="human",DataBase=c('KEGG'),return_enrich_rs=T,
            col=cols[i],savedir="./",savename=paste0("C",i,"_enrichKEGG_Radarchart.pdf")))
}

for(i in 1:4){
  genes=map_df$Gene[map_df$Cluster==i]
  assign(paste0('df',i), OneGroup_Radarchart(genelt=genes,enrichFunc=c('enrich'),db="human",DataBase=c('GO'),return_enrich_rs=T,
            col=cols[i],savedir="./",savename=paste0("C",i,"_enrichGO_Radarchart.pdf")))
}

for(i in 1:4){
  genes=map_df$Gene[map_df$Cluster==i]
  assign(paste0('df',i),OneGroup_Radarchart(genelt=genes,enrichFunc=c('enrichr'),db="human",DataBase=c('GMT'),GMT.dir='~/Biotools/GSEA/h.all.v7.4.symbols.gmt',return_enrich_rs=T,
            col=cols[i],savedir="./",savename=paste0("C",i,"_enrichHall_Radarchart.pdf")))
}

order=c(3,4,1,2)

plot_df3=data.frame(
  desc=df3$Description,
  log10Pvalue= -log10(df3$pvalue)
)

plot_df4=data.frame(
  desc=df4$Description[c(2,5,7,8,11)],
  log10Pvalue= -log10(df4$pvalue[c(2,5,7,8,11)])
)

plot_df1=data.frame(
  desc=df1$Description[c(1,2,3,4,8,9)],
  log10Pvalue= -log10(df1$pvalue[c(1,2,3,4,8,9)])
)

plot_df2=data.frame(
  desc=df2$Description[c(1,3,4,9,12)],
  log10Pvalue= -log10(df2$pvalue[c(1,3,4,9,12)])
)

plot_df=rbind(plot_df3,plot_df4,plot_df1,plot_df2)
plot_df$group=rep(rep(paste0('Clu',order)),c(2,5,6,5))

plot_lists=lapply(paste0('Clu',order),function(x){
  tmp=plot_df[plot_df$group==x,]
  tmp$desc=factor(tmp$desc,levels=rev(tmp$desc))
  ggplot(tmp,aes(x=log10Pvalue,y=desc))+
    geom_bar(fill=cols[as.numeric(substr(x,4,4))],stat="identity",width=0.4)+
    theme_bw()+labs(x='-log10PValue',y='',title=x)
})
pdf('./Monocle_clu5_enrichKEGG_barplot.pdf',width=6,height=12)
cowplot::plot_grid(plotlist=plot_lists,ncol=1)
dev.off()

grep('^CD',map_df$Gene[map_df$Cluster==1],value=T)

grep('^IL',map_df$Gene[map_df$Cluster==1],value=T)

grep('^CXCR|^CCL|^CCR',map_df$Gene[map_df$Cluster==1],value=T)

dat2=dat
dat2@assayData$exprs=Treg@assays$RNA@data

pdf('Monocle_C1_CD_Genes_branched_pseudotime.pdf',width=5,height=20)
plot_genes_branched_pseudotime(dat2[grep('^CD',map_df$Gene[map_df$Cluster==1],value=T),],
                       branch_point = 2,
                       color_by = "Cluster",cell_size=0.5,
                       ncol = 1)+scale_color_manual(values=clucol)
dev.off()

pdf('Monocle_C1_IL_Genes_branched_pseudotime.pdf',width=5,height=20)
plot_genes_branched_pseudotime(dat2[grep('^IL',map_df$Gene[map_df$Cluster==1],value=T),],
                       branch_point = 2,
                       color_by = "Cluster",cell_size=0.5,
                       ncol = 1)+scale_color_manual(values=clucol)
dev.off()

pdf('Monocle_C1_Cytokines_Genes_branched_pseudotime.pdf',width=5,height=20)
plot_genes_branched_pseudotime(dat2[grep('^CXCR|^CCL|^CCR',map_df$Gene[map_df$Cluster==1],value=T),],
                       branch_point = 2,
                       color_by = "Cluster",cell_size=0.5,
                       ncol = 1)+scale_color_manual(values=clucol)
dev.off()

Inhibitory=c('ENTPD1','TIGIT','IL10','LAG3','NT5E','HAVCR2','TGFB1','CTLA4','GZMA','TNFSF14','GZMM')
Resident<-c('ITGAE','ITGA1','CXCR6',"CX3CR1",'CRTAM','ZNF683')
Naive=c('TCF7','CCR7','SELL','LEF1')
Activited<-c('IL7R','rna_IL2RA','CD69')
Proliferation=c('STMN1','TOP2A','MKI67','CCNB1','TYMS','HMGB2','PCNA','DEK')

pdf('Monocle_C1_Function_Genes_branched_pseudotime.pdf',width=5,height=20)
plot_genes_branched_pseudotime(dat2[intersect(map_df$Gene[map_df$Cluster==1],c(Resident,Inhibitory)),],
                       branch_point = 2,
                       color_by = "Cluster",cell_size=0.5,
                       ncol = 1)+scale_color_manual(values=clucol)
dev.off()

C1_selected_genes=c('CD177','LAYN','CCR8','CD274','CD83','IL1R1','IL1R2','CXCR6','CX3CR1','CCL20','CCL22','CTLA4','TIGIT','HAVCR2')
pdf('Monocle_C1_Selected_Genes_branched_pseudotime.pdf',width=5,height=1.5*length(C1_selected_genes))
plot_genes_branched_pseudotime(dat2[C1_selected_genes,],
                       branch_point = 2,
                       color_by = "Cluster",cell_size=0.5,
                       ncol = 1)+scale_color_manual(values=clucol)
dev.off()

#C2
pdf('Monocle_C2_CD_Genes_branched_pseudotime.pdf',width=5,height=20)
plot_genes_branched_pseudotime(dat2[grep('^CD',map_df$Gene[map_df$Cluster==2],value=T),],
                       branch_point = 2,
                       color_by = "Cluster",cell_size=0.5,
                       ncol = 1)+scale_color_manual(values=clucol)
dev.off()

pdf('Monocle_C2_IL_Genes_branched_pseudotime.pdf',width=5,height=20)
plot_genes_branched_pseudotime(dat2[grep('^IL',map_df$Gene[map_df$Cluster==2],value=T),],
                       branch_point = 2,
                       color_by = "Cluster",cell_size=0.5,
                       ncol = 1)+scale_color_manual(values=clucol)
dev.off()

pdf('Monocle_C2_Cytokines_Genes_branched_pseudotime.pdf',width=5,height=20)
plot_genes_branched_pseudotime(dat2[grep('^CXCR|^CCL|^CCR',map_df$Gene[map_df$Cluster==2],value=T),],
                       branch_point = 2,
                       color_by = "Cluster",cell_size=0.5,
                       ncol = 1)+scale_color_manual(values=clucol)
dev.off()

pdf('Monocle_C2_Function_Genes_branched_pseudotime.pdf',width=5,height=20)
plot_genes_branched_pseudotime(dat2[intersect(map_df$Gene[map_df$Cluster==2],c(Naive,Proliferation)),],
                       branch_point = 2,
                       color_by = "Cluster",cell_size=0.5,
                       ncol = 1)+scale_color_manual(values=clucol)
dev.off()

C2_selected_genes=c('CD38','CDK1','CCL3','CCL4','CCL5','CCNB1','MKI67','PCNA','TOP2A')
pdf('Monocle_C2_Selected_Genes_branched_pseudotime.pdf',width=5,height=1.5*length(C2_selected_genes))
plot_genes_branched_pseudotime(dat2[C2_selected_genes,],
                       branch_point = 2,
                       color_by = "Cluster",cell_size=0.5,
                       ncol = 1)+scale_color_manual(values=clucol)
dev.off()

#scMetabolism
gene.lt=qusage::read.gmt('~/R/library/4.2.0/scMetabolism/data/KEGG_metabolism_nc.gmt')
gene.lt=gene.lt[-which(names(gene.lt)=='D-Arginine and D-ornithine metabolism')]
Treg=pochi::RunAUCell(Treg,assay="RNA",slot="counts",genesets=gene.lt, verbose=T)

meta_rs=Treg@assays$AUC@data

Treg$type<-paste0(Treg$Cluster,"_",Treg$loc)
order_type<-c()
for(i in paste0('C',0:7)){
    order_type<-c(order_type,paste0(i,"_",c("P","N",'T')))
}
mat<-Treg@assays$AUC@data
a<-sapply(rownames(meta_rs),function(x){
    sapply(order_type,function(y){
        mean(as.numeric(as.vector(mat[x,which(Treg$type==y)])))
    })
})
t(a)->a

annotation_col = data.frame(
  Sample = factor(rep(c("P","N",'T'), 8),levels=c("P","N",'T')), 
  CellType=factor(rep(paste0('C',0:7),rep(3,8)),levels=paste0('C',0:7))
)
rownames(annotation_col) = colnames(a)

samplecol=ggsci::pal_futurama("planetexpress")(3)[c(3,1,2)]
names(samplecol)=levels(factor(Treg$loc))

ann_colors = list(
  Sample = samplecol,
  CellType=clucol
)

library(pheatmap)
pdf("Treg_Cluster_sampletype_metabolism_pheatmap.pdf",width=15,height=20)
pheatmap(a,scale="row",cellwidth=10,cellheight=10,border_color=NA,annotation_col=annotation_col,annotation_colors = ann_colors,angle_col=315,cluster_rows=T,cluster_cols=F,show_rownames=T,show_colnames=T,color=rev(heatmaply::RdBu(100)))
dev.off()

a2<-sapply(rownames(meta_rs),function(x){
    sapply(paste0('C',0:7),function(y){
        mean(as.numeric(as.vector(mat[x,which(Treg$Cluster==y)])))
    })
})
t(a2)->a2
tmprs=apply(a2,1,function(x){
    names(x)[which.max(x)]
})

meta_genes=intersect(rownames(dat2),gene.lt[['Riboflavin metabolism']])
pdf('Monocle_VB2_metabolism_Genes_branched_pseudotime.pdf',width=5,height=1.5*length(meta_genes))
plot_genes_branched_pseudotime(dat2[meta_genes,],
                       branch_point = 2,
                       color_by = "Cluster",cell_size=0.5,
                       ncol = 1)+scale_color_manual(values=clucol)
dev.off()

pdf('Monocle_HIF1A_branched_pseudotime.pdf',width=5,height=1.5*2)
plot_genes_branched_pseudotime(dat2[c('HIF1A','EPAS1'),],
                       branch_point = 2,
                       color_by = "Cluster",cell_size=0.5,
                       ncol = 1)+scale_color_manual(values=clucol)
dev.off()

apoptosis_genes=intersect(map_df$Gene[map_df$Cluster==1],KEGG$symbol[KEGG$term=='Apoptosis'])
pdf('Monocle_apoptosis_genes_branched_pseudotime.pdf',width=5,height=1.5*length(apoptosis_genes))
plot_genes_branched_pseudotime(dat2[apoptosis_genes,],
                       branch_point = 2,
                       color_by = "Cluster",cell_size=0.5,
                       ncol = 1)+scale_color_manual(values=clucol)
dev.off()

Treg$Pseudotime=dat$Pseudotime
Treg$VB2=Treg@assays$AUC@data['Riboflavin metabolism',]
Treg$State=dat$State

df2=data.frame(
  pseudotime=c(as.numeric(as.vector(Treg$Pseudotime[Treg$State %in% c(1,5,2,3)])),as.numeric(as.vector(Treg$Pseudotime[Treg$State %in% c(1,5,2,4)]))),
  VB2=c(as.numeric(as.vector(Treg$VB2[Treg$State %in% c(1,5,2,3)])),as.numeric(as.vector(Treg$VB2[Treg$State %in% c(1,5,2,4)]))),
  group=rep(c('line1','line2'),c(sum(table(Treg$State)[c(1,5,2,3)]),sum(table(Treg$State)[c(1,5,2,4)])))
)
df=data.frame(
  pseudotime=as.numeric(as.vector(Treg$Pseudotime)),
  VB2=as.numeric(as.vector(Treg$VB2)),
  celltype=Treg$Cluster
)
p=ggplot()+
  ggrastr::rasterise(geom_point(data=df,aes(x=pseudotime,y=VB2,color=celltype),size=1),dpi=300)+
  geom_smooth(data=df2,aes(x=pseudotime,y=VB2,group=group,linetype=group),method='loess',color='black',size=2)+
  scale_color_manual(values=clucol)+
  xlab("Pseudotime")+ylab("VB2")+
  theme_bw()+theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0.5),panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
ggsave("Treg_Pseudotime_VB2.pdf",p,width=7,height=3.5)

gene.lt2=list(
    'Apoptosis'=KEGG$symbol[KEGG$term=='Apoptosis']
)
Treg=pochi::RunAUCell(Treg,assay="RNA",slot="counts",auc_assay_name='KEGG',genesets=gene.lt2, verbose=T)

Treg$Apoptosis=Treg@assays$KEGG@data['Apoptosis',]
df2=data.frame(
  pseudotime=c(as.numeric(as.vector(Treg$Pseudotime[Treg$State %in% c(1,5,2,3)])),as.numeric(as.vector(Treg$Pseudotime[Treg$State %in% c(1,5,2,4)]))),
  Apoptosis=c(as.numeric(as.vector(Treg$Apoptosis[Treg$State %in% c(1,5,2,3)])),as.numeric(as.vector(Treg$Apoptosis[Treg$State %in% c(1,5,2,4)]))),
  group=rep(c('line1','line2'),c(sum(table(Treg$State)[c(1,5,2,3)]),sum(table(Treg$State)[c(1,5,2,4)])))
)
df=data.frame(
  pseudotime=as.numeric(as.vector(Treg$Pseudotime)),
  Apoptosis=as.numeric(as.vector(Treg$Apoptosis)),
  celltype=Treg$Cluster
)
p=ggplot()+
  ggrastr::rasterise(geom_point(data=df,aes(x=pseudotime,y=Apoptosis,color=celltype),size=1),dpi=300)+
  geom_smooth(data=df2,aes(x=pseudotime,y=Apoptosis,group=group,linetype=group),method='loess',color='black',size=2)+
  scale_color_manual(values=clucol)+
  xlab("Pseudotime")+ylab("Apoptosis")+
  theme_bw()+theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0.5),panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
ggsave("Treg_Pseudotime_Apoptosis.pdf",p,width=7,height=3.5)

#Scenic
dir.create("./scenic/")
setwd("scenic/")

badterm <- c(
  "^RPS","^RPL","^MT-","^MTMR","^MTND","^MTRN","^MTCO","^MRPL","^MRPS","^HBA","^HBB","^MTATP",
  "^MALAT1$", "^XIST$", "^XIST_intron$",
  "^AC[0-9]+\\.", "^AP[0-9]+\\.", "^AL[0-9]+\\.", "-AS[0-9]*$", "^MIR[0-9]", "^LINC[0-9]", "^SNOR", "^sno", "^SNR",
  "^TMSB",  "-", ":", "\\.", '^KIAA',
  "^DNAJ"
)
badgenelt<-unique(grep(pattern=paste(c(badterm),collapse="|"),x=rownames(Treg),perl=T,value=T))

ct=Treg@assays$RNA@counts[-which(rownames(Treg) %in% badgenelt),]

mat<-t(as.matrix(Treg@assays$RNA@data[rownames(ct),]))
write.table(mat,"./Treg.csv",quote=F,row.names=T,col.names=T,sep=",")

export PATH=/public/workspace/liangyuan/.conda/envs/pyscenic/lib/python3.10/site-packages/bin:$PATH

pyscenic grn -o ./expr_mat.adjacencies.tsv -m grnboost2 --seed 12345 --num_workers 20 \
  ./Treg.csv /public/workspace/liangyuan/Biotools/SCENIC/hs_hgnc_curated_tfs.txt

pyscenic ctx ./expr_mat.adjacencies.tsv \
 /public/workspace/liangyuan/Biotools/SCENIC/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
 /public/workspace/liangyuan/Biotools/SCENIC/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
 --annotations_fname /public/workspace/liangyuan/Biotools/SCENIC/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
 --expression_mtx_fname ./Treg.csv \
 --mode "custom_multiprocessing" \
 --output ./regulons.csv \
 --num_workers 20

pyscenic aucell \
 ./Treg.csv \
 ./regulons.csv \
 -o ./auc_mtx.csv \
 --num_workers 20

read.csv("./auc_mtx.csv",header=T,row.names=1,sep=",")->auc
colnames(auc)<-gsub("\\.","",colnames(auc))
t(auc)->auc
as.matrix(auc)->auc

sapply(levels(factor(Treg$Cluster)),function(x){
  tmp<-auc[,which(colnames(auc) %in% colnames(Treg)[which(Treg$Cluster==x)])]
  apply(tmp,1,mean)->mean
  mean
})->a
library(pheatmap)
pdf("./Treg_Cluster_TF_pheatmap.pdf",height=50,width=8)
pheatmap(a,scale="row",cellheight = 10,cellwidth=10,angle_col=45,show_rownames=T,show_colnames=T,cluster_rows = T,cluster_cols=T,color = rev(heatmaply::RdBu(100)))
dev.off()

#
ctMat <- lapply(names(table(Treg$Cluster)), function(i) {
  as.numeric(Treg$Cluster == i)
})
ctMat <- do.call(cbind, ctMat)
colnames(ctMat) <- names(table(Treg$Cluster))
rownames(ctMat) <- colnames(caf)
library(pbapply)
library(philentropy)
auc<-t(auc) #tfs in column
rssMat <- pblapply(colnames(auc), function(i) {
  sapply(colnames(ctMat), function(j) {
    1 - JSD(rbind(auc[, i], ctMat[, j]), unit = 'log2', est.prob = "empirical")
  })
})

rssMat <- do.call(rbind, rssMat)
rownames(rssMat) <- colnames(auc)
colnames(rssMat) <- colnames(ctMat)

source("~/Biotools/pysenic/plotRegulonRank.R")
num=20
regulonsplot=lapply(paste0('C',0:7),function(x){
  coll=clucol[x]
  PlotRegulonRank(rssMat,pt.size=2, cell.type=x,topn=num,color=coll)
})
pdf("Treg_Cluster_top_TFs_regulons.pdf",width=20,height=10)
cowplot::plot_grid(plotlist=regulonsplot, nrow = 2)
dev.off()

rownames(rssMat[order(rssMat[,3],decreasing=T),])[1:20]

DefaultAssay(Treg)='RNA'
colours = colorRampPalette(c("grey89", "#FEF8E5","#FFEABC","#FDBC52","#F68523","#DD3226","#A31B3A","#5E2C81","#382F85","#28316C"))(1000)
p=FeaturePlot(Treg,raster=T,order=T,features='REL',col=c("lightgray","red"),pt.size=5,reduction="umap",ncol=1,raster.dpi=c(2048,2048))+scale_color_gradientn(colors=colours)
ggsave('Treg_REL_RNAExp.pdf',p,width=7,height=6)
VLN=VlnPlot(Treg,pt.size=0,features =c('REL'),group.by="Cluster",sort=T,cols=clucol)+
    xlab("")+ylab("")
ggsave("./Treg_REL_StackVlnplot.pdf",VLN,width=10,height=8)

Treg$REL_auc=as.numeric(as.vector(auc[,'REL']))
colours = colorRampPalette(c("grey89",'lightgrey', "#9FC6DF","#75CAFD","#6F75FF","#3C00FF","#2A00FF"))(50)
p=FeaturePlot(Treg,raster=T,order=T,features='REL_auc',col=c("lightgray","red"),pt.size=5,reduction="umap",ncol=1,raster.dpi=c(2048,2048))+
  scale_color_gradientn(colors=colours)
ggsave('./Treg_REL_AUC_fea.pdf',p,width=7,height=6)
VLN=VlnPlot(Treg,pt.size=0,features =c('REL_auc'),group.by="Cluster",sort=T,cols=clucol)+
    xlab("")+ylab("")
ggsave("./Treg_REL_auc_StackVlnplot.pdf",VLN,width=10,height=8)

#Survival analysis
#ICGC
exp=read.table('~/HCCDB18_mRNA_level3.txt',header=T,sep="\t")
exp=exp[,-1] %>% as.data.frame() %>% group_by(Symbol) %>% summarise_each(funs(mean))
rowns=exp$Symbol
exp=as.matrix(exp[,-1])
rownames(exp)=rowns
exp=exp[-1,]
saveRDS(exp,"ICGC_JP_expr.RDS")

exp=readRDS('~/Subject/HCC_CAF_scRNA/Harmony/Survival/ICGC-JP/ICGC_JP_HCC_exp.RDS')

clic=as.data.frame(t(read.table('~/Subject/HCC_CAF_scRNA/Harmony/Survival/ICGC-JP/HCCDB18.sample.txt',header=T,row.names=1,sep="\t")))
tumor_samples=rownames(clic)[which(clic$TYPE=="HCC")]

info=as.data.frame(t(read.table('~/Subject/HCC_CAF_scRNA/Harmony/Survival/ICGC-JP/HCCDB18.patient.txt',header=T,row.names=1,sep="\t")))
rownames(info)=gsub('P','S',rownames(info))

inte_samples=intersect(tumor_samples,rownames(info))
exp=exp[,inte_samples]
info=info[inte_samples,]

future::plan("multicore", workers = 5)
library(dplyr)
Treg@active.ident<-factor(Treg$Cluster)
dt.markers <- FindAllMarkers(Treg, features=rownames(Treg)[-which(rownames(Treg) %in% badgenelt)],only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dt.markers$pct_diff<-dt.markers$pct.1-dt.markers$pct.2
dt.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) ->top

genelt=lapply(levels(top$cluster),function(x){
  top$gene[top$cluster==x]
})
names(genelt)=levels(top$cluster)

ssgsea_rs=GSVA::gsva(as.matrix(exp),gset.idx.list=genelt,method="ssgsea",min.sz=3, max.sz=100, verbose=T,parallel.sz=1)

library(survminer)
library(survival)

info$state=ifelse(info$STATUS=="Alive",0,1)
info$time=as.numeric(as.vector(info$SUR))

i="C0"
info$score<-as.numeric(as.vector(ssgsea_rs[i,]))
median_score=median(info$score)
info$type=ifelse(info$score>=median_score,"high","low")
pdf(paste0("ICGC_LIHC_",i,"_ssgsea_surv.pdf"),onefile=F)
fit<-survfit(Surv(time,state)~type,data=info)
ggsurvplot(fit,title=i,data=info,conf.int=F,pval=T,risk.table=T,palette=c("#E7B800","#2E9FDF"),legend="bottom",legend.title="Score",legend.labs=c("high","low"))
dev.off()

res.cut <- surv_cutpoint(info, time = "time", event = "state",variables = c("score"))
res.cat <- surv_categorize(res.cut)

pdf(paste0("ICGC_LIHC_",i,"_best_surv.pdf"),onefile=F)
fit<-survfit(Surv(time,state)~score,data=res.cat)
ggsurvplot(fit,title=i,data=res.cat,conf.int=F,pval=T,risk.table=T,palette=c("#E7B800","#2E9FDF"),legend="bottom",legend.title="Score",legend.labs=c("high","low"))
dev.off()

#TCGA
exp=readRDS('~/Subject/HCC_CAF_scRNA/Harmony/Bulk_T_vs_N/TCGA/TCGA_HCC_Expr.RDS')
colnames(exp)=gsub('-','\\.',colnames(exp))

clic=as.data.frame(t(read.table('~/Subject/HCC_CAF_scRNA/Harmony/Survival/TCGA/HCCDB15.sample.txt',header=T,row.names=1,sep="\t")))
tumor_samples=rownames(clic)[which(clic$TYPE=="HCC")]

info=as.data.frame(t(read.table('~/Subject/HCC_CAF_scRNA/Harmony/Survival/TCGA/HCCDB15.patient.txt',header=T,row.names=1,sep="\t")))
rownames(info)=gsub('P','S',rownames(info))

inte_samples=intersect(tumor_samples,rownames(info))
exp=exp[,inte_samples]
info=info[inte_samples,]

ssgsea_rs=GSVA::gsva(as.matrix(exp),gset.idx.list=genelt,method="ssgsea",min.sz=3, max.sz=100, verbose=T,parallel.sz=1)

library(survminer)
library(survival)

info$state=ifelse(info$STATUS=="Alive",0,1)
info$time=as.numeric(as.vector(round(as.numeric(info$SURVIVAL_TIME)/30,2)))

i="C0"
info$score<-as.numeric(as.vector(ssgsea_rs[i,]))
median_score=median(info$score)
info$type=ifelse(info$score>=median_score,"high","low")
pdf(paste0("TCGA_LIHC_",i,"_ssgsea_surv.pdf"),onefile=F)
fit<-survfit(Surv(time,state)~type,data=info)
ggsurvplot(fit,title=i,data=info,conf.int=F,pval=T,risk.table=T,palette=c("#E7B800","#2E9FDF"),legend="bottom",legend.title="Score",legend.labs=c("high","low"))
dev.off()

res.cut <- surv_cutpoint(info, time = "time", event = "state",variables = c("score"))
res.cat <- surv_categorize(res.cut)

pdf(paste0("TCGA_LIHC_",i,"_best_surv.pdf"),onefile=F)
fit<-survfit(Surv(time,state)~score,data=res.cat)
ggsurvplot(fit,title=i,data=res.cat,conf.int=F,pval=T,risk.table=T,palette=c("#E7B800","#2E9FDF"),legend="bottom",legend.title="Score",legend.labs=c("high","low"))
dev.off()

#FU-HCC
tmp=readRDS('~/Subject/HCC_CAF_scRNA/Harmony/Survival/FU-HCC/Fanjia.RDS')

exp=tmp[['expr']]
colnames(exp)=paste0("P",colnames(exp))

info=tmp[['survival']]
rownames(info)=paste0("P",info$ID)

dt.markers %>% group_by(cluster) %>% top_n(n = 12, wt = avg_log2FC) ->top

genelt=lapply(levels(top$cluster),function(x){
  top$gene[top$cluster==x]
})
names(genelt)=levels(top$cluster)

ssgsea_rs=GSVA::gsva(as.matrix(exp),gset.idx.list=genelt,method="ssgsea",min.sz=3, max.sz=100, verbose=T,parallel.sz=1)

library(survminer)
library(survival)

info$state=as.numeric(as.vector(info$OS))
info$time=as.numeric(as.vector(info$OS.time))

i="C3"
info$score<-as.numeric(as.vector(ssgsea_rs[i,]))
median_score=median(info$score)
info$type=ifelse(info$score>=median_score,"high","low")
pdf(paste0("FU-HCC_",i,"_ssgsea_surv.pdf"),onefile=F)
fit<-survfit(Surv(time,state)~type,data=info)
ggsurvplot(fit,title=i,data=info,conf.int=F,pval=T,risk.table=T,palette=c("#E7B800","#2E9FDF"),legend="bottom",legend.title="Score",legend.labs=c("high","low"))
dev.off()

res.cut <- surv_cutpoint(info, time = "time", event = "state",variables = c("score"))
res.cat <- surv_categorize(res.cut)

pdf(paste0("FU-HCC_",i,"_best_surv.pdf"),onefile=F)
fit<-survfit(Surv(time,state)~score,data=res.cat)
ggsurvplot(fit,title=i,data=res.cat,conf.int=F,pval=T,risk.table=T,palette=c("#E7B800","#2E9FDF"),legend="bottom",legend.title="Score",legend.labs=c("high","low"))
dev.off()

#HCC PD1 耐药
exp=read.table('~/Subject/CQY_PSME3/data/FU-PD/FUPD_expr.txt',header=T,row.names=1,sep='\t')

group=readxl::read_xlsx('~/Subject/CQY_PSME3/data/AntiPD1_Bulk/count.xlsx',sheet=2)
groups=group$Effect

dt.markers %>% group_by(cluster) %>% top_n(n = 12, wt = avg_log2FC) ->top

genelt=lapply(levels(top$cluster),function(x){
  top$gene[top$cluster==x]
})
names(genelt)=levels(top$cluster)

ssgsea_rs=GSVA::gsva(as.matrix(exp),gset.idx.list=genelt,method="ssgsea",min.sz=3, max.sz=100, verbose=T,parallel.sz=1)

plot_list=lapply(rownames(ssgsea_rs),function(x){
    df=data.frame(
        value=as.numeric(ssgsea_rs[x,]),
        group=groups
    )
    pvalue=round(wilcox.test(value~group,df,alternative="greater")$p.value,8)
    ggplot(df,aes(x=group,y=value,color=group))+ 
        stat_boxplot(geom="errorbar", width = 0.4)+
        geom_boxplot(width=0.4,outlier.shape=NA)+
        geom_point(aes(color=group),position = position_jitter(width = 0.3, height = 0),size=2)+
        #ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test",comparisons = list(c("PD","PR"),alternative="greater"),size=8)+
        scale_color_manual(values=ggsci::pal_futurama("planetexpress")(5)[c(2,3)])+
        theme_bw()+theme(text = element_text(size = 8),axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0),panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"), legend.position="bottom") +
        ggtitle(paste0(x," P=",pvalue))+labs(x="", y = "ssGSEA score")
})
pdf('FUPD_ssgsea_PD_vs_PR.pdf',width=15,height=12)
cowplot::plot_grid(plotlist=plot_list,nrow=3)
dev.off()

