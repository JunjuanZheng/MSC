library(dplyr)  
library(ggplot2)  
library(patchwork)  
library(Seurat)  
library(clustree)  # 用于评估聚类稳定性  
library(DoubletFinder)
library(clustree)  
library(igraph)  
library(cluster)  
library(bluster)
library(SeuratWrappers)
# 设置工作目录  
projectPath <- "/mnt2/wanggd_group/zjj/BGCscRNA/MSC_Duan"  
setwd(projectPath)  

##### 1.数据导入并创建seurat对象 #####  

# 获取数据文件夹下的所有样本文件列表  
samples <- list.files(file.path(projectPath, "Data/QuantifyRawdata/"))  
seurat_list <- list()  
sample_ids <- c()  # 改为向量而不是列表  

# 读取每个样本的10x数据并创建Seurat对象  
for (sample in samples) {  
  tryCatch({  
    # 拼接文件路径  
    data.path <- file.path(projectPath, "Data/QuantifyRawdata", sample,   
                         paste0(sample, "_outs"), "filtered_cell_gene_matrix")  
    
    # 检查路径是否存在  
    if (!dir.exists(data.path)) {  
      warning(paste("Directory not found:", data.path))  
      next  
    }  
    # 读取数据  
    seurat_data<-Read10X(data.dir = data.path)  
    # 提取样本ID  
    sample_id <- sample  
    
    # 创建Seurat对象  
    seurat_obj <- CreateSeuratObject(counts = seurat_data,  
                                    project = sample_id,  
                                    min.features = 0,  
                                    min.cells = 0) 
								
	seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

    vln_plot_nCount <- VlnPlot(seurat_obj, features = "nCount_RNA", pt.size = 0) + 
    ggtitle("nCount")

    vln_plot_nFeature <- VlnPlot(seurat_obj, features = "nFeature_RNA", pt.size = 0) + 
    ggtitle("nFeature")

    vln_plot_percent_mt <- VlnPlot(seurat_obj, features = "percent.mt", pt.size = 0) + 
    ggtitle("percent.mt")

    # 将小提琴图组合起来显示
    combined_vln_plot <- vln_plot_nCount | vln_plot_nFeature | vln_plot_percent_mt
	# 保存小提琴图到文件
    output_file <- file.path(projectPath, "Output", paste0(sample,".nFeature_nCount_percent_mt.ViolinPlots.png"))
    ggsave(output_file, plot = combined_vln_plot, width = 12, height = 12, dpi = 300)
    # 提示保存成功
    message(paste("Violin plots saved to:", output_file))
	
    # 将Seurat对象添加到列表中  
    seurat_list[[sample_id]] <- seurat_obj  
    sample_ids <- c(sample_ids, sample_id)  
    
  }, error = function(e) {  
    message(paste("Error processing sample:", sample))  
    message(e)  
  })  
}  

# 检查是否成功创建了Seurat对象  
if (length(seurat_list) == 0) {  
  stop("No Seurat objects were created successfully")  
}    

#根据VlunPlot subset object	
for (sample in samples) {
    seurat_obj = seurat_list[[sample]]  
    seurat_obj <- subset(seurat_obj, subset = nCount_RNA < 50000 & 
                                              nFeature_RNA > 200 &
                                              percent.mt < 25)
	seurat_list[[sample]] = seurat_obj										  
}
# 保存初始对象  
save(seurat_list,   
    file = file.path(projectPath, "Data", "20250522_3SamplesMSC_ObjectOriginal.Rdata"))

### 过滤双细胞
############################
load(file.path(projectPath, "Data", "20250522_3SamplesMSC_ObjectOriginal.Rdata"))
# 创建双细胞检测结果目录
doublet_dir <- file.path(projectPath, "Output/DoubletFinder")
dir.create(doublet_dir, showWarnings = FALSE, recursive = TRUE)
pdf(file.path(doublet_dir, paste0("20250526_3SamplesMSC_ObjectDoublets.pdf")), width = 10, height = 8)

cc_genes <- cc.genes.updated.2019 # 使用 Seurat 包中的细胞周期基因	

# 处理每个样本 
for (i in seq_along(seurat_list)) {  
  obj_name <- names(seurat_list)[i]  
  current_obj <- seurat_list[[obj_name]]  
  
  message(sprintf("Processing %d/%d: %s", i, length(seurat_list), obj_name))  
  
  # 基础预处理  
# 基础预处理（到PCA阶段）
  current_obj <- current_obj %>%  
    NormalizeData() %>%  
    FindVariableFeatures() %>%  
    ScaleData() %>%  
    RunPCA(verbose = FALSE)
  
  # 自动参数优化
  sweep.res <- paramSweep(current_obj, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # 选择最佳pK值
  pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  # 计算预期双细胞数量（按细胞数的2.5%估算，可根据实验调整）
  nExp_poi <- round(0.025 * ncol(current_obj)) 
  
  # 运行DoubletFinder
  current_obj <- doubletFinder(
    current_obj, 
    PCs = 1:10, 
    pN = 0.25, 
    pK = pK, 
    nExp = nExp_poi, 
    reuse.pANN = FALSE, 
    sct = FALSE
  )
  
  # 提取双细胞分类结果
  df.class <- paste("DF.classifications_0.25", pK, nExp_poi, sep = "_")
  current_obj$Doublet_Classification <- current_obj@meta.data[[df.class]]
  
  # 添加双胞率计算和标注  
  doublet_count <- sum(current_obj$Doublet_Classification == "Doublet")  
  total_cells <- ncol(current_obj)  
  doublet_rate <- round(doublet_count / total_cells * 100, 2)  
  # 生成带统计信息的绘图  
  doublet_plot <- DimPlot(current_obj,   
                         group.by = "Doublet_Classification",  
                         reduction = "pca",  
                         cols = c("gray90", "red")) +  # 设置颜色  
    labs(title = paste0(obj_name, " Doublet Detection")) +  
    annotate("text",  
             x = Inf, y = Inf,  # 右上角定位  
             label = paste0("Doublet Rate: ", doublet_rate, "%\n",  
                          "Predicted: ", nExp_poi, " (", doublet_count, " found)"),  
             hjust = 1.1, vjust = 1.1,  # 微调位置  
             size = 5,  
             color = "darkred",  
             fontface = "bold") +  
    theme(  
      plot.title = element_text(size=14, face="bold"),  
      legend.position = c(0.8, 0.2)  # 调整图例位置  
    )  
  print(doublet_plot)
  # 过滤双细胞
  singlet_cells <- colnames(current_obj)[current_obj$Doublet_Classification == "Singlet"]
  current_obj <- subset(current_obj, cells = singlet_cells) 
  current_obj <- CellCycleScoring(current_obj, s.features = cc_genes$s.genes, 
                                g2m.features = cc_genes$g2m.genes, 
                                set.ident = TRUE)
  current_obj <-NormalizeData(current_obj, 
                              normalization.method = "LogNormalize", 
                              scale.factor = 10000) %>%
              ScaleData(vars.to.regress = c("S.Score", "G2M.Score", "percent.mt")) %>%   #能减少因细胞周期状态或其他技术因素导致的误差，从而更清晰地看到生物学上有意义的信号
              FindVariableFeatures() %>%
              RunPCA(verbose = FALSE)
  

  # 更新对象  
  seurat_list[[obj_name]] <- current_obj  
  
  # 保存当前进度  
  save(current_obj,  
       file = file.path(projectPath, "Data",   
                       #paste0("20250522_3SamplesMSC_ObjectAfterDoublets_", i, ".Rdata")))
                        paste0("20250526_3SamplesMSC_ObjectAfterDoublets_ScaleCellCircle_", i, ".Rdata")))					   
	
}
#save(seurat_list,  
#     file = file.path(projectPath, "Data", "20250522_3SamplesMSC_ObjectOriginal_AfterDoublets.Rdata"))  
save(seurat_list,  
     file = file.path(projectPath, "Data", "20250526_3SamplesMSC_ObjectOriginal_AfterDoubletsScaleCellCircle.Rdata"))  

#############################################################################
#################################
######### Integreation ##########
#################################
#integration_methods <- c(
  #"Seurat", "fastMNN", "Harmony",  "BBKNN",  这些方法上面用到过
  #"Uncorrected", 
  #"scVI", "MNN", "Scanorama","CSS", "LIGER", "Conos", "ComBat"
#)

library(future)
options(future.globals.maxSize = 10000 * 1024^2)  # 设置为 10 GB
counts_list <- list(
  AD   = seurat_list[[1]]@assays$RNA@layers$counts,
  BM   = seurat_list[[2]]@assays$RNA@layers$counts,
  UC   = seurat_list[[3]]@assays$RNA@layers$counts
)

unintegrated <- merge(seurat_list[[1]], seurat_list[2:3]) 
unintegrated <- NormalizeData(unintegrated)
unintegrated <- FindVariableFeatures(unintegrated)
unintegrated <- ScaleData(unintegrated)
unintegrated <- RunPCA(unintegrated)
unintegrated <- FindNeighbors(unintegrated, dims = 1:30)
unintegrated <- RunUMAP(unintegrated, dims = 1:30)
pdf(file.path(projectPath, "Output", paste0("20250522_3SamplesMSC_unIntegrationUmap.pdf")), width = 10, height = 8)
DimPlot(unintegrated, reduction = "umap", group.by = "orig.ident")
dev.off()
obj = unintegrated
############## ------------------------ 在R里面可以实现的整合方法
integration_info <- list(  
  list(reduction = "integrated.cca", cluster.name = "cca_clusters",method.use ='CCAIntegration',umap.reduction = 'umap.cca'),  
  list(reduction = "integrated.rpca", cluster.name = "rpca_clusters",method.use ='RPCAIntegration',umap.reduction = 'umap.rpca'),  
  list(reduction = "harmony", cluster.name = "harmony_clusters",method.use ='HarmonyIntegration',umap.reduction = 'umap.harmony'),  
  list(reduction = "integrated.jointPCA", cluster.name = "jointPCA_clusters",method.use = 'JointPCAIntegration',umap.reduction = 'umap.jointPCA'),
  list(reduction = "fastMNN", cluster.name = "fastMNN_clusters",method.use ='FastMNNIntegration',umap.reduction = 'umap.fastMNN'),  
  list(reduction = "integrated.scvi", cluster.name = "scVI_clusters",method.use = 'scVIIntegration')    
)

objIntegrated_List = list() 
for (i in (1:length(integration_info))) {
  #i=5
  red.use <- integration_info[[i]]$reduction  
  clst.name <- integration_info[[i]]$cluster.name  
  method.use <- integration_info[[i]]$method.use
  umap.reduction <- integration_info[[i]]$umap.reduction

  objIntegrated<- IntegrateLayers(
		object = obj,   
		method = method.use,  
		orig.reduction = "pca",       # 原始降维的名称  
		new.reduction = red.use,  
		verbose = FALSE  
    )
  # re-join layers after integration
  objIntegrated[["RNA"]] <- JoinLayers(objIntegrated[["RNA"]])

  objIntegrated <- FindNeighbors(objIntegrated, reduction = red.use, dims = 1:30)
  objIntegrated <- FindClusters(objIntegrated, resolution = 0.25,cluster.name = clst.name)
  
  objIntegrated <- RunUMAP(objIntegrated, dims = 1:30, reduction = red.use,reduction.name = umap.reduction)
  
  objIntegrated_List[[method.use]] = objIntegrated
  save(objIntegrated,file = file.path(projectPath, "Data",   
       paste0("20250522_3SamplesMSC_integratedObjDetails_",method.use,"00.Rdata")))
}
save(objIntegrated_List,file = file.path(projectPath, "Data",   
                       paste0("20250522_3SamplesMSC_integratedObjDetails_ListSeuratMethods.Rdata")))

pdf(file.path(projectPath, "Output", paste0("20250522_3SamplesMSC_IntegrationUmap.pdf")), width = 10, height = 8)
for (i in (1:length(integration_info))) {
  #i=5
  red.use <- integration_info[[i]]$reduction  
  clst.name <- integration_info[[i]]$cluster.name  
  method.use <- integration_info[[i]]$method.use
  umap.reduction <- integration_info[[i]]$umap.reduction

  objIntegrated = objIntegrated_List[[method.use]]
  p1 <- DimPlot(  
    object = objIntegrated,
    reduction = umap.reduction,  
    group.by = c("orig.ident"),  
    label.size = 2  
  ) + ggtitle("orig.ident")
  p2 <- DimPlot(  
    object = objIntegrated,
    reduction = umap.reduction,  
    group.by = c(clst.name),  
    label.size = 2  
  ) + ggtitle(red.use)
  combinedUmap = p1|p2
  print(combinedUmap) 
}  
dev.off()

#######################################
load(file = file.path(projectPath, "Data",   
                       paste0("20250522_3SamplesMSC_integratedObjDetails_ListSeuratMethods.Rdata")))
					   
resolutions <- c(0.025, 0.05, 0.1, 0.15, 0.2, 0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6)  
for (ii in names(objIntegrated_List)) {
  red.use <- integration_info[[ii]]$reduction  
  clst.name <- integration_info[[ii]]$cluster.name  
  method.use <- integration_info[[ii]]$method.use
  umap.reduction <- integration_info[[ii]]$umap.reduction

  IntegrationObj = objIntegrated_List[[method.use]]					   

  for (i in resolutions){
    IntegrationObj = FindClusters(  
      object     = IntegrationObj,  
      resolution = i,  
      algorithm  = 1 
      # cluster.name = clst.name,
    )
  }
  objIntegrated_List[[method.use]] = IntegrationObj
  objIntegrated =IntegrationObj
  save(objIntegrated,file = file.path(projectPath, "Data/Integration",paste0("20250617_3SamplesMSC_",method.use,"_FindClusters_objIntegrated.Rdata")))

}
objIntegrated = FindClusters(  
      object     = objIntegrated,  
      resolution = 0.16,  
      algorithm  = 1 
      # cluster.name = clst.name,
    )
IntegrationObj=	objIntegrated

library(clustree)    # 提供 clustree()
library(ggplot2)  
umap_pdf <- file.path(projectPath, "Output/Integration",paste0("IntegrationObj_",method.use,"_AllClusterTreeUMAP_Plots.pdf"))  
pdf(umap_pdf, width = 10, height = 8) 
#for(n in (1:length(seurat_list))) {
  name = method.use #names(seurat_list[n])
  current_obj =IntegrationObj
  clustree_plot <- clustree(current_obj, prefix = "RNA_snn_res.") +  
    ggtitle(paste0(name, " Clustering Tree")) +  
    theme(plot.title = element_text(size = 14, face = "bold"))  
  print(clustree_plot)
  best_resolution = 'RNA_snn_res.0.4'
  Idents(current_obj) = current_obj$RNA_snn_res.0.4
  umap_plot <- DimPlot(  
    current_obj,  
    reduction = umap.reduction,  
    label = TRUE  
  ) +  
    ggtitle(paste0(name, " (Final resolution = ", best_resolution, ")"))  
  print(umap_plot) 

  umap_plot01 <- DimPlot(  
    current_obj,  
    reduction = umap.reduction,
    group.by = 'RNA_snn_res.0.4',	
    label = TRUE  
  ) +  
    ggtitle(paste0(name, " RNA_snn_res.0.4")) 
  print(umap_plot01) 
  
  umap_plot01 <- DimPlot(  
    current_obj,  
    reduction = umap.reduction,
    group.by = 'orig.ident',	
    label = TRUE  
  )
  print(umap_plot01)
  
#}
dev.off()
IntegrationObj = current_obj
objIntegrated =current_obj
save(objIntegrated,file = file.path(projectPath, "Data/Integration",paste0("20250617_3SamplesMSC_",method.use,"_FindClusters_objIntegrated.Rdata")))


##################################################################
#########--------整合后 IntegrationObj 分析-------- ##############
##################################################################
load(file.path(projectPath, "Data/Integration",'20250617_3SamplesMSC_FastMNNIntegration_FindClusters_objIntegrated.Rdata') )
#objIntegrated = objIntegrated_List[['FastMNNIntegration']]
#objIntegrated = JoinLayers(objIntegrated)
#DefaultAssay(objIntegrated)='RNA'
#objIntegrated <- NormalizeData(objIntegrated)
#objIntegrated <- FindVariableFeatures(objIntegrated)
#objIntegrated <- ScaleData(objIntegrated,verbose = FALSE,features = rownames(objIntegrated),vars.to.regress = 'Phase')
load(file.path(projectPath, "Data/Integration20250829",paste0("20250829_3SamplesMSC_FastMNNIntegration_FindClusters_objIntegrated.Rdata")) )



library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)
library(ComplexHeatmap)
library(circlize)
Idents(objIntegrated) = objIntegrated$RNA_snn_res.0.18   #0.4
clusters <- objIntegrated$RNA_snn_res.0.18   #0.4
# Create a new vector with "C" prefixed to the cluster numbers
new_clusters <- paste0("C", clusters)
# Add the new clusters as a new metadata column
objIntegrated$MSC_Clusters <- new_clusters
# 确保是因子
objIntegrated$MSC_Clusters <- factor(objIntegrated$MSC_Clusters)
# 旧水平
old <- levels(objIntegrated$MSC_Clusters)
# 新水平：把数字部分 +1
new <- paste0("C", as.integer(sub("^C", "", old)) + 1)
# 替换水平
levels(objIntegrated$MSC_Clusters) <- new
objIntegrated$celltype=objIntegrated$MSC_Clusters
Idents(objIntegrated) = 'MSC_Clusters'
# Find markers for all clusters
all.markers <- FindAllMarkers(objIntegrated, 
                             only.pos = TRUE, 
                             min.pct = 0.25, 
                             logfc.threshold = 0.25)
write.csv(all.markers,file.path(projectPath, "Data/Integration20250829",'20250829_3SamplesMSC_FastMNNIntegration_FindAllMarkers_Res0.18.csv'),quote=F)

# Select top markers per cluster (e.g., top 5)
# 获取唯一基因列表并清理掉不需要的基因类型
all.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.58) %>%  # log2FC > 0.58相当于FC > 1.5
    slice_head(n = 30) %>%
    ungroup() -> topGs
uedeGenes <- unique(topGs$gene)
# 一步完成所有过滤（更高效）
TopGenes <- uedeGenes[
    !grepl("^MT-", uedeGenes) &        # 过滤线粒体基因
    !grepl("^RPL\\d+", uedeGenes) &    # 过滤核糖体大亚基蛋白基因
    !grepl("^RPS\\d+", uedeGenes) &    # 过滤核糖体小亚基蛋白基因
    !grepl("^ENSCAFG", uedeGenes)      # 过滤Ensembl犬基因ID
]
# 检查结果
length(uedeGenes)   # 过滤前基因数
length(TopGenes)    # 过滤后基因数
#DoHeatmap中默认slot = "scale.data" ，也就是用的数据是ScaleData后的结果数据集。
#而FindAllMarkers中默认slot = "data" ，也就是用的数据是NormalizeData后的数据集。
#因此就可能会造成DoHeatmap画热图时，不在前2000个高变基因中的marker不出现在热图中。
#可以基于top基因重新scale，并且抽样展示，不展示全部的细胞
objIntegrated$MSC_Clusters<-factor(x=objIntegrated$MSC_Clusters,c("C1","C2","C3","C4","C5","C6"))
#objIntegrated$Clusters<-factor(x=objIntegrated$Clusters,c("C0","C1","C2","C3","C4","C5","C6","C7","C8","C9"))
objIntegrated$Phase <- factor(objIntegrated$Phase, levels = c('G1', 'S', 'G2M'))

############### ClusterGVis 
#get top 20genes
library(Seurat)
library(ClusterGVis)
library(org.Cf.eg.db)
all.markers <- FindAllMarkers(objIntegrated, 
                             only.pos = TRUE, 
                             min.pct = 0.25, 
                             logfc.threshold = 0.25)
markers<-all.markers%>%
  dplyr::group_by(cluster)%>%  
  dplyr::filter(avg_log2FC > 0.58 & p_val_adj <= 0.05)%>%  
  dplyr::top_n(n=50,wt=avg_log2FC)
  
TopGenes <- markers$gene[
    !grepl("^MT-", markers$gene) &        # 过滤线粒体基因
    !grepl("^RPL\\d+", markers$gene) &    # 过滤核糖体大亚基蛋白基因
    !grepl("^RPS\\d+", markers$gene) &    # 过滤核糖体小亚基蛋白基因
    !grepl("^ENSCAFG", markers$gene)      # 过滤Ensembl犬基因ID
]
markersUsed=markers[which(markers$gene %in% TopGenes),]

data<-prepareDataFromscRNA(object=scRNA,
                                diffData=markersUsed,
                                showAverage=TRUE)
str(data)
enrich01<-enrichCluster(object=data,
                        OrgDb=org.Cf.eg.db,
                        type="BP",  # BP
                        organism="cfa",
                        pvalueCutoff=0.05,
                        topn=6,
                        seed=123)
enrich02<-enrichCluster(object=data,
                        OrgDb=org.Cf.eg.db,
                        type="KEGG",  # BP
                        organism="cfa",
                        pvalueCutoff=0.05,
                        topn=6,
                        seed=123)
resultsEnr = rbind(enrich01,enrich02)
write.csv(resultsEnr,file.path(projectPath, "Output/Integration20250829",'20250829_3SamplesMSC_FastMNNIntegration_FindAllMarkers_Res0.18_enrichCluster_YuanTop50.csv'))

enrich = read.csv(file.path(projectPath, "Output/Integration20250829",'20250829_3SamplesMSC_FastMNNIntegration_FindAllMarkers_Res0.18_enrichCluster_Yuan.csv'),row.names=1)
#check
head(enrich)
#绘制细胞平均表达热图
#add gene name 可以手动添加
markGenes=unique(markersUsed$gene)[sample(1:length(unique(markersUsed$gene)),40,
                                             replace=F)]
#line plot
visCluster(object=data,
           plot.type="line")
###热图
visCluster(
  object = data, 
  plot.type  = "heatmap", 
  column_names_rot = 45,
  markGenes = markGenes,
  cluster.order  = c(1:6)
)

pdf(file.path(projectPath, "Output/Integration20250829",'20250829_3SamplesMSC_FastMNNIntegration_FindAllMarkers_Res0.18visCluster_YuanTop50.pdf'),width=16,height=12)
mycol=c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#A6761D")   # "#E6AB02",  , "#666666"
heatmap01 = visCluster(
  object          = data,
  plot.type        = "both",
  column_names_rot = 60,
  show_row_dend   = FALSE,
  markGenes       = markGenes,
  markGenes.side   = "left",   ##将表达趋势折线图放到左侧
  genes.gp = c('italic',fontsize = 13),
  annoTerm.data    = enrich01,
  line.side        = "left",    #表达趋势图
  cluster.order    = c(1:6),
  sample.col=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FDBF6F","#FF7F00"),
  ctAnno.col = mycol,
  go.col = c(rep(mycol, each = 6)),#将词云颜色修改为跟注释条一致
  go.size = "pval", #按p值显示注释大小
  add.bar = F   #不要 注释的条形图
)
heatmap001 = visCluster(
  object          = data,
  plot.type        = "both",
  column_names_rot = 60,
  show_row_dend   = FALSE,
  markGenes       = markGenes,
  markGenes.side   = "left",   ##将表达趋势折线图放到左侧
  genes.gp = c('italic',fontsize = 13),
  annoTerm.data    = enrich02,
  line.side        = "left",    #表达趋势图
  cluster.order    = c(1:6),
  sample.col=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FDBF6F","#FF7F00"),
  ctAnno.col = mycol,
  go.col = c(rep(mycol, each = 6)),#将词云颜色修改为跟注释条一致
  go.size = "pval", #按p值显示注释大小
  add.bar = F   #不要 注释的条形图
)

heatmap02 = visCluster(
  object          = data,
  plot.type        = "both",
  column_names_rot = 60,
  show_row_dend   = FALSE,
  markGenes       = markGenes,
  markGenes.side   = "left",   ##将表达趋势折线图放到左侧
  genes.gp = c('italic',fontsize = 13),
  annoTerm.data    = enrich02,
  line.side        = "left",    #表达趋势图
  cluster.order    = c(1:6),
  sample.col=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FDBF6F","#FF7F00"),
  ctAnno.col = mycol,
  go.col = c(rep('#1B9E77', each = 5),rep('#D95F02', each = 2),rep('#7570B3', each = 6),
           rep('#E7298A', each = 6), rep('#66A61E', each = 4),rep('#A6761D', each = 6)),#将词云颜色修改为跟注释条一致
  go.size = "pval", #按p值显示注释大小
  add.bar = F   #不要 注释的条形图
)
heatmap03 = visCluster(
  object          = data,
  plot.type        = "both",
  column_names_rot = 60,
  show_row_dend   = FALSE,
  markGenes       = markGenes,
  markGenes.side   = "left",   ##将表达趋势折线图放到左侧
  genes.gp = c('italic',fontsize = 13),
  annoTerm.data    = enrich,
  line.side        = "left",    #表达趋势图
  cluster.order    = c(1:6),
  sample.col=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FDBF6F","#FF7F00"),
  ctAnno.col = mycol,
  go.col = c(rep(mycol, each = 6)),#将词云颜色修改为跟注释条一致
  go.size = "pval", #按p值显示注释大小
  add.bar = F   #不要 注释的条形图
)

dev.off()

###————— Conserved & specific
###############
library(SCP)
library(BiocParallel)
library(grid)

register(MulticoreParam(workers = 8, progressbar = TRUE))
# 重命名 S.Score 为 S_Score
names(objIntegrated@meta.data)[names(objIntegrated@meta.data) == "S.Score"] <- "S_Score"
# 重命名 G2M.Score 为 G2M_Score
names(objIntegrated@meta.data)[names(objIntegrated@meta.data) == "G2M.Score"] <- "G2M_Score"
save(objIntegrated,file = file.path(projectPath, "Data/Integration",paste0("20250617_3SamplesMSC_",method.use,"_FindClusters_objIntegrated.Rdata")))

load('/mnt2/wanggd_group/zjj/BGCscRNA/MSC_Duan/Data/Integration20250829/20250829_3SamplesMSC_FastMNNIntegration_FindClusters_objIntegrated.Rdata')
# 验证是否成功重命名
head(objIntegrated@meta.data)
objIntegrated$orig.ident[objIntegrated$orig.ident == "AD-MSC"] <- "AT-MSC"
# 验证修改结果
unique(objIntegrated$orig.ident)
save(objIntegrated,file = file.path(projectPath, "Data/Integration20250829",paste0("20250829_3SamplesMSC_FastMNNIntegration_FindClusters_objIntegrated.Rdata")))

objIntegrated[["RNA"]]<-as(object=objIntegrated[["RNA"]],Class="Assay")

#transition_data <- table(objIntegrated$Clusters, objIntegrated$Phase)
# 指定要展示的基因———— MSC specific genes  
##[ITGB1 (CD29), CD44, THY1 (CD90), ENG (CD105) and B2M,]  [NT5E (CD73), CD79A, PTPRC (CD45) and DRA, hardly detected ]
#genes_of_interest <- c("ITGB1", "CD44", "THY1", "ENG", "B2M", "NT5E", "CD79A", "DRA", "PTPRC")

#全的 参考  https://pmc.ncbi.nlm.nih.gov/articles/PMC11094970/ 		和 马的那篇文章	   

genes_of_interest <- c(
  "ITGB1", "CD44", 'PTPRC', 'B2M','CD79A','DRA', "GD2", "STRO-3", "FZD9", "PDGFRB", "SSEA-4", "TNFRSF10D", 
  "KIT", "LGR5", "ABCG2", "ENTPD1", "CD200", "LGR6", "EPHB2", "ROR2", 
  "MX1", "CMKLR1", "CD9", "SUSD2", "SDC2", "LY6A", "ITGAV", "ENG", 
  "HMMR", "PODXL", "CD274", "TNFAIP6", "THY1", "ITGA6", "F3", "CXCR4", 
  "GLI1", "STRO-1", "MCAM", "STRO-4", "EPHA7", "LEPR", "CD34", "ALCAM", 
  "CSPG4", "EPHA2", "NT5E", "ITGA1", "NCAM1", "NGFR", "PDGFRA", "NES", 
  "VCAM1", "BST2", "ISLR", "TLX1", "ALDH1A1", "SERPINF1", "S100A9", 
  "LRRC75A"
)

#典型标记物 CD73/NT5E、CD90/THY1、CD105/ENG 和 CD44
# 删掉 表达不高的   RUNX2   SOX9   PPARG
genes_of_interest <- c(
  "ITGB1", "CD44", "B2M", "CD79A",  "ROR2", "SDC2", "ITGAV", "TNFAIP6", 
  "THY1", "ITGA6", "F3", "ALCAM", "EPHA2", 'ITGA1',"NCAM1", 'SERPINF1'
)
					   
valid_genes <- intersect(genes_of_interest, rownames(objIntegrated))
pdf(file.path(projectPath, "Output/Integration20250829",'DotPlo_MSC_SpecificMarkers_in_AD-BM-UC_Added_ExpressedMarkers.pdf'),width = 6, height = 10)
#pdf(file.path(projectPath, "Output/Integration20250829",'DotPlo_MSC_SpecificMarkers_in_AD-BM-UC_AllMSCMarkers.pdf'),width = 6, height = 13)
ht <- GroupHeatmap(
  srt = objIntegrated,
  features = genes_of_interest,#genes,
  group.by = c("MSC_Clusters",'orig.ident'),
  heatmap_palette = 'YlOrRd',#"RdBu",
  group_palette = c("Paired",'Set3'),   
  cell_annotation = c("Phase",'G2M_Score','S_Score'),
  cell_annotation_palette = c("Dark2","Paired", "Paired"),  #Dark2
  show_row_names = TRUE, row_names_side = "left", 
  add_dot = TRUE, add_reticle = TRUE
)
print(ht$plot)
dev.off()
#png(file.path(projectPath, "Output/Integration20250829",'DotPlot_MSC_SpecificMarkers_in_AT-BM-UC_AllMSCMarkers_300dpi.png'),
png(file.path(projectPath, "Output/Integration20250829",'DotPlot_MSC_SpecificMarkers_in_AT-BM-UC_Added_ExpressedMarkers_300dpi.png'),
    width = 8, height = 14, units = "in", res = 300, pointsize = 12)
print(ht$plot)
dev.off()


# "ITGB1" "CD44"  "ENG"   "B2M"   "NT5E"  "CD79A" "PTPRC"
data("pancreas_sub")
print(pancreas_sub)
names(objIntegrated@assays)
srt = objIntegrated
# convert a v5 assay to a v3 assay
srt[["RNA"]]<-as(object=srt[["RNA"]],Class="Assay")
# convert a v3 assay to a v5 assay
#pbmc3k[["RNA"]] <- as(object = pbmc3k[["RNA"]], Class = "Assay5")

########--------########
library(patchwork)  # 确保加载patchwork包
pdf(file.path(projectPath, "Output/Integration20250829", 
                                  paste0("VlnPlo_MSC_SpecificMarkers_in_AD-BM-UC_Added_AllMarkers.pdf")),height =15,width=10)    #NoLegend()+
pdf(file.path(projectPath, "Output/Integration20250829", 
              "VlnPlo_MSC_SpecificMarkers_in_AD-BM-UC_AddedSeparately.pdf"),
    width = 10, height = 8)  # 调整宽高，宽度10高8更宽裕
n_genes_per_page <- 2
n_genes <- length(valid_genes)

for(i in seq(1, n_genes, by = n_genes_per_page)) {
  plots_list <- list()
  
  for(j in 0:(n_genes_per_page - 1)) {
    idx <- i + j
    if(idx > n_genes) break
    
    g <- valid_genes[idx]
    
    p01 <- FeatureStatPlot(
      srt,
      stat.by = g,
      assay = "RNA",
      slot = "data",
      group.by = 'orig.ident',
      plot.by = "group",
      fill.by = "group",
      legend.position = "none",
      ylab = "Expression level"
    ) + 
      xlab(NULL) + 
      theme(plot.title = element_text(hjust = 0.5))
    
    p1 <- VlnPlot(srt, 
                  cols = c("limegreen","navy","coral"),
                  pt.size = 0,
                  group.by = "orig.ident",
                  features = g, 
                  log = FALSE) + 
      theme(axis.title.x = element_blank(),
            panel.border = element_rect(fill = NA, color = "black"),
            axis.title.y = element_text(color = "black", size = 12),
            legend.position = 'none') + 
      geom_boxplot(width=0.1, fill="white", outlier.shape = NA)
    
    p_combined <- p01 / p1   # 上下排列 FeatureStatPlot 在上，VlnPlot 在下
    plots_list[[length(plots_list) + 1]] <- p_combined
  }
  
  # 把每个基因两图组合成一页的左右两列排列
  page_plot <- wrap_plots(plots_list, ncol = n_genes_per_page, byrow = TRUE) +
               plot_layout(guides = 'collect') &
               theme(plot.margin = margin(10, 10, 10, 10))
  
  print(page_plot)
}
dev.off()

####################### Clusters specific marker genes
##################
### 特定marker分子 osteoblast (成骨), chondrocyte (成软骨), and adipocyte(成脂)   
load(file.path(projectPath, "Data/Integration20250829",paste0("20250829_3SamplesMSC_FastMNNIntegration_FindClusters_objIntegrated.Rdata")) )


objIntegrated$celltype=objIntegrated$MSC_Clusters
Idents(objIntegrated) = 'MSC_Clusters'
objMarkers <- FindAllMarkers(objIntegrated, only.pos = TRUE)  # only.pos = TRUE 仅保留正向标记
# 提取每个簇 log2FC >1 的前 10 个标记基因
objMarkers = all.markers
top10 <- objMarkers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 15) %>%  # 取前 10 个
  ungroup()
library(dplyr)
library(stringr)
# 删掉 gene 以 "ENSCAFG" 开头的行，并提取为向量
# 过滤后保留的 features 向量
features <- top10 %>%
  dplyr::filter(!grepl("^ENSCAFG", gene),
                !grepl("^MT(-|\\.)", gene)) %>%
  dplyr::pull(gene)
markers = unique(features)  
pdf(file.path(projectPath, "Output/Integration20250829",paste0("20250913_3SamplesMSC_FastMNNIntegration_FindClusters_objIntegrated.TopMarkers.Dotplot.pdf")),width = 13)
DotPlot(objIntegrated, features = markers,
        cols = c("lightgrey", "#E64B35FF"),  # 颜色梯度：灰→红
        dot.scale = 6,                       # 点大小缩放因子
        scale = TRUE) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),  # X轴标签垂直
    axis.text.y = element_text(face = "italic")   +                  # Y轴斜体 
    scale_size_continuous(range = c(0.5, 5))  # 点大小范围
)  
DotPlot(objIntegrated, features = markers,
        cols = c("#F0F0F0", "#2166AC"),
        dot.scale = 5,
        scale = TRUE) +
  theme_minimal(base_size = 12) +   # 可选：整体基础字号
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    legend.position = "bottom",
    legend.text = element_text(size = 10)  # 可选：图例字号
  ) +
  labs(title = "Cell Cluster Marker Expression") +
  scale_size_continuous(range = c(0.5, 5))
dev.off()

###  
genesSpecific <- c(
  'ACADVL','RHOB', 'ACTB','NEB',"NES", "CSPG4", "MCAM", "DNMT1", "EZH2", "TOP2A",  "CCNA2", #"E2F1", "MYEL2", "MKN67",
  "ID4", "SCX", "COL11A1", "PPARG", #"CEBPD", 
  "EBF2", "STMN1", "RRM2", "RACGAP1", "SPC25", "SMC2", "KIF20B", "CENPE", "KIF23", "CENPF", "TPX2",   #"HMGA2",  "DIAP1",
  "CXCL12", "COMP", "CLU", "CHI3L1", "LUM",  "ACTA2",   #"CTSL",
 "TPM2"   # "MYL6", 
)
colors2 <- c("#96C3D8", "#F5B375", "#C0937E", "#67A59B", "#A5D38F", "#8D75AF", "#F19294", "#E45D61", "#BDA7CB")
pdf(file.path(projectPath, "Output/Integration20250829",paste0("20250913_3SamplesMSC_FastMNNIntegration_FindClusters_objIntegrated.SpecificMarkers.Dotplot.pdf")),width = 13)
DotP01=DotPlot(objIntegrated, features = genesSpecific,
        cols = colors2,  # 颜色梯度：灰→红
        dot.scale = 6,                       # 点大小缩放因子
        scale = TRUE) +RotatedAxis()
  
DotP02=DotPlot(objIntegrated, features = genesSpecific,
        cols = c("#F0F0F0", "#2166AC"),
        dot.scale = 5,
        scale = TRUE) +
  theme_minimal(base_size = 12) +   # 可选：整体基础字号
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    legend.position = "bottom",
    legend.text = element_text(size = 10)  # 可选：图例字号
  ) +
  labs(title = "Cell Cluster Marker Expression") +
  scale_size_continuous(range = c(0.5, 5))
dev.off()
  
levels(objIntegrated$MSC_Clusters)
# [1] "C1" "C2" "C3" "C4" "C5" "C6"
