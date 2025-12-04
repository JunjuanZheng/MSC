#################################################################
############ ------  dorothea TF分析
#################################################################

library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(dorothea)
library(Seurat)

projectPath='/mnt2/wanggd_group/zjj/BGCscRNA/MSC_Duan'

load('/mnt2/wanggd_group/zjj/BGCscRNA/MSC_Duan/Data/Integration20250829/20250829_3SamplesMSC_FastMNNIntegration_FindClusters_objIntegrated.Rdata')
Idents(objIntegrated) = 'origIdent_Clusters'
objIntegratedSub <- subset(objIntegrated, idents = c('AD-MSC_C7', 'UC-MSC_C7'), invert = TRUE)   

inputs_dir <- system.file("extdata", package = "decoupleR")
## We read Dorothea Regulons for Human:
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))

##1. We obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C"))


#2.计算细胞的TF活性----
## We compute Viper Scores 
Idents(objIntegratedSub) = 'MSC_Clusters'
obj=objIntegratedSub
obj <- run_viper(obj, regulon,
                options = list(method = "scale", minsize = 4, 
                               eset.filter = FALSE, cores = 4, 
                               verbose = FALSE))
obj@assays$dorothea
save(obj, file = file.path(projectPath, "Data/Integration20250829",paste0("20251202_3SamplesMSC_FastMNNIntegration_DoRothEA.Results_obj.Rdata")))

load(file.path(projectPath, "Data/Integration20250829",paste0("20251202_3SamplesMSC_FastMNNIntegration_DoRothEA.Results_obj.Rdata")))
#3 对dorothea矩阵进行降维聚类分群-----
## We compute the Nearest Neighbours to perform cluster

DefaultAssay(object = obj) <- "dorothea"
obj <- ScaleData(obj)
obj <- RunPCA(obj, features = rownames(obj), verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:10, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, reduction = 'fastMNN',reduction.name = 'umap.fastMNN')

obj <- RunUMAP(obj, dims = 1:10, umap.method = "uwot", metric = "cosine")
Idents(obj) = 'MSC_Clusters'
obj.markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, 
                               logfc.threshold = 0.25, verbose = FALSE)

top10 <- obj.markers %>% group_by(cluster) %>% top_n(15, avg_log2FC);top10
obj@assays$dorothea@data[1:4,1:4]
top10 =top10[top10$avg_log2FC > 0.25,] 
TopGenes=unique(top10$gene)

TopGenesFiltered = TopGenes[
  !grepl("^MT-", TopGenes) &
  !grepl("^RPL\\d+", TopGenes) &
  !grepl("^RPS\\d+", TopGenes) &
  !grepl("^ENSCAFG", TopGenes)  # 或者可以用 gene_short_name 看你 ID 在哪一列
]
sub_obj <- ScaleData(subset(obj,downsample=100),features=TopGenes)

# 可视化 一
pdf(file.path(projectPath, "Output/Integration20250829",paste0('20251202_3SamplesMSC_FastMNNIntegration.DoHeatmapTop15.pdf')),width = 9,height =12)
DoHeatmap(sub_obj ,TopGenesFiltered,size=3,slot = 'scale.data')
dev.off()
# 可视化 二

pdf(file.path(projectPath, "Output/Integration20250829",paste0('20251202_3SamplesMSC_FastMNNIntegration.DoRothEAResult.pdf')),width = 9,height =12)
DoHeatmap(sub_obj ,TopGenesFiltered,size=3,slot = 'scale.data')

d.all=sub_obj
head(d.all@meta.data)
Idents(d.all)= 'MSC_Clusters'
a=AverageExpression(d.all,return.seurat = TRUE)
  
a$orig.ident=rownames(a@meta.data)
head(a@meta.data);a
markers=obj.markers ; head(markers)
markers$cluster=factor(markers$cluster,levels = unique(markers$cluster) )
DoHeatmap(a,draw.lines = FALSE, slot = 'scale.data',assay = 'dorothea',
            features = markers %>%group_by(cluster) %>%dplyr::slice_max(avg_log2FC,n = 5) %>% .$gene )

DoHeatmap(a,draw.lines = FALSE, slot = 'scale.data',assay = 'dorothea',
            features = markers %>%group_by(cluster) %>%dplyr::slice_max(avg_log2FC,n = 25) %>% .$gene )

#可视化 三
#上面进行可视化的TF都是通过findmarkers找到的，我们还可以根据std值找TF
#获取Viper得分矩阵，评价细胞群的TF活性----
## We transform Viper scores, scaled by seurat, into a data frame to better 
## handling the results
viper_scores_df <- GetAssayData(obj, slot = "scale.data", 
                                assay = "dorothea") %>%
  data.frame(check.names = F) %>%
  t();viper_scores_df[1:4,1:4]
## We create a data frame containing the cells and their clusters
CellsClusters <- data.frame(cell = names(Idents(obj)), 
                            cell_type = as.character(Idents(obj)),
                            check.names = F) #也可以使用其他的分类信息
CellsClusters[1:4, ]
## We create a data frame with the Viper score per cell and its clusters
viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters);viper_scores_clusters[1:4,]
## We summarize the Viper scores by cellpopulation
summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity));summarized_viper_scores
var(1,3  )

#5细胞群间变化最大的20个TFs进行可视化----
## We select the 20 most variable TFs. (150*10 populations = )  9个细胞亚群
highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>% # 计算变异度 calculates the variance of the "avg" column. So, for each transcription factor, it computes the variance of its average score.
  ungroup() %>%
  top_n(900, var) %>%
  distinct(tf);highly_variable_tfs
## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) ;summarized_viper_scores_df

palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length, 
                   max(summarized_viper_scores_df), 
                   length.out=floor(palette_length/2)))

viper_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=14, 
                       fontsize_row = 10, 
                       color=my_color, breaks = my_breaks, 
                       #main = "DoRothEA", 
                       angle_col = 45,
                       treeheight_col = 0,  border_color = NA)
pdf(file.path(projectPath, "Output/Integration20250829",paste0('20251202_3SamplesMSC_FastMNNIntegration.DoRothEAResult_pheatmap.pdf')),width = 9,height =18)

print(viper_hmap)

dev.off()



#################################################################
############ ------PROGENy--单细胞通路活性评分
#################################################################
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(Seurat)
library(progeny)
#BiocManager::install("OmnipathR")
library(decoupleR)
library(viridisLite)
load('/mnt2/wanggd_group/zjj/BGCscRNA/MSC_Duan/Data/Integration20250829/20250829_3SamplesMSC_FastMNNIntegration_FindClusters_objIntegrated.Rdata')
Idents(objIntegrated) = 'origIdent_Clusters'
objIntegratedSub <- subset(objIntegrated, idents = c('AD-MSC_C7', 'UC-MSC_C7'), invert = TRUE)   
obj=objIntegratedSub

Idents(obj) = 'MSC_Clusters'
# We create a data frame with the specification of the cells that belong to 
## each cluster to match with the Progeny scores. 
CellsClusters <- data.frame(Cell = names(Idents(obj)), 
                            CellType = as.character(Idents(obj)), 
                            stringsAsFactors = FALSE) 
head(CellsClusters)
## We compute the Progeny activity scores and add them to our Seurat object 
## as a new assay called Progeny. 
obj <- progeny(obj, scale=FALSE, organism="Human", top=500, perm=1, return_assay = TRUE)  #"Human" Mouse
obj@assays$progeny  ## 查看 Progeny 分析结果
obj@assays$progeny %>%dim()
obj@assays$progeny@data[,1:19]
save(obj, file = file.path(projectPath, "Data/Integration20250829",paste0("20251203_3SamplesMSC_FastMNNIntegration_PROGENy.Results_obj.Rdata")))

# Assay data with 14 features for 2638 cells # First 10 features: 
# Androgen, EGFR, Estrogen, Hypoxia, JAK-STAT, MAPK, NFkB, p53, PI3K, TGFb
#https://github.com/ixxmu/mp_duty/issues/5597
# 对 Progeny assay 数据进行标准化
obj <- Seurat::ScaleData(obj, assay = "progeny")
# 提取 Progeny scores 并转换为数据框
progeny_scores_df <- as.data.frame(
  t(GetAssayData(obj, slot = "scale.data", assay = "progeny"))
) %>% 
  rownames_to_column("Cell") %>%                          # 将行名转换为列
  gather(Pathway, Activity, -Cell)                        # 对数据进行变形，将 pathway 和值放在一起
# 对应细胞类型数据框
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)
# 总结 Progeny scores 根据细胞类型进行分组
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%                         # 按 Pathway 和 CellType 进行分组
  summarise(avg = mean(Activity),                         # 计算活动的平均值
            std = sd(Activity))                           # 计算活动的标准差
# 处理总结的数据框，以便绘图
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%                                 # 删除标准差列
  spread(Pathway, avg) %>%                               # 将数据展平，使每个 Pathway 为一个列
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) # 转换为数据框

#Plot 设置热图所需参数
paletteLength <- 100  # 调色板的长度
# 创建调色板
myColor <- colorRampPalette(c("#5e8ac2", "white", "#c04343"))(paletteLength)

# 生成用于热图的 breaks 值
progenyBreaks <- c(
  seq(min(summarized_progeny_scores_df), 0, length.out = ceiling(paletteLength / 2) + 1),
  seq(max(summarized_progeny_scores_df) / paletteLength, max(summarized_progeny_scores_df), 
      length.out = floor(paletteLength / 2))
)
# 保存热图到 PDF 文件
dataPlotTemp = t(summarized_progeny_scores_df[,-1]) #转置数据框，去掉第一列（通常是细胞名或路径名）
#colnameOrdered = c('AD-MSC_C0', 'BM-MSC_C0', 'UC-MSC_C0', 
#             'AD-MSC_C1', 'BM-MSC_C1', 'UC-MSC_C1', 
#             'AD-MSC_C2', 'BM-MSC_C2', 'UC-MSC_C2', 
#             'AD-MSC_C3', 'BM-MSC_C3', 'UC-MSC_C3', 
#             'AD-MSC_C4', 'BM-MSC_C4', 'UC-MSC_C4', 
#             'AD-MSC_C5', 'BM-MSC_C5', 'UC-MSC_C5', 
#             'AD-MSC_C6', 'BM-MSC_C6', 'UC-MSC_C6', 
#			 'BM-MSC_C7',
#             'AD-MSC_C8', 'BM-MSC_C8', 'UC-MSC_C8', 
#             'AD-MSC_C9', 'BM-MSC_C9', 'UC-MSC_C9')
colnameOrdered = c('C1','C2','C3','C4','C5','C6')
dataPlot=dataPlotTemp[,colnameOrdered]
pdf(file.path(projectPath, "Output/Integration20250829",paste0('20251203_3SamplesMSC_FastMNNIntegration.PROGENyResult.pdf')),width = 9)
# 创建热图
progeny_hmap01 <- pheatmap(
  dataPlot,     
  fontsize = 14,                             # 主字体大小
  fontsize_row = 10,                         # 行字体大小
  color = myColor,                           # 使用自定义的颜色调色板
  breaks = progenyBreaks,                   # 使用生成的 breaks
  main = "PROGENy",                   # 热图标题
  angle_col = 45,                           # 列标签旋转角度
  treeheight_col = 0,                       # 列树的高度（设置为0表示不显示树状图）
  border_color = NA,                         # 不显示边框颜色
  cluster_cols = FALSE              # 禁用列聚类，按原始顺序排列
)
print(progeny_hmap01)

progeny_hmap02 = pheatmap(dataPlot,
                        fontsize=12, 
                        fontsize_row = 10, color=turbo(90), #"inferno" "magma" "cividis" "viridis"
                        #  breaks = progenyBreaks,
                        main = "PROGENy",
                        cluster_cols = FALSE,             				
                        angle_col =90, treeheight_col = 0, border_color = NA)
print(progeny_hmap02)
						
library(viridis)
library(scop)
DefaultAssay(obj) <- 'progeny'   #	'JAK-STAT' 'NFkB' 'TNFa'
Idents(obj) = obj$Clusters
FeaturePlot(obj,features = "NFkB", coord.fixed = T,label = TRUE,
                order = T, cols = viridis(10)) 
FeaturePlot(obj,features = "JAK-STAT", coord.fixed = T, order = T, 
				cols =viridis::viridis(256)   #viridis::cividis(256)  viridis::plasma(256)      viridis::turbo(10)
				) 
FeaturePlot(obj,features = "TNFa", coord.fixed = T,label = TRUE,
                order = T, cols = viridis(10)) 


FeatureDimPlot(obj, features = c('JAK-STAT','NFkB','TNFa','p53'), palette = "Spectral",
				reduction = "UMAP.FastMNN", theme_use = "theme_blank")

FeatureDimPlot(obj, features = c('JAK-STAT'), palette = "Spectral",
				reduction = "UMAP.FastMNN", theme_use = "theme_blank")
FeatureDimPlot(obj, features = c('NFkB'), palette = "Spectral",
				reduction = "UMAP.FastMNN", theme_use = "theme_blank")
FeatureDimPlot(obj, features = c('TNFa'), palette = "Spectral",
				reduction = "UMAP.FastMNN", theme_use = "theme_blank")
FeatureDimPlot(obj, features = c('p53'), palette = "Spectral",
				reduction = "UMAP.FastMNN", theme_use = "theme_blank")


dev.off()


##########################################################################
############ ------Cytosig--The cytokine signaling activity of immune cells
############ https://doi.org/10.1038/s41592-021-01274-5    https://cytosig.ccr.cancer.gov/ 线上分析
#################################################################
### The cytokine signaling activity of immune cells from scRNA-sequencing data was predicted using CytoSig (33). The scRNA-sequencing data from all cells in a cluster were averaged to create pseudo-bulk data. This was then used with the CytoSig Python package to return the predicted cytokine signaling activity scores for different immune cells or CD8 T cell subtypes.
### https://mp.weixin.qq.com/s/A2IjaBmkerBfSEyGhRPwjQ
library(Seurat)
load('/mnt2/wanggd_group/zjj/BGCscRNA/MSC_Duan/Data/Integration20250829/20250829_3SamplesMSC_FastMNNIntegration_FindClusters_objIntegrated.Rdata')
Idents(objIntegrated) = 'origIdent_Clusters'
objIntegratedSub <- subset(objIntegrated, idents = c('AD-MSC_C7', 'UC-MSC_C7'), invert = TRUE)   
obj=objIntegratedSub
sce_aver = AggregateExpression(object = obj, group.by = c('MSC_Clusters'))$RNA # 使用seurat自带求不同状态，细胞类型的平均值。
results =log2(as.matrix(sce_aver) +1)
outputFile <- file.path(projectPath, "Data/Integration20250829", paste0('20251204_3SamplesMSC_FastMNNIntegration.Cytosig.txt'))
filteredResults <- results[!startsWith(rownames(results), "ENSCAFG"), ]
write.table(filteredResults,
            file = outputFile,
            sep = "\t",          # 使用制表符分隔
            quote = FALSE,       # 不使用引号
            row.names = TRUE,    # 写入行名
            col.names = TRUE)    # 写入列名
write.table(sce_aver,file.path(projectPath, "Data/Integration20250829",paste0('20251204_3SamplesMSC_FastMNNIntegration.Cytosig.txt')))




