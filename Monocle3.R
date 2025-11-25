library(paletteer)
library(Seurat)
library(monocle3)
library(dplyr)
library(BiocParallel)
library(scop)
PrepareEnv()
projectPath='/mnt2/wanggd_group/zjj/BGCscRNA/MSC_Duan'

load('/mnt2/wanggd_group/zjj/BGCscRNA/MSC_Duan/Data/Integration20250829/20250829_3SamplesMSC_FastMNNIntegration_FindClusters_objIntegrated.Rdata')
scRNA=objIntegrated
scRNA$celltype=scRNA$MSC_Clusters
scRNA[["RNA"]]<-as(object=scRNA[["RNA"]],Class="Assay")  #Assay5转成Assay3

register(MulticoreParam(workers=4,progressbar=TRUE))
pdf(file.path(projectPath, "Output/Integration20250829","objIntegrated_Monocle3_Result.pdf"), height = 10, width = 10)
expression_matrix <- GetAssayData(scRNA, assay = 'RNA',slot = 'counts')
cell_metadata <- scRNA@meta.data 
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)
##构建细胞轨迹
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

###monocle3分析前预处理
# 归一化/预处理数据
# 这一步使用的PCA分析，dim数代表纳入的PCA数量
cds <- preprocess_cds(cds, num_dim = 25,norm_method = c("none"))
# 这个函数用于确认设定的dim数是否足够代表主要变异
plot_pc_variance_explained(cds)
#选择一个能够保证大部分Variance都纳入的PCA componets值。
# 可选(去批次处理)
cds <- align_cds(cds, alignment_group = "orig.ident")
#.降维并可视化
# 降维聚类，可选择UMAP、PCA或者TSNE
coloursBy = 'celltype'
cds <- reduce_dimension(cds,reduction_method='UMAP',preprocess_method = 'PCA')
###基因/轨迹共定位
#这里可以看基因和轨迹的共定位情况
#modelgenes <- "KRAS"
#plot_cells(cds,
#           cell_size=1.5,group_label_size=4,
#           genes=modelgenes,
#           label_cell_groups=FALSE,
#           show_trajectory_graph=FALSE)

#细胞聚类
cds <- cluster_cells(cds) #cluster your cells
plot_cells(cds, color_cells_by = 'partition')#coloursBy)

cds0 = cds

###Seurat的umap分布数据替换monocle3的umap数据
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(scRNA, reduction = "umap.fastMNN")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed


###.轨迹分析
# 轨迹推断
cds <- learn_graph(cds,
				   verbose=T,
                   use_partition=F,  #默认是T，T时是顾及全局的情况
                   learn_graph_control=list(minimal_branch_len=20,#在图修剪过程中要保留的分支直径路径的最小长度。默认值是10。
                                            euclidean_distance_ratio=10#生成树中两个末端节点的欧氏距离与生成树上允许连接的任何连接点之间的最大距离之比。默认值为1。
                      ))
				   
p1 = plot_cells(cds, 
           color_cells_by = coloursBy,
           label_groups_by_cluster=FALSE,
           cell_size=1,group_label_size=4,
           label_leaves= F, # 是否显示不同细胞结局
           label_branch_points=F, # 是否显示不同的分支节点
           trajectory_graph_color='#023858',
           trajectory_graph_segment_size = 1)
print(p1)
###定义起点-轨迹可视化
#接下来按照CytoTRACE2中的结果设定细胞发育的起点。
#定义root cell, 推断拟时方向
# 网页自定
#cds <- order_cells(cds)  
# 代码定
# a helper function to identify the root principal points:
## 指定初始细胞群和细胞亚群的列名
time_bin="C3"
get_earliest_principal_node <- function(cds, time_bin=time_bin, coloursBy) {
  # 确保 `coloursBy` 变量存在并正确
  if (!coloursBy %in% colnames(colData(cds))) {
    stop(paste("The column", coloursBy, "does not exist in colData(cds)"))
  }
  # 找到符合条件的细胞 ID
  cell_ids <- which(colData(cds)[, coloursBy] == time_bin)
  if(length(cell_ids) == 0) {
    stop(paste("No cells found for time_bin:", time_bin))
  }
  # 计算最近的主成分节点
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
					as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))
					]
  return(root_pr_nodes)
}
# 获取最近的主成分节点
root_pr_nodes <- get_earliest_principal_node(cds, time_bin=time_bin, coloursBy=coloursBy)
# 使用获得的节点顺序细胞
cds <- order_cells(cds, root_pr_nodes=root_pr_nodes)

#2D可视化
p2 = plot_cells(cds, label_cell_groups = F, 
           color_cells_by = "pseudotime", 
           label_branch_points = F, 
           label_roots =F, # 显示根
           label_leaves =F,
           graph_label_size = 0, 
           cell_size=2, 
           trajectory_graph_color='black',
           trajectory_graph_segment_size = 2)
print(p2)

p002 = plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                  label_leaves = FALSE, trajectory_graph_color = "purple", cell_size = 0.8,
                  label_branch_points = FALSE, 
                  label_roots = F) + 
  theme(legend.position="right", 
        # legend.key.size =  unit(0.2, "inches"), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 16), 
        text = element_text(size = 16), axis.text = element_text(size = 16)) + 
  # theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))  + 
  tidydr::theme_dr()+
  scale_color_gradientn(colours = c("#0088C3","#6FBCDB","#BCE5F1","#F2EFDA","#FFC482","#FF7B4C","#FF5831"))
print(p002)

p111=plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
              label_leaves = FALSE,  label_branch_points = FALSE,
              trajectory_graph_color = "white",
              cell_size=0.3) +
  theme_classic()+
  theme(axis.text.x = element_text(size=10,color="black"),
        axis.text.y = element_text(size=10,color="black"))+
  scale_color_manual(values  = c("red","blue"))
  
pdata<-p111$data

pp01 = ggplot(pdata,aes(x=data_dim_1,data_dim_2,color=cell_color))+
  geom_point(size=0.5)+
  labs(x="umap_1",y="umap_2",color="pseudotime",title="pseudotime")+
  scale_color_gradientn(colours  = c("#96b6e2","#dbe5f5","#f6c7c0","#eda69d","#ce564f"))+
  theme_classic()+
  theme(axis.text.x = element_text(size=10,color="black"),
        axis.text.y = element_text(size=10,color="black"),
        plot.title = element_text(hjust = 0.5)  # This centers the title
  )	
print(pp01)
  
pp02 = ggplot(pdata,aes(x=celltype,y=cell_color,fill=celltype))+
  geom_violin(scale=T,width=1)+
  scale_fill_manual(values =c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FDBF6F","#FF7F00"))+
  theme_classic()+
  labs(x="",y="pseudotime")+
  theme(axis.text.x = element_text(size=10,color="black"),
        axis.text.y = element_text(size=10,color="black"))
print(pp02)
#计算每个origin下的pse组成
pdata <- pdata %>%
  mutate(cell_color_category = case_when(
    cell_color >= 0 & cell_color < 10 ~ "pse 0-10",
    cell_color >= 10 & cell_color < 20 ~ "pse 10-20",
    cell_color >= 20 & cell_color <= 27 ~ "pse 20-30",
    TRUE ~ NA_character_  # 对于不在上述范围内的值，标记为 NA
  ))
# 计算每个 origin 分组下的 subtype 组成百分比
percent_data <- pdata %>%
  group_by(celltype, cell_color_category) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)
# 使用 ggplot2 绘制柱形图
ppp1=ggplot(percent_data, aes(x = celltype, y = percentage, fill = cell_color_category)) +
  geom_bar(stat = "identity", position = "stack",color="black") +
  labs(title = "          Pseudotime",
       x = "",
       y = "Percentage",
       fill = "Subtype") +
  scale_fill_manual(values =c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FDBF6F","#FF7F00"))+
  theme_classic()+theme(axis.text.x = element_text(size=10,color="black"),
                        axis.text.y = element_text(size=10,color="black"))
print(ppp1)
pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

pppp = ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(celltype, monocle3_pseudotime, median), fill = celltype)) +
  geom_boxplot() +
  scale_fill_manual(values =c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FDBF6F","#FF7F00"))
print(pppp)


###monocle3常见图谱
genes <- graph_test(cds, neighbor_graph="principal_graph", cores=10)

save(list = c("cds", "genes"), 
     file = file.path(projectPath, "Data/Integration20250829", 
                      "objIntegrated_Monocle3_cds_.Rdata"))

A=load(file.path(projectPath, "Data/Integration20250829","objIntegrated_Monocle3_cds_.Rdata"))

#挑选top10画图展示
uedeGenes <- genes %>% top_n(n=80, morans_I) %>%
  pull(gene_short_name) %>% as.character()
  
uedeGenes <- genes %>% subset(q_value < 0.05 & morans_I > 0.005) %>%
  pull(gene_short_name) %>% as.character()
  
TopGenes <- uedeGenes[
    !grepl("^MT-", uedeGenes) &        # 过滤线粒体基因
    !grepl("^RPL\\d+", uedeGenes) &    # 过滤核糖体大亚基蛋白基因
    !grepl("^RPS\\d+", uedeGenes) &    # 过滤核糖体小亚基蛋白基因
    !grepl("^ENSCAFG", uedeGenes)      # 过滤Ensembl犬基因ID
]  
  
#基因表达趋势图
pseuPlot = plot_genes_in_pseudotime(cds[TopGenes[1:10],], 
                         min_expr=0.5, ncol = 2)
print(pseuPlot)
						 
pseuPlot00 = plot_genes_in_pseudotime(cds[TopGenes[1:10],], color_cells_by="celltype", 
                         min_expr=0.5, ncol = 2)	
print(pseuPlot00)
						 
#FeaturePlot图
feaPlot=plot_cells(cds, genes=TopGenes, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
print(feaPlot)
dev.off()

###
# Filter significant genes
pdf(file.path(projectPath, "Output/Integration20250829","objIntegrated_Monocle3_Result02_TopOrderedGenes_.pdf"), height = 10, width = 10)
deg_genes <- subset(genes, q_value < 0.05)
deg_genes_Ordered = deg_genes[order(deg_genes$morans_I, decreasing = TRUE), ]
#head()
# Plot top dynamic gene
#top_gene <- rownames(deg_genes[which.max(deg_genes$morans_I), ])
top_gene <- rownames(deg_genes_Ordered)[1:20]

TopGenes <- top_gene[
    !grepl("^MT-", top_gene) &        # 过滤线粒体基因
    !grepl("^RPL\\d+", top_gene) &    # 过滤核糖体大亚基蛋白基因
    !grepl("^RPS\\d+", top_gene) &    # 过滤核糖体小亚基蛋白基因
    !grepl("^ENSCAFG", top_gene)      # 过滤Ensembl犬基因ID
]  
rowData(cds)$gene_short_name <- rownames(rowData(cds))

pseuPlot001 = plot_genes_in_pseudotime(cds[TopGenes,], 
                         min_expr=0.5, ncol = 2,color_cells_by = "pseudotime")
print(pseuPlot001)
dev.off()

library(RColorBrewer)
library(ggplot2)

# 获取 Paired 调色板
paired_colors <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FDBF6F","#FF7F00")#brewer.pal(length(unique(colData(cds)$celltype)), "Paired")
color_map <- setNames(paired_colors, levels(unique(colData(cds)$celltype)) )

temp= c('COL1A1','COL1A2','COL11A1','SERPINE1','ACTB','ACTA2')
#temp= c('RNASEK1','CRYAB1','COL1A1','COL1A2','COL11A1','SERPINE1')
pseuPlot0001 = plot_genes_in_pseudotime(cds[temp,], 
                         min_expr=0.5, ncol = 2,color_cells_by = "pseudotime")
pseuPlot00002 = plot_genes_in_pseudotime(cds[temp,], 
                         min_expr=0.5, ncol = 2,color_cells_by = "celltype") +
               scale_color_manual(values = color_map)						 
print(pseuPlot0001)
print(pseuPlot0002)
dev.off()

###########
TopGenes_use <- TopGenes[1:100]
pdf(file.path(projectPath, "Output/Integration20250829",
              "objIntegrated_Monocle3_Result02_TopOrderedGenes20251121.pdf"),
    height = 10, width = 10)
# 每页 6 个基因
genes_per_page <- 6
n_pages <- ceiling(length(TopGenes_use) / genes_per_page)
for (i in seq_len(n_pages)) {
  idx_start <- (i - 1) * genes_per_page + 1
  idx_end   <- min(i * genes_per_page, length(TopGenes_use))
  genes_this_page <- TopGenes_use[idx_start:idx_end]
  # 一次传入多个基因，ncol = 3 => 每页最多 2 行 × 3 列 = 6 个
  p <- plot_genes_in_pseudotime(
    cds[genes_this_page, ],
    min_expr = 0.5,
    color_cells_by = "celltype",
    ncol = 3
  ) + scale_color_manual(values = color_map)
  print(p)
}
dev.off()


#########################################################################################
library(patchwork)  # 或 library(cowplot)
TopGenes_use <- TopGenes[1:100]
						 
library(ClusterGVis)
mat=ClusterGVis::pre_pseudotime_matrix(cds_obj = cds, gene_list = TopGenes_use)
# kmeans
ck <- clusterData(as.matrix(mat),
                  cluster.method = "kmeans",
                  cluster.num = 3)



plots_list <- vector("list", length(TopGenes_use))
names(plots_list) <- TopGenes_use
# 先把每个基因的图生成出来
for (i in seq_along(TopGenes_use)) {
  tp <- TopGenes_use[i]
  plots_list[[i]] <- plot_genes_in_pseudotime(
    cds[tp, ],
    min_expr = 0.5,
    color_cells_by = "celltype"
  ) + scale_color_manual(values = color_map) +
    ggtitle(tp)  # 可以加上基因名
}
pdf(file.path(projectPath, "Output/Integration20250829",
              "objIntegrated_Monocle3_Result02_TopOrderedGenes20251121_9YiYe.pdf"),
    height = 10, width = 10)
genes_per_page <- 9
n_pages <- ceiling(length(plots_list) / genes_per_page)
for (i in seq_len(n_pages)) {
  idx_start <- (i - 1) * genes_per_page + 1
  idx_end   <- min(i * genes_per_page, length(plots_list))
  plots_this_page <- plots_list[idx_start:idx_end]
  
  # patchwork 把 list 的图拼成 3 列
  p_page <- wrap_plots(plots_this_page, ncol = 3)
  print(p_page)
}
dev.off()


#################################################################################
pr_graph_res <- graph_test(cds, neighbor_graph = "principal_graph")
# 按 morans_I 大小排序，选显著的动态基因
dyn_genes <- subset(pr_graph_res, q_value < 0.05)
dyn_genes <- dyn_genes[order(dyn_genes$morans_I, decreasing = TRUE), ]
# 行名一般是 gene_id，如果你有 gene_short_name，可以合并一下方便后面用
dyn_genes$gene_id <- rownames(dyn_genes)
write.csv(dyn_genes,file.path(projectPath, "Output/Integration20250829",
              "objIntegrated_Monocle3_morans_I_OrderedGenes.csv"))
#拟时序基因热图
top50 <- genes %>%
    top_n(n = 10, morans_I) %>%
    pull(gene_short_name) %>%
    as.character()



#  C4的Top marker genes
degs=read.csv(file.path(projectPath, "Data/Integration20250829",
			"20250829_3SamplesMSC_FastMNNIntegration_FindAllMarkers_Res0.18.csv"))
C4_markers = degs[degs$cluster == 'C4',]
C4_dyn_markers <- dyn_genes[dyn_genes$gene_short_name %in% C4_markers$gene, ]
write.csv(C4_dyn_markers,file.path(projectPath, "Output/Integration20250829",
              "objIntegrated_Monocle3_morans_I_OrderedGenes_C4TopMarkers.csv"))


#########################################################################################
library(patchwork)  # 或 library(cowplot)
uedeGenes <- C4_dyn_markers %>% subset(q_value < 0.05 & morans_I > 0.005) %>%
  pull(gene_short_name) %>% as.character()  
TopGenes <- uedeGenes[
    !grepl("^MT-", uedeGenes) &        # 过滤线粒体基因
    !grepl("^RPL\\d+", uedeGenes) &    # 过滤核糖体大亚基蛋白基因
    !grepl("^RPS\\d+", uedeGenes) &    # 过滤核糖体小亚基蛋白基因
    !grepl("^ENSCAFG", uedeGenes)      # 过滤Ensembl犬基因ID
]  
TopGenes_use <- TopGenes[1:100]
plots_list <- vector("list", length(TopGenes_use))
names(plots_list) <- TopGenes_use
# 先把每个基因的图生成出来
for (i in seq_along(TopGenes_use)) {
  tp <- TopGenes_use[i]
  plots_list[[i]] <- plot_genes_in_pseudotime(
    cds[tp, ],
    min_expr = 0.5,
    color_cells_by = "celltype"
  ) + scale_color_manual(values = color_map) +
    ggtitle(tp)  # 可以加上基因名
}
pdf(file.path(projectPath, "Output/Integration20250829",
              "objIntegrated_Monocle3_C4Degs_TopOrderedGenes20251124.pdf"),
    height = 10, width = 10)
genes_per_page <- 9
n_pages <- ceiling(length(plots_list) / genes_per_page)
for (i in seq_len(n_pages)) {
  idx_start <- (i - 1) * genes_per_page + 1
  idx_end   <- min(i * genes_per_page, length(plots_list))
  plots_this_page <- plots_list[idx_start:idx_end]
  
  # patchwork 把 list 的图拼成 3 列
  p_page <- wrap_plots(plots_this_page, ncol = 3)
  print(p_page)
}
dev.off()


ciliated_genes <- c("COL1A1",
                    "COL1A2",
                    "COL11A1",
                    "SERPINE1",
                    "ATF4",
                    "CCN1")  #SERPINE1 ATF4  ACTB  ACTA2
cds_subset <- cds[rowData(cds)$gene_short_name %in% ciliated_genes,]
plot_genes_violin(cds_subset, group_cells_by="celltype", ncol=2) +
      theme(axis.text.x=element_text(angle=45, hjust=1))

plot_genes_in_pseudotime(
    cds['ATF4', ],
    min_expr = 0.5,
    color_cells_by = "celltype"
  ) + scale_color_manual(values = color_map) +
    ggtitle('ATF4')  # 可以加上基因名




C4_cells <- colnames(cds)[colData(cds)$celltype == "C4"]
cds_C4   <- cds[, C4_cells]

# （可选）如果在 C4 中轨迹结构很简单，重新 learn_graph / order_cells
# cds_C4 <- learn_graph(cds_C4)
# cds_C4 <- order_cells(cds_C4)

gt_C4 <- graph_test(cds_C4, neighbor_graph = "principal_graph")
gt_C4_dyn <- subset(gt_C4, q_value < 0.05)






####################################
# 在 C4 中找“特异表达”基因
colData(cds)$C4_vs_others <- ifelse(colData(cds)$celltype == "C4", "C4", "Others")

# 拟合模型，formula 里面用这个分组变量
deg_model <- fit_models(
    cds,
    model_formula_str = "~ C4_vs_others"
)
deg_table <- coefficient_table(deg_model)
deg_table %>% filter(term != "(Intercept)") %>%
      select(gene_short_name, term, q_value, estimate)
#检查模型质量
evaluate_fits(deg_model)  #用来过滤掉表达细胞数太少的基因；拟合失败或不可靠的基因（status != "OK"）
#num_cells_expressed：该基因在多少个细胞里有表达
#status：模型是否正常（"OK" 才可信）
#null_deviance / deviance：空模型 vs 含分组变量模型的偏差
#AIC / BIC：模型优劣指标
library(dplyr)
## 1. 差异系数表（去掉截距）
coef_tbl <- coefficient_table(deg_model) %>%
  filter(term != "(Intercept)")
## 2. 拟合质量表
fit_info <- evaluate_fits(deg_model)
## 3. 只保留模型正常且表达细胞数足够的基因
fit_ok <- fit_info %>%
  filter(
    status == "OK",
    num_cells_expressed >= 10     # 阈值可以根据数据改，比如 >= 20
  ) %>%
  select(gene_id, gene_short_name, num_cells_expressed)
## 4. 合并两张表
deg_merged <- coef_tbl %>%
  inner_join(fit_ok, by = c("gene_id", "gene_short_name"))
## C4 高表达 marker（C4 > Others）
C4_markers <- deg_merged %>%
  filter(
    term == "C4_vs_othersOthers",
    q_value < 0.05,          # 显著
    estimate < -0.5          # 绝对差异足够大，阈值可调
  ) %>%
  arrange(estimate)  
## Others 高表达 marker（Others > C4）
Others_markers <- deg_merged %>%
  filter(
    term == "C4_vs_othersOthers",
    q_value < 0.05,
    estimate > 0.5
  ) %>%
  arrange(desc(estimate))
 
library(dplyr)
drop_list_cols <- function(df) {
  df %>% select(where(~ !("list" %in% class(.x))))
}
C4_markers_export     <- drop_list_cols(C4_markers)
C4_markers_filtered <- C4_markers_export %>%
  filter(
    !grepl("^MT-",      gene_short_name) &   # 过滤线粒体基因
    !grepl("^RPL\\d+",  gene_short_name) &   # 过滤核糖体大亚基蛋白基因
    !grepl("^RPS\\d+",  gene_short_name) &   # 过滤核糖体小亚基蛋白基因
    !grepl("^ENSCAFG",  gene_short_name)     # 过滤 Ensembl 犬基因 ID（没注释名字的）
  )
write.csv(
  C4_markers_filtered,
  file = file.path(
    projectPath, "Output/Integration20250829",
    "objIntegrated_Monocle3_fit_models_C4VS.others_C4markers.csv"
  ),
  row.names = FALSE
)
Others_markers_export <- drop_list_cols(Others_markers)
Others_markers_filtered <- Others_markers_export %>%
  filter(
    !grepl("^MT-",      gene_short_name) &
    !grepl("^RPL\\d+",  gene_short_name) &
    !grepl("^RPS\\d+",  gene_short_name) &
    !grepl("^ENSCAFG",  gene_id)          # 或 gene_short_name，看你 ID 在哪一列
  )
write.csv(
  Others_markers_filtered,
  file = file.path(
    projectPath, "Output/Integration20250829",
    "objIntegrated_Monocle3_fit_models_C4VS.others_Othersmarkers.csv"
  ),
  row.names = FALSE
)
#和拟时序动态基因（graph_test）结合
pr_graph_res <- graph_test(cds, neighbor_graph = "principal_graph")
dyn_genes <- pr_graph_res %>%
  filter(q_value < 0.05) %>%
  arrange(desc(morans_I)) # %>%
dyn_genes$gene_id = rownames(dyn_genes)   # 视你的对象结构调整

# C4_markers 做交集, C4 高表达 marker，又在拟时序上动态（Moran’s I 高）的基因
Temp = Others_markers_filtered #Others_markers_filtered #C4_markers_filtered
C4_dyn_markers <- Temp %>%
  inner_join(dyn_genes, by = "gene_id") %>%
  arrange(desc(morans_I))
###先已经用 q_value 过滤：例如只保留 q_value < 0.05
###在此基础上按 Moran’s I（动态程度） 排序
###再加一个绝对效应大小：比如 abs(estimate) 越大越靠前
C4_dyn_markers_top50 <- C4_dyn_markers %>%
  filter(q_value.y < 0.05) %>%                       # 保险起见再卡一次显著性
  arrange(desc(morans_I), desc(abs(estimate))) %>% # 先按动态程度，再按差异强度
  slice(1:50)
top50_genes <- unique(C4_dyn_markers_top50$gene_short_name.y)
length(top50_genes)  

p_violin <- plot_genes_violin(
  cds[top50_genesUsed[1:6], ],
  group_cells_by = "celltype", ncol=2) +
      theme(axis.text.x=element_text(angle=45, hjust=1))
p_violin

plots_list <- vector("list", length(top50_genesUsed))
names(plots_list) <- top50_genesUsed
# 先把每个基因的图生成出来
for (i in seq_along(TopGenes_use)) {
  tp <- TopGenes_use[i]
  plots_list[[i]] <- plot_genes_in_pseudotime(
    cds[tp, ],
    min_expr = 0.5,
    color_cells_by = "celltype"
  ) + scale_color_manual(values = color_map) +
    ggtitle(tp)  # 可以加上基因名
}
pdf(file.path(projectPath, "Output/Integration20250829",
              "objIntegrated_Monocle3_fit_models_C4VS.others_Othersmarkers20251125.pdf"),
    height = 10, width = 10)
genes_per_page <- 9
n_pages <- ceiling(length(plots_list) / genes_per_page)
for (i in seq_len(n_pages)) {
  idx_start <- (i - 1) * genes_per_page + 1
  idx_end   <- min(i * genes_per_page, length(plots_list))
  plots_this_page <- plots_list[idx_start:idx_end]
  
  # patchwork 把 list 的图拼成 3 列
  p_page <- wrap_plots(plots_this_page, ncol = 3)
  print(p_page)
}
dev.off()















# 提取 C4_vs_othersC4 这一项（monocle3 里会这样命名系数）
deg_C4 <- subset(
    deg_table,
    term == "C4_vs_othersC4"
)
# 取 log2 fold change 大于一定阈值的基因，作为更“特异”的基因
deg_C4_specific <- subset(deg_C4, estimate > 0.5)  # 可以根据需要调整阈值

C4_genes <- unique(deg_C4_specific$gene_short_name)
length(C4_genes)
head(C4_genes)








######################################
scRNA <- RunSlingshot(
  srt = scRNA,
  group.by = "celltype",
  reduction = "umap.fastMNN"
)
FeatureDimPlot(
  scRNA,
  features = paste0("Lineage", 1:2),
  reduction = "umap.fastMNN",
  theme_use = "theme_blank"
)
scRNA <- RunDynamicFeatures(
  srt = scRNA,
  lineages = c("Lineage1", "Lineage2"),
  n_candidates = 200
)

scRNA <- AnnotateFeatures(
  scRNA,
  species = "Homo_sapiens",
  db = c("TF", "CSPA")
)
ht <- DynamicHeatmap(
  srt = scRNA,
  lineages = c("Lineage1", "Lineage2"),
  use_fitted = TRUE,
  n_split = 6,
  reverse_ht = "Lineage1",
  species = "Homo_sapiens", #"Mus_musculus",
  db = "GO_BP",
  anno_terms = TRUE,
  anno_keys = TRUE,
  anno_features = TRUE,
  heatmap_palette = "viridis",
  cell_annotation = "celltype",
  separate_annotation = list(
    "celltype", c("S_Score", "G2M_Score")
  ),
  separate_annotation_palette = c("Paired", "Set1"),
  feature_annotation = c("TF", "CSPA"),
  feature_annotation_palcolor = list(
    c("gold", "steelblue"), c("forestgreen")
  ),
  pseudotime_label = 25,
  pseudotime_label_color = "red",
  height = 5,
  width = 2
)
print(ht$plot)

'NDC80','SPC25','CDK1','UBE2C','SGO1','CCNA2','HMGB2','CEP55' 'RYAB'
DynamicPlot(
  srt = scRNA,
  lineages = c("Lineage1", "Lineage2"),
  group.by = "celltype",
  features = c(
    'RRM2','SMC2','NUF2','UBE2C','NDC80','CENPF'
  ),
  compare_lineages = TRUE,
  compare_features = FALSE
)


