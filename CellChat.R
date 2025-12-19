#devtools::install_github("jinworks/CellChat",force =  T)
library(Seurat)
library(ggsignif)
library(ggpubr)
library(cowplot)
library(dplyr)
library(NMF)
library(ggalluvial)
library(CellChat)
library(pheatmap)
library(qs)
library(patchwork)
library(future)
options(future.globals.maxSize = 5* 1024^5)  # 增加全局大小限制


# --- 1. 设置路径和加载数据 ---
projectPath <- '/mnt2/wanggd_group/zjj/BGCscRNA/MSC_Duan'

load('/mnt2/wanggd_group/zjj/BGCscRNA/MSC_Duan/Data/Integration20250829/20250829_3SamplesMSC_FastMNNIntegration_FindClusters_objIntegrated.Rdata')

# 数据预处理
Idents(objIntegrated) <- 'origIdent_Clusters'
objIntegratedSub <- subset(objIntegrated, idents = c('AD-MSC_C7', 'UC-MSC_C7'), invert = TRUE)
Idents(objIntegratedSub) <- 'MSC_Clusters'
obj <- objIntegratedSub
# 检查细胞类型
print(table(Idents(obj)))

# --- 2. 创建 CellChat 对象 ---
# 提取数据矩阵和元数据
data.input <- GetAssayData(obj, assay = "RNA", layer = "data")  # 标准化后的数据
meta <- data.frame(labels = Idents(obj), row.names = colnames(obj))
colnames(meta) <- "labels"
# 创建 CellChat 对象
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

# --- 3. 设置配体-受体数据库 ---
# 使用人类数据库（如果是小鼠，改为 CellChatDB.mouse）
CellChatDB <- CellChatDB.human
# 可以选择使用全部数据库或子集
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")  # 只用分泌信号cytokines b and growth factors
CellChatDB.use <- CellChatDB  # 使用全部数据库
cellchat@DB <- CellChatDB.use

# --- 4. 预处理表达数据 ---
# 识别过表达的配体和受体
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4)  # 并行计算，可根据你的CPU核心数调整
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# --- 5. 推断细胞通讯网络 ---
# 计算通讯概率
cellchat <- computeCommunProb(cellchat, type = "triMean")
# 过滤掉细胞数太少的通讯
cellchat <- filterCommunication(cellchat, min.cells = 10)

# 在信号通路水平推断通讯
cellchat <- computeCommunProbPathway(cellchat)

# 聚合细胞通讯网络
cellchat <- aggregateNet(cellchat)
# --- 6. 提取并筛选与 C4 相关的通讯 ---
# 提取所有配体-受体对的通讯数据
df.net <- subsetCommunication(cellchat)
# 查看数据结构
head(df.net)
# =============================================================================
# CellChat 分析：生成其他核心可视化图表
# =============================================================================
# --- 1. 通讯数量与强度的整体视图 ---
# 这张图可以快速了解哪些细胞亚群是主要的信号发出者或接收者。
# 确保 projectPath 已定义
# projectPath <- '/mnt2/wanggd_group/zjj/BGCscRNA/MSC_Duan'
pdf(file = file.path(projectPath, "./Output/Integration20250829/CellChat_Overall_Interactions.pdf"), width = 12, height = 6)

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

print("1. 整体通讯网络图已保存。")

# --- 2. 基于信号通路的热图 ---
# 这张图展示了每个细胞亚群在发送和接收信号时的主要信号通路。
# 可以通过 vertex.receiver 参数调整顺序。
library(ComplexHeatmap)
p_heatmap <- netVisual_heatmap(cellchat, 
                               measure = "weight", 
                               color.heatmap = "Reds")

# b. 使用 pdf() 和 draw() 保存为 PDF 文件
pdf(file = file.path(projectPath, "./Output/Integration20250829/CellChat_Signaling_Pathway_Heatmap.pdf"), width = 10, height = 8)
# ComplexHeatmap 对象需要被显式“绘制”出来
draw(p_heatmap) 
# 关闭图形设备，完成保存
dev.off()
print("2. 信号通路热图已保存。")

###
cellchat@netP$pathways
# [1] "COLLAGEN"  "LAMININ"   "FN1"       "PTPRM"     "THBS"      "CDH"
# [7] "SEMA3"     "CADM"      "NEGR"      "PTN"       "NCAM"      "ADGRG"
# [13] "PROS"      "MPZ"       "PERIOSTIN" "VEGF"      "SLIT"      "PCDH"
# [19] "ADGRL"     "Glutamate" "JAM"       "EGF"       "NOTCH"     "EPHA"
# [25] "PDGF"
#pathways.show <- c("SEMA3") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

#### plotGeneExpression
pathway_list <- cellchat@netP$pathways
pdf(file = file.path(projectPath, "./Output/Integration20250829/CellChat_AllSigPathwayGenes.pdf"), width = 10, height = 8)
for (pathway_name in pathway_list) {
  # 打印当前正在处理的通路，方便追踪进度
  message(paste("正在为通路繪製图形:", pathway_name))
  # 使用 tryCatch 来捕获并跳过可能出错的通路
  # enriched.only = TRUE 时，如果一个通路没有富集的配体-受体对，plotGeneExpression会报错
  tryCatch({
    # 生成小提琴图
    # enriched.only = TRUE: 只显示在通讯推断中被认为是显著的（过表达）的基因
    # type = "violin": 可以换成 "dot" 来生成点图
    p <- plotGeneExpression(cellchat, 
                            signaling = pathway_name, 
                            enriched.only = TRUE, 
                            type = "violin")
    # 将生成的ggplot对象打印到PDF文件中
    print(p)
  }, error = function(e) {
    # 如果发生错误，打印一条消息并跳过这个通路
    message(paste("警告：无法为通路", pathway_name, "生成图形。可能是因为没有富集的基因。错误信息:", e$message))
  })
}
dev.off()

#######All Sig Pathways   https://theislab.github.io/interaction-tools/14-CellChat.html#5_network_analysis  
# 1. 通路总览热图
pathways.show <-cellchat@netP$pathways
pdf(file = file.path(projectPath, "./Output/Integration20250829/CellChat_IdentifyCommunicationPatterns.pdf"), width = 10, height = 8)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 3)
p_outgoing <- netAnalysis_dot(cellchat, pattern = "outgoing")
print(p_outgoing)
netAnalysis_river(cellchat, pattern = "outgoing",cutoff = 0.4)

cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = 3)
p_incoming <- netAnalysis_dot(cellchat, pattern = "incoming")
print(p_incoming)
netAnalysis_river(cellchat, pattern = "incoming",cutoff = 0.5)
dev.off()




######################################################################################################
# 筛选与 C4 相关的通讯（C4 作为 source 或 target）
df.C4 <- df.net %>%
  filter(source == "C4" | target == "C4")

# 进一步筛选：只保留 C4 与 C1,C2,C3,C5,C6 之间的通讯
other_clusters <- c("C1", "C2", "C3", "C5", "C6")
df.C4.filtered0 <- df.C4 %>%
  filter(
    (source == "C4" & target %in% other_clusters) |
    (target == "C4" & source %in% other_clusters)
  )

sigPathways=c('THBS','CADM','MPZ','VEGF','EGF','EPHA','SLIT','EGF','CADM')
df.C4.filtered <- df.C4.filtered0[df.C4.filtered0$pathway_name %in% sigPathways, ]

方法1：气泡图（Bubble Plot）— 最常用
r
# 使用 ggplot2 绑定配体-受体对进行气泡图可视化
library(ggplot2)
#pdf(file.path(projectPath, "./Output/Integration20250829/CellChat_C4_plots.pdf"),
#    width = 16, height = 10)
# 创建 interaction_name_2 列（如果没有的话）
df.C4.filtered$interaction <- paste(df.C4.filtered$ligand, df.C4.filtered$receptor, sep = " - ")
df.C4.filtered$communication <- paste(df.C4.filtered$source, df.C4.filtered$target, sep = " -> ")

# 气泡图
p1 <- ggplot(df.C4.filtered, aes(x = communication, y = interaction)) +
  geom_point(aes(size = prob, color = prob)) +
  scale_color_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = median(df.C4.filtered$prob)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 8)) +
  labs(x = "Communication (Source -> Target)", 
       y = "Ligand - Receptor", 
       color = "Comm. Prob.",
       size = "Comm. Prob.",
       title = "C4-related Cell-Cell Communications")

print(p1)
# 保存图片
ggsave(file.path(projectPath, "./Output/Integration20250829/C4_communication_bubble.pdf"), p1, width = 10, height = 12)
方法2：分开画 C4 发出和接收的信号
r
# C4 作为信号发送方
df.C4.source <- df.C4.filtered %>% filter(source == "C4")

# C4 作为信号接收方
df.C4.target <- df.C4.filtered %>% filter(target == "C4")

# --- C4 发出的信号 ---
p2 <- ggplot(df.C4.source, aes(x = target, y = paste(ligand, receptor, sep = " - "))) +
  geom_point(aes(size = prob, color = prob)) +
  scale_color_gradient(low = "lightblue", high = "darkred") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Target Cell", y = "Ligand - Receptor", 
       title = "Signals FROM C4 to Other Clusters")

# --- C4 接收的信号 ---
p3 <- ggplot(df.C4.target, aes(x = source, y = paste(ligand, receptor, sep = " - "))) +
  geom_point(aes(size = prob, color = prob)) +
  scale_color_gradient(low = "lightblue", high = "darkred") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Source Cell", y = "Ligand - Receptor", 
       title = "Signals TO C4 from Other Clusters")
# 组合图
p23 = p2 + p3
print(p23)
ggsave(file.path(projectPath, "./Output/Integration20250829/C4_outgoing_incoming_signals.pdf"), p2 + p3, width = 14, height = 10)
方法3：热图展示通讯强度
r
library(tidyr)
library(pheatmap)

# 创建宽格式矩阵
df.C4.filtered$pair <- paste(df.C4.filtered$ligand, df.C4.filtered$receptor, sep = "_")
df.C4.filtered$comm <- paste(df.C4.filtered$source, df.C4.filtered$target, sep = "_")

mat <- df.C4.filtered %>%
  select(pair, comm, prob) %>%
  pivot_wider(names_from = comm, values_from = prob, values_fill = 0) %>%
  column_to_rownames("pair") %>%
  as.matrix()

# 绘制热图
pheatmap(mat, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         color = colorRampPalette(c("white", "yellow", "red"))(100),
         main = "C4-related Communication Probability",
         fontsize_row = 8,
         fontsize_col = 10)
#方法4：使用 CellChat 内置函数（针对特定LR对）

# 获取 C4 相关的配体-受体对名称
LR.pairs <- unique(df.C4.filtered$interaction_name)

# 使用 CellChat 内置的气泡图函数
netVisual_bubble(cellchat, 
                 sources.use = "C4",  # C4 作为发送方
                 targets.use = c("C1", "C2", "C3", "C5", "C6"),
                 remove.isolate = TRUE)

netVisual_bubble(cellchat, 
                 sources.use = c("C1", "C2", "C3", "C5", "C6"),
                 targets.use = "C4",  # C4 作为接收方
                 remove.isolate = TRUE)

# 如果只想画特定的信号通路
netVisual_bubble(cellchat, 
                 sources.use = "C4",
                 targets.use = c("C1", "C2", "C3", "C5", "C6"),
                 signaling = c("CCL", "CXCL", "FGF"),  # 指定通路
                 remove.isolate = TRUE)
#方法5：弦图（Chord Diagram）

# 针对特定信号通路画弦图
netVisual_chord_gene(cellchat, 
                     sources.use = c("C4"), 
                     targets.use = c("C1", "C2", "C3", "C5", "C6"),
                     lab.cex = 0.8,
                     legend.pos.y = 30)





########################################



# 按照 prob（通讯概率）排序，筛选最显著的通讯
df.C4.top <- df.C4.filtered %>%
  arrange(desc(prob)) %>%
  head(50)  # 取前50个最显著的，可以根据需要调整

print(df.C4.top)

# --- 7. 绘制气泡图（类似上传的图片风格）---

# 方法1：使用 CellChat 内置函数绘制气泡图
# 筛选特定的配体-受体对
selected_LR <- data.frame(interaction_name = unique(df.C4.top$interaction_name))

# 绘制气泡图：C4 作为 source
p1 <- netVisual_bubble(cellchat, 
                       sources.use = "C4", 
                       targets.use = other_clusters,
                       signaling = NULL,
                       pairLR.use = selected_LR,
                       remove.isolate = TRUE,
                       sort.by.target = TRUE,
                       angle.x = 45) +
  ggtitle("C4 -> Other Clusters") +
  theme(plot.title = element_text(hjust = 0.5))

# 绘制气泡图：C4 作为 target
p2 <- netVisual_bubble(cellchat, 
                       sources.use = other_clusters, 
                       targets.use = "C4",
                       signaling = NULL,
                       pairLR.use = selected_LR,
                       remove.isolate = TRUE,
                       sort.by.source = TRUE,
                       angle.x = 45) +
  ggtitle("Other Clusters -> C4") +
  theme(plot.title = element_text(hjust = 0.5))

# 合并图片   类似于https://www.nature.com/articles/s41467-020-18916-5/figures/6
p_combined <- p1 + p2

# 保存图片
ggsave(filename = file.path(projectPath, "./Output/Integration20250829/CellChat_C4_bubble_plot.pdf"),
       plot = p_combined, width = 16, height = 10)
ggsave(filename = file.path(projectPath, "./Output/Integration20250829/CellChat_C4_bubble_plot.png"),
       plot = p_combined, width = 16, height = 10, dpi = 300)



###################################
# 查看数据库中的信号通路名称
unique(CellChatDB.human$interaction$pathway_name)
# 或者查看你cellchat对象中富集的信号通路
unique(cellchat@LR$LRsig$pathway_name)

# Cytokines (趋化因子通路)
# Growth factors (生长因子通路)

# 提取富集的配体-受体对（使用通路名称）
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL", "CXCL", 
                                                         "VEGF", "FGF", 
                                                         "HGF", "IGF", 
                                                         "PDGF", "PGF"))))
pairLR.use
netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1,2,3,4,6,7,8), remove.isolate = FALSE,pairLR.use = pairLR.use)



###################################




# --- 8. 自定义气泡图---
# =============================================================================
# 合并气泡图：修正版
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# --- 准备数据（这部分与你原来的一样）---
df.plot.combined <- df.C4.filtered %>%
  mutate(
    LR_pair = gsub("_", " ", interaction_name_2),
    neg_log10_pval = -log10(pval + 1e-10),
    log2_prob = log2(prob + 1e-10),
    direction = ifelse(source == "C4", "C4_source", "C4_target"),
    x_label = ifelse(source == "C4", target, source)
  ) %>%
  filter(pval < 0.05)

df.plot.combined <- df.plot.combined %>%
  mutate(x_combined = paste0(direction, "_", x_label))

top_LR <- df.plot.combined %>%
  group_by(LR_pair) %>%
  summarise(max_prob = max(prob)) %>%
  arrange(desc(max_prob)) %>%
  head(30) %>% # 增加到30个，确保有足够的数据点
  pull(LR_pair)

df.plot.final <- df.plot.combined %>%
  filter(LR_pair %in% top_LR)

# --- 动态设置 X 轴顺序（核心修正部分）---
# 1. 获取实际存在于数据中的 source 和 target 细胞类型
source_clusters_present <- unique(df.plot.final$x_label[df.plot.final$direction == "C4_target"])
target_clusters_present <- unique(df.plot.final$x_label[df.plot.final$direction == "C4_source"])

# 2. 动态生成 X 轴顺序
x_order <- c(
  paste0("C4_source_", sort(target_clusters_present)),
  paste0("C4_target_", sort(source_clusters_present))
)

# 3. 将 x_combined 设置为因子，并指定顺序
df.plot.final$x_combined <- factor(df.plot.final$x_combined, levels = x_order)

# 设置 Y 轴顺序
df.plot.final$LR_pair <- factor(df.plot.final$LR_pair, 
                                 levels = rev(sort(unique(df.plot.final$LR_pair))))

# --- 创建 X 轴标签 ---
x_labels <- paste0('C',gsub("C4_source_|C4_target_", "", levels(df.plot.final$x_combined)))

# --- 绘制合并的气泡图（使用动态坐标）---
# 计算分组标注的位置
n_source_group <- length(target_clusters_present)
n_target_group <- length(source_clusters_present)

p_merged <- ggplot(df.plot.final, aes(x = x_combined, y = LR_pair)) +
  geom_point(aes(size = neg_log10_pval, color = log2_prob)) +
  scale_color_gradient2(
    low = "blue", mid = "yellow", high = "red",
    midpoint = median(df.plot.final$log2_prob, na.rm = TRUE),
    name = "log2(mean\n(ligand, receptor))"
  ) +
  scale_size_continuous(
    range = c(1, 8), name = "-log10(P value)",
    breaks = c(2, 4, 6, 8, 10) # 调整刻度
  ) +
  scale_x_discrete(labels = x_labels, drop = FALSE) + # drop=FALSE 确保所有x轴位置都显示
  
  # --- 动态添加分组标注 ---
  # 第一组 (C4 -> Others)
  annotate("segment", x = 0.5, xend = n_source_group + 0.5, y = -0.5, yend = -0.5, color = "black", linewidth = 0.8) +
  annotate("text", x = (1 + n_source_group) / 2, y = -1.2, label = "C4 -> Others", size = 4) +
  
  # 第二组 (Others -> C4)
  annotate("segment", x = n_source_group + 0.5, xend = n_source_group + n_target_group + 0.5, y = -0.5, yend = -0.5, color = "black", linewidth = 0.8) +
  annotate("text", x = n_source_group + (1 + n_target_group) / 2, y = -1.2, label = "Others -> C4", size = 4) +
  
  # 添加垂直分隔线
  {if (n_source_group > 0 && n_target_group > 0) 
    geom_vline(xintercept = n_source_group + 0.5, linetype = "dashed", color = "grey50", linewidth = 0.5)
  } +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"), # 调整角度和颜色
    axis.text.y = element_text(color = "black"),
    axis.title = element_blank(),
    panel.grid.major = element_line(color = "grey95"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.margin = margin(t = 10, r = 10, b = 40, l = 10)
  ) +
  coord_cartesian(clip = "off")

# 显示图形
print(p_merged)

# 保存
ggsave(filename = file.path(projectPath, "./Output/Integration20250829/CellChat_C4_merged_bubble_plot_FIXED.pdf"),
       plot = p_merged, width = 14, height = 12)
ggsave(filename = file.path(projectPath, "./Output/Integration20250829/CellChat_C4_merged_bubble_plot_FIXED.pdf"),
       plot = p_merged, width = 14, height = 12, dpi = 300)

# --- 9. 其他可视化（可选）---

# 圆形图：显示 C4 与其他 cluster 的通讯强度
pdf(file.path(projectPath, "./Output/Integration20250829/CellChat_C4_circle_plot.pdf"), width = 10, height = 10)
netVisual_circle(cellchat@net$weight, 
                 sources.use = which(levels(cellchat@idents) == "C4"),
                 targets.use = which(levels(cellchat@idents) %in% other_clusters),
                 vertex.weight = as.numeric(table(cellchat@idents)),
                 weight.scale = TRUE, 
                 label.edge = FALSE,
                 title.name = "Interaction weights/strength: C4")
dev.off()



# --- 10. 保存 CellChat 对象 ---
saveRDS(cellchat, file = file.path(projectPath, "./Output/Integration20250829/CellChat_MSC_object.rds"))

# 保存筛选后的通讯结果
write.csv(df.C4.filtered, 
          file = file.path(projectPath, "CellChat_C4_communications.csv"),
          row.names = FALSE)

print("分析完成！")
print(paste("结果已保存到:", projectPath))









###############################################
library(celltalker)
library(Seurat)
library(SeuratData)
library(tidyverse)
#运行celltalker
objInteractions <- celltalk(input_object=obj,
                                metadata_grouping="MSC_Clusters",
                                ligand_receptor_pairs=ramilowski_pairs,
                                number_cells_required=100,
                                min_expression=1000,
                                max_expression=20000,
                                scramble_times=10)
## Identify top statistically significant interactions
top_stats <- objInteractions %>%
  mutate(fdr=p.adjust(p_val,method="fdr")) %>%
  filter(fdr<0.05) %>%
  group_by(cell_type1) %>%
  top_n(3,interact_ratio) %>%
  ungroup()
#绘制circos图
pdf(file = file.path(projectPath, "./Output/Integration20250829/Celltalk_Overall_Interactions.pdf"), width = 6, height = 6)
colors_use <- RColorBrewer::brewer.pal(n=length(unique(obj$MSC_Clusters)),"Set2")
circos_plot(ligand_receptor_frame=top_stats,
            cell_group_colors=colors_use,
            ligand_color="blue",
            receptor_color="red",
            cex_outer=0.5,
            cex_inner=0.4)
dev.off()
