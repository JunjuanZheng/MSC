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
options(future.globals.maxSize = 5 * 1024^3)  # 2 GiB

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
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")  # 只用分泌信号
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

# 筛选与 C4 相关的通讯（C4 作为 source 或 target）
df.C4 <- df.net %>%
  filter(source == "C4" | target == "C4")

# 进一步筛选：只保留 C4 与 C1,C2,C3,C5,C6 之间的通讯
other_clusters <- c("C1", "C2", "C3", "C5", "C6")
df.C4.filtered <- df.C4 %>%
  filter(
    (source == "C4" & target %in% other_clusters) |
    (target == "C4" & source %in% other_clusters)
  )

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

# --- 8. 自定义气泡图（更接近上传图片的风格）---
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


# =============================================================================
# CellChat 分析：生成其他核心可视化图表
# =============================================================================

# 确保 projectPath 已定义
# projectPath <- '/mnt2/wanggd_group/zjj/BGCscRNA/MSC_Duan'

# --- 1. 通讯数量与强度的整体视图 ---
# 这张图可以快速了解哪些细胞亚群是主要的信号发出者或接收者。
pdf(file = file.path(projectPath, "./Output/Integration20250829/CellChat_Overall_Interactions.pdf"), width = 12, height = 6)
# 设置布局：1行2列，这样两个图会并排排列
par(mfrow = c(1, 2)) 

# b. 绘制第一个图：通讯数量
p001 = netVisual_circle(cellchat@net$count, 
                 vertex.weight = as.numeric(table(cellchat@idents)), 
                 weight.scale = TRUE, 
                 label.edge = FALSE, 
                 title.name = "Number of interactions")

# c. 绘制第二个图：通讯强度
p002 = netVisual_circle(cellchat@net$weight, 
                 vertex.weight = as.numeric(table(cellchat@idents)), 
                 weight.scale = TRUE, 
                 label.edge = FALSE, 
                 title.name = "Interaction strength")
p0012 = p001|p002
print(p0012)
# d. 关闭图形设备，完成保存
dev.off()
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

# --- 3. 主要信号通路的Chord图（弦图） ---
# Chord图非常适合可视化特定信号通路中复杂的相互作用。
# 我们以COLLAGEN通路为例，因为它在你的数据中似乎很重要。
# 你可以查看 cellchat@netP$pathways 来获取所有可用的通路名称。
# 查看所有通路
# print(cellchat@netP$pathways) 
pathways_to_plot <- c("COLLAGEN", "FN1", "LAMININ") # 选择你感兴趣的通路
for (pathway in pathways_to_plot) {
  if (pathway %in% cellchat@netP$pathways) {
    p_chord <- netVisual_aggregate(cellchat, 
                                   signaling = pathway, 
                                   layout = "chord")
    
    ggsave(filename = file.path(projectPath,'./Output/Integration20250829/', paste0("CellChat_Chord_", pathway, ".pdf")), 
           plot = p_chord, width = 8, height = 8)
    
    print(paste("3. ", pathway, "通路的Chord图已保存。"))
  }
}

# --- 4. 信号通路的角色模式识别 ---
# CellChat可以将信号通路分为几类（如INCOMING, OUTGOING, BOTH），并识别每个细胞亚群在其中的角色。
# 这对于理解细胞的功能定位非常有帮助。
# a. 计算并可视化传出（outgoing）模式
# selectK(cellchat, pattern = "outgoing") # 运行这行来帮助确定最佳的模式数量
nPatterns_outgoing <- 3 # 假设我们选择3种模式
cellchat <- identifyCommunicationPatterns(cellchat, 
                                          pattern = "outgoing", 
                                          k = nPatterns_outgoing)
p_pattern_out <- netAnalysis_river(cellchat, pattern = "outgoing")
ggsave(filename = file.path(projectPath, "CellChat_Pattern_Outgoing.pdf"), 
       plot = p_pattern_out, width = 10, height = 8)

# b. 计算并可视化传入（incoming）模式
# selectK(cellchat, pattern = "incoming") # 运行这行来帮助确定最佳的模式数量
nPatterns_incoming <- 4 # 假设我们选择4种模式
cellchat <- identifyCommunicationPatterns(cellchat, 
                                          pattern = "incoming", 
                                          k = nPatterns_incoming)
p_pattern_in <- netAnalysis_river(cellchat, pattern = "incoming")

ggsave(filename = file.path(projectPath, "CellChat_Pattern_Incoming.pdf"), 
       plot = p_pattern_in, width = 10, height = 8)

print("4. 信号通路角色模式图已保存。")


# --- 5. 信号网络中心性分析（可选，但非常有用）---
# 这可以帮助识别在整个通讯网络中最重要的信号通路。
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
p_centrality <- netAnalysis_signalingRole_network(cellchat, 
                                                  signaling = pathways_to_plot, 
                                                  width = 12, height = 5, 
                                                  font.size = 10)

ggsave(filename = file.path(projectPath, "CellChat_Signaling_Centrality.pdf"), 
       plot = p_centrality, width = 10, height = 5)

print("5. 信号网络中心性分析图已保存。")

print("所有额外的 CellChat 图表已生成并保存！")
