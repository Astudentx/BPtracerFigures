
package <- c('readxl','magrittr',
             'ggtree','reshape2','ggplot2','dplyr','tidyr',
             'readr','purrr','tibble','stringr','forcats',
             'ggplot2','ggsci','ggridges','RColorBrewer','pheatmap',
             'patchwork','ggpubr','zyzPackage','cowplot','ggrepel','networkD3',"ggrepel","caret","vegan")
lapply(package, function(x){library(x, character.only = T)})


library(bptracer)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)


# Read and save phylogenetic tree
# tree_raw <- read.tree(file = "01.rawdata/Pangenomes_Tree/zyz-18640.tree.nwk")
# save(tree_raw,file = "02.datasave/03.Tree.Rdata")

load(file = "02.datasave/02.HiIndex.Rdata")
load(file = "02.datasave/00.RefseqDatabase.Rdata")
load(file = "01.rawdata/UniPathogenDB/UniPathogenDB.Rdata")
load(file = "02.datasave/03.Tree.Rdata")
attach(NewPA_Mini)
attach(final_bio_risk)


tree <- tree_raw  # Pipeline: load tree and metadata, group by phylum, find clade roots, highlight, annotate, export
tree$tip.label %>% head()
tree$edge.length %>% head()
tree$node.label %>% head()

# Phylogenetic analysis
## Genome info per phylogenetic node  


# Companion data processing
df1 <- final_bio_risk %>% select(GCFID,Species) %>% unique()  # GCF to species

# Get Refseq metadata
df2 <- refseq_240519_usedGCFtax %>% select(assembly_accession,Species.ID,5:11) # Taxonomy and pathogen flag
colnames(df2)[1] <- c("GCFID") 
attach(NewPA_Mini)
NewPA_Mini$ID <- as.numeric(NewPA_Mini$ID)
df2$Species.ID <- as.numeric(df2$Species.ID)
df2_PA <- left_join(df2,NewPA_Mini,by=c("Species.ID"="ID")) %>% select(1:9,13)
uniq_c(df2_PA$GCFID) %>% head()
df2_PA$Type[is.na(df2_PA$Type)] <- "NonPathogens"


# Species-by-gene counts
df3 <-  final_bio_risk %>% select(Species,GeneKind,Type) %>% uniq_c_df()
df4 <- dcast(
  data = df3,
  formula = Species+Type ~ GeneKind,
  value.var = "Count",
  fill = 0
) %>% as.data.frame()
df4 <- df4 %>% select("Species","ARGs","MGEs","MRGs","SGs","VFs")


GCF_SpeciesGene <- left_join(df2_PA,df4,by = join_by(Species))
GCF_SpeciesGene[is.na(GCF_SpeciesGene)] <- 0
colnames(GCF_SpeciesGene)[1] <- "label"
rownames(GCF_SpeciesGene) <- GCF_SpeciesGene$GCFID



GCF_SpeciesGene <-  profile_top_col(inputtable = GCF_SpeciesGene,n = 9,measure.vars = "Phylum",replaceName = "Others")
GCF_SpeciesGene$group <- GCF_SpeciesGene$Phylum
GCF_SpeciesGene_melt <- melt(GCF_SpeciesGene,measure.vars = 11:15)




## Tree filtering
library(ape)
head(tree$tip.label)
tree$tip.label %>% length()
uniqID_sel <- tree$tip.label         # 最终使用版本
tree_bac <- keep.tip(tree, uniqID_sel)
tree_tips <- tibble(label = tree_bac$tip.label)
tree_sel1 <- left_join(tree_bac,GCF_SpeciesGene,by="label")
# tree_sel1 <- left_join(tree_bac,GCF_SpeciesGene_melt,by="label")
tree_sel1@extraInfo
tree_sel1@extraInfo %>% nrow()
tree_sel1@extraInfo$Type %>% unique()

temp <- tree_sel1@extraInfo
df1 <- as_tibble(tree_sel1) %>% select(Phylum,node) %>% na.omit()
cls <- split(df1$node, df1$Phylum)
uniq_c(df1$Phylum)

tree_sel2 <- groupOTU(tree_sel1, cls, group_name = "PhylumGroup")
tree_df2     <- as_tibble(tree_sel2)
# groupOTU sets mixed nodes to "0"


parent_group <- tree_df2 %>% select(node, group_parent = PhylumGroup)
uniq_c(parent_group$group_parent)


find_group_roots <- function(tree_df, group_name, group_col = "PhylumGroup") {
  # Nodes whose parent left the group mark clade roots
  idx <- which(tree_df[[group_col]] == group_name)
  
  if (length(idx) == 0) return(integer(0))
  
  df_g <- tree_df[idx, ]
  
  roots <- df_g$node[match(df_g$parent, tree_df$node) %>% { tree_df[[group_col]][.] != group_name }]
  
  roots <- roots[!is.na(roots)]
  unique(roots)
}


find_root_with_mrca <- function(tree_obj, tree_df, group_name) {
  # Prefer pure clade roots; fallback to MRCA when dispersed
  roots <- find_group_roots(tree_df, group_name)
  
  if (length(roots) > 0) {
    return(roots)
  }
  
  message(sprintf("No pure clade root found for %s; fallback to MRCA.", group_name))
  
  labels <- tree_df$label[tree_df$Phylum == group_name]
  labels <- labels[!is.na(labels)]
  
  if (length(labels) <= 1) return(integer(0))
  
  mrca <- getMRCA(tree_obj@phylo, labels)
  
  if (is.null(mrca)) return(integer(0))
  
  return(mrca)
}


Phylum_list <- setdiff(unique(tree_df2$Phylum), c("0", "Others", NA))
df_clade <- data.frame()
for (g in Phylum_list) {  # Locate roots per phylum for geom_hilight
  roots <- find_root_with_mrca(tree_sel2, tree_df2, g)
  if (length(roots) == 0) next
  df_clade <- rbind(df_clade,
                    data.frame(node = roots,
                               Phylum = g))
}

df_clade <- df_clade[!is.na(df_clade$node), ]



# Plotting
phylum_levels <- sort(unique(tree_sel2@extraInfo$group))
color_mapping_all <- color_scheme_bp("color_zyz6")[seq_along(phylum_levels)]
names(color_mapping_all) <- phylum_levels



color_mapping_line <- color_mapping_all
color_mapping_fill <- color_mapping_all[df_clade$Phylum]


phylum_levels_PA <- sort(unique(na.omit(tree_sel2@extraInfo$Type)))
color_mapping_PA <- c("#33395d","white","#c5d8f1","#8c2e2f")
names(color_mapping_PA) <- phylum_levels_PA


p4 <- ggtree(tree_sel2,size=0.2,aes(color=PhylumGroup), layout = "fan")+
  scale_color_manual(values =color_mapping_line)+
  new_scale_fill()+
  geom_hilight(data = df_clade,aes(node=node,fill=Phylum) ,show.legend = F,alpha=0.15,size=0.35) +
  scale_fill_manual(values = color_mapping_fill)+
  theme(legend.title=element_blank())
p4



p4_1 <- ggtree(tree_sel2,size=0.2,aes(color=PhylumGroup), layout = "slanted")+
  scale_color_manual(values =color_mapping_line)+
  new_scale_fill()+
  geom_hilight(data = df_clade,aes(node=node,fill=Phylum) ,show.legend = F,alpha=0.3,size=0.35) +
  scale_fill_manual(values = color_mapping_fill)+
  theme(legend.title=element_blank())
p4_1


p5 <- p4_1 +
  new_scale_color() +
  geom_tippoint(
    data    = function(x) x[!is.na(x$Type), ],
    aes(color = Type),
    size   = .2,
    alpha  = .6
  ) +
  scale_color_manual(values = color_mapping_PA)
p5



p6 <- p5+ new_scale_fill() +
  geom_fruit(
    data    = function(x) x[!is.na(x$Type), ],
    geom    = geom_tile,
    mapping = aes(fill = Type),
    pwidth  = 0.4,
    offset  = 0.02
  ) +
  scale_fill_manual(
    values = color_mapping_PA,
    name   = "Type"
  )
p6


p7 <- p6 + 
  new_scale_color() +
  geom_fruit(geom = geom_col,color="#d74a49",mapping = aes(x = ARGs),pwidth = 0.1,offset = 0.1)+
  geom_fruit(geom = geom_col,color="#5e62a9",mapping = aes(x = MGEs),pwidth = 0.1)+
  geom_fruit(geom = geom_col,color="#1b4552",mapping = aes(x = MRGs),pwidth = 0.1)+
  geom_fruit(geom = geom_col,color="#8ba0a4",mapping = aes(x = VFs),pwidth = 0.1)+
  geom_fruit(geom = geom_col,color="black",mapping = aes(x = SGs),pwidth = 0.1)+
  scale_color_manual(values = c("ARGs" = "#d74a49", "MGEs" = "#5e62a9", "MRGs" = "#1b4552","VFs" = "#8ba0a4", "SGs" = "black"),
                     name = "Gene Categories")  # 设置图例标题
p7


# 从 p6 中取出 tip 的数据（带 label 和 5 个基因列）
df_gene <- p6$data %>%
  dplyr::filter(isTip) %>%
  dplyr::select(label, ARGs, MGEs, MRGs, VFs, SGs) %>%
  tidyr::pivot_longer(
    cols      = c(ARGs, MGEs, MRGs, VFs, SGs),
    names_to  = "GeneCategory",
    values_to = "value"
  )



p8 <- p7 +
  new_scale_color() +   # 开一套新的 fill 映射，避免跟前面的 Type 冲突
  geom_fruit(
    data    = df_gene,
    geom    = geom_col,
    mapping = aes(
      y    = label,          # 按 tip 对齐
      x    = value,          # 柱子的长度
      color = GeneCategory    # 用填充色区分基因类别
    ),
    pwidth      = 0.5,       # 整个堆叠柱的宽度
    offset      = 0.1,      # 跟 Type 色带之间的距离
    # position    = "stack"
    # orientation = "y"        # 建议显式指定
  ) +
  scale_color_manual(
    values = c(
      "ARGs" = "#d74a49",
      "MGEs" = "#5e62a9",
      "MRGs" = "#1b4552",
      "VFs"  = "#8ba0a4",
      "SGs"  = "black"
    ),
    name = "Gene Categories"
  )
p8

# save(tree_sel1,tree_sel2,tree_df,tree_df2,file = "02.datasave/01.Tree_bak.Rdata")
# load(file = "02.datasave/01.Tree_bak.Rdata")
ggsave(filename = "03.figureUP/M3_Tree.pdf",width = 12,height = 8)
# ggsave(filename = "03.figureUP/M3_Tree2.pdf",width = 14,height = 10)
