library(tidyverse)
library(Seurat)
library(infercnv)
library(cowplot)
library(ggdendro)
library(dendextend)

## load data ------------------------------------------

## load preprocessed seurat object  
seu_epi_pri <- read_rds("_data/computed/final/seu_epi_final.rds")

## global vars
cell_type_colors <- setNames(unique(seu_epi_pri$cell_type_epi_color), unique(seu_epi_pri$cell_type_epi_custom))

cols.use <- c(
  p001 = '#E69F00', #palette_OkabeIto[1],
  p007 = '#E69F00', #palette_OkabeIto[1],
  p008 = '#56B4E9', #palette_OkabeIto[2],
  p009 = '#009E73', #palette_OkabeIto[3],
  p012 = "#ff0198",
  p013 = '#0072B2', #palette_OkabeIto[5],
  p014 = '#D55E00', #palette_OkabeIto[6],
  p016 = '#CC79A7', #palette_OkabeIto[7],
  p017 = "#451077",
  p020 = "#cccc80",
  p021 = "#FF0000", 
  p025 = "#f6cbcc",
  p026 = "#F0E442",
  CMS1 = "#eba83a", CMS2 = "#027eb5", CMS3 = "#d684ae", CMS4 = "#00a881",
  `CMS1,CMS2` = "#779378", `CMS1,CMS3` = "#E19674", `CMS1,CMS4` = "#76A85E", 
  `CMS2,CMS3` = "#6C81B2", `CMS2,CMS4` = "#01939B",  
  G1 = "#94b6d2",
  S = "#dc4040",
  G2M = "#7aa77f",
  `NA` = "grey",
  Normal = "steelblue",
  Tumor = "red",
  PDO = "steelblue",
  CNA = "red",
  CNN = "grey"
)

pids <- sort(unique(seu_epi_pri$source_id))

## prepare infercnv inputs ---------------------------------

## 1. count matrix: cells x genes
rcm <- seu_epi_pri@assays$RNA@counts
cell_subset <- seu_epi_pri$cell_id[!(str_detect(seu_epi_pri$cell_type_epi_custom, "^TC") & seu_epi_pri$sample_origin == "Normal")]
rcm <- rcm[,cell_subset]
colnames(rcm) <- str_replace_all(colnames(rcm), ":", "_")

## 2. two column annotation table: cell_id, cell_type 
af <- FetchData(seu_epi_pri, c("cell_type_epi_simple", "sample_origin", "source_id", "cell_id", "sample_id")) %>%
  filter(cell_id %in% cell_subset) %>%
  mutate(cell_type_infercnv = ifelse(sample_origin == "Tumor", paste0("malignant_", source_id), as.character(cell_type_epi_simple)),
         cell_id = str_replace_all(cell_id, ":", "_")) %>% 
  select(cell_id, cell_type_infercnv) %>% 
  as_tibble
write.table(af, file = "_data/infercnv/annotation.txt", sep = "\t", col.names = F, row.names = F)

## 3. four column gene location table: gene, chr, start, end
## used the one from infercnv website
gof <- read_tsv("_data/infercnv/gencode_v21_gen_pos.complete.txt", col_names = c("gene", "chr", "start", "end"))
common_genes <- intersect(gof$gene, rownames(rcm))

## 4. defining cell types for baseline reference
baseline_cell_types <- sort(unique(af$cell_type_infercnv)[!str_detect(unique(af$cell_type_infercnv), "malignant")])

## run infer cnv --------------------------------------------
# 
# ## run time: ~12h
# 
# ## create infercnv object
# cnv_obj <- infercnv::CreateInfercnvObject(
#   raw_counts_matrix = rcm[common_genes,],
#   annotations_file = "_data/infercnv/annotation.txt",
#   gene_order_file = "_data/infercnv/gencode_v21_gen_pos.complete.txt",
#   ref_group_names = baseline_cell_types
# ) 
# 
# ## run infercnv
# cnv_obj <- infercnv::run(
#   cnv_obj, cutoff=0.1, out_dir="_data/infercnv/out", 
#   cluster_by_groups=T, denoise=F, HMM=F, 
#   num_threads=parallel::detectCores()
# )
# write_rds(cnv_obj, "_data/infercnv/infercnv_results.rds")

cnv_obj <- read_rds("_data/infercnv/infercnv_results.rds")


## custom plotting --------------------------------------------

gof <- filter(gof, chr %in% paste0("chr", 1:22),
              gene %in% common_genes) %>% 
  mutate(chr = ordered(paste0(" ", chr), levels = paste0(" chr", 1:22)))

## infercnv expression results
expr_all <- cnv_obj@expr.data
expr_x <- cnv_obj@expr.data[,str_detect(colnames(cnv_obj@expr.data), "^p0[0-9][0-9]t")]

## infercnv dendrogram
dendro_lines <- read_lines("_data/infercnv/out/infercnv.observations_dendrogram.txt")
dendro_list <- lapply(dendro_lines, function(x) phylogram::read.dendrogram(text = x)) %>% setNames(pids)
cell_order <- lapply(dendro_lines, function(x) str_split(x, "\\(|:|,|\\)") %>% unlist %>% .[str_detect(., "p0")]) %>% unlist

## find clones
clone_tbl <- lapply(dendro_list, function(x) enframe(cutree(x, k = 2), "cell_id", "clone")) %>% 
  bind_rows(.id = "source_id") %>% 
  mutate(clone = as.character(clone)) %>% 
  bind_rows(
    FetchData(seu_epi_pri, c("source_id", "cell_id")) %>%
      mutate(cell_id = str_replace_all(cell_id, ":", "_"),
             clone = "Normal") %>%
      filter(!str_detect(seu_epi_pri$cell_type_epi_custom, "^TC") & seu_epi_pri$sample_origin == "Normal")
  )

clone_scores <- clone_tbl %>% 
  group_by(source_id, clone) %>%
  do(cna_score = mean(colSums(abs(expr_all[,.$cell_id]))/nrow(expr_all))) %>%
  ungroup() %>%
  unnest(cna_score) %>%
  spread(clone, cna_score) %>%
  mutate(mean_n = mean(Normal, na.rm = T)) %>%
  gather(clone, cna_score, -mean_n, -source_id) %>%
  mutate(cna_score = cna_score/mean_n) %>%
  mutate(max_n = max(cna_score[.$clone == "Normal"], na.rm = T)) %>%
  mutate(cna_clone = ifelse(cna_score > max_n, "CNA", "CNN")) %>%
  select(source_id, cna_score, cna_clone, clone) %>%
  #filter(clone != "Normal") %>% 
  left_join(clone_tbl, by = c("clone", "source_id")) %>%
  mutate(cell_id = ordered(cell_id, levels = cell_order)) %>%
  mutate(helper_var_c = " CN")

write_tsv(select(clone_scores, -helper_var_c), "_data/_tab/infercnv_clone_scores.tsv")
                    
segments_tbl <- lapply(dendro_list, function(x) segment(dendro_data(x)) %>% mutate(helper_var_d = "Dendrogram")) %>% 
  bind_rows(.id = "source_id") 

## cell annotation
anno_data <- as_tibble(FetchData(seu_epi_pri, c("cell_type_epi_custom", "source_id", "sample_origin", "cell_id"))) %>%
  mutate(cell_id = str_replace_all(cell_id, ":", "_"),
         cell_id = ordered(cell_id, levels = cell_order)) %>%
  filter(sample_origin == "Tumor") %>%
  mutate(helper_var_p = " Patient",
         helper_var_c = " Cell type")

## joining data 
plot_data <- as.data.frame(expr_x) %>% 
  as_tibble(rownames = "gene") %>%
  gather(cell_id, value, -gene) %>%
  left_join(gof, by = "gene") %>% 
  mutate(cell_id = ordered(cell_id, levels = cell_order)) %>%
  #filter(chr == 13) %>% 
  left_join(anno_data, by = "cell_id") %>% 
  mutate(cell_id = ordered(cell_id, levels = cell_order)) %>% 
  mutate(value_cut = ifelse(value > 1.15, 1.15, ifelse(value < 0.85, 0.85, value))) %>% 
  group_by(chr) %>% 
  mutate(rank = rank(start)) %>%
  ungroup

## sub plot objects -----------------------------------

plot_cnv <- ggplot(plot_data) +
  geom_tile(aes(rank, cell_id, fill=value_cut)) +
  theme_void() +
  facet_grid(source_id~chr, scales = "free", space = "free") +
  scale_fill_gradient2(midpoint = 1, low = scales::muted("blue"), high = scales::muted("red")) +
  theme(panel.spacing.x = unit(0, "npc"),
        panel.spacing.y = unit(0.003, "npc"),
        strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        strip.text.y = element_blank()) +
  labs(fill = "inferCNV\nexpression")

plot_anno_cell <- ggplot(anno_data) +
  geom_tile(aes(helper_var_c, cell_id, fill = cell_type_epi_custom)) +
  theme_void() +
  facet_grid(source_id~helper_var_c, scales = "free", space = "free") +
  scale_fill_manual(values = cell_type_colors) +
  guides(fill = F) +
  theme(panel.spacing.x = unit(0, "npc"),
        panel.spacing.y = unit(0.003, "npc"),
        strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        strip.text.y = element_blank()) 

plot_anno_patient <- ggplot(anno_data) +
  geom_tile(aes(helper_var_p, cell_id, fill = source_id)) +
  theme_void() +
  facet_grid(source_id~helper_var_p, scales = "free", space = "free") +
  scale_fill_manual(values = cols.use) +
  guides(fill = F) +
  theme(panel.spacing.x = unit(0, "npc"),
        panel.spacing.y = unit(0.003, "npc"),
        strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        strip.text.y = element_blank())

plot_anno_clone <- ggplot(filter(clone_scores, clone != "Normal")) +
  geom_tile(aes(helper_var_c, cell_id, fill = cna_clone)) +
  theme_void() +
  facet_grid(source_id~helper_var_c, scales = "free", space = "free") +
  scale_fill_manual(values = cols.use) +
  guides(fill = F) +
  theme(panel.spacing.x = unit(0, "npc"),
        panel.spacing.y = unit(0.003, "npc"),
        strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        strip.text.y = element_blank())
                       
plot_dendro <- ggplot(segments_tbl) +
  geom_segment(aes(x=-y,y=x,xend=-yend,yend=xend),size=0.5) +
  theme_void() +
  facet_grid(source_id~helper_var_d, scales = "free", space = "free") +
  theme(panel.spacing = unit(0, "npc"),
        panel.spacing.x = unit(0, "npc"),
        panel.spacing.y = unit(0.003, "npc"),
        strip.text = element_blank()) +
  scale_y_continuous(expand=c(0,0))

pg_anno <- plot_grid(plot_dendro, plot_anno_patient, plot_anno_clone, plot_anno_cell, nrow = 1, align = "h", rel_widths = c(0.625, 0.125, 0.125, 0.125))

ggsave("_fig/infercnv_heat_anno.pdf", pg_anno, width = 2.4, height = 20)
                       
pg <- plot_grid(plot_dendro, plot_anno_patient, plot_anno_clone, plot_anno_cell, plot_cnv, nrow = 1, align = "h", rel_widths = c(0.1, 0.02, 0.02, 0.02, 0.84))

ggsave("_fig/infercnv_heat.png", pg, width = 15, height = 20)

                       
## reference cell plot --------------------------------

expr_n <- cnv_obj@expr.data[,str_detect(colnames(cnv_obj@expr.data), "^p0[0-9][0-9]n")]

## cell annotation
anno_data_ref <- as_tibble(FetchData(seu_epi_pri, c("cell_type_epi_custom", "source_id", "sample_origin", "cell_id"))) %>%
  mutate(cell_id = str_replace_all(cell_id, ":", "_")) %>%
  filter(sample_origin == "Normal", !str_detect(cell_type_epi_custom, "^TC")) %>%
  arrange(cell_type_epi_custom, desc(source_id)) %>% 
  mutate(helper_var_p = " Patient",
         helper_var_c = " Cell type",
         cell_id = ordered(cell_id, levels = .$cell_id))

cell_order_ref <- anno_data_ref$cell_id

## joining data 
plot_data_ref <- as.data.frame(expr_n) %>% 
  as_tibble(rownames = "gene") %>%
  gather(cell_id, value, -gene) %>%
  left_join(gof, by = "gene") %>% 
  mutate(cell_id = ordered(cell_id, levels = cell_order_ref)) %>% 
  left_join(anno_data_ref, by = "cell_id") %>% 
  mutate(value_cut = ifelse(value > 1.15, 1.15, ifelse(value < 0.85, 0.85, value))) %>% 
  group_by(chr) %>% 
  mutate(rank = rank(start)) %>%
  ungroup

## reference sub plot objects -------------------------------

plot_cnv_ref <- ggplot(plot_data_ref) +
  geom_tile(aes(rank, cell_id, fill=value_cut)) +
  theme_void() +
  facet_grid(cell_type_epi_custom~chr, scales = "free", space = "free") +
  scale_fill_gradient2(midpoint = 1, low = scales::muted("blue"), high = scales::muted("red")) +
  theme(panel.spacing.x = unit(0, "npc"),
        panel.spacing.y = unit(0.003, "npc"),
        strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        strip.text.y = element_blank()) +
  labs(fill = "inferCNV\nexpression")

plot_anno_cell_ref <- ggplot(anno_data_ref) +
  geom_tile(aes(helper_var_c, cell_id, fill = cell_type_epi_custom)) +
  theme_void() +
  facet_grid(cell_type_epi_custom~helper_var_c, scales = "free", space = "free") +
  scale_fill_manual(values = cell_type_colors) +
  guides(fill = F) +
  theme(panel.spacing.x = unit(0, "npc"),
        panel.spacing.y = unit(0.003, "npc"),
        strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        strip.text.y = element_blank()) 

plot_anno_patient_ref <- ggplot(anno_data_ref) +
  geom_tile(aes(helper_var_p, cell_id, fill = source_id)) +
  theme_void() +
  facet_grid(cell_type_epi_custom~helper_var_p, scales = "free", space = "free") +
  scale_fill_manual(values = cols.use) +
  guides(fill = F) +
  theme(panel.spacing.x = unit(0, "npc"),
        panel.spacing.y = unit(0.003, "npc"),
        strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        strip.text.y = element_blank())

plot_anno_empty <- ggplot(anno_data_ref) + 
  geom_tile(aes(helper_var_p, cell_id), fill = "white") +
  theme_void() +
  facet_grid(cell_type_epi_custom~helper_var_p, scales = "free", space = "free") +
  scale_fill_manual(values = cols.use) +
  guides(fill = F) +
  theme(panel.spacing.x = unit(0, "npc"),
        panel.spacing.y = unit(0.003, "npc"),
        strip.text = element_blank())

pg_anno_ref <- plot_grid(plot_anno_patient_ref, plot_anno_cell_ref, nrow = 1, align = "h", rel_widths = c(1/2, 1/2))

ggsave("_fig/infercnv_heat_anno_ref.pdf", pg_anno_ref, width = 0.6, height = 10)
                       
pg_ref <- plot_grid(plot_anno_empty, plot_anno_patient_ref, plot_anno_cell_ref, plot_cnv_ref, nrow = 1, align = "h", rel_widths = c(0.12, 0.02, 0.02, 0.84))

ggsave("_fig/infercnv_heat_ref.png", pg_ref, width = 15, height = 10)



## inferCNV UMAP ------------------------

umap_data <- as_tibble(FetchData(seu_epi_pri, c("UMAP_1", "UMAP_2", "source_id", "cell_id", "sample_origin", "cell_type_epi_custom"))) %>%
  filter(sample_origin == "Tumor") %>% 
  mutate(cell_id = str_replace_all(cell_id, ":", "_")) %>% 
  left_join(clone_scores, by = c("cell_id", "source_id"))

ggplot() + 
  geom_point(aes(UMAP_1, UMAP_2, color = cna_clone), 
             size = 0.01, data = umap_data) +
  geom_point(aes(UMAP_1, UMAP_2), size = 0.1, color = "red", 
             data = filter(umap_data, cna_clone == "CNA")) +
  scale_color_manual(values = cols.use) +
  coord_fixed() + 
  theme_void()

ggsave("_fig/infercnv_umap.pdf", width = 4, height = 3)

ggplot() + 
  geom_point(aes(UMAP_1, UMAP_2, color = cna_clone), 
             size = 0.01, data = umap_data) +
  geom_point(aes(UMAP_1, UMAP_2), size = 0.1, color = "red", 
             data = filter(umap_data, cna_clone == "CNA")) +
  scale_color_manual(values = cols.use) +
  coord_fixed() + 
  facet_wrap(~source_id, ncol = 4) +
  theme_void()

ggsave("_fig/infercnv_umap_per_patient.pdf", width = 10, height = 6)

comp_tbl <- umap_data %>% 
  group_by(cna_clone, cell_type_epi_custom) %>% 
  tally() %>%
  group_by(cna_clone) %>%
  mutate(nrel = n/sum(n)*100)

comp_tbl_purity <- umap_data %>% 
  group_by(cna_clone, source_id) %>% 
  tally() %>%
  group_by(source_id) %>%
  mutate(nrel = n/sum(n)*100)%>%
  ungroup %>%
  arrange(desc(cna_clone), nrel) %>%
  mutate(source_id = ordered(source_id, levels = unique(source_id)))

comp_tbl_patient <- umap_data %>% 
  group_by(source_id, cna_clone, cell_type_epi_custom) %>% 
  tally() %>%
  group_by(source_id, cna_clone) %>%
  mutate(nrel = n/sum(n)*100)

ggplot(comp_tbl, aes(cna_clone, nrel, fill = fct_rev(cell_type_epi_custom))) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  scale_fill_manual(values = cell_type_colors) +
  theme_cowplot() +
  theme(legend.position = "bottom", 
        panel.border = element_rect(linetype = 1, color = "black", size = 1),
        axis.line = element_blank(),
        strip.text = element_blank()) +
  facet_grid(cna_clone~., scales = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(y = "Fraction of cells [%]", x = "", fill = "") 

ggsave("_fig/infercnv_composition.pdf", width = 6, height = 4)
  
ggplot(comp_tbl_patient, aes(cna_clone, nrel, fill = fct_rev(cell_type_epi_custom))) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  scale_fill_manual(values = cell_type_colors) +
  theme_cowplot() +
  theme(strip.background = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5), 
        legend.position = "bottom", 
        panel.border = element_rect(linetype = 1, color = "black", size = 1),
        axis.line = element_blank()) +
  facet_grid(source_id~., scales = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(y = "Fraction of cells [%]", x = "", fill = "") 

ggsave("_fig/infercnv_composition_per_patient.pdf", width = 4, height = 8)
       

ggplot(comp_tbl_purity, aes(source_id, nrel, fill = fct_rev(cna_clone))) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  scale_fill_manual(values = cols.use) +
  theme_cowplot() +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(), 
        legend.position = "bottom", 
        panel.border = element_rect(linetype = 1, color = "black", size = 1),
        axis.line = element_blank()) +
  facet_grid(source_id~., scales = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(y = "Fraction of cells [%]", x = "", fill = "") 

ggsave("_fig/infercnv_tumor_purity.pdf", width = 4, height = 5)

write_tsv(comp_tbl, "_data/_tab/infercnv_composition.tsv")
write_tsv(comp_tbl_patient, "_data/_tab/infercnv_composition_per_patient.tsv")
write_tsv(comp_tbl_purity, "_data/_tab/infercnv_composition_purity.tsv")
                       

