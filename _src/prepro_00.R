## example run:
## Rscript prepro_00.R --run_mode anchoring --sample_names epi_p007t epi_p008t epi_p009t1 epi_p009t2 epi_p012t epi_p013t epi_p007n epi_p008n epi_p009n1 epi_p009n2 epi_p012n epi_p013n --ref_samples epi_p007t epi_p008t epi_p009t1 epi_p009t2 epi_p012t epi_p013t
## Rscript prepro_00.R --run_mode annotation --sample_names epi_p007t epi_p008t epi_p009t1 epi_p009t2 epi_p012t epi_p013t epi_p007n epi_p008n epi_p009n1 epi_p009n2 epi_p012n epi_p013n
## Rscript ./prepro_00.R --run_mode anchoring --anchoring FALSE --sample_names epi_p007t epi_p008t epi_p009t1 epi_p009t2 epi_p012t epi_p013t epi_p014t epi_p016t epi_p017t epi_p007n epi_p008n epi_p009n1 epi_p009n2 epi_p012n epi_p013n epi_p014n epi_p016n epi_p017n imm_p007t imm_p008t imm_p009t1 imm_p009t2 imm_p012t imm_p013t imm_p014t imm_p016t imm_p017t imm_p007n imm_p008n imm_p009n1 imm_p009n2 imm_p012n imm_p013n imm_p014n imm_p016n imm_p017n str_p007t str_p008t str_p009t1 str_p009t2 str_p012t str_p013t str_p014t str_p016t str_p017t str_p007n str_p009n1 str_p009n2 str_p012n str_p013n str_p014n str_p016n str_p017n p009ot0 p009ot2 p013ot0 p013ot3 p014ot0

library(tidyverse)


## parameter parsing ----------------------------------------

cargs <- list()
cargs$run_mode <- NULL
cargs$sample_names <- NULL
cargs$ref_samples <- NULL
cargs$seu_obj <- NULL
cargs$rel_comp <- 10
cargs$max_cells <- 2000
cargs$nFeature_lower <- 500
cargs$nFeature_upper <- 5000
cargs$nCount_lower <- 1000
cargs$nCount_upper <- 50000
cargs$nFrac_s_lower <- 0.3
cargs$nFrac_s_upper <- 0.9
cargs$nFrac_u_lower <- 0.1
cargs$nFrac_u_upper <- 0.7
cargs$nFrac_a_lower <- 0
cargs$nFrac_a_upper <- 0.2
cargs$pMT_lower <- 0
cargs$pMT_upper <- 0.8
cargs$pHB_lower <- 0
cargs$pHB_upper <- 0.1
cargs$cluster_res <- 2
cargs$nbin <- 24
cargs$anchoring <- TRUE
cargs$demux <- FALSE
cargs$root_dir <- "/extra/flo/sc/sc18"
cargs$grouping_var <- "sample_origin"

cargs_user <- commandArgs(trailingOnly = T) %>% 
  paste0(collapse = " ") %>% 
  str_split("--") %>% 
  lapply(str_split, " ") %>% 
  unlist(recursive = F) %>% 
  lapply(function(x) set_names(list(x[-1]), x[1])) %>% 
  unlist(recursive = F) %>% 
  lapply(function(x) x[x!=""]) %>% 
  .[-1]

for (i in 1:length(cargs_user)) {
  cargs[[names(cargs_user)[i]]] <- cargs_user[[names(cargs_user)[i]]]
}

print(cargs)

knit_root_dir <- cargs$root_dir

setwd(knit_root_dir)

         
## functions ----------------------------------------------
         
run_qc <- function(sample_name = cargs$sample_names, 
                   loom_file = paste0("_data/velocyto/", sample_name, ".loom"), 
                   nFeature_lower = cargs$nFeature_lower,
                   nFeature_upper = cargs$nFeature_upper,
                   nCount_lower = cargs$nCount_lower,
                   nCount_upper = cargs$nCount_upper,
                   nFrac_s_lower = cargs$nFrac_s_lower,
                   nFrac_s_upper = cargs$nFrac_s_upper,
                   nFrac_u_lower = cargs$nFrac_u_lower,
                   nFrac_u_upper = cargs$nFrac_u_upper,
                   nFrac_a_lower = cargs$nFrac_a_lower,
                   nFrac_a_upper = cargs$nFrac_a_upper,
                   pMT_lower = cargs$pMT_lower,
                   pMT_upper = cargs$pMT_upper,
                   pHB_lower = cargs$pHB_lower,
                   pHB_upper = cargs$pHB_upper,
                   rel_comp = cargs$rel_comp,
                   demux = cargs$demux) {
  rmarkdown::render(
    "_src/prepro_01_qc.Rmd", params = list(
      sample_name = sample_name,
      loom_file = loom_file,
      nFeature_lower = nFeature_lower,
      nFeature_upper = nFeature_upper,
      nCount_lower = nCount_lower,
      nCount_upper = nCount_upper,
      nFrac_s_lower = nFrac_s_lower,
      nFrac_s_upper = nFrac_s_upper,
      nFrac_u_lower = nFrac_u_lower,
      nFrac_u_upper = nFrac_u_upper,
      nFrac_a_lower = nFrac_a_lower,
      nFrac_a_upper = nFrac_a_upper,
      pMT_lower = pMT_lower,
      pMT_upper = pMT_upper,
      pHB_lower = pHB_lower,
      pHB_upper = pHB_upper,
      rel_comp = rel_comp,
      demux = demux
    ),
    output_file = paste0("01_qc_", sample_name, "_strict.html"),
    output_dir = "_html",
    knit_root_dir = knit_root_dir
  )
}

run_sub <- function(seu_obj = cargs$seu_obj,
                    max_cells = cargs$max_cells) {
  rmarkdown::render(
    "_src/prepro_02_sub.Rmd", params = list(
      seu_obj = seu_obj,
      max_cells = max_cells
    ),
    output_file = paste0("02_sub_", substr(basename(cargs$seu_obj), 1, 100), ".html"),
    output_dir = "_html",
    knit_root_dir = knit_root_dir
  )
}

run_anchoring <- function(sample_names = cargs$sample_names,
                          ref_samples = cargs$ref_samples,
                          rel_comp = cargs$rel_comp,
                          anchoring = cargs$anchoring) {
  rmarkdown::render(
    "_src/prepro_03_anchoring.Rmd", params = list(
      sample_names = sample_names,
      ref_samples = ref_samples,
      rel_comp = rel_comp,
      anchoring = anchoring
    ),
    output_file = paste0("03_anchoring_", substr(paste0(sample_names, collapse = "_"), 1, 100), ".html"),
    output_dir = "_html",
    knit_root_dir = knit_root_dir
  )
}

run_annotation <- function(seu_obj = cargs$seu_obj,
                           cluster_res = cargs$cluster_res,
                           nbin = cargs$nbin,
                           grouping_var = cargs$grouping_var) {
  rmarkdown::render(
    "_src/prepro_04_annotation.Rmd", params = list(
      seu_obj = seu_obj,
      cluster_res = cluster_res,
      nbin = nbin,
      grouping_var = grouping_var
    ),
    output_file = paste0("04_anno_", substr(basename(cargs$seu_obj), 1, 100), ".html"),
    output_dir = "_html",
    knit_root_dir = knit_root_dir
  )
}

if (cargs$run_mode == "qc") { run_qc()
} else if (cargs$run_mode == "subset") { run_sub()                         
} else if (cargs$run_mode == "anchoring") { run_anchoring()
} else if (cargs$run_mode == "annotation") { run_annotation()
} else {stop("invalid run mode")}

