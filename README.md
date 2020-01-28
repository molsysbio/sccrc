# Analysis pipeline for CRC patient samples 

The preprocessing pipeline for scRNA-seq has four modules consisting of parameterized `Rmarkdown` scripts: 

1. QC+filtering: `prepro_01_qc.Rmd`
    * takes velocyto.loom file and creates seu_raw object (assumes path: "_data/_patients/velocyto/p007n.loom")
    * needs filtering parameters for qc (e.g. mito reads)
2. Subsampling: `prepro_02_sub.Rmd`
    * takes a Seurat object as input file and outputs a Seurat object as output with N cells
2. Anchoring: `prepro_03_anchoring.Rmd`
    * takes multiple qc-filtered Seurat objects and merges (and anchors) them into one merged (and one anchored) Seurat object
3. Annotation: `prepro_04_ann.Rmd`
    * takes a seurat object and annotates it with preliminary cell type information and gene expression signatures

Wrapper functions for module calls are described in `prepro_00.R`. All calls are listed in `prepro.sh`. R library dependencies are listed within each module file.

After preprocessing, figures were generated as described in `scCRC_paper_figures.Rmd`.
