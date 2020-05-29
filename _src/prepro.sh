#!/bin/bash
ROOT_DIR=$1

## QCs
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p007n
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p007t
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p008n
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p008t
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p009n1
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p009n2
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p009t1
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p009t2
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p012n
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p012t
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p013n
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p013t
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p014n
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p014t
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p016n
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p016t
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p017n
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p017t
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p009ot0
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p009ot2
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p013ot0
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p013ot3
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p020n --demux TRUE
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p020t --demux TRUE
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p021n --demux TRUE
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p021t --demux TRUE
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p025n --demux TRUE
Rscript prepro_00.R --run_mode qc --root_dir $ROOT_DIR --sample_names p025t --demux TRUE

## merging and anchoring 
Rscript prepro_00.R --run_mode anchoring --anchoring FALSE --sample_names epi_p007n_strict epi_p007t_strict epi_p008n_strict epi_p008t_strict epi_p009n1_strict epi_p009n2_strict epi_p009t1_strict epi_p009t2_strict epi_p012n_strict epi_p012t_strict epi_p013n_strict epi_p013t_strict epi_p014n_strict epi_p014t_strict epi_p016n_strict epi_p016t_strict epi_p017n_strict epi_p017t_strict epi_p020n_strict epi_p020t_strict epi_p021n_strict epi_p021t_strict epi_p025n_strict epi_p025t_strict epi_p026t_strict --root_dir $ROOT_DIR --filter_warm TRUE

Rscript prepro_00.R --run_mode anchoring --anchoring FALSE --sample_names imm_p007n_strict imm_p007t_strict imm_p008n_strict imm_p008t_strict imm_p009n1_strict imm_p009n2_strict imm_p009t1_strict imm_p009t2_strict imm_p012n_strict imm_p012t_strict imm_p013n_strict imm_p013t_strict imm_p014n_strict imm_p014t_strict imm_p016n_strict imm_p016t_strict imm_p017n_strict imm_p017t_strict imm_p020n_strict imm_p020t_strict imm_p021n_strict imm_p021t_strict imm_p025n_strict imm_p025t_strict imm_p026t_strict --root_dir $ROOT_DIR --filter_warm TRUE

Rscript prepro_00.R --run_mode anchoring --anchoring FALSE --sample_names str_p007n_strict str_p007t_strict str_p008t_strict str_p009n1_strict str_p009n2_strict str_p009t1_strict str_p009t2_strict str_p012n_strict str_p012t_strict str_p013n_strict str_p013t_strict str_p014n_strict str_p014t_strict str_p016n_strict str_p016t_strict str_p017n_strict str_p017t_strict str_p021t_strict str_p025n_strict str_p025t_strict str_p026t_strict --root_dir $ROOT_DIR --filter_warm TRUE

Rscript prepro_00.R --run_mode annotation --seu_obj _data/computed/anchored/seu_merge_epi_p007n_strict_epi_p007t_strict_epi_p008n_strict_epi_p008t_strict_epi_p009n1_strict_epi_p009n2_strict_epi_p009t1_strict_epi_p009t2_strict_epi_p012n_strict_epi_p012t_strict_epi_p013n_strict_epi_p013t.rds --root_dir $ROOT_DIR --cluster_res 0.8

Rscript prepro_00.R --run_mode annotation --seu_obj _data/computed/anchored/seu_merge_imm_p007n_strict_imm_p007t_strict_imm_p008n_strict_imm_p008t_strict_imm_p009n1_strict_imm_p009n2_strict_imm_p009t1_strict_imm_p009t2_strict_imm_p012n_strict_imm_p012t_strict_imm_p013n_strict_imm_p013t.rds --root_dir $ROOT_DIR --cluster_res 0.8

Rscript prepro_00.R --run_mode annotation --seu_obj _data/computed/anchored/seu_merge_str_p007n_strict_str_p007t_strict_str_p008t_strict_str_p009n1_strict_str_p009n2_strict_str_p009t1_strict_str_p009t2_strict_str_p012n_strict_str_p012t_strict_str_p013n_strict_str_p013t_strict_str_p014n.rds --root_dir $ROOT_DIR --cluster_res 0.4
