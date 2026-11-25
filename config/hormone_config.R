MERGED_RDS <- "/data01/liuxy/Project/CK2_added_v3_Results/seurat_obj_with_CK2_v3.rds"
ANALYSIS_MODE <- "annotated"
OUT_DIR <- "/data01/liuxy/Project/激素扰动评分"
RESUME <- TRUE
CELLTYPE_COL <- NULL
CONDITION_COL <- "treatment"
SAMPLE_COL <- "sample_id"
EXPECTED_REPLICATES <- NULL
MIN_REPLICATES_PER_CONDITION <- 2L
PAIR_COL <- NULL
COVARIATE_COLS <- character()
MIX_NEIGHBOR_BLOCK_COL <- NULL
CELLTYPE_KEEP <- NULL
CONDITION_MAP <- list(
 CK=c("CK","Control","CTRL","REF","Reference","对照"),
 ABA=c("ABA","Abscisic_acid","AbscisicAcid"),
 eBL=c("eBL","EBL","epiBL","Epibrassinolide"),
 IAA=c("IAA","Indole_3_acetic_acid","IndoleAceticAcid","Auxin"),
 MEJA=c("MEJA","MeJA","Methyl_jasmonate","MethylJasmonate","MJ"),
 SA=c("SA","Salicylic_acid","SalicylicAcid"))
