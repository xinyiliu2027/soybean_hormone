# Soybean hormone single-cell atlas

Target repository: xinyiliu2027/soybean_hormone. Target Pages path: /soybean_hormone/#atlas.

The local website contains real CK, ABA, eBL, IAA, MEJA and SA cells from the previously exported merged RDS. It includes five confirmed marker expression vectors. It does not yet contain all-gene hormone summaries or downstream computed results. Do not reuse abiotic stress DE results as hormone contrasts.

The shared configuration is config/hormone_config.R and reflects the user's exact treatment aliases, treatment/sample_id metadata columns, variable group replicate counts, and minimum two biological replicates. Do not overwrite the source RDS.

scripts/annotation preserves the six annotation scripts from soybean_abiotic. These start from already loaded Seurat samples/reference objects; the source repository does not supply a complete soybean raw-data QC pipeline. A strawberry DoubletFinder example must not be relabeled as soybean QC.

source_extracted contains the four DOCX code texts without scientific changes. scripts/hormone contains reviewed working copies; R syntax parsing passes. Non-executable comments and repeated blank lines were removed from 12 scripts (693 lines), with identical R parse trees before and after that cleanup. Two document transcription issues were fixed: a bare Chinese heading in Monocle3, and a trailing ===== after mt$Celltype in hdWGCNA. The unguarded Monocle3 main_resume() call is disabled to prevent accidental execution of an old saved run. The original TF code mixes hormone and stress results; the original PAGA code includes nine treatments and embeds Python. The hdWGCNA code reads a module-scored object ; its reviewed copy accepts unequal replicate counts with a minimum of two per condition. Existing cached results and biological trajectory roots must be audited before running any of these scripts.

The repository has been created. See the deployment record for publication status. Analysis scripts have not been run on the hormone RDS in this session. Raw QC uses reviewed thresholds supplied through config/samples.csv, not invented thresholds. Organelle filtering and DoubletFinder require a further sample-specific audit.
