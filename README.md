# Tumor_immune_rejection

Based on S. Naulaerts' code from https://github.com/BVDElab/Tumour_immune_rejection.git .
I was asked to produce violin plots checking Stabilin-1 expression in the different conditions and cell types.

### Copy/pasted from the original repository:

2 samples were processed: clonidine and untreated for MC38-OVA mice with tumour. Tumour tissue was extracted and sent for sequencing. Samples were processed with Scanpy, using the external algorithm PhenoGraph for subpopulation detection. Clusters were manually annotated, differential testing was performed with wilcoxon rank sum tests implemented in Diffxpy. Enrichment analysis has been conducted with GSEApy.

To run, execute the Python script Main.py. Extract the two ZIP files from the GEO project and place the corresponding folders in "Data". You will also need to place the set enrichment definitions in a folder named "Bader_GSEA_GMTs". These you can find here: http://baderlab.org/GeneSets
