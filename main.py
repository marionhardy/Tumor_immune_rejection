# -*- coding: utf-8 -*-

#%% Package import

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy
import random
from collections import Counter
import diffxpy.api as de
import os
from gseapy.plot import gseaplot


#%% Data Loading and set up

plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['svg.fonttype'] = 'none'

# Make sure environment is always the same
random.seed(770)
np.random.seed(770)

# Load data
untreated = sc.read_10x_mtx("Data/untreated")
sc.pl.highest_expr_genes(untreated, n_top=20, )
sns.despine()
plt.tight_layout()
plt.savefig("Top20-most-expressed-genes_untreated.svg")
plt.close()

treated = sc.read_10x_mtx("Data/treated")
sc.pl.highest_expr_genes(treated, n_top=20, )
sns.despine()
plt.tight_layout()
plt.savefig("Top20-most-expressed-genes_treated.svg")
plt.close()

#%% Set up regression

# Mark mitochondrial content
untreated.var['mt'] = untreated.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
treated.var['mt'] = treated.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'

# QC
sc.pp.calculate_qc_metrics(untreated, qc_vars=['mt'], inplace=True)
sc.pp.calculate_qc_metrics(treated, qc_vars=['mt'], inplace=True)

sc.pl.violin(untreated, ['n_genes_by_counts', 'total_counts', 'pct_counts_in_top_100_genes', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
plt.savefig("untreated_qc.svg")
plt.close()
print("untreated", "n_genes_by_counts", np.mean(untreated.obs["n_genes_by_counts"]))
print("untreated", "total_counts", np.mean(untreated.obs["total_counts"]))
print("untreated", "pct_counts_mt", np.mean(untreated.obs["pct_counts_mt"]))

sc.pl.violin(treated, ['n_genes_by_counts', 'total_counts', 'pct_counts_in_top_100_genes', 'pct_counts_mt'], jitter=0.4,
             multi_panel=True)
plt.savefig("treated_qc.svg")
plt.close()
print("treated", "n_genes_by_counts", np.mean(treated.obs["n_genes_by_counts"]))
print("treated", "total_counts", np.mean(treated.obs["total_counts"]))
print("treated", "pct_counts_mt", np.mean(untreated.obs["pct_counts_mt"]))

sc.pl.scatter(untreated, x='total_counts', y='n_genes_by_counts')
plt.savefig("totalcounts-ngenes_untreated.svg")
plt.close()

sc.pl.scatter(treated, x='total_counts', y='n_genes_by_counts')
plt.savefig("totalcounts-ngenes_treated.svg")
plt.close()

# Concatenate
adata = untreated.concatenate(treated)

# Filter genes and cells out
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

#%% Normalization

# Normalize per 10k reads + logtransform, reduce variability
sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed=False, inplace=True)
sc.pp.log1p(adata, base=2)
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
adata = adata[:, [x for x in adata.var_names if not x.startswith("mt-")]]


# Variable gene detection
adata.obs["treatment"] = adata.obs["batch"].map({"0": "untreated", "1": "treated"})
sc.pp.highly_variable_genes(adata, flavor="seurat")
sc.pl.highly_variable_genes(adata)
plt.savefig("highly_variable_prop.svg")
plt.close()

#%% PCA

# PCA
highly_variable = adata.var.loc[adata.var["highly_variable"]].index
sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True, random_state=770)
sc.pl.pca(adata, color=["treatment"], cmap="YlOrRd")
plt.tight_layout()
sns.despine()
plt.savefig("batch_correction.svg")
plt.close()

# Cumulative variance explained:
cml_var_explained = np.cumsum(adata.uns['pca']['variance_ratio'])
x = range(len(adata.uns['pca']['variance_ratio']))
y = cml_var_explained
plt.scatter(x, y, s=4)
plt.xlabel('PC')
plt.ylabel('Cumulative variance explained')
plt.title('Cumulative variance explained by PCs')
plt.savefig("components_needed_pervariance.svg")
plt.close()


#%% Do PCA on highly variable genes

subdata = adata[:, highly_variable]
sc.external.tl.phenograph(subdata, clustering_algo="leiden", seed=770, k=25, n_iterations=5000)
adata.obs["pheno_leiden"] = subdata.obs["pheno_leiden"]
sc.tl.rank_genes_groups(adata, 'pheno_leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
plt.savefig("rank_genes_groups.svg")
plt.close()


#%%  Data projections

sc.tl.pca(adata, n_comps=26, use_highly_variable=True)
sc.pp.neighbors(adata, n_neighbors=25, n_pcs=26)
sc.tl.umap(adata)
sc.pl.umap(adata, color=["pheno_leiden"], legend_loc="on data")
sns.despine()
plt.tight_layout()
plt.savefig("leiden_umap.svg")
plt.close()

# Plot per treatment

sc.tl.umap(adata)
sc.pl.umap(adata[adata.obs["treatment"].isin(["treated"])], color=["pheno_leiden"], legend_loc="on data")
sns.despine()
plt.tight_layout()
plt.show()
plt.savefig("leiden_umap_clonidine.svg")
plt.close()

sc.tl.umap(adata)
sc.pl.umap(adata[adata.obs["treatment"].isin(["untreated"])], color=["pheno_leiden"], legend_loc="on data")
sns.despine()
plt.tight_layout()
plt.savefig("leiden_umap_untreated.svg")
plt.close()

# Mappings we identified
mapdic = {0: "Macrophages_2", 1: "Macrophages_1", 2: "Macrophages_3", 3: "Macrophages_1",
          4: "Tumour cells", 5: "Macrophages_5", 6: "Macrophages_4", 7: "CD8+ T cells", 
          8: "Mfap5+ fibroblasts", 9: "Erythrocytes", 10: "CD4+ T cells", 
          11: "Neutrophils", 12: "Tumour cells",
          13: "Dendritic cells", 14: "Macrophages_3", 15: "Endothelial cells",
          16: "CD8+ T cells", 17: "Endothelial cells"}

adata.obs["Pheno_name"] = adata.obs["pheno_leiden"].replace(mapdic)
adata.obs["Pheno_name"].fillna("Unknown", inplace=True)

mycmap = {'CD4+ T cells': sns.color_palette()[0], 'CD8+ T cells': sns.color_palette()[1],
          'Macrophages_1': "goldenrod",
          'Macrophages_2': "salmon",
          'Macrophages_3': "red",
          'Macrophages_4': "darkcyan",
          'Macrophages_5': "hotpink",
          "Tumour cells": "sienna",
          "Mfap5- fibroblasts": sns.color_palette()[7],
          "Mfap5+ fibroblasts": sns.color_palette()[8],
          "Endothelial cells": "pink", "Erythrocytes": "purple",
          "Neutrophils": "limegreen", "Dendritic cells": "green"}

sc.pl.umap(adata, color="Pheno_name", legend_loc="on data", palette=mycmap, legend_fontweight="normal")
plt.tight_layout()
sns.despine()
plt.savefig("UMAP_celltype.svg")
plt.close()

for tment in set(adata.obs["treatment"]):
    tindexes = adata.obs.loc[adata.obs["treatment"] == tment].index
    subdata = adata[tindexes, :]
    sc.pl.umap(subdata, color="Pheno_name", legend_loc="on data", palette=mycmap, legend_fontweight="normal")
    plt.tight_layout()
    sns.despine()
    plt.savefig(tment + ".svg")
    plt.close()
    
    
#%% Stabilin-1 expression

#adata = sc.read_h5ad('/Volumes/BVDE_MHY/Python/scRNAseq_Stefan/adata')

sc.pl.umap(adata[adata.obs["treatment"].isin(["untreated"])], color="Stab1")

sc.pl.umap(adata[adata.obs["treatment"].isin(["treated"])], color="Stab1")


adata = adata[adata.obs.treatment.isin(['untreated','treated'])].copy()
df = sc.get.obs_df(adata, ['Stab1', 'treatment', 'Pheno_name'])
df = df.set_index(['Stab1','treatment']).stack().reset_index()
df.columns = ['AvgExpr', 'treatment','whatever','Pheno_name']


plt.figure(figsize=(20,8))
sns.violinplot(data=df, x='Pheno_name', y='AvgExpr', hue="treatment",
                split=False, linewidth=2, width=0.8)
plt.xticks(rotation=45)
plt.tight_layout()
plt.title('Stab1 expression')
plt.savefig('Stab1_expr_violin.svg', dpi = 100, bbox_inches='tight')
plt.close()

plt.figure(figsize=(20,8))
sns.violinplot(data=df, x='Pheno_name', y='AvgExpr', hue="treatment",
                split=False, linewidth=2, width=0.8)
sns.stripplot(data=df, x='Pheno_name', y='AvgExpr', hue="treatment", jitter=True,
              dodge=True, zorder = 1, size = 4, alpha = 0.5,
              linewidth=1)
plt.xticks(rotation=45)
plt.tight_layout()
plt.title('Stab1 expression')
plt.savefig('Stab1_points_expr_violin.svg', dpi = 100, bbox_inches='tight')


#%%  Piecharts

countlist = []
for agroup in set(adata.obs["Pheno_name"]):
    cdic = Counter(adata.obs["treatment"].loc[adata.obs["Pheno_name"] == agroup])
    countlist.append([agroup, "untreated", cdic["untreated"]])
    countlist.append([agroup, "treated", cdic["treated"]])
countlist = pd.DataFrame(countlist, columns=["Cell_type", "Treatment", "Count"])
countlist.to_excel("Cell_Counts.xlsx")

# Differential testing
for x in mapdic.values():
    if not os.path.exists(x):
        os.mkdir(x)

for acol in set(adata.obs["Pheno_name"]):
    subindexes = list(adata.obs.loc[adata.obs["Pheno_name"] == acol].index)
    subdata = adata[subindexes]

    test = de.test.rank_test(data=subdata, grouping="treatment", is_logged=True)
    tester = test.summary().loc[
        test.summary()["gene"].isin([x for x in test.summary()["gene"] if not x.startswith("mt")])]
    tester.to_csv(os.path.join(acol, "Diffxpy_StatisticsTreatedvsUntreated.tsv"), sep="\t")

    tester = tester.loc[(tester["gene"].isin(highly_variable))]
    rnk = tester[["gene", "log2fc"]]

    pre_res = gseapy.prerank(rnk=rnk,
                             gene_sets="Bader_GSEA_GMTs/Mouse_Human_MSigdb_February_05_2021_symbol.gmt",
                             processes=4, permutation_num=1000,
                             outdir=os.path.join(acol, 'GSEA_treatment', 'prerank_report_msigdb'), no_plot=True,
                             seed=770)
    terms = pre_res.res2d.index
    for i in range(0, len(terms)):
        gseaplot(rank_metric=pre_res.ranking, term=terms[i], **pre_res.results[terms[i]],
                 ofname=os.path.join(acol, "GSEA_treatment", "prerank_report_msigdb",
                                     terms[i].split("%")[0].replace("/", "-") + ".svg"))
        plt.close()
