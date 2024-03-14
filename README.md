# scPAFA

<img src="https://github.com/ZhuoliHuang/scPAFA/assets/61071877/cde987e2-bfc7-4ece-b718-494a0318cb67" align="right" width="192" height="160">



**Single Cell Pathway Activity Factor Analysis**

A Python library designed for large-scale single-cell datasets allowing rapid PAS computation and uncovering biologically interpretable disease-related multicellular pathway modules, which are low-dimensional representations of disease-related PAS variance in multiple cell types.



# Workflow

<img src="https://github.com/ZhuoliHuang/scPAFA/assets/61071877/b8bdee9e-b98f-467a-b345-7ffb5acfbfd9" width="800" height="400">

# Installation

We recommend using scPAFA in a virtual environment.
```
conda create -n scPAFA_env python=3.10
```

### Install from PyPi

```
conda activate scPAFA_env
pip install scPAFA
```

### Install from GitHub

In your workdir
```
conda activate scPAFA_env
git clone https://github.com/ZhuoliHuang/scPAFA
cd ./scPAFA
python setup.py install
```
# Tutorial

**Pathway input:** Download pathway information and generate a pathway dictionary
  
The pathway input of scPAFA is a Python dictionary, each item with a pathway name as a key and a list of genes as values.

(1) Download pathway collection

Pathway collection can be downloaded from [MsigDB](https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H) ('JSON bundle' is recommended), or [NCATS bioplanet](https://tripod.nih.gov/bioplanet/download/pathway.csv). Users can also use a custom pathway collection.

(2) Generate pathway dictionary

We provided [examples](https://github.com/ZhuoliHuang/scPAFA/blob/main/tutorial/generate_pathway_input.ipynb) of constructing pathway dictionary from the MsigDB and NCATS bioplanet databases.

**Step1:** Calculate Pathway Activity Score

In step1, single-cell gene expression matrix and collection of pathways are used to compute PAS by ‘fast_Ucell’([example](https://github.com/ZhuoliHuang/scPAFA/blob/main/tutorial/step1_generate_PAS_UCell.ipynb)) or ‘fast_score_genes’([example](https://github.com/ZhuoliHuang/scPAFA/blob/main/tutorial/step1_generate_PAS_score_genes.ipynb)). These functions are more computationally efficient implementation of UCell and AddModuleScore (also known as ‘score_genes’ in Scanpy)

**Step2~3:** Pseudobulk processing and MOFA model training

In step 2, the single-cell PAS matrix is reformatted into a suitable input ([a long-table-like pandas dataframe](https://github.com/bioFAM/mofapy2/blob/master/mofapy2/notebooks/getting_started_python.ipynb)) for [Multi-Omics Factor Analysis (MOFA)](https://biofam.github.io/MOFA2/index.html) along with cell-level metadata including sample/donor, cell type, and technical batch information. In step 3, MOFA model is trained to capture variance in PAS among different samples. Notably, MOFA contains general framework (single-group framework) and multi-group framework,  [the aim of multi-group framework is to find out which sources of variability are shared between the different groups](https://biofam.github.io/MOFA2/faq.html). In the presence of clearly known batch effects, we recommend using multi-group MOFA+ framework for correction. We provided [examples](https://github.com/ZhuoliHuang/scPAFA/blob/main/tutorial/steps2%263_Pseudobulk_processing_and_MOFA_model_training.ipynb) of steps 2 and 3.

**Step4:** Downstream analysis of the MOFA Model

In step 4, together with sample-level clinical metadata, disease-related multicellular pathway modules (latent factor and its corresponding weights of pathways across cell types) can be identified by statistical analysis. Downstream analyses include characterization and interpretation of multicellular pathway modules, sample/donor stratification, and visualization of high-weighted pathways ([example](https://github.com/ZhuoliHuang/scPAFA/blob/main/tutorial/step4_Downstream_analysis_of_the_MOFA_Model.ipynb)).

# Supported Systems
Including but not limited to：

Ubuntu 20.04

Windows 10/11

# Citation
scPAFA is now reported on biorxiv.

Uncovering disease-related multicellular pathway modules on large-scale single-cell transcriptomes with scPAFA

Zhuoli Huang, Yuhui Zheng, Weikai Wang, Wenwen Zhou, Chen Wei, Xiuqing Zhang, Xin Jin, Jianhua Yin

bioRxiv 2024.03.11.584023; doi: https://doi.org/10.1101/2024.03.11.584023
