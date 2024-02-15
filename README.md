<img src="https://github.com/ZhuoliHuang/scPAFA/assets/61071877/3b3de70c-0bb6-438b-84e8-aaa705897390" align="right" width="240" height="200">

# scPAFA

**Single Cell Pathway Activity Factor Analysis**

A Python library designed for large-scale single-cell datasets allowing rapid PAS computation and uncovering biologically interpretable disease-related multicellular pathway modules, which are low-dimensional representations of disease-related PAS variance in multiple cell types.


# Workflow

<img src="https://github.com/ZhuoliHuang/scPAFA/assets/61071877/b8bdee9e-b98f-467a-b345-7ffb5acfbfd9" width="800" height="400">


**Pathway input:** Download pathway information and generate pathway dictionary
  
The pathway input of scPAFA is a Python dictionary, each item with a pathway name as a key and a list of genes as values.

(1) Download pathway collection

Pathway collection can be downloaded from MsigDB (https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H, 'JSON bundle' is recommended), or NCATS bioplanet (https://tripod.nih.gov/bioplanet/download/pathway.csv). Users can also use a custom pathway collection.

(2) Generate pathway dictionary

We provided examples of constructing pathway dictionary from the MsigDB and NCATS bioplanet databases.(https://github.com/ZhuoliHuang/scPAFA/blob/main/tutorial/generate_pathway_input.ipynb)

**Step1:** Calculate Pathway Activity Score

In step1, single-cell gene expression matrix and collection of pathways are used to compute PAS by ‘fast_Ucell’ or ‘fast_score_genes’. 

These functions are more computationally efficient implementation of UCell and AddModuleScore (also known as ‘score_genes’ in Scanpy)

(1) Method based on UCell 

(2) Methods based on AddMouduleScore(Seurat)/scanpy.tl.score_genes(Scanpy)

**Step2~3:** Pseudobulk processing and MOFA model training
In step 2, single-cell PAS matrix is reformatted into a suitable input for MOFA along with cell-level metadata including sample/donor, cell type, and technical batch information. 


**step4:** Downstream analysis of the MOFA Model
