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

(1) Method based on UCell 

(2) Methods based on AddMouduleScore(Seurat)/scanpy.tl.score_genes(Scanpy)

**Step2:** Pseudobulk processing and MOFA model training

**Step3:** MOFA model training

(1)  single-group MOFA model training

(2)  multi-group MOFA model training

**step4:** Downstream analysis of the MOFA Model
