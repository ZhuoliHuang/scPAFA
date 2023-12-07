from setuptools import setup, find_packages

setup(
    name='scPAFA',
    version='0.1.0',
    author='Zhuoli Huang',
    description='Single Cell Pathway Activity Factor Analysis',
    packages=find_packages(),
    install_requires=[
        'scanpy>=1.9.0',
        'anndata>=0.9.0',
        'mofax == 0.3.6',
        'pandas>=2.0.3',
        'plotnine',
        'statsmodels',
        'pathos',
        'tqdm',
        'nheatmap'
    ]
)