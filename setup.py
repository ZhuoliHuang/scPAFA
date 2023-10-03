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
        'pandas>=2',
        'dask',
        'pandarallel'
    ]
)