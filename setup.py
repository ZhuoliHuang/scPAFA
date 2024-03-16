from setuptools import setup, find_packages

setup(
    name='scPAFA',
    version='0.1.2',
    author='Zhuoli Huang',
    author_email='bioagr_huangzl@163.com',
    url='https://github.com/ZhuoliHuang/scPAFA',
    license='GPL-3.0',
    description='Single Cell Pathway Activity Factor Analysis',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    packages=find_packages(),
    install_requires=[
        'scanpy>=1.9.0',
        'anndata>=0.9.0',
        'mofax == 0.3.6',
        'pandas>=2.0.3',
        'scipy == 1.11.3',
        'plotnine',
        'statsmodels',
        'pathos',
        'tqdm',
        'nheatmap',
        'seaborn'
    ]
)
