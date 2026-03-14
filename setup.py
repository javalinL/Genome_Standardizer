from setuptools import setup, find_packages

setup(
    name='genome-standardizer',
    version='3.0.0',
    description='A robust pipeline for standardizing complex polyploid genomes.',
    author='L. javalin',
    author_email='beiningjia412@gmail.com',
    packages=find_packages(),
    install_requires=[
        'biopython>=1.78',
        'tqdm>=4.60.0',
    ],
    entry_points={
        'console_scripts': [
            'gstd=genome_standardizer.main:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
