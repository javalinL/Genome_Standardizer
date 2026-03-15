from setuptools import setup, find_packages

setup(
    name='genome-standardizer',
    version='4.1.8',
    description='A robust pipeline for standardizing complex polyploid genomes, '
                'with an optional Terminal User Interface (TUI).',
    author='L. javalin',
    author_email='beiningjia412@gmail.com',
    packages=find_packages(),
    install_requires=[
        'biopython>=1.78',
        'tqdm>=4.60.0',
    ],
    extras_require={
        'tui': [
            'textual>=0.47.0',
        ],
    },
    package_data={
        'genome_standardizer': [
            'locales/*.json',
            'tui/*.tcss',
        ],
    },
    entry_points={
        'console_scripts': [
            'gstd=genome_standardizer.main:main',
            'gstd-tui=genome_standardizer.tui.app:run_tui',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires='>=3.8',
)
