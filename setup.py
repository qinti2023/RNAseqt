from setuptools import setup, find_packages

setup(
    name="RNAseqt",
    version="1.1.0",
    description="A Python package to generate Snakemake workflow files",
    author="qinti",
    author_email="qinti@zju.edu.cn",
    url="https://github.com/qinti2023/RNAseqt",
    packages=find_packages(),
    install_requires=[
    ],
    entry_points={
        "console_scripts": [
            "RNAseqt = RNAseqt.main:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
