from setuptools import setup, find_packages

setup(
    name="scvalue",
    version="0.0.1",
    author="Li Huang",
    author_email="hl@ism.cams.cn",
    description="scValue: value-based subsampling of large-scale single-cell transcriptomic data",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/LHBCB/scvalue",
    packages=find_packages(),
    install_requires=[
        "scikit-learn>=1.5.2",
        "scipy>=1.14.0",
        "numpy>=1.26.4",
        "pandas>=1.5.3",
        "joblib>=1.4.2"
    ],
    license='BSD-3-Clause',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
)

