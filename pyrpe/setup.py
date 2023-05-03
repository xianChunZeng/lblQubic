"""Python tools for Jaqal"""

from setuptools import setup, find_packages

name = "pyRPE"
description = "Python implementation of Robust Phase Estimation"
version = "1.0"

setup(
    name=name,
    description=description,
    long_description=open("README.md", "r").read(),
    long_description_content_type="text/markdown",
    version=version,
    author="Antonio Russo, Kenneth Rudinger",
    author_email="arusso@sandia.gov",
    packages=find_packages(include=["quapack.*"], where="src"),
    package_dir={"": "src"},
    package_data={"jaqalpaq.parser": ["jaqal_grammar.lark"]},
    install_requires=["numpy"],
    python_requires=">=3.6.5",
    platforms=["any"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "License :: OSI Approved :: BSD Software License",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Physics",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Unix",
    ],
)
