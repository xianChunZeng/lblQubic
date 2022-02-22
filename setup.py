import setuptools
from setuptools.command.build_py import build_py
import sys
import subprocess



with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="chipcalibration",
    version="0.0.1",
    description="QuBiC calibration repo",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://gitlab.com/LBL-QubiC/software",
	#    project_urls={
	#    "Bug Tracker": "https://github.com/pypa/sampleproject/issues",
	#},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: LICENSE",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    python_requires=">=3.6",
	include_package_data=True,
	package_data={
		"": ["*.json"],
		},
	extra_requires={
		"numpy": ['numpy>1.20'],
		"scipy": ['scipy>1.7'],
		"sklearn": ['scipy>=0.0'],
		},
	install_requires=[
        #"qubic",
		"numpy",
		"scipy",
		"sklearn"
		],
)
