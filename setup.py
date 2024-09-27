###
# setup.py for package substratemine
#
# @author Robert QIAO (School of Biological Sciences, Flinders University, SA, Australia)
###

import setuptools

with open("README.md", "r") as readme_file:
    long_description = readme_file.read()

setuptools.setup(
    name="substrateminer",
    version="0.9.0rc1",
    author="Robert QIAO",
    author_email="robert.qiao@flinders.edu.au",
    description="A python package to discover enzyme substrates based on sequence consensus",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/DPP4ResearchGroup/substrateminer", # point to github.io upon release
    packages=setuptools.find_packages(),
    classifiers=[
        "Development Status :: 1 - Planning",
        "Environment :: Console",
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Education",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='>=3.10.13',
)
