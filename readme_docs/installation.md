Due to complex dependency requirements of `substrateminer`, conda is recommended here. Please ensure that you have conda installed on your system. If you do not have conda installed, please refer to the [Miniconda installation guide](https://docs.conda.io/en/latest/miniconda.html).

Firstly, download a copy of the latest release of `substrateminer` from the [GitHub Releases page](https://github.com/DPP4ResearchGroup/substrateminer/releases/latest) to a chosen local path before setup the required `conda` environment as instructed below. `substrateminer` supports both Linux and MacOS platforms, please choose the appropriate environment file based on your platform as follows:
```
# For Linux users
$ cd substrateminer
$ conda env create -f environment-Linux.yml # Linux support
$ conda activate substrateminer
```
OR
```
# For MacOS users
$ cd substrateminer
$ conda env create -f environment-macOS.yml # MacOS support
$ conda activate substrateminer
```

Then install `substrateminer` package as below:
```
$ pip install .      # utility mode
$ pip install -e .   # debug mode
```
Testing the installation with the help command
```
$ substrateminer --help
```