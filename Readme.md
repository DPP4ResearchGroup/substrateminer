# substrateminer

![Header Image](docs/assets/png/logo-2.jpg)

## Overview

[![status](https://joss.theoj.org/papers/09b68390d13060894366a1626b45aa05/status.svg)](https://joss.theoj.org/papers/09b68390d13060894366a1626b45aa05)

`substrateminer` is a python package that offer a suite of discovery tools to investigate enzyme substrate repertorie based on sequence cleavage consensus.

## CI/CD Status

### UnitTest Status

| Branch | `main`                                                                                                                                                                                                                                     | `develop`                                                                                                                                                                                                                                    | `features`                                                                                                                                                                                                                                         |
|:-------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Linux  | ![substrateminer-main](https://github.com/DPP4ResearchGroup/substrateminer/actions/workflows/python-ci.yml/badge.svg?branch=main)                                                                                                          | ![substrateminer-dev](https://github.com/DPP4ResearchGroup/substrateminer/actions/workflows/python-ci.yml/badge.svg?branch=develop)                                                                                                          | ![substrateminer-features](https://github.com/DPP4ResearchGroup/substrateminer/actions/workflows/python-ci.yml/badge.svg?branch=features)                                                                                                          |
| macOS  | [![substrateminer-main](https://github.com/DPP4ResearchGroup/substrateminer/actions/workflows/substrateminer-mac.yml/badge.svg?branch=main)](https://github.com/DPP4ResearchGroup/substrateminer/actions/workflows/substrateminer-mac.yml) | [![substrateminer-dev](https://github.com/DPP4ResearchGroup/substrateminer/actions/workflows/substrateminer-mac.yml/badge.svg?branch=develop)](https://github.com/DPP4ResearchGroup/substrateminer/actions/workflows/substrateminer-mac.yml) | [![substrateminer-features](https://github.com/DPP4ResearchGroup/substrateminer/actions/workflows/substrateminer-mac.yml/badge.svg?branch=features)](https://github.com/DPP4ResearchGroup/substrateminer/actions/workflows/substrateminer-mac.yml) |

### Documentation Status

| Page           | Status                                                                                                                                                                                                                                        |
|:---------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| substrateminer | [![pages-build-deployment](https://github.com/DPP4ResearchGroup/substrateminer/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/DPP4ResearchGroup/substrateminer/actions/workflows/pages/pages-build-deployment) |

## TL;DR

### Installation Guide

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

#### Requirements

`substrateminer` requires the following dependencies:
- Python >= 3.10.11
- BioPython >= 1.84
- Numpy
- SciPy
- Pandas
- matplotlib <= 3.6.0
- weblogo
- requests
- click
- PyYAML
- pillow
  Optional binary dependencies for multiple sequence alignment:
- Clustal Omega >= 1.2.4
- MUSCLE >= 5.1
- MAFFT >= 7.475

### Quick Start Guide

`substrateminer` provides three main categories of functionalities namely `motif`, `miner`, and `pathfinder`. `substrateminer` also integrates multi-sequence alignment tools to facilitate the analysis.

```
Usage: substrateminer [OPTIONS] COMMAND [ARGS]...

  A suite to tools to discover enzyme substrates based on sequence consensus.

  Main entry point for substrateminer CLI.

Options:
  --help  Show this message and exit.

Commands:
  miner       Filter amino acid sequences from a reference file.
  motif       Interface for consensus sequence determination.
  msa         Interface for multi-sequence alignments
  pathfinder  Find the pathological/molecular path for a substrate.
```

#### Usage Examples

<a id="motif"></a> **Motif**

Consensus can be derived from a collection of sequences using the `consensus` subcommand.

```
$ substrateminer motif consensus -i unittests/data/msa_align.fas -O . 
```

The conservation of the conseqnsus can be visualised using the `weblogo` subcommand.

```
substrateminer motif weblogo -i unittests/data/weblogo_align.fas -o weblogo_output.png
```

<a id="miner"></a> **Miner**

To identify potential substrates (degradome) from a collection of sequences (this is commonly proteom of a species), the `miner` subcommand can be used.

```
$ substrateminer miner --referencefile unittests/data/test-uniprot.txt --config unittests/test-config.yml --filtermode size --outmode inline
```

<a id="pathfinder"></a> **Pathfinder**

To identify the molecular path for a substrate, the `pathfinder` subcommand can be used.

```
$ substrateminer pathfinder -i unittests/data/uniprot_id_short.txt -o path.txt -a
```

## Construct a customised workflow

`substrateminer` is designed to provide a suite of methods to investigate enzyme substrate repertorie based on sequence cleavage consensus. The package is modular and extensible and can be used to design custom workflows. The following demonstrates a typical workflow:

![Design Workflow](png/workflow.png)

### Methods and Functions Overview

##### Multiple Sequence Alignment (MSA)

```
usage: msa.py [-h] -i INPUT -o OUTPUT -m METHOD

Perform multiple sequence alignment

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input file path
  -o OUTPUT, --output OUTPUT
                        Output file path
  -m METHOD, --method METHOD
                        Alignment method (clustalomega, mafft, muscle)
```

##### Motif

```
usage: consensus.py [-h] {consensus,weblogo} ...

Determine consensus sequence from a multiple sequence alignment (MSA) and draw summative plots.

positional arguments:
  {consensus,weblogo}
    consensus          Determine consensus sequence from a multiple sequence alignment (MSA) and
                       draw sequence entropy and gap frequency plots.
    weblogo            Generate a weblogo image from an input file.

options:
  -h, --help           show this help message and exit
```

###### Consensus

```
usage: consensus.py consensus [-h] -i Input alignment file in FASTA format.
                              [-o Output gap stripped FASTA file name] [-O Output directory]
                              [-c Method for removing insertions] [-t Gap frequency threshold]
                              [-f]

options:
  -h, --help            show this help message and exit
  -i Input alignment file in FASTA format.
                        Filename for FASTA alignment
  -o Output gap stripped FASTA file name
                        Output FASTA filename. If not given will use name of input FASTA file as
                        template to name output files.
  -O Output directory   Output directory for all output files. If not given will use directory of
                        input FASTA file.
  -c Method for removing insertions
                        Desired method for removing insertions. 1 = Positions with gap frequencies
                        < threshold (0.5 default, change with -t flag). 2 = Positions with residue
                        as most frequent character. 3 = Positions with residues in a specific
                        sequence. If not given will ask for user input upon running script. See
                        README for further explantion of methods.
  -t Gap frequency threshold
                        Gap frequecy threshold to define a consensus positions. Only valid for
                        Option 1 for removing insertions. Must be a value between 0 and 1
                        (default: 0.5)
  -f                    Include flag to prevent saving images of MSA data analysis.
```

###### Weblogo

```
usage: consensus.py weblogo [-h] -i INPUTFILE -o FILENAME [-s RESOLUTION] [-F FILETYPE]

options:
  -h, --help     show this help message and exit
  -i INPUTFILE   Input alignment file/self-aligned file in FASTA/text format.
  -o FILENAME    Output filename for weblogo image.
  -s RESOLUTION  Resolution of the weblogo image.
  -F FILETYPE    File type of the output image.
```

##### Miner

```
Usage: substrateminer miner [OPTIONS]

  Filter amino acid sequences from a reference file.

Options:
  --referencefile TEXT            The reference file containing sequences.
                                  [required]
  --referencetype [swiss|genbank|embl]
                                  The type of reference file. Default is
                                  swiss.
  --filtermode [size|motif|loc]   The mode of filtering. Default is size.
                                  [required]
  --config TEXT                   The path to the configuration file.
  --stats                         Generate statistics for the filtered
                                  sequences.
  --outmode [all|file|inline]     The output mode for the filtered sequences.
  --outputfilename TEXT           The output file name for the filtered
                                  sequences.
  --outputfiletype [fasta|text|txt|genbank|swiss]
                                  The output file type for the filtered
                                  sequences.
  --help                          Show this message and exit.
```

##### Pathfinder

```
Usage: substrateminer pathfinder [OPTIONS]

  Find the pathological/molecular path for a substrate.

Options:
  -i, --input PATH     Input file path
  -o, --output TEXT    Output file path
  -a, --api            Use KEGG API to retrieve pathways and diseases
  -u, --uniprots TEXT  UniProt ID for a protein, comma-separated for multiple
                       IDs (e.g., P12345,Q67890) or space-separated for
                       multiple IDs (e.g., "P12345 Q67890")
  -g, --orgs TEXT      Organism code for the KEGG API (default: hsa)
  --help               Show this message and exit.
```

## GitHub Actions CI manual

### UnitTests Sequence

CI/CD is carried out with GitHub Actions workflow and consists following steps:

- Checks out the repository.
- Sets up Python.
- Caches the conda environment.
- Installs Miniconda and creates the conda environment.
- Runs CLI tests with pip.
- Runs unit tests with unittest.

## License

This project is licensed under the MIT License. See the [LICENSE](readme_docs/LICENSE) file for details.

## Issue/Bug Reporting

Any issues you encounter with `substrateminer`, please report by [open a bug issue](https://github.com/manifestoso/substrateminer/issues) and provide as much details as possible, including examples, error messages and environment setup will be highly appreciated.

## Contributing

We welcome contributions to `substrateminer`. To contribute, please follow the steps below:
1. Fork the repository to your designated location.
2. Create a new branch with a descriptive name for your proposed feature and/or bugfix.
3. Make your changes and commit them with clear and concise commit messages.
4. Push your changes to your forked repository.
5. Submit a pull request.

> ****Major changes:**** please [open an issue](https://github.com/manifestoso/substrateminer/issues) first to discuss what you would like to change.
