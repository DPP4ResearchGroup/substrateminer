# substrateminer
`substrateminer` is a python package that offer a suite of discovery tools to investigate enzyme substrate repertorie based on sequence cleavage consensus.

## CI/CD Status
### UnitTest Status
| Branch | `pd` | `develop` | `features` |
|:--------|:--------|:--------|:--------|
|Linux|![substrateminer-main](https://github.com/manifestoso/substrateminer/actions/workflows/python-ci.yml/badge.svg?branch=main)|![substrateminer-dev](https://github.com/manifestoso/substrateminer/actions/workflows/python-ci.yml/badge.svg?branch=develop)|![substrateminer-features](https://github.com/manifestoso/substrateminer/actions/workflows/python-ci.yml/badge.svg?branch=features)|
|macOS| | | |

### GitHub Pages Status
| Page | Status |
|:--------|:--------|
|substrateminer|[![pages-build-deployment](https://github.com/manifestoso/substrateminer/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/manifestoso/substrateminer/actions/workflows/pages/pages-build-deployment)

## TL;DR
`substrateminer` provides three main categories of functionalities namely `motif`, `miner`, and `pathfinder`, `substrateminer` also integrates multi-sequence alignment tools to facilitate the analysis.
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
### Installation/Setup
Download and setup conda environment befefore install `substrateminer` package with pip.
```
$ cd substrateminer
$ conda env create -f environment-Linux.yml
$ conda activate substrateminer
```
Install `substrateminer` package
```
$ pip install .      # utility mode
$ pip install -e .   # debug mode
```
Testing the installation with the help command
```
$ substrateminer --help
```

### Examples:
#### Motif
Consensus can be derived from a collection of sequences using the `consensus` subcommand.
```
$ substrateminer motif consensus -i unittests/data/msa_align.fas -O . 
```
The conservation of the conseqnsus can be visualised using the `weblogo` subcommand.
```
substrateminer motif weblogo -i unittests/data/weblogo_align.fas -o weblogo_output.png
```

#### Miner
To identify potential substrates (degradome) from a collection of sequences (this is commonly proteom of a species), the `miner` subcommand can be used.
```
$ substrateminer miner --referencefile unittests/data/test-uniprot.txt --config unittests/test-config.yml --filtermode size --outmode inline
```

#### Pathfinder
To identify the molecular path for a substrate, the `pathfinder` subcommand can be used.
```
```

## Requirements
`substrateminer` requires the following dependencies:
- Python 3.10.11 or later
- BioPython 1.84 or later
- Numpy
- SciPy
- Pandas
- matplotlib 3.6.0 or older
- weblogo
- requests
- click
- PyYAML
- pillow


Optional binary dependencies for multiple sequence alignment:
- Clustal Omega 1.2.4 or later
- MUSCLE 5.1 or later
- MAFFT 7.475 or later

## Installation
Due to complex dependency requirements of `substrateminer`, conda is recomand to here. 
```
```

## Methods and Functions

### Motif
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
#### Consensus
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
#### Weblogo
```
usage: consensus.py weblogo [-h] -i INPUTFILE -o FILENAME [-s RESOLUTION] [-F FILETYPE]

options:
  -h, --help     show this help message and exit
  -i INPUTFILE   Input alignment file/self-aligned file in FASTA/text format.
  -o FILENAME    Output filename for weblogo image.
  -s RESOLUTION  Resolution of the weblogo image.
  -F FILETYPE    File type of the output image.
```
### Miner
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

### Pathfinder
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
## GitHub Actions CI - UnitTests sequence

CI/CD is carried out with GitHub Actions workflow and consists following steps:

- Checks out the repository.
- Sets up Python.
- Caches the conda environment.
- Installs Miniconda and creates the conda environment.
- Runs CLI tests with pip.
- Runs unit tests with unittest.

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Issue/Bug Reports and Contributions
We welcome contributions to `substrateminer`. If you would like to contribute, please fork the repository and submit a pull request. For major changes, please [open an issue](https://github.com/manifestoso/substrateminer/issues) first to discuss what you would like to change.
