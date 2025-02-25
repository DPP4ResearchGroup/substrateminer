# substrateminer
`substrateminer` is a python package that offer a suite of discovery tools to investigate enzyme substrate repertorie based on sequence cleavage consensus.
`substrateminer` is a python package that offer a suite of discovery tools to investigate enzyme substrate repertorie based on sequence cleavage consensus.

## CI/CD Status
<<<<<<< HEAD
|Branch|`main`|`develop`|`features`|
|---|---|---|
|Status|
[![substrateminer releasing](https://github.com/manifestoso/substrateminer/actions/workflows/python-ci.yml/badge.svg?branch=main)](https://github.com/manifestoso/substrateminer/actions/workflows/python-ci.yml)|[![substrateminer releasing](https://github.com/manifestoso/substrateminer/actions/workflows/python-ci.yml/badge.svg?branch=develop)](https://github.com/manifestoso/substrateminer/actions/workflows/python-ci.yml)|[![substrateminer releasing](https://github.com/manifestoso/substrateminer/actions/workflows/python-ci.yml/badge.svg?branch=features)](https://github.com/manifestoso/substrateminer/actions/workflows/python-ci.yml)|
=======
### UnitTest Status
| Branch | `pd` | `develop` | `features` |
|:--------|:--------|:--------|:--------|
|Linux|![substrateminer-main](https://github.com/manifestoso/substrateminer/actions/workflows/python-ci.yml/badge.svg?branch=main)|![substrateminer-dev](https://github.com/manifestoso/substrateminer/actions/workflows/python-ci.yml/badge.svg?branch=develop)|![substrateminer-features](https://github.com/manifestoso/substrateminer/actions/workflows/python-ci.yml/badge.svg?branch=features)|
|macOS| [![substrateminer-main](https://github.com/manifestoso/substrateminer/actions/workflows/substrateminer-mac.yml/badge.svg?branch=main)](https://github.com/manifestoso/substrateminer/actions/workflows/substrateminer-mac.yml)| [![substrateminer-dev](https://github.com/manifestoso/substrateminer/actions/workflows/substrateminer-mac.yml/badge.svg?branch=develop)](https://github.com/manifestoso/substrateminer/actions/workflows/substrateminer-mac.yml) | [![substrateminer-features](https://github.com/manifestoso/substrateminer/actions/workflows/substrateminer-mac.yml/badge.svg?branch=features)](https://github.com/manifestoso/substrateminer/actions/workflows/substrateminer-mac.yml) |

### GitHub Pages Status
| Page | Status |
|:--------|:--------|
|substrateminer|[![pages-build-deployment](https://github.com/manifestoso/substrateminer/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/manifestoso/substrateminer/actions/workflows/pages/pages-build-deployment)
>>>>>>> a2e63c9 (fix #36 - doc prerelease cr)

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
Examples:
### Motif
Consensus can be derived from a collection of sequences using the `consensus` subcommand.
```
$ substrateminer motif consensus -i input.fasta -O output/results/
```
The conservation of the conseqnsus can be visualised using the `weblogo` subcommand.
```
substrateminer motif weblogo -i input.fasta -o weblogo_output.png
```

### Miner
To identify potential substrates (degradome) from a collection of sequences (this is commonly proteom of a species), the `miner` subcommand can be used.
```
$ substrateminer miner -i input.fasta -o output.fasta
$ substrateminer msa -i input.fasta -o output.fasta
$ substrateminer pathfinder -i input.fasta -o output.fasta
```

```

## Requirements
`substrateminer` requires the following dependencies:
- Python 3.10.11 or later
- BioPython 1.84 or later

Optional binary dependencies for multiple sequence alignment:
- Clustal Omega 1.2.4 or later
- MUSCLE 5.1 or later
- MAFFT 7.475 or later

## Installation
Dependencies
```


## Methods and Functions
### Multiple Sequence Alignment (MSA)
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

### Consensus

## Contributing
### Feedback
Feedback is greatly appreciated. If you have any questions, comments, or concerns, please feel free to open an issue on the project's issue tracker.
## License
`substrateminer` is distributed under the MIT license.
