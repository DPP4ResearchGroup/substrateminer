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