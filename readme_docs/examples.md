#### Motif
Consensus can be derived from a collection of sequences using the `consensus` subcommand.
```
$ substrateminer motif consensus -i unittests/data/msa_align.fas -O . 
```
The conservation of the consensus can be visualised using the `weblogo` subcommand.
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
$ substrateminer pathfinder -i unittests/data/uniprot_id_short.txt -o path.txt -a
```
