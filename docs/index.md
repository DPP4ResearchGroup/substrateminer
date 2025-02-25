---
layout: default
---


```sh
substrateminer msa -i <input_file> -o <output_file> -m <method>
```

- `<input_file>`: Path to the input file containing sequences in FASTA format.
- `<output_file>`: Path to the output file where the alignment results will be saved.
- `<method>`: Alignment method to use. Options are `clustalomega`, `mafft`, or `muscle`.

### Example

```sh
substrateminer msa -i sequences.fasta -o alignment.aln -m clustalomega
```

## Running Tests

To run the unit tests, use the following command:

```sh
conda run -n substrate-miner python -m unittest discover -s unittests -p "test_*.py"
```

## GitHub Actions CI

The project includes a GitHub Actions workflow for continuous integration. The workflow is defined in `.github/workflows/python-ci.yml` and performs the following steps:

- Checks out the repository.
- Sets up Python.
- Caches the conda environment.
- Installs Miniconda and creates the conda environment.
- Runs CLI tests with pip.
- Runs unit tests with unittest.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contributing

We welcome contributions to SubstrateMiner. If you would like to contribute, please fork the repository and submit a pull request. For major changes, please open an issue first to discuss what you would like to change.

## Contact

If you have any questions or need further assistance, please contact us at [your-email@example.com](mailto:your-email@example.com).