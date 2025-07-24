CI/CD is carried out with GitHub Actions workflow and consists following steps:

- Checks out the repository.
- Sets up Python.
- Caches the conda environment.
- Installs Miniconda and creates the conda environment.
- Runs CLI tests with pip.
- Runs unit tests with unittest.