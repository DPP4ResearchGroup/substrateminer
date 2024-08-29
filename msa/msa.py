from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Align.Applications import MafftCommandline
import argparse
from Bio import SeqIO

def perform_msa(input_file, output_file, method):
    """
    Perform multiple sequence alignment using ClustalOmega or MAFFT.

    Parameters:
    - input_file (str): Path to the input file containing sequences in FASTA format.
    - output_file (str): Path to the output file where the alignment results will be saved.
    - method (str): Alignment method to use. Options are "clustalomega" or "mafft".

    Raises:
    - ValueError: If an invalid alignment method is specified.

    Returns:
    - None

    Other Bio.Align.Applications options:
    - EmbossCommandline: This class provides a command line wrapper for the EMBOSS suite of sequence analysis tools. It can be used to perform multiple sequence alignment using EMBOSS tools like "needle" or "water".

    Note:
    - The input file should contain sequences in FASTA format.
    - The alignment results will be printed to the console and saved to the output file.
    """
    
    try:
        # Read input sequences from a file
        sequences = SeqIO.parse(input_file, "fasta")

        # Perform multiple sequence alignment using the specified method
        if method == "clustalomega":
            clustalomega_cline = ClustalOmegaCommandline(infile=input_file, outfile=output_file, verbose=True, auto=True)
            stdout, stderr = clustalomega_cline()
        elif method == "mafft":
            mafft_cline = MafftCommandline(input=input_file)
            stdout, stderr = mafft_cline()
        else:
            raise ValueError("Invalid alignment method specified")

        # Print the alignment results
        print(stdout)

        # Print the error message, if any
        print(stderr)

        # Save the alignment results to a file
        with open(output_file, "w") as f:
            f.write(stdout)
    except Exception as e:
        print(f"An error occurred during multiple sequence alignment: {str(e)}")

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Perform multiple sequence alignment")
    parser.add_argument("-i", "--input", help="Input file path", required=True)
    parser.add_argument("-o", "--output", help="Output file path", required=True)
    parser.add_argument("-m", "--method", help="Alignment method (clustalomega, mafft)", required=True)
    args = parser.parse_args()

    # Perform multiple sequence alignment
    perform_msa(args.input, args.output, args.method)