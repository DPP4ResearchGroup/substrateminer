from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline, MafftCommandline, MuscleCommandline

import argparse
import click
import sys
import shutil
import os


def perform_msa(input_file, output_file, method):
    """
    Perform multiple sequence alignment using ClustalOmega, MAFFT and MUSCLE.

    Parameters:
    - input_file (str): Path to the input file containing sequences in FASTA format.
    - output_file (str): Path to the output file where the alignment results will be saved.
    - method (str): Alignment method to use. Options are "clustalomega" or "mafft".

    Raises:
    - ValueError: If an invalid alignment method is specified.

    Returns:
    - None

    Other Bio.Align.Applications options:
    - EmbossCommandline: This class provides a command line wrapper for the EMBOSS suite of
        sequence analysis tools. It can be used to perform multiple sequence alignment using
        EMBOSS tools like "needle" or "water".

    Note:
    - The input file should contain sequences in FASTA format.
    - The alignment results will be printed to the console and saved to the output file.
    """
    
    try:
        # Read input sequences from a file
        sequences = SeqIO.parse(input_file, "fasta")

        # Check if the output file exists and remove it if it does
        if os.path.exists(output_file):
            os.remove(output_file)

        # Perform multiple sequence alignment using the specified method
        if method == "clustalomega":
            clustalomega_cline = ClustalOmegaCommandline(infile=input_file,\
                outfile=output_file, verbose=True, auto=True)
            stdout, stderr = clustalomega_cline()
        elif method == "mafft":
            mafft_exe = shutil.which("mafft")
            if mafft_exe is None:
                raise FileNotFoundError("MAFFT executable not found.\
                    Please ensure it is installed and in your PATH.")
            mafft_cline = MafftCommandline(mafft_exe, input=input_file)
            stdout, stderr = mafft_cline()
            with open(output_file, "w") as mafft_w_handle:
                mafft_w_handle.write(stdout)
        elif method == "muscle":
            muscle_exe = shutil.which("muscle")
            if muscle_exe is None:
                raise FileNotFoundError("MUSCLE executable not found.\
                    Please ensure it is installed and in your PATH.")
            muscle_cline = MuscleCommandline(muscle_exe, input=input_file,\
                out=output_file)
            stdout, stderr = muscle_cline()
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


@click.command('msa', context_settings=dict(ignore_unknown_options=True))
@click.argument('args', nargs=-1, type=click.UNPROCESSED)
def msa_switch(args):
    """
    Interface for multi-sequence alignments
    """
    
    msa_click(sys.argv[2:])

def msa_click(args):
    # debug block
    #print("%%%%%%%%%%%%%%%%%%%%%%%  DEBUG  %%%%%%%%%%%%%%%%%%%%%%%")
    #print("click pathway")
    #print(args)
    #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    
    if not any(args):
        sys.argv = ['msa.py', '-h']
        
    sys.argv = ['msa.py'] + list(args)
    # Perform multiple sequence alignment
    msa_main()
    
def msa_main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Perform multiple sequence alignment")
    parser.add_argument("-i", "--input", help="Input file path", required=True)
    parser.add_argument("-o", "--output", help="Output file path", required=True)
    parser.add_argument("-m", "--method", help="Alignment method (clustalomega, mafft, muscle)", required=True)
    args = parser.parse_args()

    # Perform multiple sequence alignment
    perform_msa(args.input, args.output, args.method)

if __name__ == "__main__":
    if sys.argv[1:] and not any(item in ['-h', '--help'] for item in sys.argv[1:]):
        # debug block
        print("%%%%%%%%%%%%%%%%%%%%%%%  DEBUG  %%%%%%%%%%%%%%%%%%%%%%%")
        print("click pathway")
        #print(sys.argv[1:])
        #print([item for item in sys.argv[1:]])
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        msa_click(sys.argv[1:])
    else:
        msa_main()