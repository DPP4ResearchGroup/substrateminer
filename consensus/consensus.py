import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
import csv
import subprocess
import consensus_funcs as modules

def consensus(inputfile: str, outputfile: str, outputdir: str, consensus_option: int, threshold: float, save_figures):
    """
    Determines the consensus sequence from a multiple sequence alignment (MSA) file.
    Args:
        inputfile (str): Path to the input MSA file.
        outputfile (str): Path to the output file for the consensus sequence.
        outputdir (str): Path to the output directory for the generated files.
        consensus_option (int): Method for determining consensus positions:
            1: Positions with gap frequencies < threshold (0.5 default, change with -t flag)
            2: Positions with residue as most frequent character
            3: Positions with residues in a specific sequence (you will give sequence ID)
        threshold (float): Gap frequency threshold for consensus_option 1. Must be between 0 and 1.
        save_figures (bool): Flag indicating whether to save figures.
    Returns:
        None
    Raises:
        FileNotFoundError: If the input file is not found.
        ValueError: If the consensus_option is not 1, 2, or 3.
        ValueError: If the threshold is not between 0 and 1.
    Notes:
        - The input file must be in FASTA format.
        - The sequences in the input MSA must be aligned.
        - The consensus sequence is determined based on the specified consensus_option.
        - The generated files include the consensus sequence, residue frequencies, sequence entropies, and gap stripped alignment.
        - Figures of sequence entropies and gap frequencies can be saved if save_figures is True.
    """
    
    #FIGURE SETTINGS
    mpl.rcParams['axes.titlesize'] = 18
    mpl.rcParams['axes.labelsize'] = 18
    mpl.rcParams['xtick.labelsize'] = 14
    mpl.rcParams['ytick.labelsize'] = 14
    mpl.rcParams['axes.facecolor'] = 'FFFFFF'
    mpl.rcParams['axes.edgecolor'] = '000000'
    mpl.rcParams['axes.linewidth'] = 1.0
    mpl.rcParams['axes.labelweight'] = 'regular'
    mpl.rcParams['xtick.major.pad'] = 3
    mpl.rcParams['ytick.major.pad'] = 3
    plt.rcParams['font.family'] = 'sans-serif'
    
    print()
    print(f'Reading file: {inputfile}')

    #read FASTA and populate lists for sequences and IDs
    seqs = list()
    names = list()
    
    #flag to strip gaps when cleaning sequence
    strip_gaps = False
    
    try:
        with open(inputfile, 'r') as n:
            seqs, ids = modules.read_fasta(n, strip_gaps)
    except FileNotFoundError:
        print(f'Could not find file: {inputfile}')
        sys.exit(1)
        
    # If no sequences were added to list, file was not in FASTA format
    if len(seqs) < 2:
        print('Provided file is not in FASTA format.')
        sys.exit(1)
        
    # Ensures sequences are aligned (all sequences have same length)
    if not modules.is_fasta_aligned(seqs):
        print('Sequences in provided FASTA are not aligned')
        sys.exit(1)

    num_seqs = len(seqs)

    # standard 20 amino acid alphabet
    res_list = list('ACDEFGHIKLMNPQRSTVWY-')

    # Matrix for residue frequencies at each position
    marginals = modules.marginal_frequencies(seqs, res_list)

    # Calculating consensus sequences for method based on user input
    consensus_positions = list()
    consensus_sequence = list()

    # Get user input for how to determine positions to include from MSA (handling gaps)
    # Users must eneter '1', '2', or '3'
    # Exits script if they do not enter valid input in five tries
    consensus_choice = consensus_option
    n_tries = 5
    while str(consensus_choice) not in ['1', '2', '3'] and n_tries > 0:
        consensus_choice = str(input('\nWhich method for determining consensus positions do you want?\n\
            ***Note: Option 1 recommended.***\n\
                1: Positions with gap frequencies < threshold (0.5 default, change with -t flag)\n\
                2: Positions with residue as most frequent character\n\
                3: Positions with residues in a specific sequence (you will give sequence ID)\n'))
        n_tries -= 1
        if n_tries == 0:
            print('Invalid responses. Must enter 1, 2, or 3. Exiting program...')
            sys.exit(1)
    print()

    # Includes all positions for which gap frequency < 0.5
    # In other words, the sum of frequencies of 20 residues is > 0.5
    # This option can include positions where a gap is the most frequent
    if consensus_choice == '1':
        threshold = threshold
        if threshold < 0 or threshold > 1:
            print('Gap frequency threshold must be between 0 and 1')
            sys.exit(1)
    
        for i, j in enumerate(marginals):
            if j[-1] < threshold:
                max_res = np.argmax(j[:-1])
                consensus_positions.append(i)
                consensus_sequence.append(res_list[max_res])
                
    # Includes all positions for which a residue is the most frequent
    # i.e., eliminates positions for which the most frequent occurrence is a gap
    elif consensus_choice == '2':
        for i, j in enumerate(marginals):
            max_res = np.argmax(j)
            if max_res != len(j) - 1:
                consensus_positions.append(i)
                consensus_sequence.append(res_list[max_res])
                
    # Includes all positions occupied by residues in a user-defined reference sequence
    else:
        #Take user input for sequence ID of reference sequence
        ref_seq = input('What is the ID of your reference sequence? ')
        
        if len(ref_seq) == 0:
            print('\nNo sequence ID was given...')
            sys.exit(1)
        
        # Users may not give sequence ID with leading '>' character
        # IDs from FASTA all have leading '>' character, so must add to reference ID
        if ref_seq[0] != '>':
            ref_seq = ''.join(['>', ref_seq])
            
        try:
            ref_index = ids.index(ref_seq)
        except ValueError:
            print('\nCould not find sequence ID in set. Check your MSA for the correct ID...')
            sys.exit(1)
        
        for i, (j, k) in enumerate(zip(seqs[ref_index], marginals)):
            if j != '-':
                consensus_positions.append(i)
                max_res = np.argmax(k[:-1])
                consensus_sequence.append(res_list[max_res])
                
    consensus_sequence = ''.join(consensus_sequence)
    print(f'Consensus sequence: {consensus_sequence}')
    print()

    # Removes positions with high gap frequencies from all sequences in alignment
    # Calculates residue frequencies for gap stripped alignment
    seqs_gap_stripped = list()
    for i in seqs:
        seqs_gap_stripped.append(''.join([i[j] for j in consensus_positions]))
    marginals_gap_stripped = modules.marginal_frequencies(seqs_gap_stripped, res_list)

    # Calculates sequence entropies for all poisitions in the gap stripped alignment
    seq_entropies = modules.seq_entropy(marginals_gap_stripped)

    # Creates matrix of residue frequencies at all positions in gap stripped alignment
    # in a format of CSV
    out_marginals = []
    for i, j in enumerate(marginals_gap_stripped):
        if i == 0:
            out_marginals.append([''] + res_list)
            out_marginals.append([f'Position {i+1}'] + list(j))
        else:
            out_marginals.append([f'Position {i+1}'] + list(j))
            
    # Calculating gap frequencies
    # Gap is last element in each row of marginal frequencies
    gap_frequencies = [i[-1] for i in marginals_gap_stripped]

    # Setting path for file outputs
    #set default path to inputfile directory
    path = os.path.dirname(os.path.abspath(inputfile)) 

    out_file_prefix = outputfile
    # Intermediate output stack
    if not out_file_prefix:
        out_file_prefix = os.path.splitext(os.path.basename(inputfile))[0]
    out_file_dir = outputdir
    if not out_file_dir:
        out_file_dir = os.path.dirname(os.path.abspath(inputfile))
    
    print (out_file_prefix)
    
    # Save figures if -f flag is not given
    if not save_figures:
        #Save figure of alignment sequence entropies
        fig, ax = plt.subplots()
        ax.stem(np.arange(1, len(seq_entropies) + 1), seq_entropies, use_line_collection = True)
        ax.set_xlabel('Sequence position')
        ax.set_ylabel('Sequence entropy (bits)')
        ax.set_ylim(0, 4.5)
        fig.savefig(f'{out_file_dir}/{out_file_prefix}_sequenceEntropies.png', bbox_inches = 'tight', dpi = 300)
        plt.close()
        
        #Save figure of gap stripped alignment gap frequencies
        fig, ax = plt.subplots()
        ax.stem(np.arange(1, len(gap_frequencies) + 1), gap_frequencies, use_line_collection = True)
        ax.set_xlabel('Sequence position')
        ax.set_ylabel('Gap frequency')
        ax.set_ylim(0, 1)
        fig.savefig(f'{out_file_dir}/{out_file_prefix}_gapFrequencies.png', bbox_inches = 'tight', dpi = 300)
        plt.close()
                
    #Create a FASTA of the gap stripped alignment
    out_fasta = modules.make_fasta(ids, seqs_gap_stripped)

    #Save files
    with open(f'{out_file_dir}/{out_file_prefix}_residueFrequencies.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(out_marginals)
        
    #Write summary of analysis to file
    with open(f'{out_file_dir}/{out_file_prefix}_consensus_output.txt', 'w', newline='') as f:
        f.write(f'Determined consensus sequence for: {inputfile}\n\n')
        f.write('Parameters used:\n')
        if consensus_choice ==  1:
            f.write(f'Method for removing insertions: {consensus_choice}\n')
            f.write(f'Gap frequency threshold: {threshold}\n\n')
        else:
            f.write(f'Method for removing insertions: {consensus_choice}\n\n')
        f.write(f'>Consensus_sequence\n{consensus_sequence}\n\n')
        f.write(f'MSA sequence entropy per residue: {np.mean(seq_entropies):.3f}\n\n')
        f.write(f'Wrote gap stripped alignment to: {out_file_prefix}_gapStrip.txt\n\n')
        f.write(f'Wrote CSV of residue frequencies to: {out_file_prefix}_residueFrequencies.csv\n\n')
        if save_figures:
            f.write(f'Wrote plot of sequence entropies to: {out_file_prefix}_sequenceEntropies.png\n\n')
            f.write(f'Wrote plot of gap frequencies to: {out_file_prefix}_gapFrequencies.png\n\n')

    with open(f'{out_file_dir}/{out_file_prefix}_gapStrip.txt', 'w') as f:
        f.writelines(out_fasta)
        
    print(f'Wrote summary of output to: {out_file_dir}/{out_file_prefix}_consensus_output.txt')

def weblogo(inputfile: str, outputfile: str, resolution: int, filetype: str):
    """
    Generates a weblogo image using the specified input file, output file, resolution, and file type.

    Parameters:
    - inputfile (str): The path to the input file.
    - outputfile (str): The path to save the generated weblogo image.
    - resolution (int): The resolution of the weblogo image.
    - filetype (str): The file type of the output image.

    Returns:
    None
    """
    
    p_subscript = 'Protein_Sequence_Consensus'
    weblogo_cline = f'weblogo -f {inputfile} -o {outputfile} -F {filetype} --resolution {resolution} -A protein -P {p_subscript}'
    try:
        subprocess.run(weblogo_cline, shell=True)
    except Exception as e:
        print(f'Error occurred while running weblogo: {e}')
    

if __name__ == '__main__':
        
    parser = argparse.ArgumentParser(description='Determine consensus sequence from a multiple sequence alignment (MSA) and draw summative plots.')
    subparsers = parser.add_subparsers(dest='command', help='sub-command help')
    subparsers.required = True

    # Create the parser for subcommand 'consensus'
    parser_consensus = subparsers.add_parser('consensus', help='Determine consensus sequence from a multiple sequence alignment (MSA) and draw sequence entropy and gap frequency plots.')
    parser_consensus.add_argument('-i', required=True, metavar = 'Input alignment file in FASTA format.', help='Filename for FASTA alignment')
    parser_consensus.add_argument('-o', metavar = 'Output gap stripped FASTA file name', help='Output FASTA filename. If not given will use name of input FASTA file as template to name output files.')
    parser_consensus.add_argument('-O', metavar = 'Output directory', help='Output directory for all output files. If not given will use directory of input FASTA file.')
    parser_consensus.add_argument('-c', default='0', metavar = 'Method for removing insertions',\
        help='Desired method for removing insertions.\
                1 = Positions with gap frequencies < threshold (0.5 default, change with -t flag).\
                2 =  Positions with residue as most frequent character.\
                3 = Positions with residues in a specific sequence. If not given will ask for user input upon running script. See README for further explantion of methods.'\
    )
    parser_consensus.add_argument('-t', type=float, default = 0.5, metavar = 'Gap frequency threshold', help='Gap frequecy threshold to define a consensus positions. Only valid for Option 1 for removing insertions. Must be a value between 0 and 1 (default: 0.5)')
    parser_consensus.add_argument('-f', action='store_true', help='Include flag to prevent saving images of MSA data analysis.')
    parser_consensus.set_defaults(func=lambda args:\
        consensus(inputfile = args.i, outputfile = args.o, outputdir = args.O, consensus_option = args.c, threshold = args.t, save_figures = args.f))
    
    # Create the parser for subcommand 'weblogo'
    parser_weblogo = subparsers.add_parser('weblogo', help='Generate a weblogo image from an input file.')
    parser_weblogo.add_argument('-i', required=True, help = 'Input alignment file/self-aligned file in FASTA/text format.', metavar='INPUTFILE')
    parser_weblogo.add_argument('-o', required=True, help = 'Output filename for weblogo image.', metavar='FILENAME')
    parser_weblogo.add_argument('-s', default=300, help = 'Resolution of the weblogo image.', metavar='RESOLUTION')
    parser_weblogo.add_argument('-F', default='png', choices=['png', 'jpeg', 'svg', 'pdf'], help = 'File type of the output image.', metavar='FILETYPE')
    parser_weblogo.set_defaults(func=lambda args:\
        weblogo(inputfile = args.i, outputfile = args.o, resolution = args.s, filetype = args.F))
    
    
    # Parse the arguments and call the appropriate function
    args = parser.parse_args()
    args.func(args)
    
    #Show help if no arguments passed
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)