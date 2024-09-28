#!/usr/bin/env python
###########################################################################################################
# Project: substrateminer
#        Amino Acid Sequence Filter 
#        This script is used to filter amino acid sequences from a reference file.
#           by taking a reference file and amino acid targets as input.
#           Then filters the amino acid sequences from the reference file based on the target files.
###########################################################################################################  

# Import required libraries
import click
import yaml
import sys
import os
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from miner import miner_funcs as modulePepProcess
from records.record_class import motifFeatures, filterOutput

from Bio import SwissProt
from Bio import SeqIO

## Aux Functions
def convert_bytes(size_in_bytes):
    if size_in_bytes >= 1024 ** 3:  # If size is in GB
        size_in_gb = size_in_bytes / (1024 ** 3)
        return f"{size_in_gb:.2f} GB"
    elif size_in_bytes >= 1024 ** 2:  # If size is in MB
        size_in_mb = size_in_bytes / (1024 ** 2)
        return f"{size_in_mb:.2f} MB"
    elif size_in_bytes >= 1024:  # If size is in KB
        size_in_kb = size_in_bytes / 1024
        return f"{size_in_kb:.2f} KB"
    else:
        return f"{size_in_bytes} Bytes"
    
## Record Handle Switch
def record_handle_switch(handle, referencetype, referencefile):
    record_switch = {
        "embl": case_embl,
        "swiss": case_swiss,
        "genbank": case_genbank,   
    }
    
    records = record_switch.get(referencetype)(handle, referencefile)
    return records

def case_swiss(handle, referencefile):
    records = list(SwissProt.parse(handle))
    record_size_handle = convert_bytes(sys.getsizeof(records))
    
    print("++++++++++++++++++++++++++++++++++++++++++++")
    print("Done Big Reference List Parsing with SwissProt")
    print(f"Size of variable records: {record_size_handle}")
    print("++++++++++++++++++++++++++++++++++++++++++++")
    
    record_swiss = SeqIO.index(referencefile, "swiss")
    return records, record_swiss

def case_genbank(handle, referencefile):
    records = list(SeqIO.parse(handle, "genbank"))
    record_size_handle = convert_bytes(sys.getsizeof(records))

    print("++++++++++++++++++++++++++++++++++++++++++++")
    print("Done Big Reference List Parsing with GenBank")
    print(f"Size of variable records: {record_size_handle}")
    print("++++++++++++++++++++++++++++++++++++++++++++")
    
    return records

def case_embl(handle, referencefile):
    records = list(SeqIO.parse(handle, "embl"))
    record_size_handle = convert_bytes(sys.getsizeof(records))
    
    print("++++++++++++++++++++++++++++++++++++++++++++")
    print("Done Big Reference List Parsing with EMBL")
    print(f"Size of variable records: {record_size_handle}")
    print("++++++++++++++++++++++++++++++++++++++++++++")

def validate_output(ctx, param, value):
    if ctx.params.get('outmode') != 'inline' and not value:
        raise click.BadParameter('outputfilename must be provided if outmode is not inline')
    return value
def validate_outputfiletype(ctx, param, value):
    if ctx.params['outputfilename'] and not value:
        return 'fasta'
    return value

@click.command('miner')
@click.option('--referencefile', required=True, help='The reference file containing sequences.')
@click.option('--referencetype', default="swiss", type=click.Choice(["swiss", "genbank", "embl"]),\
    help='The type of reference file. Default is swiss.')
@click.option('--filtermode', required=True, type=click.Choice(["size", "motif", "loc"]),\
    help='The mode of filtering. Default is size.')
@click.option("--config", default="config.yml",\
    help="The path to the configuration file.")
@click.option("--stats", is_flag=True,\
    help="Generate statistics for the filtered sequences.")
@click.option("--outmode", default="all", type=click.Choice(["all", "file", "inline"]),\
    help="The output mode for the filtered sequences.")
@click.option("--outputfilename", callback=validate_output,\
    help="The output file name for the filtered sequences.")
@click.option("--outputfiletype", callback=validate_outputfiletype,\
    type=click.Choice(["fasta", "text", "txt", "genbank", "swiss"]),\
    help="The output file type for the filtered sequences.")
def miner_main(referencefile, referencetype, filtermode, config, stats,\
    outmode, outputfilename, outputfiletype) -> list[filterOutput]:
    """
    Filter amino acid sequences from a reference file.
    
    referencefile (str): Path to the reference file containing sequences.
    referencetype (str): Type of the reference file (e.g., fasta, fastq).
    filtermode (str): Mode of filtering. Can be 'size', 'motif', or 'loc'.
    config (str): Path to the configuration file in YAML format.
    outmode (str): Mode for output generation.
    outputfilename (str): Name of the output file.
    outputfiletype (str): Type of the output file.
    
    ValueError: If an invalid filter mode is provided.
    SystemExit: If there is an invalid configuration for the specified filter mode.
    """
    
    try:
        out_stack = seq_filter(referencefile, referencetype, filtermode, config, stats, outmode, outputfilename, outputfiletype)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except SystemExit as e:
        sys.exit(1)
        
    return out_stack
    
def seq_filter(referencefile, referencetype, filtermode, config, stats,\
    outmode, outputfilename, outputfiletype) -> list[filterOutput]:
    """
    Filter amino acid sequences from a reference file.
    
    Args:
        referencefile (str): Path to the reference file containing sequences.
        referencetype (str): Type of the reference file (e.g., fasta, fastq).
        filtermode (str): Mode of filtering. Can be 'size', 'motif', or 'loc'.
        config (str): Path to the configuration file in YAML format.
        outmode (str): Mode for output generation.
        outputfilename (str): Name of the output file.
        outputfiletype (str): Type of the output file.
    
    Raises:
        ValueError: If an invalid filter mode is provided.
        SystemExit: If there is an invalid configuration for the specified filter mode.
    
    Returns:
        list[filterOutput]: A list of filterOutput objects containing the filtered sequences.
    """
    
    # configure the module
    # read configure file
    with open(config, 'r') as cf:
        config = yaml.safe_load(cf)
    
    # configurations
    size_type_stack = []
    size_value_stack = []
    
    output_base = config['ops']['output_dir']
    # Check if output_base path exists, if not create it
    try:
        if not os.path.exists(output_base):
            os.makedirs(output_base)
    except OSError as e:
        print(f"Error creating directory {output_base}: {e}")
        sys.exit(1)
    
    # output storage
    out_size_stack = {}
    out_loc_stack = {}
    out_motif_stack = {}
    
    # intrim storage stack
    out_step_stack = []
    
    # filter mode switch
    if filtermode == "size":
        try:
            size_stack = config['size']
            for size_stack_item in size_stack:
                size_type = size_stack_item['type']
                size_value = size_stack_item['value']
                if size_type not in ["max", "min", "equal"]:
                    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                    print("Invalid size type for 'size', please check the configuration file.")
                    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                    sys.exit(1)
                size_type_stack.append(size_type)
                size_value_stack.append(size_value)
        except KeyError:
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            print("Invalid size configuration for 'size', please check the configuration file.")
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            sys.exit(1)

    elif filtermode == "motif":
        try:
            motif_type = config['motif']['type']
            motif_consensus = config['motif']['consensus']
            motif_target_type = config['motif']['target_motif_type']
            if motif_type == "exo":
                motif_terminal = config['motif']['terminal']
            elif motif_type == "endo":
                motif_terminal = "X"
                
        except KeyError:
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            print("Invalid motif configuration for 'motif', please check the configuration file.")
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            sys.exit(1)
        
        # Sanity - Ensure the motif contains only valid amino acids
        valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWYX")
        for motif in motif_consensus:
            if not all(aa in valid_amino_acids for aa in motif):
                raise ValueError("The target motif contains invalid amino acids. For unknown amino acids, use 'X'.")

        valid_motif_type = set(["exo", "endo"])
        valid_motif_terminal = set(["N", "C", "X"])
        if motif_type not in valid_motif_type:
            raise ValueError("Invalid motif type, please check the configuration file.")
        if motif_terminal not in valid_motif_terminal:
            raise ValueError("Invalid motif terminal, please check the configuration file.")
        
    elif filtermode == "loc":
        try:
            cellular_loc = config['cellular_location']
        except KeyError:
            print("Invalid cellular location configuration 'cellular_loc', please check the configuration file.")
            sys.exit(1)
    else:
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        raise ValueError("Invalid filter mode, please check the configuration file.")
        
    # setup the full records
    handle = open(referencefile, "r")
    if referencetype == "swiss":
        ## TODO: homo records extraction routine
        records, record_handle = record_handle_switch(handle, referencetype, referencefile)
    else:
        records = record_handle_switch(handle, referencetype, referencefile)   
    handle.close()

    ## debug/testing
    # print (size_type)
          
    # Generating File Cluster for each DPP Consensus search:
    full_records = records
    
    # Generate Target Sequence List
    target_seq = []
    
    outmode_switch = {
        "all": "all",
        "file": "file",
        "inline": "inline"
    }

    if filtermode == "size":
        intrim_records = full_records
        size_stack_counter = 0
        
        if stats:
            size_stats = True
        else:
            size_stats = False

        for size_stack_item in size_stack:
            local_min = size_value_stack[size_stack_counter]
            target_size_min_switch = {
                "max": 0,
                "min": local_min,
                "equal": local_min,
            }
            local_max = size_value_stack[size_stack_counter]
            target_size_max_switch = {
                "max": local_max,
                "min": 0,
                "equal": local_max,
            }
            
            out_size_stack = modulePepProcess.seperate_target_size(records = intrim_records,\
                target_size_max = target_size_max_switch[size_stack_item['type']], target_size_min = target_size_min_switch[size_stack_item['type']],\
                output_base = output_base, output_format = outmode_switch[outmode], output_file_name = outputfilename, output_file_format = outputfiletype,\
                stats = size_stats)
            intrim_records = out_size_stack['target']
            out_step_stack.append(out_size_stack)
            
            # stage prompts
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            print(f"Filter: {filtermode}, Stage: {size_stack_counter + 1}")
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            
            # debug block
            #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            #print(size_stack_item['type'])
            #print(target_size_max_switch[size_stack_item['type']])
            #print(out_size_stack['nontarget'][0].features[0].location.end)
            #print(len(out_size_stack['target']))
            #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            
            size_stack_counter += 1

        # stage prompts
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        print(f"Filter: {filtermode} completed. Filtered records: {len(out_size_stack['target'])}")
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        
        # Debug Block
        #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        #print(len(out_size_stack['target']))
        #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        out_stack = out_step_stack
        
    elif filtermode == "loc":
        intrim_records = full_records
    
        if stats:
            loc_stats = True
        else:
            loc_stats = False
            
        for location in cellular_loc:
            print(location, outmode_switch[outmode])
            out_loc_stack = modulePepProcess.seperate_cellular_location(records = intrim_records, target_orgin = location,\
                output_base = output_base, output_format = outmode_switch[outmode], output_file_name = outputfilename, output_file_format = outputfiletype,\
                stats = loc_stats)
            
            out_step_stack.append(out_loc_stack)
            # stage prompts
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            print(f"Filter: {filtermode}, Stage: {location}")
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    
            # debug block   
            #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            #print(len(out_loc_stack['target']))
            #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    
        # stage prompts
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        print(f"Filter: {filtermode} completed. Filtered records: {len(out_loc_stack['target'])}")
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        
        out_stack = out_step_stack
        
    elif filtermode == "motif":
        intrim_records = full_records
        
        if stats:
            motif_stats = True
        else:
            motif_stats = False
        
        for motif in motif_consensus:
            motif_target_features = motifFeatures(type = motif_type, motif = motif, terminal = motif_terminal)

            out_motif_stack = modulePepProcess.seperate_target_motif(records = intrim_records,\
                target_motif_features = motif_target_features,\
                output_base = output_base, output_format = outmode_switch[outmode], output_file_name = outputfilename, output_file_format = outputfiletype,\
                stats = motif_stats)

            out_step_stack.append(out_motif_stack)
            
            # stage prompts
            print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
            print(f"Filter: {filtermode} completed. Filtered records: {len(out_motif_stack['target'])}")
            print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
            
            # debug block
            #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            #print(len(out_motif_stack['target']))
            #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

        out_stack = out_step_stack
        
    else:
        raise ValueError("Invalid filter mode, please check the configuration file.")

    # return output conformation
    return_out_stack = []
    for each_stack in out_step_stack:
        return_out = filterOutput(each_stack['target'], each_stack['nontarget'])
        return_out_stack.append(return_out)
    
    return return_out_stack


if __name__ == "__main__":
    miner_main()