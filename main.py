#############################################
# This is the main entry point for the CLI. #
#############################################

import click

import consensus.consensus as consensus
import msa.msa as msa
import miner.miner as miner
import pathfinder.pathfinder as pathfinder

@click.group()
def substrateminer():
    """
    A suite to tools to discover enzyme substrates based on sequence consensus.
    
    Main entry point for substrateminer CLI.
    """
    pass

substrateminer.add_command(miner.miner_main)
substrateminer.add_command(msa.msa_switch)
substrateminer.add_command(consensus.consensus_switch)
substrateminer.add_command(pathfinder.pathfinder_main)


if __name__ == '__main__':
    substrateminer()
