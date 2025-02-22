##################################################
#            Test Description                    #
##################################################
# This file contains the unit tests for the      #
# MSA CLI response.                              #
##################################################

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../')))

import unittest
from unittest.mock import patch, MagicMock
from substrate_miner.msa import msa

class TestMSACli(unittest.TestCase):

    @patch('substrate_miner.msa.msa.msa_click')
    @patch('sys.argv', ['-i', '/data/msa.fasta', '-o', '/output/msa_output.aln', '-m', 'clustalomega'])
    def test_msa_click(self, mock_msa_click):
        msa.msa_click(sys.argv)
        mock_msa_click.assert_called_once()

    @patch('substrate_miner.msa.msa.msa_main')
    @patch('sys.argv', ['msa'])
    def test_msa_main(self, mock_msa_main):
        msa.msa_main()
        mock_msa_main.assert_called_once()

    @patch('substrate_miner.msa.msa.msa_main')
    @patch('sys.argv', ['msa', '-h'])
    def test_msa_main_help(self, mock_msa_main):
        msa.msa_main()
        mock_msa_main.assert_called_once()

    @patch('substrate_miner.msa.msa.perform_msa')
    @patch('argparse.ArgumentParser.parse_args',\
        return_value=MagicMock(input="/data/msa.fasta", output="/output/msa_output.aln", method="clustalomega"))
    def test_msa_main_args(self, mock_parse_args, mock_perform_msa):
        msa.msa_main()
        mock_perform_msa.assert_called_once()
        
if __name__ == '__main__':
    unittest.main()