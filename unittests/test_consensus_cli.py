##################################################
#            Test Description                    #
##################################################
# This file contains the unit tests for the      #
# consensus CLI response.                        #
##################################################

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../')))

import unittest
from unittest.mock import patch, mock_open
from click.testing import CliRunner

from substrate_miner.consensus import consensus 
import sys
import os
    
class TestConsensusFunctions(unittest.TestCase):

    def setUp(self):
        self.runner = CliRunner()

    @patch('substrate_miner.consensus.consensus.consensus_click')
    @patch('sys.argv', ['substrate_miner', 'consensus', 'unittests/data/test-consensus.fas', 'outputs', '-c', '1'])
    def test_consensus_cli(self, mock_consensus_click):
        consensus.consensus_click(sys.argv)
        mock_consensus_click.assert_called_once()
        
    @patch('substrate_miner.consensus.consensus.consensus_main')
    @patch('sys.argv', ['consensus'])
    def test_consensus_main(self, mock_consensus_main):
        consensus.consensus_main()
        mock_consensus_main.assert_called_once()
  
    @patch('substrate_miner.consensus.consensus.consensus_main')
    @patch('sys.argv', ['weblogo'])
    def test_consensus_main_logo(self, mock_consensus_main):
        consensus.consensus_main()
        mock_consensus_main.assert_called_once()  

    def test_consensus(self):

        consensus.consensus('unittests/data/test-consensus.fas', '', 'unittests', '1', 0.5, False)    
        
        assert os.path.exists('unittests/test-consensus_sequenceEntropies.png')
        assert os.path.exists('unittests/test-consensus_gapFrequencies.png')
        
        with open('unittests/test-consensus_gapStrip.txt', 'r') as outfile:
            output_content = outfile.read()
        with open('unittests/results/test-consensus_gapStrip.txt', 'r') as expoutfile:
            expected_content = expoutfile.read()
        assert output_content == expected_content

if __name__ == '__main__':
    unittest.main()