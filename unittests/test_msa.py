from substrate_miner.msa import msa
import unittest
from unittest.mock import patch, mock_open
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../')))

class TestMSA(unittest.TestCase):

    def test_perform_msa(self):
        with patch('subprocess.Popen') as mock_popen:
            mock_popen.return_value.communicate.return_value = (b'output', b'error')
            
            msa.perform_msa(os.path.abspath('unittests/data/msa.fasta'), os.path.abspath('outputs/output.aln'), 'clustalomega')
            mock_popen.assert_called_once()
    
    @patch('subprocess.Popen')
    def test_perform_msa_exception(self, mock_popen):
        mock_popen.side_effect = Exception()
        msa.perform_msa('unittests/data/msa.fasta', 'outputs/output.aln', 'clustalomega')
        
    @patch('subprocess.Popen')
    @patch('builtins.open', new_callable=mock_open)
    def test_perform_msa_output(self, mock_open, mock_popen):
        mock_popen.return_value.communicate.return_value = (b'output', b'error')
        msa.perform_msa(os.path.abspath('unittests/data/msa.fasta'), os.path.abspath('outputs/output.aln'), 'clustalomega')
        mock_open.assert_called_once()
        
    def test_perform_msa_mafft(self):
        with patch('subprocess.Popen') as mock_popen:
            mock_popen.return_value.communicate.return_value = (b'output', b'error')
            msa.perform_msa(os.path.abspath('unittests/data/msa.fasta'), os.path.abspath('outputs/output.aln'), 'mafft')
            mock_popen.assert_called_once()
    
    @patch('subprocess.Popen')
    def test_perform_msa_mafft_exception(self, mock_popen):
        mock_popen.side_effect = Exception()
        msa.perform_msa(os.path.abspath('unittests/data/msa.fasta'), os.path.abspath('outputs/output.aln'), 'clustalomega')
        
    @patch('subprocess.Popen')
    @patch('builtins.open', new_callable=mock_open)
    def test_perform_msa_mafft_output(self, mock_open, mock_popen):
        mock_popen.return_value.communicate.return_value = (b'output', b'error')
        msa.perform_msa(os.path.abspath('unittests/data/msa.fasta'), os.path.abspath('unittests/output.aln'), 'mafft')
        mock_open.assert_called_once()

    def test_perform_msa_performance(self):
        msa.perform_msa(os.path.abspath('unittests/data/msa.fasta'), os.path.abspath('unittests/msa_output.aln'), 'clustalomega')
        assert os.path.exists('unittests/msa_output.aln')

    def test_perform_msa_output_comparison(self):
        msa.perform_msa(os.path.abspath('unittests/data/msa.fasta'), os.path.abspath('unittests/msa_output.aln'), 'mafft')
        with open('unittests/msa_output.aln', 'r') as output_file:
            output_content = output_file.read()
        with open('unittests/results/msa_mafft.aln', 'r') as expected_file:
            expected_content = expected_file.read()
        self.assertEqual(output_content, expected_content)
        
if __name__ == '__main__':
    unittest.main()