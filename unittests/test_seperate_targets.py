##################################################
#            Test Description                    #
##################################################
# This file contains the unit tests for the      #
# separate_target_size function.                 #
##################################################

import sys
import yaml
import os
import unittest
from io import StringIO
from unittest.mock import patch

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../')))
from substrate_miner.miner import miner


class TestSeparateTargets(unittest.TestCase):

    # Example test case with valid input
    referencefile = os.path.join('unittests/data/test-uniprot.txt')
    referencetype = "swiss"
    filtermode = "size"
    config = os.path.abspath(os.path.join(os.path.dirname(__file__), 'test-config.yml'))
    outmode = "inline"
    outputfilename = "unitest-out"
    outputfiletype = 'fasta'
    stats = False

    with open(config, 'r') as file:
        test_config = yaml.safe_load(file)

    def test_separate_target_size_valid_input(self):
        # example of a valid input
        out_stack = miner.seq_filter(self.referencefile, self.referencetype, self.filtermode, self.config, self.stats, self.outmode, self.outputfilename, self.outputfiletype)
        result_seqs = [record.sequence for record in out_stack[0].target]
        result = result_seqs[0:5]

        # debug block
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        print("Test mode: ", self.filtermode)
        print(result)
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

        expected_result = [
            'MFRRLTFAQLLFATVLGIAGGVYIFQPVFEQYAKDQKELKEKMQLVQESEEKKS',
            'MRWQEMGYIFYPRKLR',
            'MTQRAGAAMLPSALLLLCVPGCLTVSGPSTVMGAVGESLSVQCRYEEKYKTFNKYWCRQPCLPIWHEMVETGGSEGVVRSDQVIITDHPGDLTFTVTLENLTADDAGKYRCGIATILQEDGLSGFLPDPFFQVQVLVSSASSTENSVKTPASPTRPSQCQGSLPSSTCFLLLPLLKVPLLLSILGAILWVNRPWRTPWTES',
            'MGDQPCASGRSTLPPGNAREAKPPKKRCLLAPRWDYPEGTPNGGSTTLPSAPPPASAGLKSHPPPPEK',
            'MLLLLLLLLLLPPLVLRVAASRCLHDETQKSVSLLRPPFSQLPSKSRSSSLTLPSSRDPQPLRIQSCYLGDHISDGAWDPEGEGMRGGSRALAAVREATQRIQAVLAVQGPLLLSRDPAQYCHAVWGDPDSPNYHRCSLLNPGYKGESCLGAKIPDTHLRGYALWPEQGPPQLVQPDGPGVQNTDFLLYVRVAHTSKCHQETVSLCCPGWSTAAQSQLTAALTSWAQRRGFVMLPRLCLKLLGSSNLPTLASQSIRITGPSVIAYAACCQLDSEDRPLAGTIVYCAQHLTSPSLSHSDIVMATLHELLHALGFSGQLFKKWRDCPSGFSVRENCSTRQLVTRQDEWGQLLLTTPAVSLSLAKHLGVSGASLGVPLEEEEGLLSSHWEARLLQGSLMTATFDGAQRTRLDPITLAAFKDSGWYQVNHSAAEELLWGQGSGPEFGLVTTCGTGSSDFFCTGSGLGCHYLHLDKGSCSSDPMLEGCRMYKPLANGSECWKKENGFPAGVDNPHGEIYHPQSRCFFANLTSQLLPGDKPRHPSLTPHLKEAELMGRCYLHQCTGRGAYKVQVEGSPWVPCLPGKVIQIPGYYGLLFCPRGRLCQTNEDINAVTSPPVSLSTPDPLFQLSLELAGPPGHSLGKEQQEGLAEAVLEALASKGGTGRCYFHGPSITTSLVFTVHMWKSPGCQGPSVATLHKALTLTLQKKPLEVYHGGANFTTQPSKLLVTSDHNPSMTHLRLSMGLCLMLLILVGVMGTTAYQKRATLPVRPSASYHSPELHSTRVPVRGIREV'
        ]

        self.assertEqual(result, expected_result)

    def test_separate_target_loc_valid_input(self):
        # Test case with size
        filtermode = "loc"

        # example of a valid input
        out_stack = miner.seq_filter(self.referencefile, self.referencetype, filtermode, self.config, self.stats, self.outmode, self.outputfilename, self.outputfiletype)
        result_seqs = [record.sequence for record in out_stack[0].target]
        result = result_seqs[1:3]

        # debug block
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        print("Test mode: ", filtermode)
        print(result)
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

        expected_result = [
            'MGDQPCASGRSTLPPGNAREAKPPKKRCLLAPRWDYPEGTPNGGSTTLPSAPPPASAGLKSHPPPPEK',
            'MTAEDSTAAMSSDSAAGSSAKVPEGVAGAPNEAALLALMERTGYSMVQENGQRKYGGPPPGWEGPHPQRGCEVFVGKIPRDVYEDELVPVFEAVGRIYELRLMMDFDGKNRGYAFVMYCHKHEAKRAVRELNNYEIRPGRLLGVCCSVDNCRLFIGGIPKMKKREEILEEIAKVTEGVLDVIVYASAADKMKNRGFAFVEYESHRAAAMARRKLMPGRIQLWGHQIAVDWAEPEIDVDEDVMETVKILYVRNLMIETTEDTIKKSFGQFNPGCVERVKKIRDYAFVHFTSREDAVHAMNNLNGTELEGSCLEVTLAKPVDKEQYSRYQKAARGGGAAEAAQQPSYVYSCDPYTLAYYGYPYNALIGPNRDYFVKAGSIRGRGRGAAGNRAPGPRGSYLGGYSAGRGIYSRYHEGKGKQQEKGYELVPNLEIPTVNPVAIKPGTVAIPAIGAQYSMFPAAPAPKMIEDGKIHTVEHMISPIAVQPDPASAAAAAAAAAAAAAAVIPTVSTPPPFQGRPITPVYTVAPNVQRIPTAGIYGASYVPFAAPATATIATLQKNAAAAAAMYGGYAGYIPQAFPAAAIQVPIPDVYQTY'
        ]

        self.assertEqual(result, expected_result)

    def test_separate_target_size_no_matching_records(self):
        # Test case with size
        filtermode = "motif"

        # example of a valid input
        out_stack = miner.seq_filter(self.referencefile, self.referencetype, filtermode, self.config, self.stats, self.outmode, self.outputfilename, self.outputfiletype)
        result_seqs = [record.sequence for record in out_stack[0].target]
        result = result_seqs[0:2]

        # debug block
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        print("Test mode: ", filtermode)
        print(result)
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

        expected_result = [
            'MSGQLERCEREWHELEGEFQELQETHRIYKQKLEELAALQTLCSSSISKQKKHLKDLKLTLQRCKRHASREEAELVQQMAANIKERQDVFFDMEAYLPKKNGLYLNLVLGNVNVTLLSNQAKFAYKDEYEKFKLYLTIILLLGAVACRFVLHYRVTDEVFNFLLVWYYCTLTIRESILISNGSRIKGWWVSHHYVSTFLSGVMLTWPNGPIYQKFRNQFLAFSIFQSCVQFLQYYYQRGCLYRLRALGERNHLDLTVEGFQSWMWRGLTFLLPFLFCGHFWQLYNAVTLFELSSHEECREWQVFVLAFTFLILFLGNFLTTLKVVHAKLQKNRGKTKQP',
            'MSSGNYQQSEALSKPTFSEEQASALVESVFGLKVSKVRPLPSYDDQNFHVYVSKTKDGPTEYVLKISNTKASKNPDLIEVQNHIIMFLKAAGFPTASVCHTKGDNTASLVSVDSGSEIKSYLVRLLTYLPGRPIAELPVSPQLLYEIGKLAAKLDKTLQRFHHPKLSSLHRENFIWNLKNVPLLEKYLYALGQNRNREIVEHVIHLFKEEVMTKLSHFRECINHGDLNDHNILIESSKSASGNAEYQVSGILDFGDMSYGYYVFEVAITIMYMMIESKSPIQVGGHVLAGFESITPLTAVEKGALFLLVCSRFCQSLVMAAYSCQLYPENKDYLMVTAKTGWKHLQQMFDMGQKAVEEIWFETAKSYESGISM'
        ]

        self.assertEqual(result, expected_result)

    # Test to check the output from the terminal
    @patch('sys.stdout', new_callable=StringIO)  # Mock sys.stdout
    def test_separate_target_loc_cli(self, mock_stdout):
        # Test case with location
        self.test_separate_target_loc_valid_input()
        output = mock_stdout.getvalue()
        
        expected_result_handle = []
        expected_result_handle_stack = self.test_config['tests']
        for each_handle in expected_result_handle_stack:
            if each_handle['name'] == 'test_target_loc_cli':
                expected_result_handle.append(each_handle['expect_results'])

        if len(expected_result_handle) != 1:
            raise ValueError("Expected result handle not found ...")

        # debug block
        # sys.stdout.write("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        # for item in expected_result_handle:
        #     sys.stdout.write(item + "\n")
        # sys.stdout.write("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")

        with open(expected_result_handle[0], 'r') as file:
            expected_result = file.readlines()

        for expected in expected_result:
            with self.subTest(expected=expected): 
                self.assertIn(expected, output)


if __name__ == '__main__':
    unittest.main()