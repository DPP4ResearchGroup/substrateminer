##################################################
#            Test Description                    #
##################################################
# This file contains the unit tests for the      #
# pathfinder function.                           #
##################################################

import unittest
import click
from click.testing import CliRunner
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../')))
from substrate_miner.pathfinder.pathfinder\
    import pathfinder_main, validate_extension, validate_orgs


class TestPathfinder(unittest.TestCase):

    def setUp(self):
        self.runner = CliRunner()

    def test_validate_extension_valid(self):
        ctx = None
        param = None
        value = 'test.txt'
        self.assertEqual(validate_extension(ctx, param, value), value)

    def test_validate_extension_invalid(self):
        ctx = None
        param = None
        value = 'test.pdf'
        with self.assertRaises(click.BadParameter):
            validate_extension(ctx, param, value)

    def test_validate_orgs_with_api(self):
        ctx = type('', (), {})()  # Create a dummy object
        ctx.params = {'api': True}
        param = None
        value = None
        self.assertEqual(validate_orgs(ctx, param, value), value)

    def test_validate_orgs_without_api(self):
        ctx = type('', (), {})()  # Create a dummy object
        ctx.params = {'api': False}
        param = None
        value = 'hsa'
        self.assertEqual(validate_orgs(ctx, param, value), value)

    def test_validate_orgs_without_api_no_org(self):
        ctx = type('', (), {})()  # Create a dummy object
        ctx.params = {'api': False}
        param = None
        value = None
        with self.assertRaises(click.BadParameter):
            validate_orgs(ctx, param, value)

    def test_pathfinder_main_no_input(self):
        result = self.runner.invoke(pathfinder_main,\
            ['--uniprots', 'Q6XP50', '--api'])
        self.assertEqual(result.exit_code, 0)

    def test_pathfinder_main_no_uniprots(self):
        result = self.runner.invoke(pathfinder_main,\
            ['--input', 'unittests/data/uniprot_id_test.txt', '--api'])
        self.assertEqual(result.exit_code, 0)

    def test_pathfinder_main_no_input_no_uniprots(self):
        result = self.runner.invoke(pathfinder_main)
        self.assertEqual(result.exit_code, 1)
        self.assertIn("Error: Please provide either an input file or a UniProt ID.", result.output)

    def test_pathfinder_main_with_api(self):
        result = self.runner.invoke(pathfinder_main,\
            ['--uniprots', 'Q6XP50', '--api'])
        self.assertEqual(result.exit_code, 0)

    def test_pathfinder_main_with_output(self):
        result = self.runner.invoke(pathfinder_main,\
            ['--uniprots', 'Q6XP50', '--output', 'unittests/pathfinder_output.txt', '--api'])
        self.assertEqual(result.exit_code, 0)

    def test_pathfinder_main_with_output_no_uniprots(self):
        result = self.runner.invoke(pathfinder_main,\
            ['--input', 'unittests/data/uniprot_id_test.txt',\
            '--output', 'unittests/pathfinder_output.txt', '--api'])
        self.assertEqual(result.exit_code, 0)

if __name__ == '__main__':
    unittest.main()