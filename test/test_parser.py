import unittest
from rnannot.parser import parse_args


class ParserTestCase(unittest.TestCase):
    def test(self):
        self.assertRaises(Exception, parse_args, ['a', '--platform'])
        self.assertRaises(Exception, parse_args, ['a', '--name'])
        self.assertRaises(Exception, parse_args, ['a', 'b', 'c'])  # only one file or two files should be provided
        self.assertRaises(Exception, parse_args, ['a', 'b', '--sra'])  # only one SRA run file should be provided
