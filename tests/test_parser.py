import unittest
from rnannot.parser import parse_args


class ParserTestCase(unittest.TestCase):
    def test(self):
        try:
            parse_args(['-i', './example/197043.tsv', '-g' ,'../GCA_000696855.1_Hvit_1.0_genomic.fna.gz', '-d'])
        except:
            print('Parser erroneously parses the correct command.')
        try:
            parse_args(['-i', './example/197043.tsv', '-g' ,'../GCA_000696855.1_Hvit_1.0_genomic.fna.gz'])
        except:
            print('Parser erroneously parses the correct command.')
